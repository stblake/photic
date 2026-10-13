"""Decoupling the water from the bottom and depth (EXPERIMENTS.md, step 1-2).

photic's search (a profile over common depths, then a refinement of the two best; LM through
invert_lab) on tier 3's OSOAA pixels in three ways:

  joint     everything fitted per pixel window, as photic does
  true      the water (ln P, ln G, ln X, Y per scene) held at the truth: the ceiling
  deep      the water of each scene fitted first to its deep-water Rrs (photic's own model,
            no bottom, depth at H_MAX), then held fixed for the shallow pixels

Needs an invert_lab built with Y_MODEL=3 (Y as a parameter) and BOTTOM_BOUNDS=1, given by
$PHOTIC_INVERT_LAB. The deep-water Rrs of a scene is the OSOAA lookup table with no bottom at
its deepest node (32 m).

  decouple.py [--realistic] [--noise 0.0005] [water ...]
"""

import argparse
import json
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import error_budget as EB
import phaseA as A
import tier3_scenes as t3

N_WATER = 4
H_MAX = 60.0


def deep_rrs(water, n):
    lut_dir = os.path.join(A.DATA, "runs", "tier3", t3.LUT_DIR)
    out = []
    for k in range(n):
        lut = json.load(open(os.path.join(lut_dir, t3.scene_key(water, k) + ".json")))
        out.append([lut["%g" % l]["R0"][-1] for l in t3.BANDS])
    return out


def deep_water(water, k, R):
    """Water parameters of scene k fitted to its deep-water Rrs (one region, no bottom)."""
    P, G, X, Y, theta_air, tide = t3.SCENE_SETS[water][k]
    lines = ["SETUP 1 %d 1 0 3" % len(t3.BANDS),
             "SCENE 0 %.10g 0 %.10g %s" % (t3.theta_w(theta_air), tide, " ".join("%g" % l for l in t3.BANDS)),
             "RRS 0 0 %s" % " ".join("%.12g" % v for v in R), "END"]
    setup = "\n".join(lines) + "\n"
    s = A.lab(setup, ["START 30 0.3"])[0]
    v = np.array(s["START"])
    f = np.zeros(len(v))
    v[0] = H_MAX
    f[-N_WATER:] = 1.0  # the water only
    r = A.lab(setup, ["LM %s %s" % (A.vec(v), A.vec(f))])[0]
    return np.array(r["params"][-N_WATER:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--realistic", action="store_true")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("water", nargs="*", default=["clear", "moderate", "turbid"])
    a = ap.parse_args()
    if a.realistic:
        t3.ALBEDO_BY_BOTTOM = EB.REALISTIC
    for w in a.water:
        H, b, g = A.scene_grids(w, a.scenes, True, a.noise)
        deep = np.concatenate([deep_water(w, k, R) for k, R in enumerate(deep_rrs(w, a.scenes))])
        truth_w = np.concatenate([[np.log(s[0]), np.log(s[1]), np.log(s[2]), 1.0] for s in t3.SCENE_SETS[w][:a.scenes]])
        pix = [(i, j) for i in range(0, 32, 2) for j in (8, 20, 28)]

        def one(ij):
            i, j = ij
            setup, truth, nr, origin, phys = A.pixel_setup(w, a.scenes, i, j, H, b, g, 2, noise=a.noise)
            n_par = len(truth)
            free_all = np.zeros(n_par, dtype=bool)
            water_fixed = free_all.copy()
            water_fixed[n_par - N_WATER*a.scenes:] = True
            out = [H[i, j]]
            out.append(EB.search(setup, truth, nr, origin, free_all)[0])          # joint
            out.append(EB.search(setup, truth, nr, origin, water_fixed)[0])       # water true
            t2 = truth.copy()
            t2[n_par - N_WATER*a.scenes:] = deep
            out.append(EB.search(setup, t2, nr, origin, water_fixed)[0])          # water from deep
            return out

        with ThreadPoolExecutor(12) as ex:
            rows = np.array(list(ex.map(one, pix)))
        dw = deep.reshape(a.scenes, N_WATER)
        tw = truth_w.reshape(a.scenes, N_WATER)
        print("%s, %d scenes, OSOAA%s, noise %g" % (w, a.scenes, ", realistic bottoms" if a.realistic else "", a.noise))
        print("   deep-water fit vs truth, per scene: P %s  G %s  X %s  Y %s" % tuple(
            "/".join("%.3g:%.3g" % (np.exp(dw[k, c]) if c < 3 else dw[k, c], np.exp(tw[k, c]) if c < 3 else tw[k, c])
                     for k in range(a.scenes)) for c in range(4)))
        print("   %-8s" % "" + "".join("%17s" % b for b in ("0-2", "2-5", "5-10", "10-15", "15-20", "20-30")))
        for k, name in enumerate(("joint", "true", "deep")):
            cells = []
            for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 15), (15, 20), (20, 30)):
                m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
                e = rows[m, k + 1] - rows[m, 0]
                cells.append("%+6.2f (%5.2f)" % (e.mean(), np.sqrt((e**2).mean())))
            print("   %-8s" % name + "".join("%17s" % c for c in cells), flush=True)


if __name__ == "__main__":
    main()
