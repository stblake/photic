"""Error budget for the residual depth bias on OSOAA scenes (EXPERIMENTS.md, step 3).

photic's per-pixel search (B5/B6: a profile over common depths, then a full refinement, by
samodel_lm through drivers/invert_lab's LM command) is run on tier 3's OSOAA pixels with one
component held at its true value at a time:

  none          everything fitted (Y fitted per scene)
  Y             the backscattering slope Y (true 1)
  water         P, G, X and Y
  bottom        the bottom albedos B_k (true: 0.6 of the pixel's bottom type)
  Y+bottom      Y and the bottom

The component whose truth removes the bias is the one the model gets wrong. Needs an
invert_lab built with Y_MODEL=3 (Y as a parameter), given by $PHOTIC_INVERT_LAB.

  error_budget.py [--noise 0] [water ...]
"""

import argparse
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import phaseA as A
import tier3_scenes as t3

H_GRID = [0.5, 1.0, 2.0, 3.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0]
N_WATER = 4  # ln P, ln G, ln X, Y (Y_MODEL 3)
NB = 3
CONFIGS = ["none", "Y", "water", "bottom", "Y+bottom"]
# Medians of the in situ classes (Heron Reef, Moreton Bay) at 560 nm, in photic's units.
REALISTIC = {"sand": 0.5, "seagrass": 0.07, "coral": 0.10}


def masks(nr, n, config, n_par):
    """Indices held at the truth (n_par from the packed truth, so any bottom layout works)."""
    fixed = np.zeros(n_par, dtype=bool)
    off = n_par - N_WATER*n
    if config in ("Y", "Y+bottom"):
        fixed[[off + N_WATER*k + 3 for k in range(n)]] = True
    if config == "water":
        fixed[off:] = True
    if config in ("bottom", "Y+bottom"):
        fixed[nr:off] = True
    return fixed


def search(setup, truth, nr, origin, fixed):
    """photic's search with the fixed parameters at the truth."""
    starts = A.lab(setup, ["START %g 0.3" % h for h in H_GRID])
    cmds = []
    for s in starts:
        v = np.array(s["START"])
        v[fixed] = truth[fixed]
        f = (~fixed).astype(float)
        f[:nr] = 0.0  # depth held on the grid
        cmds.append("LM %s %s" % (A.vec(v), A.vec(f)))
    prof = A.lab(setup, cmds)
    order = np.argsort([p["rrs_error"] for p in prof])[:2]
    cmds = []
    for k in order:
        v = np.array(prof[k]["params"])
        f = (~fixed).astype(float)
        cmds.append("LM %s %s" % (A.vec(v), A.vec(f)))
    ref = A.lab(setup, cmds)
    best = min(ref, key=lambda r: r["rrs_error"])
    return A.depth_of(best, nr, origin), best["rrs_error"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("--realistic", action="store_true", help="realistic albedo per bottom (B8)")
    ap.add_argument("water", nargs="*", default=["clear", "moderate", "turbid"])
    a = ap.parse_args()
    if a.realistic:
        t3.ALBEDO_BY_BOTTOM = REALISTIC
    for w in a.water:
        H, b, g = A.scene_grids(w, a.scenes, True, a.noise)
        pix = [(i, j) for i in range(0, 32, 2) for j in (8, 20, 28)]

        def one(ij):
            i, j = ij
            setup, truth, nr, origin, phys = A.pixel_setup(w, a.scenes, i, j, H, b, g, 2, noise=a.noise)
            return [H[i, j]] + [search(setup, truth, nr, origin, masks(nr, a.scenes, c, len(truth)))[0] for c in CONFIGS]

        with ThreadPoolExecutor(12) as ex:
            rows = np.array(list(ex.map(one, pix)))
        print("%s, %d scenes, OSOAA%s, noise %g: depth bias (RMSE) by depth bin, m" % (w, a.scenes, ", realistic bottoms" if a.realistic else "", a.noise))
        print("  %-10s" % "held true" + "".join("%15s" % b for b in ("0-5", "5-10", "10-15", "15-20", "20-30")))
        for k, c in enumerate(CONFIGS):
            cells = []
            for lo, hi in ((0, 5), (5, 10), (10, 15), (15, 20), (20, 30)):
                m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
                e = rows[m, k + 1] - rows[m, 0]
                cells.append("%+6.2f (%5.2f)" % (e.mean(), np.sqrt((e**2).mean())))
            print("  %-10s" % c + "".join("%15s" % x for x in cells), flush=True)


if __name__ == "__main__":
    main()
