"""B7 of EXPERIMENTS.md: does a red-edge band (Sentinel-2 B5, 705 nm) anchor backscattering?

Where water absorbs strongly (a_w(705) ~ 0.7 m^-1) the bottom is invisible below a few metres
and rrs fixes b_b. This compares photic on the same OSOAA scenes with the four usual bands
(443-665 nm) and with 705 nm added: same pixels, same noise realisation (paired).

The scenes are tier 3's, but with realistic albedos per bottom (sand 0.6, seagrass 0.2,
coral 0.25), since tier 3's 0.6 for seagrass would make its reflectance at 705 nm about 2.
The OSOAA lookup tables are extended to 705 nm and cached with tier 3's.

  b7_nir.py [--water moderate] [--noise 0.0005]
"""

import argparse
import json
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import phaseA as A
import tier2_osoaa as t2
import tier3_scenes as t3

BANDS5 = [443.0, 490.0, 560.0, 665.0, 705.0]
ALBEDO = {"sand": 0.6, "seagrass": 0.2, "coral": 0.25}


def osoaa_config(work):
    data = A.DATA
    osoaa = dict(root=os.path.join(data, "tools", "OSOAA_scalar"),
                 exe=os.path.join(data, "tools", "OSOAA_scalar", "exe", "OSOAA_MAIN.exe"),
                 black_ocean=os.path.join(data, "tools", "OSOAA_blackocean", "exe", "OSOAA_MAIN.exe"),
                 db=os.path.join(work, "db"), phase=t2.phase_table())
    for d in ["MIE_AER", "MIE_HYD", "SURF"]:
        os.makedirs(os.path.join(osoaa["db"], d), exist_ok=True)
    return osoaa


def ensure_band(water, k, lam, tables, osoaa, work, jobs=12):
    """The cached tier 3 lookup table for a scene, extended to band lam if missing."""
    lut = t3.build_lut(water, k, tables, osoaa, work, jobs)
    key = "%g" % lam
    if key in lut:
        return lut
    P, G, X, Y, theta_air, tide = t3.SCENE_SETS[water][k]
    name = t3.scene_key(water, k)
    t2.WATER_TYPES[name] = (P, G, X, Y)
    cases = [(H, alb) for H in t3.DEPTH_NODES for alb in t3.ALBEDO_NODES]

    def one(args):
        H, alb = args
        c = dict(water=name, H=H, bottom="sand", theta_air=theta_air, lam=lam)
        ph, _ = t2.photic_side(c, tables)
        wd = os.path.join(work, "osoaa", name, "%g_%g_%g" % (lam, H, alb))
        black = osoaa["black_ocean"] if (H == t3.DEPTH_NODES[0] and alb == 0.0) else None
        return args, t2.run_osoaa(c, ph, osoaa["exe"], osoaa["root"], osoaa["phase"], wd, osoaa["db"],
                                  exe_black_ocean=black, albedo=alb)

    results = [one(cases[0])]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results += list(pool.map(one, cases[1:]))
    surface = results[0][1]["Rrs_surface_osoaa"]
    R0, Aa, C = [], [], []
    for H in t3.DEPTH_NODES:
        v = {a: r["Rrs_osoaa"] + r["Rrs_surface_osoaa"] - surface for (h, a), r in results if h == H}
        a, b = v[0.5] - v[0.0], v[1.0] - v[0.0]
        c = (b - 2.0*a)/(b - a) if abs(b - a) > 1e-15 else 0.0
        R0.append(v[0.0])
        Aa.append(b*(1.0 - c))
        C.append(c)
    lut[key] = dict(R0=R0, A=Aa, C=C)
    json.dump(lut, open(os.path.join(work, "lut", name + ".json"), "w"))
    return lut


def scenes(water, n, tables, luts, noise, seed=12345):
    H, bottom = t3.truth()
    rng = np.random.default_rng(seed)
    grids = []
    for k in range(n):
        P, G, X, Y, theta_air, tide = t3.SCENE_SETS[water][k]
        bands = []
        for lam in BANDS5:
            rho = np.zeros(H.shape)
            for name, b in t2.BOTTOMS.items():
                rho[bottom == name] = ALBEDO[name]*tables[(lam, b)]["rho"]
            g = t3.lut_rrs(luts[k], lam, H + tide, rho)
            bands.append(g + (rng.normal(0.0, noise, H.shape) if noise > 0 else 0.0))
        grids.append(bands)
    return H, bottom, grids


def run(water, n, H, bottom, grids, noise, nbands):
    t3.BANDS = BANDS5[:nbands]
    sub = [[g[m] for m in range(nbands)] for g in grids]

    def one(a):
        i, j = a
        setup, truth, nr, org, phys = A.pixel_setup(water, n, i, j, H, bottom, sub, 2, noise=noise)
        r = A.lab(setup, ["PHOTIC -1"])[0]
        return H[i, j], r["depth"], r["deep"], r["albedo"]

    pix = [(i, j) for i in range(0, 32, 2) for j in (8, 20, 28)]
    with ThreadPoolExecutor(12) as ex:
        return np.array(list(ex.map(one, pix)))


def report(label, rows):
    ok = rows[:, 2] < 0.5
    print("%s: depth given at %d/%d" % (label, ok.sum(), len(rows)))
    for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 15), (15, 20), (20, 30)):
        m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
        mk = m & ok
        e = rows[mk, 1] - rows[mk, 0]
        s = ("bias %6.2f rmse %6.2f tol %.2f" % (e.mean(), np.sqrt((e**2).mean()),
             np.mean(np.abs(e) <= 0.5 + 0.1*rows[mk, 0]))) if mk.any() else "-"
        print("   %2d-%2d coverage %2d/%2d  %s" % (lo, hi, mk.sum(), m.sum(), s))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--water", default="moderate")
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("--noise", type=float, default=0.0005)
    a = ap.parse_args()
    work = os.path.join(A.DATA, "runs", "tier3")
    tables = t2.photic_tables(A.DRIVER, BANDS5)
    osoaa = osoaa_config(work)
    luts = [ensure_band(a.water, k, 705.0, tables, osoaa, work) for k in range(a.scenes)]
    H, bottom, grids = scenes(a.water, a.scenes, tables, luts, a.noise)
    for nb in (4, 5):
        report("%s, %d bands, noise %g" % (a.water, nb, a.noise), run(a.water, a.scenes, H, bottom, grids, a.noise, nb))


if __name__ == "__main__":
    main()
