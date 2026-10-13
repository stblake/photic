"""B13 of EXPERIMENTS.md: how much can sun-angle (and tide) differences between scenes tell
us about depth?

Cramer-Rao lower bound on the depth of one pixel from photic's forward model
(photic_ref.forward, Y fixed at 1), with Gaussian Rrs noise sigma in every band and scene.
The parameters are depth H, the bottom albedo B (one bottom type) and the water (ln P, ln G,
ln X), either free per scene (as photic fits it) or shared by all scenes (the limit of tying
the water, B12). The Fisher matrix is J'J/sigma^2 with J by central differences; the bound
on H is sqrt of the H diagonal of its inverse.

  b13_crb.py
"""

import numpy as np

import photic_ref
import tier2_osoaa as t2
import tier3_scenes as t3
import phaseA as A

BANDS = [443.0, 490.0, 560.0, 665.0]
SIGMA = 5.0e-4
B = 0.6


VIEWS = None  # subsurface view zenith per scene (deg); None means nadir


def rrs_vector(theta, H, B_alb, water, tides, thetas, tables, bottom=0):
    """Rrs for every scene and band. water: list of (P, G, X) per scene."""
    out = []
    views = VIEWS or [0.0]*len(water)
    for (P, G, X), tide, th, tv in zip(water, tides, thetas, views):
        for lam in BANDS:
            t = tables[(lam, bottom)]
            out.append(photic_ref.forward(lam, P, G, X, B_alb, t["rho"], H, tide, th, tv, 0.01, 0.01, 0.0,
                                          t["a0"], t["a1"], t["aw"], t["bbw"], Y=1.0)["Rrs"])
    return np.array(out, dtype=float)


def crb(H, water, tides, thetas, tables, shared):
    n = len(water)
    if shared:
        x0 = np.array([H, B] + list(np.log(water[0])))
        def f(x):
            w = [tuple(np.exp(x[2:5]))]*n
            return rrs_vector(None, x[0], x[1], w, tides, thetas, tables)
    else:
        x0 = np.array([H, B] + [v for w in water for v in np.log(w)])
        def f(x):
            w = [tuple(np.exp(x[2 + 3*k:5 + 3*k])) for k in range(n)]
            return rrs_vector(None, x[0], x[1], w, tides, thetas, tables)
    J = []
    for i in range(len(x0)):
        h = 1e-5*max(abs(x0[i]), 1.0)
        xp, xm = x0.copy(), x0.copy()
        xp[i] += h
        xm[i] -= h
        J.append((f(xp) - f(xm))/(2*h))
    J = np.array(J).T
    F = J.T @ J/SIGMA**2
    # Depth is unidentifiable if the Fisher matrix is (numerically) singular.
    if np.linalg.cond(F) > 1e14:
        return np.inf
    return float(np.sqrt(max(np.linalg.inv(F)[0, 0], 0.0)))


def main():
    tables = t2.photic_tables(A.DRIVER, BANDS)
    waters = {"clear": (0.02, 0.02, 0.004), "moderate": (0.08, 0.15, 0.02), "turbid": (0.30, 0.50, 0.08)}
    configs = [
        ("2 scenes, same sun 30, view 0/12 deg, no tide", [30.0, 30.0], [0.0, 0.0], [0.0, 12.0]),
        ("2 scenes, same sun 30, view 0/0, no tide", [30.0, 30.0], [0.0, 0.0], [0.0, 0.0]),
        ("2 scenes, tier 3 angles (30, 45 deg), tides 0/0.4", [30.0, 45.0], [0.0, 0.4]),
        ("2 scenes, same angle (30, 30), tides 0/0.4", [30.0, 30.0], [0.0, 0.4]),
        ("2 scenes, tier 3 angles, no tide difference", [30.0, 45.0], [0.0, 0.0]),
        ("2 scenes, wide angles (20, 60 deg), tides 0/0.4", [20.0, 60.0], [0.0, 0.4]),
        ("4 scenes, tier 3 angles (25-50), tides -0.3..0.8", [30.0, 45.0, 25.0, 50.0], [0.0, 0.4, -0.3, 0.8]),
    ]
    depths = [2.0, 5.0, 10.0, 15.0, 20.0]
    print("Cramer-Rao bound on depth (m), noise %g per band; water free per scene / shared" % SIGMA)
    for wname, w in waters.items():
        print(wname)
        global VIEWS
        for cfg in configs:
            label, angles, tides = cfg[:3]
            VIEWS = [t3.theta_w(v) for v in cfg[3]] if len(cfg) > 3 else None
            thetas = [t3.theta_w(a) for a in angles]
            water = [w]*len(angles)
            cells = []
            for H in depths:
                free = crb(H, water, tides, thetas, tables, shared=False)
                tied = crb(H, water, tides, thetas, tables, shared=True)
                cells.append("%5.2f/%5.2f" % (min(free, 99), min(tied, 99)))
            print("  %-48s %s" % (label, "  ".join("H=%g:%s" % (H, c) for H, c in zip(depths, cells))))


if __name__ == "__main__":
    main()
