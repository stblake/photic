"""A physics-based first estimate of depth, H0, per pixel (EXPERIMENTS.md, initial depth).

From photic's forward model, the bottom signal in band i is
    D_i = rrs_i - rrs_dp,i = (rho_i/pi - rrs_dp,i) exp(-M kappa_i H),   M = sec theta_w + D_B.
For a red band (665 nm), where water absorption dominates and is known, and a green band
(560 nm), the bottom brightness cancels in the ratio:
    H0 = [ln(D_560/D_665) - ln(r_560/r_665)] / (M (kappa_665 - kappa_560)),
with r_i = rho_i/pi - rrs_dp,i ~ rho_i/pi and the bottom's spectral ratio rho_560/rho_665 from
photic's table (averaged over its types). kappa and rrs_dp come from photic's own starting
water (P, G, X from the scene's Rrs ratios). Where the red bottom signal is lost in the noise
(D_665 below 3 sigma), H0 falls back to the depth at which the brightest plausible bottom
(albedo 0.75) would still give D_560: H0 = ln((0.75 rho_560/pi - rrs_dp)/D_560)/(M kappa_560).
Each scene gives H0 + tide; the scenes are averaged after removing the tides.

  h0_estimate.py [--realistic] [--noise 0.0005]     accuracy of H0 against the truth by depth
"""

import argparse

import numpy as np

import phaseA as A
import photic_ref
import tier2_osoaa as t2
import tier3_scenes as t3

REALISTIC = {"sand": 0.5, "seagrass": 0.07, "coral": 0.10}
RHO_MAX = 0.75


def water_start(R):
    """photic's starting P, G, X from one scene's Rrs at the four bands (samodel_initial)."""
    lam = np.array(t3.BANDS)
    r490, r550, r640 = (np.interp(l, lam, R) for l in (490.0, 550.0, 640.0))
    P = 0.072*(r490/r550)**-1.7
    G = 1.5*P
    aw640 = 0.329  # Pope & Fry at 640 nm (photic's table interpolated)
    X = 30.0*aw640*r640
    return np.clip(P, 5e-4, 1.0), np.clip(G, 5e-4, 2.0), np.clip(X, 1e-4, 0.5)


def water(P, G, X, th, tables):
    out = {}
    for lam in (560.0, 665.0):
        t = tables[(lam, 0)]
        r = photic_ref.forward(lam, P, G, X, 0.0, t["rho"], 1.0e3, 0.0, th, 0.0, 0.01, 0.01, 0.0,
                               t["a0"], t["a1"], t["aw"], t["bbw"], Y=1.0)
        out[lam] = (float(r["rrs_dp"]), float(r["K"]), float(r["u"]))
    return out


def h0(R_scenes, tides, thetas, tables, sigma, shape_ratio, waters=None):
    """H0 from each scene's Rrs at the bands; returns (H0, which estimator). waters: per scene,
    {lam: (rrs_dp, kappa, u)} for 560 and 665 nm; by default from photic's starting water."""
    est, kinds = [], []
    for k, (R, tide, th) in enumerate(zip(R_scenes, tides, thetas)):
        rrs = R/(0.5 + 1.5*R)
        out = waters[k] if waters else water(*water_start(R), th, tables)
        i560, i665 = t3.BANDS.index(560.0), t3.BANDS.index(665.0)
        D560 = rrs[i560] - out[560.0][0]
        D665 = rrs[i665] - out[665.0][0]
        u = out[560.0][2]
        M = 1.0/np.cos(np.radians(th)) + 1.04*np.sqrt(1 + 5.4*u)
        sig = sigma/0.5 if sigma > 0 else 2e-5  # rrs noise
        if D665 > 3*sig and D560 > 0:
            dk = out[665.0][1] - out[560.0][1]
            Hs = (np.log(D560/D665) - np.log(shape_ratio))/(M*dk)
            kinds.append("ratio")
        elif D560 > 0:
            t = tables[(560.0, 0)]
            Hs = np.log(max(RHO_MAX*t["rho"]/np.pi - out[560.0][0], 1e-6)/D560)/(M*out[560.0][1])
            kinds.append("bound")
        else:
            Hs = 30.0
            kinds.append("deep")
        est.append(np.clip(Hs, 0.1, 40.0) - tide)
    return float(np.clip(np.mean(est), 0.05, 40.0)), kinds[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--realistic", action="store_true")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("--water-ref", choices=["start", "true", "image"], default="start",
                    help="where rrs_dp and kappa come from: photic's starting water, the true water, "
                         "or the deepest rows of the image (rrs_dp) with kappa from the true water")
    ap.add_argument("water", nargs="*", default=["clear", "moderate", "turbid"])
    a = ap.parse_args()
    if a.realistic:
        t3.ALBEDO_BY_BOTTOM = REALISTIC
    tables = t2.photic_tables(A.DRIVER, t3.BANDS)
    ratios = [tables[(560.0, b)]["rho"]/tables[(665.0, b)]["rho"] for b in t2.BOTTOMS.values()]
    shape_ratio = float(np.exp(np.mean(np.log(ratios))))
    for w in a.water:
        H, bottom, grids = A.scene_grids(w, a.scenes, True, a.noise)
        sets = t3.SCENE_SETS[w][:a.scenes]
        tides = [s[5] for s in sets]
        thetas = [t3.theta_w(s[4]) for s in sets]
        waters = None
        if a.water_ref in ("true", "image"):
            waters = [water(*s[:3], th, tables) for s, th in zip(sets, thetas)]
        if a.water_ref == "image":
            # rrs_dp from the deepest rows (H >= 22 m), the image's own deep-water signal.
            for k in range(a.scenes):
                for lam in (560.0, 665.0):
                    Rd = np.median(grids[k][t3.BANDS.index(lam)][28:, :])
                    waters[k][lam] = (Rd/(0.5 + 1.5*Rd),) + waters[k][lam][1:]
        rows = []
        for i in range(0, 32, 2):
            for j in (8, 20, 28):
                R = [np.array([g[m][i, j] for m in range(len(t3.BANDS))]) for g in grids]
                est, kind = h0(R, tides, thetas, tables, a.noise, shape_ratio, waters)
                rows.append((H[i, j], est, kind))
        print("%s, %d scenes, OSOAA%s, noise %g, water from %s: H0 against the truth" % (
            w, a.scenes, ", realistic bottoms" if a.realistic else "", a.noise, a.water_ref))
        for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 15), (15, 20), (20, 30)):
            m = [r for r in rows if lo <= r[0] < hi]
            e = np.array([r[1] - r[0] for r in m])
            kinds = {k: sum(1 for r in m if r[2] == k) for k in ("ratio", "bound", "deep")}
            print("   %2d-%2d bias %+6.2f rmse %6.2f within tol %.2f  (%s)" % (
                lo, hi, e.mean(), np.sqrt((e**2).mean()),
                np.mean([abs(r[1] - r[0]) <= 0.5 + 0.1*r[0] for r in m]),
                ", ".join("%s %d" % kv for kv in kinds.items() if kv[1])))


if __name__ == "__main__":
    main()
