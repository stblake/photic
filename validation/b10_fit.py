"""B10 of EXPERIMENTS.md: refit Lee's shallow-water coefficients to exact radiative transfer.

photic's forward model (Lee et al. 1998, 1999):

  rrs_dp = (g0 + g1 u) u
  rrs    = rrs_dp [1 - exp(-(sec theta_w + D_C) K H)] + rho/pi exp(-(sec theta_w + D_B) K H)
  D_C = c0 sqrt(1 + c1 u),  D_B = b0 sqrt(1 + b1 u)                 (nadir view)

with Lee's (g0, g1, c0, c1, b0, b1) = (0.084, 0.170, 1.03, 2.4, 1.04, 5.4), fitted to
Hydrolight. Here they are refitted to the tier 2 OSOAA cases (subsurface rrs), minimising
the relative error. The fit uses clear and moderate water and is tested on turbid water,
using only physical cases (bottom reflectance <= 1).

  b10_fit.py   prints Lee's and the refitted errors on both splits, and the coefficients.
"""

import glob
import json
import os

import numpy as np
from scipy.optimize import least_squares

DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))
LEE = np.array([0.084, 0.170, 1.03, 2.4, 1.04, 5.4])


def load():
    rows = [json.load(open(f)) for f in glob.glob(os.path.join(DATA, "runs", "tier2", "cases", "*", "result.json"))]
    rows = [r for r in rows if r.get("rrs_osoaa") and r["bottom_albedo"] <= 1.0]
    return rows


def model(c, r):
    g0, g1, c0, c1, b0, b1 = c
    u, K, H = r["u"], r["K"], r["H"]
    sec = 1.0/np.cos(np.radians(r["theta_w"]))
    rho = r["bottom_albedo"]
    rrs_dp = (g0 + g1*u)*u
    DC = c0*np.sqrt(1.0 + c1*u)
    DB = b0*np.sqrt(1.0 + b1*u)
    H = np.where(np.isfinite(H), H, 1.0e6)
    return rrs_dp*(1.0 - np.exp(-(sec + DC)*K*H)) + rho/np.pi*np.exp(-(sec + DB)*K*H)


def table(rows):
    keys = ["u", "K", "H", "theta_w", "bottom_albedo", "rrs_osoaa"]
    return {k: np.array([float(r[k]) if r[k] is not None else np.inf for r in rows]) for k in keys}


def errors(c, t):
    return model(c, t)/t["rrs_osoaa"] - 1.0


def summary(label, c, t):
    e = np.abs(errors(c, t))
    return "%-28s median %.2f%%  p95 %.2f%%  max %.2f%%  bias %+.2f%%" % (
        label, 100*np.median(e), 100*np.percentile(e, 95), 100*e.max(), 100*np.mean(errors(c, t)))


def main():
    rows = load()
    train = table([r for r in rows if r["water"] in ("clear", "moderate")])
    test = table([r for r in rows if r["water"] == "turbid"])
    fit = least_squares(lambda c: errors(c, train), LEE, x_scale=LEE, bounds=(LEE*0.3, LEE*3.0))
    print("cases: %d train (clear, moderate), %d test (turbid)" % (len(train["u"]), len(test["u"])))
    print(summary("Lee, train", LEE, train))
    print(summary("refit, train", fit.x, train))
    print(summary("Lee, test (turbid)", LEE, test))
    print(summary("refit, test (turbid)", fit.x, test))
    print("coefficients g0 g1 c0 c1 b0 b1:")
    print("  Lee   " + " ".join("%.4f" % v for v in LEE))
    print("  refit " + " ".join("%.4f" % v for v in fit.x))
    # Where the error is: by kappa*H.
    for lo, hi in ((0, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 8), (8, np.inf)):
        for label, c in (("Lee", LEE), ("refit", fit.x)):
            kh = train["K"]*np.where(np.isfinite(train["H"]), train["H"], 1e6)
            m = (kh >= lo) & (kh < hi)
            e = errors(c, {k: v[m] for k, v in train.items()})
            print("  kappaH %4.1f-%-4.1f %-5s bias %+.2f%% median |e| %.2f%% (n %d)" % (
                lo, hi, label, 100*e.mean(), 100*np.median(np.abs(e)), m.sum()))
    json.dump(dict(lee=LEE.tolist(), refit=fit.x.tolist()), open(os.path.join(DATA, "runs", "b10_coefficients.json"), "w"))


if __name__ == "__main__":
    main()
