"""Tier 4: an empirical depth prior fitted on the training sites' soundings only.

EXPERIMENTS.md Phase D. For each pixel, the median Rrs over a site's scenes (each scene's
inputs as photic gets them) gives the features; the model is fitted on the training sites
(Marina and the gulf entrance) and applied to every site, so at Bundegi and Tantabiddi it is
out of sample. Two forms:
  ratio   Stumpf et al. (2003): H = m1 ln(1000 Rrs_blue)/ln(1000 Rrs_green) - m0;
  linear  Lyzenga-style: H = a0 + sum_i a_i ln(Rrs_i) over the four bands.
The better form by leave-one-site-out between the two training sites is kept. Writes
derived/tier4/<sites>/<site>/prior.npy (depth, positive down, NaN where no estimate) and
prior.json (form, coefficients, scores).

  tier4_prior.py [sites]      (default sites; e.g. sites5, sites16)
"""

import glob
import json
import os
import sys

import numpy as np

DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))
TRAIN = ["marina", "gulf_entrance"]
BANDS = [443, 492, 560, 665]
MAX_FIT_DEPTH = 15.0  # beyond this the signal is mostly lost; fit where it is informative


def read_asc(path):
    with open(path) as f:
        hdr = [next(f) for _ in range(6)]
    nodata = float(hdr[5].split()[1])
    a = np.loadtxt(path, skiprows=6)[::-1]  # photic's row 0 is south
    a[a == nodata] = np.nan
    return a


def features(site_dir):
    """Median Rrs over the site's scenes, per band."""
    med = {}
    for b in BANDS:
        files = sorted(glob.glob(os.path.join(site_dir, "s*_b%d.asc" % b)))
        stack = np.array([read_asc(f) for f in files])
        med[b] = np.nanmedian(stack, axis=0)
    return med


def design(med, form):
    if form == "ratio":
        with np.errstate(invalid="ignore", divide="ignore"):
            x = np.log(1000.0*med[492])/np.log(1000.0*med[560])
        return np.stack([np.ones_like(x), x], axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.stack([np.ones_like(med[443])] + [np.log(med[b]) for b in BANDS], axis=-1)


def fit(Xs, ys):
    X = np.concatenate(Xs)
    y = np.concatenate(ys)
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y) & (y < MAX_FIT_DEPTH)
    coef, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    return coef


def predict(X, coef):
    H = X @ coef
    H[~np.all(np.isfinite(X), axis=-1)] = np.nan
    return np.clip(H, 0.1, 40.0)


def score(H, ref):
    ok = np.isfinite(H) & np.isfinite(ref)
    e = H[ok] - ref[ok]
    tol = 0.5 + 0.1*ref[ok]
    return dict(n=int(ok.sum()), rmse=float(np.sqrt(np.mean(e**2))), bias=float(e.mean()),
                within_tol=float(np.mean(np.abs(e) <= tol)))


def main(sites_set="sites"):
    root = os.path.join(DATA, "derived", "tier4", sites_set)
    sites = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    meds = {s: features(os.path.join(root, s)) for s in sites}
    refs = {s: np.load(os.path.join(root, s, "reference.npz"))["depth"] for s in sites}
    report = {}
    for form in ("ratio", "linear"):
        # Leave one training site out: fit on the other, score on it.
        cv = {}
        for held in TRAIN:
            other = [s for s in TRAIN if s != held]
            coef = fit([design(meds[s], form).reshape(-1, design(meds[s], form).shape[-1]) for s in other],
                       [refs[s].ravel() for s in other])
            cv[held] = score(predict(design(meds[held], form), coef), refs[held])
        report[form] = dict(cv=cv, cv_rmse=float(np.mean([v["rmse"] for v in cv.values()])))
    form = min(report, key=lambda f: report[f]["cv_rmse"])
    coef = fit([design(meds[s], form).reshape(-1, design(meds[s], form).shape[-1]) for s in TRAIN],
               [refs[s].ravel() for s in TRAIN])
    out = dict(form=form, coef=coef.tolist(), train=TRAIN, max_fit_depth=MAX_FIT_DEPTH, forms=report, sites={})
    for s in sites:
        H = predict(design(meds[s], form), coef)
        np.save(os.path.join(root, s, "prior.npy"), H)
        out["sites"][s] = score(H, refs[s])
    json.dump(out, open(os.path.join(root, "prior.json"), "w"), indent=1)
    for f, r in report.items():
        print("%-6s leave-one-site-out RMSE %.2f  %s" % (f, r["cv_rmse"], {k: round(v["rmse"], 2) for k, v in r["cv"].items()}))
    print("chosen: %s, coefficients %s" % (form, np.round(coef, 3).tolist()))
    for s, v in out["sites"].items():
        print("  %-14s %-9s n %5d  RMSE %5.2f  bias %+5.2f  within tol %3.0f%%" % (
            s, "train" if s in TRAIN else "held out", v["n"], v["rmse"], v["bias"], 100*v["within_tol"]))


if __name__ == "__main__":
    main(*(sys.argv[1:] or []))
