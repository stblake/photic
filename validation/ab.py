"""Quick paired A/B test of the SA model on a few hundred out-of-sample points.

    ab.py A B [--n 300] [--seed 1] [--split test] [--jobs 8]

A and B are variants of photic: a git ref (commit, branch, tag) or "." for the working tree,
optionally followed by ":" and extra compiler flags, so a compile-time switch can be tested
without editing the code:

    ab.py HEAD .                     the working tree against the last commit
    ab.py . .:-DLUT=0                the lookup-table shortcut off
    ab.py b7-physics .:-DSEED_MODE=2

Each variant is built once (cached by commit, or by the working tree's contents, and flags)
under $PHOTIC_VALIDATION_DATA/ab/build/. Its photic must support MODEL ... POINTS.

The points are a fixed, depth-stratified sample (by --seed) from the out-of-sample split of
EXPERIMENTS.md's rules: --split test (the default) draws from tier 3's turbid water and tier
4's Bundegi and Tantabiddi; --split train from tier 3's clear and moderate water and tier 4's
Marina and gulf entrance. Half the points come from each tier, shared equally between its
cases, and within a case equally between depth bins as far as they allow.

Both variants run photic on the full scenes with POINTS, so only the sampled cells are
modelled: with the spatial window, scene setup and deep-water handling of a full run, but
without the lookup table and carry-over that neighbouring cells give in a full run, so results
can differ slightly (a few cm) from a full run's at the same cells.

Reports, per tier and overall: coverage (cells given a depth); on the cells both variants give
a depth, RMSE, median |error|, the share within tolerance (0.5 m + 10%), B's wins and the
paired difference in RMSE (B - A) with a 95% bootstrap interval, and a verdict (a change
counts only if the interval excludes zero and it is at least --min-effect of A's RMSE); RMSE by depth; the depth
sigma's calibration (share of |error| <= sigma, ideally 68%); and run time. Results go to
$PHOTIC_VALIDATION_DATA/ab/results/. The full benchmark (bench.py) remains the acceptance test.
"""

import argparse
import datetime
import glob
import hashlib
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tier2_osoaa as t2
import tier3_scenes as t3

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))
AB = os.path.join(DATA, "ab")
DRIVER = os.path.join(HERE, "drivers", "forward_dump")

SPLITS = {
    "test": dict(water=["turbid"], sites=["bundegi", "tantabiddi"]),
    "train": dict(water=["clear", "moderate"], sites=["marina", "gulf_entrance"]),
}
# Tier 3 scenarios per water type: (scenes, noise, realistic bottoms).
T3_SCENARIOS = [(2, 5e-4, False), (4, 5e-4, False), (2, 5e-4, True)]
DEPTH_BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 30)]
NOISE_CM = 0.01  # a win needs |e_B| < |e_A| - 1 cm


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout.strip()


# Variants.

def model_sources():
    """The working tree's model sources (tracked and untracked, not ignored)."""
    files = sh(["git", "ls-files", "-co", "--exclude-standard", "model"], cwd=REPO).split("\n")
    return sorted(f for f in files if f and os.path.isfile(os.path.join(REPO, f)))


def build(spec):
    """Build a variant; returns (label, photic path)."""
    ref, _, flags = spec.partition(":")
    if ref == ".":
        h = hashlib.sha256()
        for f in model_sources():
            h.update(f.encode())
            h.update(open(os.path.join(REPO, f), "rb").read())
        key = "wt-" + h.hexdigest()[:12]
    else:
        key = sh(["git", "rev-parse", "--short=12", ref + "^{commit}"], cwd=REPO)
    if flags:
        key += "-" + hashlib.sha256(flags.encode()).hexdigest()[:8]
    out = os.path.join(AB, "build", key)
    photic = os.path.join(out, "model", "photic")
    if not os.path.exists(photic):
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out)
        if ref == ".":
            for f in model_sources():
                os.makedirs(os.path.dirname(os.path.join(out, f)), exist_ok=True)
                shutil.copy2(os.path.join(REPO, f), os.path.join(out, f))
        else:
            archive = subprocess.run(["git", "archive", key[:12], "model"], cwd=REPO, check=True,
                                     capture_output=True).stdout
            subprocess.run(["tar", "-x", "-C", out], input=archive, check=True)
        subprocess.run(["make", "-C", os.path.join(out, "model"), "-j8", "photic",
                        "CFLAGS=-O3 " + flags], check=True, capture_output=True)
        json.dump(dict(spec=spec, key=key), open(os.path.join(out, "variant.json"), "w"))
    if "POINTS /cells to model/" not in open(os.path.join(out, "model", "bam.c")).read():
        raise SystemExit("%s: this photic has no MODEL ... POINTS, which ab.py needs" % spec)
    return key, photic


# Cases: each is a directory with the full scene inputs, run.bam and truth.npy (photic's
# orientation: row 0 south), under ab/cases/.

def tier3_cases(waters):
    tables = None
    cases = []
    H, bottom = t3.truth()
    for water in waters:
        luts = None
        for n, noise, realistic in T3_SCENARIOS:
            name = "t3_%s_n%d_noise%g%s" % (water, n, noise, "_realistic" if realistic else "")
            d = os.path.join(AB, "cases", name)
            if not os.path.exists(os.path.join(d, "truth.npy")):
                if tables is None:
                    tables = t2.photic_tables(DRIVER, t3.BANDS)
                if luts is None:
                    # The lookup tables are cached by bench.py's tier 3 run; build_lut only reads them.
                    work = os.path.join(DATA, "runs", "tier3")
                    luts = [t3.build_lut(water, k, tables, None, work, 1)
                            for k in range(max(s[0] for s in T3_SCENARIOS))]
                t3.ALBEDO_BY_BOTTOM = t3.REALISTIC_ALBEDO if realistic else None
                try:
                    t3.write_inputs(d, water, n, H, bottom, tables, luts, noise, False, False)
                finally:
                    t3.ALBEDO_BY_BOTTOM = None
                np.save(os.path.join(d, "truth.npy"), H)
            cases.append(dict(name=name, tier=3, dir=d))
    return cases


def tier4_cases(sites):
    cases = []
    for site in sites:
        src = os.path.join(DATA, "derived", "tier4", os.environ.get("PHOTIC_TIER4_SITES", "sites"), site)
        d = os.path.join(AB, "cases", "t4_" + site)
        if not os.path.exists(os.path.join(d, "truth.npy")):
            os.makedirs(d, exist_ok=True)
            for f in glob.glob(os.path.join(src, "*.asc")) + [os.path.join(src, "run.bam")]:
                dst = os.path.join(d, os.path.basename(f))
                if not os.path.exists(dst):
                    os.symlink(f, dst)
            np.save(os.path.join(d, "truth.npy"), np.load(os.path.join(src, "reference.npz"))["depth"])
        cases.append(dict(name="t4_" + site, tier=4, dir=d))
    return cases


def sample(cases, n, seed):
    """Depth-stratified points per case: {case name: [(i, j), ...]}."""
    rng = np.random.default_rng(seed)
    points = {}
    for tier in (3, 4):
        tc = [c for c in cases if c["tier"] == tier]
        for k, c in enumerate(tc):
            quota = n//2//len(tc) + (1 if k < n//2 % len(tc) else 0)
            H = np.load(os.path.join(c["dir"], "truth.npy"))
            pools = []
            for lo, hi in DEPTH_BINS:
                ii, jj = np.nonzero(np.isfinite(H) & (H >= lo) & (H < hi))
                order = rng.permutation(len(ii))
                pools.append([(int(ii[o]), int(jj[o])) for o in order])
            chosen = []
            while len(chosen) < quota and any(pools):
                for p in pools:  # round robin over the depth bins
                    if p and len(chosen) < quota:
                        chosen.append(p.pop())
            points[c["name"]] = chosen
    return points


# Running.

def asc_header(case_dir):
    first = sorted(glob.glob(os.path.join(case_dir, "s1_b*.asc")))[0]
    with open(first) as f:
        return [next(f) for _ in range(6)]


def run_case(photic, case, pts, rundir):
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    for f in glob.glob(os.path.join(case["dir"], "*.asc")):
        os.symlink(os.path.realpath(f), os.path.join(rundir, os.path.basename(f)))
    H = np.load(os.path.join(case["dir"], "truth.npy"))
    mask = np.zeros(H.shape)
    for i, j in pts:
        mask[i, j] = 1.0
    with open(os.path.join(rundir, "points.asc"), "w") as f:
        f.writelines(asc_header(case["dir"]))
        for i in range(H.shape[0] - 1, -1, -1):  # north row first; photic's row 0 is south
            f.write(" ".join("%g" % v for v in mask[i]) + "\n")
    lines = []
    for line in open(os.path.join(case["dir"], "run.bam")).read().split("\n"):
        if line.startswith("MODEL Lee_Semi_Analytical"):
            lines += ["READ GRID ab_points points.asc", line + " POINTS ab_points"]
        elif line.strip():
            lines.append(line)
    open(os.path.join(rundir, "run.bam"), "w").write("\n".join(lines) + "\n")
    t0 = time.time()
    with open(os.path.join(rundir, "photic.log"), "w") as log:
        subprocess.run([photic, "run.bam"], cwd=rundir, stdin=subprocess.DEVNULL, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    seconds = time.time() - t0
    try:
        Hm = -t3.read_nc(os.path.join(rundir, "modelled_H.nc"))
    except (OSError, IndexError):
        raise SystemExit("photic failed: see %s" % os.path.join(rundir, "photic.log"))
    try:
        S = t3.read_nc(os.path.join(rundir, "modelled_H_sigma.nc"))
    except (OSError, IndexError):
        S = np.full(H.shape, np.nan)
    ii, jj = np.array(pts).T
    return dict(truth=H[ii, jj], depth=Hm[ii, jj], sigma=S[ii, jj], seconds=seconds)


# Scoring.

def paired(truth, a, b, sa, sb, rng, n_boot=2000):
    va = np.isfinite(a) & (a > 0)
    vb = np.isfinite(b) & (b > 0)
    both = va & vb
    out = dict(n=int(truth.size), coverage_A=float(va.mean()), coverage_B=float(vb.mean()),
               n_both=int(both.sum()))
    if both.sum() < 2:
        return out
    t, ea, eb = truth[both], a[both] - truth[both], b[both] - truth[both]
    rmse = lambda e: float(np.sqrt(np.mean(e**2)))
    tol = 0.5 + 0.1*t
    idx = rng.integers(0, t.size, (n_boot, t.size))
    d = np.sqrt(np.mean(eb[idx]**2, axis=1)) - np.sqrt(np.mean(ea[idx]**2, axis=1))
    out.update(rmse_A=rmse(ea), rmse_B=rmse(eb), d_rmse=rmse(eb) - rmse(ea),
               d_rmse_ci=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
               bias_A=float(ea.mean()), bias_B=float(eb.mean()),
               median_abs_A=float(np.median(np.abs(ea))), median_abs_B=float(np.median(np.abs(eb))),
               within_tol_A=float(np.mean(np.abs(ea) <= tol)), within_tol_B=float(np.mean(np.abs(eb) <= tol)),
               wins_B=float(np.mean(np.abs(eb) < np.abs(ea) - NOISE_CM)),
               wins_A=float(np.mean(np.abs(ea) < np.abs(eb) - NOISE_CM)))
    for key, s, e in (("A", sa[both], ea), ("B", sb[both], eb)):
        if np.isfinite(s).all() and (s > 0).any():
            out["sigma_cover_" + key] = float(np.mean(np.abs(e) <= s))
    out["by_depth"] = {}
    for lo, hi in DEPTH_BINS:
        m = (t >= lo) & (t < hi)
        if m.sum():
            out["by_depth"]["%g-%g" % (lo, hi)] = dict(n=int(m.sum()), rmse_A=rmse(ea[m]), rmse_B=rmse(eb[m]))
    return out


def verdict(s, min_effect):
    """B better or worse only if the 95% interval of dRMSE excludes 0 and |dRMSE| is at least
    min_effect of A's RMSE (with 300 points, sub-centimetre changes can be 'significant')."""
    if "d_rmse_ci" not in s:
        return "too few points"
    lo, hi = s["d_rmse_ci"]
    dc = s["coverage_B"] - s["coverage_A"]
    material = abs(s["d_rmse"]) >= min_effect*s["rmse_A"]
    if hi < 0 and material:
        v = "B better"
    elif lo > 0 and material:
        v = "B worse"
    elif hi < 0 or lo > 0:
        v = "no material difference (dRMSE %+.1f%%)" % (100*s["d_rmse"]/s["rmse_A"])
    else:
        v = "no clear difference"
    return v + ("" if abs(dc) < 0.02 else ", coverage %+.0f%%" % (100*dc))


def report(res):
    out = ["A/B %s  vs  %s  (split %s, %d points, seed %d)" % (
        res["A"], res["B"], res["split"], res["n_points"], res["seed"]), ""]
    head = "%-34s %4s %9s %9s %13s %22s %9s %9s %11s" % (
        "", "n", "cover A/B", "RMSE A", "RMSE B", "dRMSE [95% CI]", "tol A/B", "wins B/A", "sigma A/B")
    out.append(head)
    pct = lambda x: "%.0f" % (100*x) if x is not None else "-"
    for name, s in res["groups"].items():
        if "d_rmse_ci" not in s:
            out.append("%-34s %4d %4s/%-4s  (too few points)" % (name, s["n"], pct(s["coverage_A"]), pct(s["coverage_B"])))
            continue
        out.append("%-34s %4d %4s/%-4s %9.2f %13.2f %+7.2f [%+.2f, %+.2f] %4s/%-4s %4s/%-4s %5s/%-5s" % (
            name, s["n"], pct(s["coverage_A"]), pct(s["coverage_B"]), s["rmse_A"], s["rmse_B"], s["d_rmse"],
            s["d_rmse_ci"][0], s["d_rmse_ci"][1], pct(s["within_tol_A"]), pct(s["within_tol_B"]),
            pct(s["wins_B"]), pct(s["wins_A"]), pct(s.get("sigma_cover_A")), pct(s.get("sigma_cover_B"))))
    out.append("")
    for name in ("tier 3", "tier 4"):
        s = res["groups"].get(name, {})
        if s.get("by_depth"):
            out.append("%s RMSE by depth (A -> B): %s" % (name, "; ".join(
                "%s m %.2f -> %.2f (n %d)" % (k, v["rmse_A"], v["rmse_B"], v["n"]) for k, v in s["by_depth"].items())))
    out.append("")
    for name in ("tier 3", "tier 4", "all"):
        if name in res["groups"]:
            out.append("%-7s %s" % (name + ":", verdict(res["groups"][name], res["min_effect"])))
    out.append("run time: A %.0f s, B %.0f s" % (res["seconds_A"], res["seconds_B"]))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("A")
    ap.add_argument("B")
    ap.add_argument("--n", type=int, default=300, help="number of points (default 300)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--split", choices=sorted(SPLITS), default="test")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--min-effect", type=float, default=0.02,
                    help="smallest change in RMSE, as a fraction of A's, that counts (default 0.02)")
    args = ap.parse_args()

    (key_a, photic_a), (key_b, photic_b) = build(args.A), build(args.B)
    split = SPLITS[args.split]
    cases = tier3_cases(split["water"]) + tier4_cases(split["sites"])
    points = sample(cases, args.n, args.seed)
    tag = "%s_n%d_s%d" % (args.split, args.n, args.seed)

    jobs = [(v, photic, c) for v, photic in ((key_a, photic_a), (key_b, photic_b)) for c in cases]
    with ThreadPoolExecutor(args.jobs) as ex:
        outs = list(ex.map(lambda job: run_case(job[1], job[2], points[job[2]["name"]],
                                                os.path.join(AB, "runs", job[0], tag, job[2]["name"])), jobs))
    runs = {(v, c["name"]): o for (v, _, c), o in zip(jobs, outs)}

    rng = np.random.default_rng(0)
    groups = {}
    sel = {c["name"]: [c["name"]] for c in cases}
    sel.update({"tier 3": [c["name"] for c in cases if c["tier"] == 3],
                "tier 4": [c["name"] for c in cases if c["tier"] == 4],
                "all": [c["name"] for c in cases]})
    for g, names in sel.items():
        cat = lambda v, k: np.concatenate([runs[(v, nm)][k] for nm in names])
        groups[g] = paired(cat(key_a, "truth"), cat(key_a, "depth"), cat(key_b, "depth"),
                           cat(key_a, "sigma"), cat(key_b, "sigma"), rng)
    res = dict(A=args.A, B=args.B, key_A=key_a, key_B=key_b, split=args.split, seed=args.seed,
               min_effect=args.min_effect,
               n_points=sum(len(p) for p in points.values()),
               date=datetime.datetime.now().isoformat(timespec="seconds"),
               seconds_A=sum(runs[(key_a, c["name"])]["seconds"] for c in cases),
               seconds_B=sum(runs[(key_b, c["name"])]["seconds"] for c in cases), groups=groups)
    os.makedirs(os.path.join(AB, "results"), exist_ok=True)
    path = os.path.join(AB, "results", "%s_vs_%s_%s.json" % (key_a, key_b, tag))
    json.dump(res, open(path, "w"), indent=1)
    print(report(res))
    print("\nresults: %s" % path)


if __name__ == "__main__":
    main()
