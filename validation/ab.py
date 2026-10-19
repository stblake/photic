"""Quick paired A/B test of the SA model on a few hundred out-of-sample points.

    ab.py A B [--n 300] [--seed 1] [--split test] [--jobs 8]

A and B are variants of photic: a git ref (commit, branch, tag) or "." for the working tree,
optionally followed by ":" and extra compiler flags, so a compile-time switch can be tested
without editing the code:

    ab.py HEAD .                     the working tree against the last commit
    ab.py . .:-DLUT=0                the lookup-table shortcut off
    ab.py b7-physics .:-DSEED_MODE=2
    ab.py HEAD HEAD@truth            an oracle start: the search refined from the true depth
    ab.py HEAD@prior .@prior         a rough depth prior (truth +-20%) as DEPTHS
    ab.py HEAD@emp .:-DPRIOR_PROFILE=1@emp   tier 4's empirical prior (tier4_prior.py, fitted on
                                     the training sites' soundings) as DEPTHS
    ab.py HEAD HEAD%sites5           tier 4 from another site set (here with 705 nm)
    ab.py HEAD HEAD#-s4              tier 4 without scene s4 (or "#s1,s3", or overpass dates
                                     "#20190428,20210914"); --sites sites16 for all 16 scenes
    ab.py HEAD .~DROP                another model: MODEL DROP in place
                                     of Lee_Semi_Analytical, with the same inputs

A variant is ref[:flags][%sites][#scenes][@start][^rinf][~model].

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
import re
import shutil
import subprocess
import threading
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
NC_LOCK = threading.Lock()  # the netCDF4 library is not thread-safe


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout.strip()


# Variants.

def model_sources():
    """The working tree's model sources (tracked and untracked, not ignored)."""
    files = sh(["git", "ls-files", "-co", "--exclude-standard", "model"], cwd=REPO).split("\n")
    return sorted(f for f in files if f and os.path.isfile(os.path.join(REPO, f)))


def parse_variant(spec):
    """ref[:flags][%sites][#scenes][@start][^rinf][~model] -> (build spec, sites set or None,
    scenes or None, start or None, rinf or None, model or None). ^image or ^lut passes each scene's
    deep-water Rrs as SCENE RINF (rinf.py), which seeds the depth search with H0; ^rsigma the tier 4
    noise as RSIGMA from the same deep pixels (tier 3 already passes its simulated noise); combine
    with +. ~model runs another MODEL (DROP) in place of Lee_Semi_Analytical."""
    spec, _, model = spec.partition("~")
    spec, _, rinf_mode = spec.partition("^")
    if rinf_mode and not set(rinf_mode.split("+")) <= {"image", "lut", "rsigma"}:
        raise SystemExit("unknown ^%s (image, lut or rsigma, combined with +)" % rinf_mode)
    start = None
    if spec.endswith(("@truth", "@prior", "@emp")):
        spec, _, start = spec.rpartition("@")
    spec, _, scenes = spec.partition("#")
    spec, _, sites = spec.partition("%")
    return spec, sites or None, scenes or None, start, rinf_mode or None, model or None


def scene_names(case, scenes):
    """The tier 4 scenes (s1, s2, ...) a "#" selection keeps: a comma-separated list of scene
    names or overpass dates (YYYYMMDD); a leading "-" keeps all but those."""
    meta = json.load(open(case["meta"]))
    dates = [d[:10].replace("-", "") for d in meta["sites"][case["site"]]["scenes"]]
    names = ["s%d" % (k + 1) for k in range(len(dates))]
    drop = scenes.startswith("-")
    tokens = scenes.lstrip("-").split(",")
    hit = [n for n, d in zip(names, dates) if n in tokens or d in tokens]
    # A date the site does not have (its coverage was too low) is skipped; anything else is an error.
    unknown = [t for t in tokens if t not in names and t not in dates and not (t.isdigit() and len(t) == 8)]
    if unknown:
        raise SystemExit("%s: no scene %s (scenes: %s)" % (case["name"], ", ".join(unknown),
                                                          ", ".join("%s=%s" % nd for nd in zip(names, dates))))
    return [n for n in names if n not in hit] if drop else hit


def build(spec):
    """Build a variant (its build spec, without # or @); returns (label, photic path)."""
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


def tier4_cases(sites, sites_set="sites"):
    cases = []
    for site in sites:
        src = os.path.join(DATA, "derived", "tier4", sites_set, site)
        name = "t4_" + site if sites_set == "sites" else "t4_%s_%s" % (sites_set, site)
        d = os.path.join(AB, "cases", name)
        if not os.path.exists(os.path.join(d, "truth.npy")):
            os.makedirs(d, exist_ok=True)
            for f in glob.glob(os.path.join(src, "*.asc")) + [os.path.join(src, "run.bam")]:
                dst = os.path.join(d, os.path.basename(f))
                if not os.path.exists(dst):
                    os.symlink(f, dst)
            np.save(os.path.join(d, "truth.npy"), np.load(os.path.join(src, "reference.npz"))["depth"])
        cases.append(dict(name=name, tier=4, dir=d, site=site,
                          meta=os.path.join(DATA, "derived", "tier4", sites_set, "meta.json")))
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


def prior_depth(case, H, start):
    """The depth passed as DEPTHS (positive down): the truth ("truth"); the truth +-20%
    ("prior", as tier 3's prior scenarios); or the empirical prior fitted on the training sites'
    soundings ("emp", tier4_prior.py; tier 3 has no soundings and uses the +-20% prior)."""
    if start == "truth":
        return H
    if start == "emp" and case["tier"] == 4:
        path = os.path.join(os.path.dirname(case["meta"]), case["site"], "prior.npy")
        if not os.path.exists(path):
            raise SystemExit("no %s: run tier4_prior.py %s" % (path, os.path.basename(os.path.dirname(case["meta"]))))
        return np.load(path)
    return H*(1.0 + np.random.default_rng(12345).normal(0.0, 0.2, H.shape))


def run_case(photic, case, pts, rundir, fresh=False, start=None, scenes=None, rinf_mode=None, model=None):
    done = os.path.join(rundir, "done.json")
    if not fresh and os.path.exists(done):
        # The same build, sample and case were run already (rundir encodes all three).
        out = {k: np.array(v, dtype=float) if isinstance(v, list) else v for k, v in json.load(open(done)).items()}
        if "misfit" not in out:  # runs cached before the misfit was recorded
            out["misfit"] = np.full(len(pts), np.nan)
        return out
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
    if start:
        # The true depth ("truth", an oracle) or a rough estimate of it ("prior": +-20%, as in
        # tier 3's prior scenarios) as the empirical depth, negative down. Cells without a truth
        # keep the normal search.
        P = prior_depth(case, H, start)
        with open(os.path.join(rundir, "truth_start.asc"), "w") as f:
            f.writelines(asc_header(case["dir"])[:5] + ["NODATA_value -9999\n"])
            for i in range(H.shape[0] - 1, -1, -1):
                f.write(" ".join("%.3f" % (-v) if np.isfinite(v) else "-9999" for v in P[i]) + "\n")
    for line in open(os.path.join(case["dir"], "run.bam")).read().split("\n"):
        if line.startswith("MODEL Lee_Semi_Analytical") and scenes and case["tier"] == 4:
            keep = scene_names(case, scenes)
            head, rest = line.split(" IN ", 1)
            line = head + " IN " + " ".join(keep) + " OUT " + rest.split(" OUT ", 1)[1]
        if line.startswith("MODEL Lee_Semi_Analytical"):
            if model:
                line = line.replace("MODEL Lee_Semi_Analytical", "MODEL " + model, 1)
            lines += ["READ GRID ab_points points.asc"]
            if start:
                lines += ["READ GRID ab_truth truth_start.asc"]
                line = line.replace(" DEPTHS prior", "") + " DEPTHS ab_truth"
            lines += [line + " POINTS ab_points"]
        elif line.strip():
            lines.append(line)
    bam = "\n".join(lines) + "\n"
    if rinf_mode:
        import rinf
        parts = set(rinf_mode.split("+"))
        if "lut" in parts and case["tier"] == 3:
            water, n = re.match(r"t3_([a-z]+)_n(\d+)_", case["name"]).groups()
            bam = rinf.add_rinf(bam, rinf.lut_rinf(water, int(n), DATA))
        elif parts & {"lut", "image"}:
            bam = rinf.add_rinf(bam, rinf.image_rinf_for_dir(case["dir"]))
        if "rsigma" in parts and case["tier"] == 4:
            bam = rinf.add_rsigma(bam, rinf.image_rinf_for_dir(case["dir"], fn=rinf.image_rsigma))
    import spectra_inputs
    open(os.path.join(rundir, "run.bam"), "w").write(spectra_inputs.with_spectra(bam, case.get("meta"), case.get("site")))
    t0 = time.time()
    with open(os.path.join(rundir, "photic.log"), "w") as log:
        subprocess.run([photic, "run.bam"], cwd=rundir, stdin=subprocess.DEVNULL, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    seconds = time.time() - t0
    with NC_LOCK:
        try:
            Hm = -t3.read_nc(os.path.join(rundir, "modelled_H.nc"))
        except (OSError, IndexError, RuntimeError):
            raise SystemExit("photic failed: see %s" % os.path.join(rundir, "photic.log"))
        try:
            S = t3.read_nc(os.path.join(rundir, "modelled_H_sigma.nc"))
        except (OSError, IndexError, RuntimeError):
            S = np.full(H.shape, np.nan)
        try:
            E = t3.read_nc(os.path.join(rundir, "modelled_error.nc"))
        except (OSError, IndexError, RuntimeError):
            E = np.full(H.shape, np.nan)
    ii, jj = np.array(pts).T
    out = dict(truth=H[ii, jj], depth=Hm[ii, jj], sigma=S[ii, jj], misfit=E[ii, jj], seconds=seconds)
    if start:
        out["prior"] = prior_depth(case, H, start)[ii, jj]
    json.dump({k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in out.items()}, open(done, "w"))
    return out


# Scoring.

def paired(truth, a, b, sa, sb, rng, fa=None, fb=None, n_boot=2000):
    va = np.isfinite(a) & (a > 0)
    vb = np.isfinite(b) & (b > 0)
    both = va & vb
    out = dict(n=int(truth.size), coverage_A=float(va.mean()), coverage_B=float(vb.mean()),
               n_both=int(both.sum()))
    # Over all points: a hit is a depth within tolerance, a miss a depth outside it, and no
    # depth is neither; net = hits - misses per point (+1, -1 or 0), so abstaining is neutral.
    tol_all = 0.5 + 0.1*truth
    net = lambda v, x: np.where(v, np.where(np.abs(np.where(v, x, 0.0) - truth) <= tol_all, 1.0, -1.0), 0.0)
    na, nb = net(va, a), net(vb, b)
    idx = rng.integers(0, truth.size, (n_boot, truth.size))
    dn = (nb[idx] - na[idx]).mean(axis=1)
    out.update(hits_A=float(np.mean(na > 0)), hits_B=float(np.mean(nb > 0)),
               misses_A=float(np.mean(na < 0)), misses_B=float(np.mean(nb < 0)),
               d_net=float(np.mean(nb - na)), d_net_ci=[float(np.percentile(dn, 2.5)), float(np.percentile(dn, 97.5))])
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
    if fa is not None and np.isfinite(fa[both]).all() and np.isfinite(fb[both]).all():
        # The Rrs misfit (%) each variant reached: with the same model, a lower misfit at a
        # worse depth means the misfit's minimum is not at the truth.
        out.update(misfit_A=float(np.median(fa[both])), misfit_B=float(np.median(fb[both])),
                   misfit_lower_B=float(np.mean(fb[both] < fa[both] - 1.0e-3)),
                   misfit_lower_A=float(np.mean(fa[both] < fb[both] - 1.0e-3)))
    for key, s, e in (("A", sa[both], ea), ("B", sb[both], eb)):
        if np.isfinite(s).all() and (s > 0).any():
            out["sigma_cover_" + key] = float(np.mean(np.abs(e) <= s))
    out["by_depth"] = {}
    for lo, hi in DEPTH_BINS:
        m = (t >= lo) & (t < hi)
        if m.sum():
            out["by_depth"]["%g-%g" % (lo, hi)] = dict(n=int(m.sum()), rmse_A=rmse(ea[m]), rmse_B=rmse(eb[m]))
    return out


def prior_score(truth, prior, depth):
    """The prior on its own, over all points (a prior always gives a depth where it exists),
    and the variant's net gain over it: photic adds skill only if it beats its prior."""
    tol = 0.5 + 0.1*truth
    vp = np.isfinite(prior) & (prior > 0)
    vd = np.isfinite(depth) & (depth > 0)
    net = lambda v, x: np.where(v, np.where(np.abs(np.where(v, x, 0.0) - truth) <= tol, 1.0, -1.0), 0.0)
    e = (prior - truth)[vp]
    return dict(hits=float(np.mean(net(vp, prior) > 0)), misses=float(np.mean(net(vp, prior) < 0)),
                rmse=float(np.sqrt(np.mean(e**2))) if e.size else float("nan"),
                d_net=float(np.mean(net(vd, depth) - net(vp, prior))))


def verdict(s, min_effect):
    """On net (hits - misses over all points, so coverage counts) and on RMSE where both give a
    depth. Each is better or worse only if its 95% interval excludes 0 and the change is
    material: net by at least min_effect (points per point), RMSE by min_effect of A's."""
    def one(d, ci, scale, sign):
        lo, hi = ci
        material = abs(d) >= min_effect*scale
        if (hi < 0 or lo > 0) and material:
            return "better" if sign*d > 0 else "worse"
        return "same" if (hi < 0 or lo > 0) else "unclear"
    v = "net %s (%+.1f pts)" % (one(s["d_net"], s["d_net_ci"], 1.0, 1), 100*s["d_net"])
    if "d_rmse_ci" in s:
        v += ", RMSE %s (%+.0f%%)" % (one(s["d_rmse"], s["d_rmse_ci"], s["rmse_A"], -1), 100*s["d_rmse"]/s["rmse_A"])
    dc = s["coverage_B"] - s["coverage_A"]
    return v + ("" if abs(dc) < 0.02 else ", coverage %+.0f%%" % (100*dc))


def report(res):
    out = ["A/B %s  vs  %s  (split %s, %d points, seed %d)" % (
        res["A"], res["B"], res["split"], res["n_points"], res["seed"]), ""]
    head = "%-34s %4s %9s %9s %9s %22s %15s %22s %9s %9s" % (
        "", "n", "cover A/B", "hits A/B", "miss A/B", "d net [95% CI] (pts)", "RMSE A -> B", "dRMSE [95% CI]",
        "wins B/A", "sigma A/B")
    out.append(head)
    pct = lambda x: "%.0f" % (100*x) if x is not None else "-"
    for name, s in res["groups"].items():
        row = "%-34s %4d %4s/%-4s %4s/%-4s %4s/%-4s %+6.1f [%+5.1f, %+5.1f]  " % (
            name, s["n"], pct(s["coverage_A"]), pct(s["coverage_B"]), pct(s["hits_A"]), pct(s["hits_B"]),
            pct(s["misses_A"]), pct(s["misses_B"]), 100*s["d_net"], 100*s["d_net_ci"][0], 100*s["d_net_ci"][1])
        if "d_rmse_ci" in s:
            row += "%6.2f -> %-6.2f %+6.2f [%+.2f, %+.2f] %4s/%-4s %4s/%-4s" % (
                s["rmse_A"], s["rmse_B"], s["d_rmse"], s["d_rmse_ci"][0], s["d_rmse_ci"][1],
                pct(s["wins_B"]), pct(s["wins_A"]), pct(s.get("sigma_cover_A")), pct(s.get("sigma_cover_B")))
        else:
            row += "(too few points with a depth from both)"
        out.append(row)
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
    for name in ("tier 3", "tier 4", "all"):
        for key in ("A", "B"):
            p = res["groups"].get(name, {}).get("prior_" + key)
            if p:
                out.append("%-7s prior alone (%s): hits %.0f%%, misses %.0f%%, RMSE %.2f; %s minus prior: net %+.1f pts" % (
                    name + ":", key, 100*p["hits"], 100*p["misses"], p["rmse"], key, 100*p["d_net"]))
    for name in ("tier 3", "tier 4", "all"):
        s = res["groups"].get(name, {})
        if "misfit_A" in s:
            out.append("%-7s Rrs misfit median A %.2f%%, B %.2f%%; lower in B at %.0f%%, in A at %.0f%% of points" % (
                name + ":", s["misfit_A"], s["misfit_B"], 100*s["misfit_lower_B"], 100*s["misfit_lower_A"]))
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
    ap.add_argument("--fresh", action="store_true", help="re-run cases already run for a build")
    ap.add_argument("--sites", default="sites", help="tier 4 site set under derived/tier4 (default sites; "
                    "sites16 has all 16 scenes)")
    ap.add_argument("--min-effect", type=float, default=0.02,
                    help="smallest change that counts: in net score (points per point) and in RMSE (as a "
                         "fraction of A's) (default 0.02)")
    args = ap.parse_args()

    variants = {}
    for spec in (args.A, args.B):
        bspec, sites_set, scenes, start, rinf_mode, model = parse_variant(spec)
        key, photic = build(bspec)
        # Site sets, scene selections, starts, RINF and models run the same build; their runs are kept apart.
        key += ("%" + sites_set if sites_set else "") + ("#" + scenes.replace(",", "+") if scenes else "") + \
            ("@" + start if start else "") + ("^" + rinf_mode if rinf_mode else "") + ("~" + model if model else "")
        variants[spec] = (key, photic, scenes, start, sites_set, rinf_mode, model)
    key_a, photic_a = variants[args.A][:2]
    key_b, photic_b = variants[args.B][:2]
    opts = {v[0]: v[2:] for v in variants.values()}
    split = SPLITS[args.split]
    cases = tier3_cases(split["water"]) + tier4_cases(split["sites"], args.sites)
    # A variant's own tier 4 site set (%sites): the same sites and survey, other inputs.
    alt = {}
    for key, (_, _, sites_set, _, _) in opts.items():
        if sites_set:
            for c in tier4_cases(split["sites"], sites_set):
                alt[(key, c["site"])] = c
    case_for = lambda key, c: alt.get((key, c.get("site")), c) if c["tier"] == 4 else c
    points = sample(cases, args.n, args.seed)
    tag = "%s_n%d_s%d" % (args.split, args.n, args.seed) + ("" if args.sites == "sites" else "_" + args.sites)

    jobs = [(v, photic, c) for v, photic in ((key_a, photic_a), (key_b, photic_b)) for c in cases]
    with ThreadPoolExecutor(args.jobs) as ex:
        outs = list(ex.map(lambda job: run_case(job[1], case_for(job[0], job[2]), points[job[2]["name"]],
                                                os.path.join(AB, "runs", job[0], tag, job[2]["name"]),
                                                args.fresh, opts[job[0]][1], opts[job[0]][0], opts[job[0]][3],
                                                opts[job[0]][4]), jobs))
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
                           cat(key_a, "sigma"), cat(key_b, "sigma"), rng,
                           cat(key_a, "misfit"), cat(key_b, "misfit"))
        for v, label in ((key_a, "A"), (key_b, "B")):
            if all("prior" in runs[(v, nm)] for nm in names):
                groups[g]["prior_" + label] = prior_score(cat(key_a, "truth"), cat(v, "prior"), cat(v, "depth"))
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
