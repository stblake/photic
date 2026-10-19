"""Experiments and figures for the DROP paper (paper/MULTISPECTRAL_SEMI_ANALYTICAL_BATHY.tex).

    paper_profile.py ablate     # each component switched off in turn: full tiers 3 and 4
    paper_profile.py scenes     # tier 4 with every subset of the scenes at each site
    paper_profile.py figures    # figures and the numbers quoted in the paper

DROP (MODEL DROP, formerly Lee_Profile; EXPERIMENTS.md Phase F) is built from the working tree by
ab.py's cached builds, with its PROFILE_* settings overridden by -D flags. The semi-analytical
reference is photic HEAD's full benchmark run (--sa, default head-d9d14d6).
"""

import argparse
import itertools
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import ab
import tier3_scenes
import tier4_score

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = ab.DATA
RUNS = os.path.join(DATA, "runs")
OUT = os.path.join(DATA, "runs", "paper")
SITES_DIR = os.path.join(DATA, "derived", "tier4", "sites")
FIGS = os.path.join(HERE, "..", "paper", "msab_figures")
PY = sys.executable
SITES = ["bundegi", "tantabiddi", "marina", "gulf_entrance"]
SITE_LABEL = {"bundegi": "Bundegi", "tantabiddi": "Tantabiddi", "marina": "Marina",
              "gulf_entrance": "Gulf entrance"}

# Each ablation switches one component off, or moves one setting, from the frozen model.
ABLATIONS = {
    "default": "",
    "no_glint": "-DPROFILE_GLINT=0",
    "site_reference": "-DPROFILE_RADIUS=0",
    "no_ambiguity": "-DPROFILE_T_AMB=0",
    "t_deep_9": "-DPROFILE_T_DEEP=9",
    "t_deep_64": "-DPROFILE_T_DEEP=64",
    "eps_0.03": "-DPROFILE_EPS=0.03",
    "no_smoothing": "-DPROFILE_SMOOTH=1",
    "red_only": "-DPROFILE_GREEN_FRACTION=1",
    "kappa_offset": "-DPROFILE_KAPPA=1",
    "kappa_jerlov": "-DPROFILE_KAPPA=2",
    "kappa_offset_local": "-DPROFILE_KAPPA=1 -DPROFILE_KAPPA_LOCAL=1",
    "kappa_jerlov_local": "-DPROFILE_KAPPA=2 -DPROFILE_KAPPA_LOCAL=1",
    "kappa_offset_om": "-DPROFILE_KAPPA=1 -DPROFILE_OFFSET_MODEL=1",
    "kappa_jerlov_om": "-DPROFILE_KAPPA=2 -DPROFILE_OFFSET_MODEL=1",
    "bright_clip": "-DPROFILE_BRIGHT=1",
    "bright_scale": "-DPROFILE_BRIGHT=2",
    "jerlov_bright_clip": "-DPROFILE_KAPPA=2 -DPROFILE_BRIGHT=1",
    "jerlov_bright_scale": "-DPROFILE_KAPPA=2 -DPROFILE_BRIGHT=2",
}


def photic_for(flags):
    _, photic = ab.build("." + (":" + flags if flags else ""))
    return photic


def cmd_ablate(a):
    def one(item):
        name, flags = item
        cid = "paper-" + name
        if os.path.exists(os.path.join(RUNS, cid, "tier4.json")) and not a.force:
            return name
        subprocess.run([PY, "bench.py", "run", "--tier", "3", "--tier", "4", "--model", "DROP",
                        "--photic", photic_for(flags), "--cid", cid], cwd=HERE, check=True,
                       stdout=subprocess.DEVNULL)
        return name
    # Build serially (make in parallel would race on ab's cache), then run in parallel.
    for flags in ABLATIONS.values():
        photic_for(flags)
    with ThreadPoolExecutor(a.jobs) as ex:
        for name in ex.map(one, ABLATIONS.items()):
            print("done:", name, flush=True)


def run_subset(photic, site, names, rundir):
    """Run DROP at a site with the scenes in names; returns the depth (positive down)."""
    site_dir = os.path.join(SITES_DIR, site)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    lines = []
    for line in open(os.path.join(site_dir, "run.bam")):
        tok = line.split()
        if line.startswith("READ GRID") and tok[2].split("_")[0] not in names:
            continue
        if line.startswith("SCENE NAME") and tok[2] not in names:
            continue
        if line.startswith("MODEL"):
            i = tok.index("OUT")
            line = "MODEL DROP IN %s %s\n" % (" ".join(names), " ".join(tok[i:]))
        lines.append(line)
    for f in os.listdir(site_dir):
        if f.endswith(".asc") and f.split("_")[0] in names:
            os.symlink(os.path.join(site_dir, f), os.path.join(rundir, f))
    open(os.path.join(rundir, "run.bam"), "w").write("".join(lines))
    subprocess.run([photic, "run.bam"], cwd=rundir, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT, check=False)
    return -tier3_scenes.read_nc(os.path.join(rundir, "modelled_H.nc"))


def cmd_scenes(a):
    photic = photic_for("")
    jobs = []
    for site in SITES:
        bam = open(os.path.join(SITES_DIR, site, "run.bam")).read()
        scenes = [l.split()[2] for l in bam.splitlines() if l.startswith("SCENE NAME")]
        for k in range(1, len(scenes) + 1):
            for sub in itertools.combinations(scenes, k):
                jobs.append((site, list(sub)))
    ref = {s: np.load(os.path.join(SITES_DIR, s, "reference.npz"))["depth"] for s in SITES}

    def one(job):
        site, names = job
        rundir = os.path.join(OUT, "scenes", site, "+".join(names))
        H = run_subset(photic, site, names, rundir)
        m = tier4_score.metrics(ref[site], H)
        shutil.rmtree(rundir, ignore_errors=True)
        return dict(site=site, scenes=names, k=len(names), **{key: m.get(key) for key in
                    ("n_reference", "n_compared", "frac_within_tol", "rmse", "bias")})
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(one, jobs))
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, "scenes.json"), "w"), indent=1)
    for site in SITES:
        for k in sorted({r["k"] for r in res if r["site"] == site}):
            rs = [r for r in res if r["site"] == site and r["k"] == k]
            hits = [100*(r["frac_within_tol"] or 0)*r["n_compared"]/r["n_reference"] for r in rs]
            print("%-14s k=%d  hits %5.1f%% (range %4.1f-%4.1f, %d subsets)" % (site, k, np.mean(hits), min(hits),
                                                                             max(hits), len(rs)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("ablate")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--force", action="store_true")
    p = sub.add_parser("scenes")
    p.add_argument("--jobs", type=int, default=6)
    p = sub.add_parser("figures")
    p.add_argument("--sa", default="head-d9d14d6-t4", help="tier 4 run of the semi-analytical model")
    p.add_argument("--sa3", default="head-d9d14d6", help="tier 3 run of the semi-analytical model")
    a = ap.parse_args()
    if a.cmd == "figures":
        import paper_figures
        paper_figures.main(a)
    else:
        {"ablate": cmd_ablate, "scenes": cmd_scenes}[a.cmd](a)
