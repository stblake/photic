#!/usr/bin/env python3
"""photic physics validation benchmark.

    bench.py run --tier 1 [--tier 2 ...] [--quick]   run tiers, store results for this commit
    bench.py score                                    assemble scorecards/<sha>.json and .md
    bench.py compare A B                              compare two scorecards (sha, path or "baseline")

Data (downloads, RT runs, imagery) live in $PHOTIC_VALIDATION_DATA, by default
~/photic_dev/validation_data; only scripts and scorecards are kept in git. See README.md.
"""

import argparse
import datetime
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.environ.get("PHOTIC_VALIDATION_DATA",
                      os.path.join(os.path.expanduser("~"), "photic_dev", "validation_data"))
SCORECARDS = os.path.join(HERE, "scorecards")


def git(*args):
    return subprocess.run(["git", "-C", REPO, *args], capture_output=True, text=True).stdout.strip()


def commit_id():
    sha = git("rev-parse", "--short", "HEAD")
    dirty = git("status", "--porcelain", "--", "model") != ""
    return sha + ("-dirty" if dirty else "")


def results_dir(cid):
    d = os.path.join(DATA, "runs", cid)
    os.makedirs(d, exist_ok=True)
    return d


def build_drivers():
    subprocess.run(["make", "-s", "-C", HERE], check=True)


def run_tier(tier, quick, photic=None, cid=None, only=None):
    photic = photic or os.path.join(REPO, "model", "photic")
    cid = cid or commit_id()
    if tier == 1:
        import tier1
        return tier1.run(os.path.join(HERE, "drivers", "forward_dump"),
                         os.path.join(REPO, "paper", "photic.tex"),
                         os.path.join(DATA, "raw", "tables", "pope97.txt"))
    if tier == 2:
        import tier2_osoaa
        work = os.path.join(DATA, "runs", "tier2" + ("_quick" if quick else ""))
        os.makedirs(work, exist_ok=True)
        return tier2_osoaa.run(os.path.join(HERE, "drivers", "forward_dump"),
                               os.path.join(DATA, "tools", "OSOAA_scalar"), work, quick=quick,
                               black_ocean_root=os.path.join(DATA, "tools", "OSOAA_blackocean"))
    if tier == 3:
        import tier3_scenes
        work = os.path.join(DATA, "runs", "tier3")
        os.makedirs(work, exist_ok=True)
        return tier3_scenes.run(os.path.join(HERE, "drivers", "forward_dump"), photic, DATA, work, quick=quick,
                                photic_work=os.path.join(work, cid, "photic"), only=only)
    if tier == 4:
        import tier4_score
        work = os.path.join(DATA, "runs", "tier4", cid)
        os.makedirs(work, exist_ok=True)
        return tier4_score.run(photic, DATA, work)
    raise SystemExit("tier %d is not implemented yet" % tier)


def cmd_run(args):
    build_drivers()
    import tier3_scenes
    tier3_scenes.MODEL = args.model
    cid = args.cid or commit_id()
    for tier in args.tier:
        print("running tier %d for %s ..." % (tier, cid), flush=True)
        result = run_tier(tier, args.quick, photic=args.photic, cid=cid, only=args.only)
        path = os.path.join(results_dir(cid), "tier%d.json" % tier)
        if args.only and os.path.exists(path):
            # Merge the selected scenarios into the existing results.
            old = json.load(open(path))
            old.setdefault("scenarios", {}).update(result.get("scenarios", {}))
            result = old
        import spectra_inputs
        result["_meta"] = dict(commit=cid, quick=args.quick, photic=args.photic, only=args.only, model=args.model,
                               spectra=spectra_inputs.describe(),
                               date=datetime.datetime.now().isoformat(timespec="seconds"))
        with open(path, "w") as f:
            json.dump(result, f, indent=1, sort_keys=True)
        print("  wrote %s" % path)


def flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = prefix + "." + k if prefix else k
        if isinstance(v, dict):
            out.update(flatten(v, key))
        elif isinstance(v, (int, float, bool)) and not k.startswith("_"):
            out[key] = v
    return out


def cmd_score(args):
    cid = args.cid or commit_id()
    rdir = results_dir(cid)
    card = dict(commit=cid, date=datetime.datetime.now().isoformat(timespec="seconds"))
    for name in sorted(os.listdir(rdir)):
        if name.startswith("tier") and name.endswith(".json"):
            card[name[:-5]] = json.load(open(os.path.join(rdir, name)))
    os.makedirs(SCORECARDS, exist_ok=True)
    path = os.path.join(SCORECARDS, cid + ".json")
    with open(path, "w") as f:
        json.dump(card, f, indent=1, sort_keys=True)
    with open(path[:-5] + ".md", "w") as f:
        f.write("# photic validation scorecard %s (%s)\n\n" % (cid, card["date"]))
        f.write("| metric | value |\n|---|---|\n")
        for k, v in flatten(card).items():
            f.write("| %s | %s |\n" % (k, ("%.6g" % v) if isinstance(v, float) else v))
    print("wrote %s and .md" % path)


def load_card(ref):
    if os.path.exists(ref):
        return json.load(open(ref))
    path = os.path.join(SCORECARDS, ref + ".json")
    if not os.path.exists(path):
        raise SystemExit("no scorecard %s" % path)
    return json.load(open(path))


def cmd_compare(args):
    a, b = flatten(load_card(args.a)), flatten(load_card(args.b))
    keys = sorted(set(a) | set(b))
    print("%-70s %14s %14s %12s" % ("metric", args.a, args.b, "change"))
    n_changed = 0
    for k in keys:
        va, vb = a.get(k), b.get(k)
        if va == vb:
            continue
        n_changed += 1
        if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and not isinstance(va, bool):
            change = "%+.3g" % (vb - va)
        else:
            change = "changed"
        print("%-70s %14s %14s %12s" % (k, va, vb, change))
    print("\n%d of %d metrics differ" % (n_changed, len(keys)))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--tier", type=int, action="append", required=True)
    r.add_argument("--quick", action="store_true")
    r.add_argument("--only", help="run only the tier 3 scenarios whose names contain this, and merge")
    r.add_argument("--photic", help="photic binary to score (default ../model/photic)")
    r.add_argument("--cid", help="name for the results (default: the commit)")
    r.add_argument("--model", default="Lee_Semi_Analytical",
                   help="the MODEL photic runs in tiers 3 and 4 (Lee_Semi_Analytical or Lee_Profile)")
    sc = sub.add_parser("score")
    sc.add_argument("--cid", help="results to score (default: the commit)")
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    args = p.parse_args()
    sys.path.insert(0, HERE)
    dict(run=cmd_run, score=cmd_score, compare=cmd_compare)[args.cmd](args)


if __name__ == "__main__":
    main()
