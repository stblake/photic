"""The spectral data the benchmark gives photic (model/data/README.md).

The benchmark names its data explicitly, so its baselines do not move when photic's defaults do:

  PHOTIC_BENCH_BOTTOM   bottom library CSV (default model/data/bottom_hydrolight.csv, the
                        library the baselines were scored with; photic itself falls back to a
                        built-in bright sand when a command file has no READ BOTTOM)
  PHOTIC_BENCH_WATER    "legacy" (default): the 10 nm Pope & Fry, Morel and Lee tables the
                        baselines were scored with; "default": photic's own defaults; or a
                        comma-separated list of water CSVs (READ WATER, in order)
  PHOTIC_BENCH_SRF      1: tier 4 scenes get their sensor's spectral response (SCENE ... SRF),
                        from the scene files in the sites' meta.json; unset or 0: band centres

The drivers (drivers/*.c) read the same bottom library and the legacy water tables.
"""

import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATADIR = os.path.join(REPO, "model", "data")
LEGACY_WATER = ["aw_pope_fry_1997.csv", "bbw_morel_1974.csv", "aphi_lee_1994.csv"]
SENSORS = {"S2A": "S2A_MSI.csv", "S2B": "S2B_MSI.csv", "S2C": "S2C_MSI.csv", "L8_": "L8_OLI.csv",
           "L9_": "L9_OLI.csv"}


def bottom_library():
    return os.path.abspath(os.environ.get("PHOTIC_BENCH_BOTTOM") or os.path.join(DATADIR, "bottom_hydrolight.csv"))


def water_files():
    w = os.environ.get("PHOTIC_BENCH_WATER") or "legacy"
    if w == "legacy":
        return [os.path.join(DATADIR, f) for f in LEGACY_WATER]
    if w == "default":
        return []
    return [os.path.abspath(f) for f in w.split(",")]


def srf_enabled():
    return os.environ.get("PHOTIC_BENCH_SRF", "0") not in ("", "0")


def describe():
    return dict(bottom=bottom_library(), water=os.environ.get("PHOTIC_BENCH_WATER") or "legacy", srf=srf_enabled())


def scene_srfs(meta_path, site):
    """{photic scene name: SRF file} for a tier 4 site, from its sites set's meta.json."""
    meta = json.load(open(meta_path))
    files = {s["isodate"]: os.path.basename(s["file"]) for s in meta["scenes"]}
    out = {}
    for k, iso in enumerate(meta["sites"][site]["scenes"]):
        out["s%d" % (k + 1)] = os.path.join(DATADIR, "srf", SENSORS[files[iso][:3]])
    return out


def with_spectra(bam, meta_path=None, site=None):
    """A command file with the benchmark's READ BOTTOM and READ WATER first and, for tier 4 with
    PHOTIC_BENCH_SRF, SRF on each SCENE line."""
    lines = bam.splitlines()
    head = []
    if not any(line.startswith("READ BOTTOM") for line in lines):
        head.append("READ BOTTOM %s" % bottom_library())
    if not any(line.startswith("READ WATER") for line in lines):
        head += ["READ WATER %s" % f for f in water_files()]
    if srf_enabled() and meta_path and site:
        srfs = scene_srfs(meta_path, site)
        for i, line in enumerate(lines):
            if line.startswith("SCENE NAME") and " SRF " not in line:
                lines[i] = line + " SRF " + srfs[line.split()[2]]
    return "\n".join(head + lines) + ("\n" if bam.endswith("\n") else "")


def read_csv(path):
    """A photic data CSV (model/data/README.md) as {column: array}, wavelength included."""
    import numpy as np
    rows = [line.strip() for line in open(path) if line.strip() and not line.startswith("#")]
    names = rows[0].split(",")
    data = [r.split(",") for r in rows[1:] if r[0].isdigit()]
    return {n: np.array([float(d[i]) for d in data]) for i, n in enumerate(names)}
