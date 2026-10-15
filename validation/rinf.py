"""Deep-water Rrs (SCENE RINF) for photic's runs, which seeds the depth search with H0
(samodel_h0, SEED_MODE; EXPERIMENTS.md, "Class-aware H0 and seeding the search").

  image  sounding-free, from each scene's own pixels: rank them by red (bands >= 660 nm)
         smoothed over 7 x 7 cells, keep the darkest 20%, then the darkest 15% of those by
         smoothed green, and take the median spectrum. Every bottom is far brighter than deep
         water in the red, so the darkest red is deeper than the red's reach, and among those the
         darkest green is the deepest. (Ranking raw single pixels selects noise in clear water.)
  lut    tier 3 only: the OSOAA water-leaving Rrs with no bottom at the deepest lookup-table
         node (32 m), i.e. the exact deep-water reflectance.

The noise of each band (SCENE RSIGMA) can be estimated from the same deep cells (image_rsigma).

The benchmark uses the PHOTIC_RINF and PHOTIC_RSIGMA (image) environment variables (unset: none,
as before); ab.py a "^image", "^lut" or "^rsigma" variant suffix, combined with "+".
"""

import os
import re

import numpy as np

SMOOTH = 7
RED_FRACTION = 0.2
GREEN_FRACTION = 0.15


def mode():
    m = os.environ.get("PHOTIC_RINF", "")
    if m not in ("", "image", "lut"):
        raise SystemExit("PHOTIC_RINF must be image or lut, not %r" % m)
    return m or None


def read_asc(path):
    with open(path) as f:
        header = {}
        for _ in range(6):
            k, v = f.readline().split()
            header[k.lower()] = float(v)
        a = np.loadtxt(f)
    a[a == header["nodata_value"]] = np.nan
    return a[::-1]  # photic's orientation: row 0 south


def smooth(field, valid, size=SMOOTH):
    from scipy import ndimage
    n = ndimage.uniform_filter(valid.astype(float), size, mode="constant")
    s = ndimage.uniform_filter(np.where(valid, field, 0.0), size, mode="constant")
    return s/np.maximum(n, 1e-9)


def deep_cells(bands, lam):
    """The deep cells of one scene (a mask) and its band grids as an array."""
    R = np.array(bands)
    lam = np.asarray(lam, dtype=float)
    valid = np.all(np.isfinite(R) & (R > 0), axis=0)
    red = smooth(R[lam >= 660.0].sum(axis=0), valid)
    green = smooth(R[int(np.argmin(np.abs(lam - 560.0)))], valid)
    cand = valid & (red <= np.percentile(red[valid], 100*RED_FRACTION))
    return cand & (green <= np.percentile(green[cand], 100*GREEN_FRACTION)), R


def image_rinf(bands, lam):
    """Deep-water Rrs (list, one per band) of one scene from its band grids (Rrs)."""
    deep, R = deep_cells(bands, lam)
    return [float(np.median(R[b][deep])) for b in range(len(lam))]


def image_rsigma(bands, lam):
    """Noise (list, one per band, in Rrs) of one scene: 1.4826 x the median absolute deviation of
    its deep cells' raw Rrs. An upper bound: it includes real variation of the deep water."""
    deep, R = deep_cells(bands, lam)
    return [float(1.4826*np.median(np.abs(R[b][deep] - np.median(R[b][deep])))) for b in range(len(lam))]


def scenes_of(run_bam):
    """[(name, band grid names, wavelengths)] of the SCENE lines, and {grid: file} of READ GRID."""
    files, scenes = {}, []
    for line in run_bam.splitlines():
        tok = line.split()
        if line.startswith("READ GRID"):
            files[tok[2]] = tok[3]
        elif line.startswith("SCENE NAME"):
            i, j = tok.index("BANDS"), tok.index("WAVELENGTHS")
            bands = tok[i + 1:j]
            scenes.append((tok[2], bands, [float(x) for x in tok[j + 1:j + 1 + len(bands)]]))
    return scenes, files


def add_rinf(run_bam, rinf):
    """run.bam with RINF inserted after each SCENE's THETASUN value; rinf: {scene: [Rrs]}."""
    out = []
    for line in run_bam.splitlines():
        if line.startswith("SCENE NAME") and " RINF " not in line:
            name = line.split()[2]
            if name in rinf:
                line = re.sub(r"(THETASUN \S+)", lambda m: m.group(1) + " RINF " +
                              " ".join("%.6g" % v for v in rinf[name]), line, count=1)
        out.append(line)
    return "\n".join(out) + ("\n" if run_bam.endswith("\n") else "")


def image_rinf_for_dir(case_dir, run_bam=None, fn=None):
    """{scene: [Rrs]} for the scenes of a run directory, by the image rule (fn: image_rinf, or
    image_rsigma for the noise)."""
    fn = fn or image_rinf
    run_bam = run_bam or open(os.path.join(case_dir, "run.bam")).read()
    scenes, files = scenes_of(run_bam)
    return {name: fn([read_asc(os.path.join(case_dir, files[b])) for b in bands], lam)
            for name, bands, lam in scenes}


def add_rsigma(run_bam, rsigma):
    """run.bam with RSIGMA appended to each SCENE line that has none; rsigma: {scene: [Rrs]}."""
    out = []
    for line in run_bam.splitlines():
        if line.startswith("SCENE NAME") and " RSIGMA " not in line and line.split()[2] in rsigma:
            line += " RSIGMA " + " ".join("%.6g" % v for v in rsigma[line.split()[2]])
        out.append(line)
    return "\n".join(out) + ("\n" if run_bam.endswith("\n") else "")


def rsigma_mode():
    m = os.environ.get("PHOTIC_RSIGMA", "")
    if m not in ("", "image"):
        raise SystemExit("PHOTIC_RSIGMA must be image, not %r" % m)
    return m or None


def lut_rinf(water, n_scenes, data, lut_dir="lut", bands=(443.0, 490.0, 560.0, 665.0)):
    """{scene: [Rrs]} for tier 3: OSOAA's no-bottom Rrs at the deepest lookup-table node, from the
    tables bench.py's tier 3 run caches (runs/tier3/<lut_dir>/<water>_s<k>.json)."""
    import json
    out = {}
    for k in range(n_scenes):
        lut = json.load(open(os.path.join(data, "runs", "tier3", lut_dir, "%s_s%d.json" % (water, k + 1))))
        out["s%d" % (k + 1)] = [float(lut["%g" % lam]["R0"][-1]) for lam in bands]
    return out
