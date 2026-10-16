"""Tier 4: run photic on the Exmouth Gulf sites and score it against DoT soundings.

Inputs come from tier4_prepare.py (derived/tier4/sites/<site>/). Each site is run in its
own directory under the results root, because photic writes its outputs to the working
directory. Depths are compared below MSL ~ AHD; the published photic Exmouth Gulf product
(README) is scored on the same cells as a "previous photic" reference.
"""

import glob
import os
import shutil
import subprocess
import time

import numpy as np

import tier3_scenes

DEPTH_BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 30)]


def metrics(H_ref, H):
    valid = np.isfinite(H_ref) & np.isfinite(H) & (H > 0)
    out = dict(n_reference=int(np.isfinite(H_ref).sum()), n_compared=int(valid.sum()))
    if not valid.any():
        return out
    e = H[valid] - H_ref[valid]
    tol = 0.5 + 0.1*H_ref[valid]
    out.update(bias=float(e.mean()), rmse=float(np.sqrt(np.mean(e**2))), mae=float(np.mean(np.abs(e))),
               median_abs=float(np.median(np.abs(e))), frac_within_tol=float(np.mean(np.abs(e) <= tol)),
               r=float(np.corrcoef(H[valid], H_ref[valid])[0, 1]) if valid.sum() > 2 else None)
    out["by_depth"] = {}
    for lo, hi in DEPTH_BINS:
        m = valid & (H_ref >= lo) & (H_ref < hi)
        if m.any():
            d = H[m] - H_ref[m]
            out["by_depth"]["%g-%g" % (lo, hi)] = dict(n=int(m.sum()), bias=float(d.mean()),
                                                      rmse=float(np.sqrt(np.mean(d**2))))
    return out


def published_on_grid(path, lon, lat):
    import rasterio
    with rasterio.open(path) as src:
        LON, LAT = np.meshgrid(lon, lat)
        vals = np.array([v[0] for v in src.sample(zip(LON.ravel(), LAT.ravel()))], dtype=float)
        vals[vals == src.nodata] = np.nan
    vals = vals.reshape(LON.shape)
    return -vals  # published depths are negative down


def run(photic, data, workroot, only=None):
    sites_dir = os.path.join(data, "derived", "tier4", os.environ.get("PHOTIC_TIER4_SITES", "sites"))
    published = glob.glob(os.path.join(data, "raw", "published", "exmouth_gulf", "*.tif"))
    cases = {}
    for site_dir in sorted(glob.glob(os.path.join(sites_dir, "*"))):
        if not os.path.isdir(site_dir):
            continue
        site = os.path.basename(site_dir)
        if only and site not in only:
            continue
        rundir = os.path.join(workroot, site)
        shutil.rmtree(rundir, ignore_errors=True)
        os.makedirs(rundir)
        for f in glob.glob(os.path.join(site_dir, "*.asc")):
            os.symlink(f, os.path.join(rundir, os.path.basename(f)))
        bam = open(os.path.join(site_dir, "run.bam")).read()
        import rinf
        if rinf.mode():
            # Deep-water Rrs as SCENE RINF, from each scene's own deep pixels (PHOTIC_RINF; rinf.py).
            bam = rinf.add_rinf(bam, rinf.image_rinf_for_dir(site_dir, bam))
        if rinf.rsigma_mode():
            # Noise per band as SCENE RSIGMA, from the same deep pixels (PHOTIC_RSIGMA; rinf.py).
            bam = rinf.add_rsigma(bam, rinf.image_rinf_for_dir(site_dir, bam, rinf.image_rsigma))
        import spectra_inputs
        open(os.path.join(rundir, "run.bam"), "w").write(
            spectra_inputs.with_spectra(bam, os.path.join(os.path.dirname(site_dir), "meta.json"), site))
        cases[site] = (site_dir, rundir)

    def execute(rundir):
        t0 = time.time()
        with open(os.path.join(rundir, "photic.log"), "w") as log:
            subprocess.run([photic, "run.bam"], cwd=rundir, stdin=subprocess.DEVNULL, stdout=log,
                           stderr=subprocess.STDOUT, check=False)
        return time.time() - t0

    # photic is single-threaded: run the sites in parallel, then score them.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=max(len(cases), 1)) as pool:
        elapsed = dict(zip(cases, pool.map(lambda c: execute(c[1]), cases.values())))
    results = {}
    for site, (site_dir, rundir) in cases.items():
        ref = np.load(os.path.join(site_dir, "reference.npz"))
        H = -tier3_scenes.read_nc(os.path.join(rundir, "modelled_H.nc"))
        res = dict(seconds=elapsed[site], photic=metrics(ref["depth"], H))
        if published:
            P = published_on_grid(published[0], ref["lon"], ref["lat"])
            res["published_photic"] = metrics(ref["depth"], P)
        results[site] = res
    return dict(sites=results)
