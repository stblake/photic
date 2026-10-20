"""An empirical baseline for the Arabian Gulf: photic's Lyzenga log-linear model fitted to
ICESat-2 on a split reference.

  prepare   the six haze-screened Landsat 8 scenes (tier4_gulf.py) as a median composite of each
            band on a 0.0003 deg grid over the site's box (.asc), and the ATL24 seabed photons
            split by pass into a calibration half (soundings for photic, .xyz) and a test half
            (derived/tier4/gulf_lyzenga/)
  cv        choose the bands, the deep-water signal and the window by 4-fold cross-validation
            by pass within the calibration half (cv.json); the test half is not used
  fit       photic's MODEL lyzenga_log_linear with the chosen settings on the whole calibration
            half, scored once on the test half, with DROP and the 2019 product on the same cells
            (score.json), and the map in the style of photic's PLOT GRID ... PALETTE gmt
  correlation  the fitted depths against ICESat-2, test and calibration halves, as density plots
            (runs/gulf_lyzenga/final/correlation.png)

The split is by ICESat-2 pass (granule), alternating by date after a seeded shuffle, so the two
halves are separate tracks rather than neighbouring photons of one track. A cell is the median of
at least 3 photons spreading at most 0.5 m (tier4_gulf.binned); a test cell that also holds
calibration photons is dropped. Depths are below MSL; the composite is at the scenes' mean tide,
which is added to the calibration depths (the depth the light saw) and removed again from the
modelled depths before scoring.
"""

import glob
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tier4_gulf as g
import tier4_prepare as tp

DATA = tp.DATA
OUT = os.path.join(DATA, "derived", "tier4", "gulf_lyzenga")
WORK = os.path.join(DATA, "runs", "gulf_lyzenga")
PHOTIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model", "photic")
PRODUCT_2019 = os.path.expanduser("~/Documents/bathymetry_data/dubai_to_abu_dhabi/H_sa_mean.nc")
CELL = tp.L8_CELL
SEED = 1
FOLDS = 4

# Candidate settings for the cross-validation: band sets (nm), the deep-water signal (the 1st
# percentile of the image, or none: ln L, Lyzenga et al. 2006) and the window radius.
BAND_SETS = [(443, 483, 561), (483, 561), (443, 483, 561, 655), (483, 561, 655), (443, 561)]
LSM_CHOICES = ["percentile", "none"]
RADII = [1, 2, 3]


def grid():
    s, w, n, e = g.LIMIT
    lon = w + CELL*np.arange(int(round((e - w)/CELL)) + 1)
    lat = s + CELL*np.arange(int(round((n - s)/CELL)) + 1)
    return lon, lat


def composite(lon, lat):
    """{band: median Rrs over the scenes}, and the scenes' mean tide (GOT4.10 at the box centre)."""
    from scipy.spatial import cKDTree
    LON, LAT = np.meshgrid(lon, lat)
    stacks = {wl: [] for wl in tp.L8_BANDS}
    tides = []
    tide_at = g.tide_for({"box": (g.LIMIT[1], g.LIMIT[3], g.LIMIT[0], g.LIMIT[2], CELL)})
    for l2w in sorted(glob.glob(os.path.join(DATA, "derived", "tier4", "acolite_gulf_l8", "*", "*L2W.nc"))):
        attrs, slon, slat, rrs = tp.read_l2w(l2w)
        ok = np.isfinite(slon) & np.isfinite(slat)
        dist, idx = cKDTree(np.column_stack([slon[ok], slat[ok]])).query(np.column_stack([LON.ravel(), LAT.ravel()]))
        far = (dist > 1.5*CELL).reshape(LON.shape)
        for wl in tp.L8_BANDS:
            src = min(rrs, key=lambda x: abs(x - wl))
            band = rrs[src][ok][idx].reshape(LON.shape)
            band[far | ~(band > 0)] = np.nan
            stacks[wl].append(band)
        tides.append(tide_at(attrs["isodate"]))
    return {wl: np.nanmedian(np.array(v), axis=0) for wl, v in stacks.items()}, float(np.mean(tides)), tides


def split_photons():
    """Calibration and test photon sets: (lon, lat, depth, granule) arrays, by pass."""
    R = np.load(g.REFERENCE)
    if "granule" not in R.files:
        raise SystemExit("atl24_gulf.npz has no granule index; rerun tier4_gulf.py reference")
    passes = np.unique(R["granule"])
    rng = np.random.default_rng(SEED)
    order = rng.permutation(passes)
    cal = np.isin(R["granule"], order[::2])
    keys = ("lon", "lat", "depth", "granule")
    return ({k: R[k][cal] for k in keys}, {k: R[k][~cal] for k in keys}, len(order[::2]), len(order[1::2]))


def bin_cells(P, lon, lat):
    j = np.rint((P["lon"] - lon[0])/CELL).astype(int)
    i = np.rint((P["lat"] - lat[0])/CELL).astype(int)
    m = (i >= 0) & (i < len(lat)) & (j >= 0) & (j < len(lon))
    key, z = i[m]*len(lon) + j[m], P["depth"][m]
    order = np.argsort(key, kind="stable")
    key, z = key[order], z[order]
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    depth = np.full(len(lat)*len(lon), np.nan)
    occupied = np.zeros(len(lat)*len(lon), dtype=bool)
    occupied[key] = True
    for a, b in zip(starts, np.r_[starts[1:], len(key)]):
        sel = z[a:b]
        if len(sel) >= 3 and np.std(sel) <= 0.5 and np.median(sel) > 0:
            depth[key[a]] = np.median(sel)
    return depth.reshape(len(lat), len(lon)), occupied.reshape(len(lat), len(lon))


def prepare():
    os.makedirs(OUT, exist_ok=True)
    lon, lat = grid()
    bands, mean_tide, tides = composite(lon, lat)
    names = []
    for wl, b in bands.items():
        name = "b%d" % wl
        tp.write_asc(os.path.join(OUT, name + ".asc"), b, lon, lat, CELL)
        names.append(name)
    water = np.all([np.isfinite(b) for b in bands.values()], axis=0)
    cal, test, n_cal, n_test = split_photons()
    cal_depth, cal_occupied = bin_cells(cal, lon, lat)
    test_depth, _ = bin_cells(test, lon, lat)
    cal_depth[~water] = np.nan
    test_depth[~water | cal_occupied] = np.nan
    i, j = np.nonzero(np.isfinite(cal_depth))
    with open(os.path.join(OUT, "calibration.xyz"), "w") as f:
        for a, b in zip(i, j):
            # photic's soundings: lon lat depth (positive down), here at the composite's mean tide.
            f.write("%.7f %.7f %.3f\n" % (lon[b], lat[a], cal_depth[a, b] + mean_tide))
    np.savez(os.path.join(OUT, "reference.npz"), lon=lon, lat=lat, calibration=cal_depth, test=test_depth,
             mean_tide=mean_tide)
    meta = dict(cell=CELL, shape=[len(lat), len(lon)], bands=names, scene_tides=tides, mean_tide=mean_tide,
                seed=SEED, passes_calibration=n_cal, passes_test=n_test,
                cells_calibration=int(np.isfinite(cal_depth).sum()), cells_test=int(np.isfinite(test_depth).sum()))
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
    print(json.dumps(meta, indent=1))


def write_xyz(path, depth, lon, lat, tide):
    i, j = np.nonzero(np.isfinite(depth))
    with open(path, "w") as f:
        for a, b in zip(i, j):
            f.write("%.7f %.7f %.3f\n" % (lon[b], lat[a], depth[a, b] + tide))


def run_photic(rundir, soundings, bands, lsm, radius, extra="", plot=None, name="H"):
    """MODEL lyzenga_log_linear on the composite; returns the depth grid (positive down, nan
    where none) and photic's printed fit."""
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    names = ["b%d" % b for b in bands]
    lines = ["READ GRID %s %s" % (n, os.path.join(OUT, n + ".asc")) for n in names]
    lines.append("READ SOUNDINGS %s" % soundings)
    words = "LSM " + " ".join("0" for _ in bands) if lsm == "none" else "PERCENTILE 1"
    lines.append("MODEL lyzenga_log_linear GRIDS %s BATHYMETRY %s %s RADIUS %d DEPTHS 0 30%s"
                 % (" ".join(names), name, words, radius, extra))
    lines.append("WRITE GRID %s H.nc" % name)
    if plot:
        lines += plot
    lines.append("EXIT")
    open(os.path.join(rundir, "run.bam"), "w").write("\n".join(lines) + "\n")
    # PLOT asks PGPLOT for a device on stdin: a PNG file in the run directory.
    with open(os.path.join(rundir, "photic.log"), "w") as log:
        subprocess.run([PHOTIC, "run.bam"], cwd=rundir, input=b"lyzenga.png/PNG\n" if plot else None,
                       stdin=None if plot else subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
    import tier3_scenes
    H = -tier3_scenes.read_nc(os.path.join(rundir, "H.nc"))
    fit = [l.rstrip() for l in open(os.path.join(rundir, "photic.log")) if l.startswith("  ")]
    return H, fit


def scores(ref, H):
    v = np.isfinite(ref)
    got = v & np.isfinite(H)
    e = H[got] - ref[got]
    out = dict(n_reference=int(v.sum()), with_depth=float(got.sum()/max(v.sum(), 1)),
               hits=float(np.sum(np.abs(e) <= 0.5 + 0.1*ref[got])/max(v.sum(), 1)),
               rmse=float(np.sqrt(np.mean(e**2))) if e.size else None, bias=float(e.mean()) if e.size else None,
               by_depth={})
    for lo, hi in [(0, 2), (2, 5), (5, 10), (10, 20), (20, 30)]:
        m = v & (ref >= lo) & (ref < hi)
        if not m.any():
            continue
        mg = m & got
        d = H[mg] - ref[mg]
        out["by_depth"]["%g-%g" % (lo, hi)] = dict(
            n=int(m.sum()), with_depth=float(mg.sum()/m.sum()), hits=float(np.sum(np.abs(d) <= 0.5 + 0.1*ref[mg])/m.sum()),
            bias=float(d.mean()) if d.size else None, rmse=float(np.sqrt(np.mean(d**2))) if d.size else None)
    return out


def cv():
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    water = np.isfinite(np.loadtxt(os.path.join(OUT, "b561.asc"), skiprows=6)[::-1])
    cal, _, _, _ = split_photons()
    passes = np.random.default_rng(SEED + 1).permutation(np.unique(cal["granule"]))
    folds = []
    for f in range(FOLDS):
        held = np.isin(cal["granule"], passes[f::FOLDS])
        train, train_occupied = bin_cells({k: v[~held] for k, v in cal.items()}, lon, lat)
        test, _ = bin_cells({k: v[held] for k, v in cal.items()}, lon, lat)
        train[~water] = np.nan
        test[~water | train_occupied] = np.nan
        path = os.path.join(OUT, "cv_train_%d.xyz" % f)
        write_xyz(path, train, lon, lat, meta["mean_tide"])
        folds.append((path, test))
    configs = [(b, l, r) for b in BAND_SETS for l in LSM_CHOICES for r in RADII]

    def one(job):
        (bands, lsm, radius), f = job
        H, _ = run_photic(os.path.join(WORK, "cv", "%s_%s_r%d_f%d" % ("-".join(map(str, bands)), lsm, radius, f)),
                          folds[f][0], bands, lsm, radius, " TIDE %.4f" % meta["mean_tide"])
        test = folds[f][1]
        m = np.isfinite(test) & np.isfinite(H)
        e = H[m] - test[m]
        return (bands, lsm, radius), dict(n=int(np.isfinite(test).sum()), got=int(m.sum()), sse=float(np.sum(e**2)),
                                          hits=int(np.sum(np.abs(e) <= 0.5 + 0.1*test[m])))
    jobs = [(c, f) for c in configs for f in range(FOLDS)]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(one, jobs))
    table = []
    for c in configs:
        r = [x for k, x in results if k == c]
        n, got = sum(x["n"] for x in r), sum(x["got"] for x in r)
        table.append(dict(bands=list(c[0]), lsm=c[1], radius=c[2], hits=sum(x["hits"] for x in r)/n,
                          with_depth=got/n, rmse=float(np.sqrt(sum(x["sse"] for x in r)/max(got, 1)))))
    table.sort(key=lambda t: -t["hits"])
    for t in table:
        print("%-22s %-10s r%d  hits %.3f  with depth %.3f  rmse %.2f" % (t["bands"], t["lsm"], t["radius"], t["hits"],
                                                                          t["with_depth"], t["rmse"]))
    json.dump(dict(folds=FOLDS, configs=table, chosen=table[0]), open(os.path.join(OUT, "cv.json"), "w"), indent=1)


def on_cells(path, lon, lat, cells):
    """A depth grid (positive down) sampled at the reference cells' centres, nearest neighbour."""
    import netCDF4
    from scipy.interpolate import RegularGridInterpolator
    with netCDF4.Dataset(path) as d:
        x = np.array(d["lons"][:], float)
        y = np.array(d["lats"][:], float)
        Z = np.ma.filled(d["bathy"][:].astype(float), np.nan)
    Z[(Z <= -43.0) | (Z >= 0.0)] = np.nan  # the 2019 product's missing value and land
    f = RegularGridInterpolator((y, x), -Z, method="nearest", bounds_error=False)
    out = np.full(cells.shape, np.nan)
    i, j = np.nonzero(np.isfinite(cells))
    out[i, j] = f(np.column_stack([lat[i], lon[j]]))
    return out


def drop_on_cells(lon, lat, cells):
    """DROP's haze-screened Landsat 8 depths (runs/tier4/gulf-heldout-l8-clear) at the cells."""
    import tier3_scenes
    out = np.full(cells.shape, np.nan)
    i, j = np.nonzero(np.isfinite(cells))
    for t in sorted(glob.glob(os.path.join(DATA, "runs", "tier4", "gulf-heldout-l8-clear", "t*"))):
        r = np.load(os.path.join(DATA, "derived", "tier4", "sites_gulf_l8", os.path.basename(t), "reference.npz"))
        H = -tier3_scenes.read_nc(os.path.join(t, "modelled_H.nc"))
        H[~(H > 0)] = np.nan
        ti = np.rint((lat[i] - r["lat"][0])/CELL).astype(int)
        tj = np.rint((lon[j] - r["lon"][0])/CELL).astype(int)
        m = (ti >= 0) & (ti < len(r["lat"])) & (tj >= 0) & (tj < len(r["lon"]))
        out[i[m], j[m]] = H[ti[m], tj[m]]
    return out


def fit():
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    c = json.load(open(os.path.join(OUT, "cv.json")))["chosen"]
    rundir = os.path.join(WORK, "final")
    # As the 2019 product's figure (dubai_to_abu_dhabi_gmt.png): photic's own plot.
    plot = ["SET BACKGROUND WHITE", "SET PAGESIZE 100", "SET LINEWIDTH 20", "PLOT GRID Dubai_to_Abu_Dhabi PALETTE gmt"]
    H, fit_log = run_photic(rundir, os.path.join(OUT, "calibration.xyz"), c["bands"], c["lsm"], c["radius"],
                            " TIDE %.4f ERROR H_error" % meta["mean_tide"], plot, name="Dubai_to_Abu_Dhabi")
    test = ref["test"]
    result = dict(settings=c, fit=fit_log, lyzenga=scores(test, H))
    if os.path.exists(PRODUCT_2019):
        result["product_2019"] = scores(test, on_cells(PRODUCT_2019, lon, lat, test))
    # DROP covers only its nine tiles; all three on the test cells inside them.
    D = drop_on_cells(lon, lat, test)
    in_tiles = np.full(test.shape, False)
    for t in glob.glob(os.path.join(DATA, "derived", "tier4", "sites_gulf_l8", "t*")):
        r = np.load(os.path.join(t, "reference.npz"))
        ii = (lat >= r["lat"][0] - CELL/2) & (lat <= r["lat"][-1] + CELL/2)
        jj = (lon >= r["lon"][0] - CELL/2) & (lon <= r["lon"][-1] + CELL/2)
        in_tiles |= np.outer(ii, jj)
    T = np.where(in_tiles, test, np.nan)
    result["tiles"] = dict(lyzenga=scores(T, H), drop=scores(T, D))
    if os.path.exists(PRODUCT_2019):
        result["tiles"]["product_2019"] = scores(T, on_cells(PRODUCT_2019, lon, lat, T))
    json.dump(result, open(os.path.join(OUT, "score.json"), "w"), indent=1)
    print("\n".join(fit_log))
    print(json.dumps({k: v for k, v in result.items() if k not in ("fit",)}, indent=1))


def correlation():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import tier3_scenes
    rundir = os.path.join(WORK, "final")
    H = -tier3_scenes.read_nc(os.path.join(rundir, "H.nc"))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    top = 26.0
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), constrained_layout=True)
    for ax, name, title in zip(axes, ("test", "calibration"), ("Test half (held out)", "Calibration half (fitted)")):
        R = ref[name]
        m = np.isfinite(R) & np.isfinite(H)
        x, y = R[m], H[m]
        e = y - x
        hits = np.mean(np.abs(e) <= 0.5 + 0.1*x)
        r = np.corrcoef(x, y)[0, 1]
        slope, icpt = np.polyfit(x, y, 1)
        hb = ax.hexbin(x, y, gridsize=70, extent=(0, top, 0, top), mincnt=1, cmap="viridis", norm=LogNorm())
        t = np.array([0.0, top])
        ax.plot(t, t, color="k", lw=1.2, label="1:1")
        ax.fill_between(t, t - (0.5 + 0.1*t), t + (0.5 + 0.1*t), color="k", alpha=0.08, lw=0,
                        label="hit tolerance (0.5 m + 10%)")
        ax.plot(t, icpt + slope*t, color="tab:red", lw=1.2, ls="--", label="least-squares line")
        ax.set_xlim(0, top)
        ax.set_ylim(0, top)
        ax.set_aspect("equal")
        ax.set_xlabel("ICESat-2 ATL24 depth below MSL (m)")
        ax.set_ylabel("SDB depth below MSL (m)")
        ax.set_title(title)
        ax.text(0.03, 0.97, "n = %d cells\nhits = %.1f%%\nRMSE = %.2f m\nbias = %+.2f m\nr = %.3f\nslope = %.2f, intercept = %.2f m"
                % (m.sum(), 100*hits, np.sqrt(np.mean(e**2)), e.mean(), r, slope, icpt),
                transform=ax.transAxes, va="top", fontsize=9, bbox=dict(facecolor="white", alpha=0.85, lw=0))
        ax.legend(loc="lower right", fontsize=8)
        fig.colorbar(hb, ax=ax, shrink=0.8, label="cells")
    fig.suptitle("Dubai to Abu Dhabi: Lyzenga SDB (Landsat 8) against ICESat-2", fontsize=12)
    out = os.path.join(rundir, "correlation.png")
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == "__main__":
    commands = {"prepare": prepare, "cv": cv, "fit": fit, "correlation": correlation}
    if len(sys.argv) == 2 and sys.argv[1] in commands:
        commands[sys.argv[1]]()
    else:
        raise SystemExit(__doc__)
