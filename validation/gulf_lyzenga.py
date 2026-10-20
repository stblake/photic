"""An empirical baseline for the Arabian Gulf: photic's Lyzenga log-linear model fitted to
ICESat-2 on a split reference.

  prepare   the six haze-screened Landsat 8 scenes (tier4_gulf.py) as a median composite of each
            band on a 0.0003 deg grid over the site's box (.asc), and the ATL24 seabed photons
            split by pass into a calibration half (soundings for photic, .xyz) and a test half
            (derived/tier4/gulf_lyzenga/)
  cv        choose the input (the composite, or all scenes jointly), the bands, the deep-water
            signal, the window and POWER by 4-fold cross-validation by pass within the calibration
            half (cv.json); the test half is not used
  fit [composite|joint /bands/ /radius/ power|nopower]
            photic's MODEL lyzenga_log_linear with the chosen settings (or those given) on the whole calibration
            half, scored once on the test half, with DROP and the 2019 product on the same cells
            (score.json), and the map in the style of photic's PLOT GRID ... PALETTE gmt
  correlation  the fitted depths against ICESat-2, test and calibration halves, as density plots
            (runs/gulf_lyzenga/final/correlation.png)
  points    the ICESat-2 reference itself as coloured points, in the same format and colour scale as
            the bathymetry map (runs/gulf_lyzenga/final/icesat2.png)

Inputs: the composite (one image at the scenes' mean tide), or every scene at once: one Lyzenga
regression with each scene's log-bands as predictors, depth = a0 + sum_s sum_b a_sb ln L_sb. The
joint model includes every weighted sum of per-scene fits as a special case, so one least-squares
solve finds the best combination; the scenes' tides (constant per scene) go into the intercept, so
the target is depth below MSL (TIDE 0). A cell needs every scene.

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
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tier4_gulf as g
import tier4_prepare as tp

DATA = tp.DATA
# $GULF_REGION (tier4_gulf.py): palm (default) or scene, the whole Landsat 8 scene. The scene's
# band grids are photic NetCDF (16-bit packed; .asc would be ~6 GB of text) on ACOLITE's projected
# grid, and it reuses the palm's cross-validated settings (a sweep at 30 M cells is impractical).
OUT = os.path.join(DATA, "derived", "tier4", "gulf_lyzenga" + g.SUFFIX)
WORK = os.path.join(DATA, "runs", "gulf_lyzenga" + g.SUFFIX)
CV_PALM = os.path.join(DATA, "derived", "tier4", "gulf_lyzenga", "cv.json")
EXT = ".asc" if g.REGION == "palm" else ".nc"
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
INPUTS = ["composite", "joint"]
POWER_CHOICES = [False, True]  # MODEL lyzenga_log_linear ... POWER: the inverted power law, deep side
COLOUR_RANGE = (-25.0, 0.0)  # m, the maps' colour scale (PLOT GRID ... MIN MAX), shared by both maps
NC_LOCK = threading.Lock()  # netCDF4 is not thread-safe
POINTS_CELL = 0.0006  # deg, the ICESat-2 points map's cells: each a visible dot


def projected_l2w():
    return sorted(glob.glob(os.path.join(DATA, "derived", "tier4", g.ACOLITE, "*", "*L2W_projected.nc")))


def grid():
    if g.REGION == "palm":
        s, w, n, e = g.LIMIT
        lon = w + CELL*np.arange(int(round((e - w)/CELL)) + 1)
        lat = s + CELL*np.arange(int(round((n - s)/CELL)) + 1)
        return lon, lat
    # ACOLITE's projected grid: x (cell centres) ascending, y descending; photic's rows run south.
    import netCDF4
    with netCDF4.Dataset(projected_l2w()[0]) as d:
        return np.array(d["x"][:], float), np.array(d["y"][:], float)[::-1]


def write_band(path, grid_, lon, lat, nodata=-9999.0):
    """A band grid as .asc, or as photic's NetCDF (nc.c: lons, lats, bathy(lats, lons) as shorts
    with add_offset, scale_factor, missing_value and decompressed_missing_value)."""
    if path.endswith(".asc"):
        tp.write_asc(path, grid_, lon, lat, CELL, nodata)
        return
    import netCDF4
    v = grid_[np.isfinite(grid_)]
    lo, hi = (float(v.min()), float(v.max())) if v.size else (0.0, 1.0)
    scale = (hi - lo)/65000.0 if hi > lo else 1.0
    offset = 0.5*(hi + lo)
    packed = np.where(np.isfinite(grid_), np.rint((grid_ - offset)/scale), -32768).astype(np.int16)
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("lons", len(lon))
        d.createDimension("lats", len(lat))
        x = d.createVariable("lons", "f4", ("lons",))
        y = d.createVariable("lats", "f4", ("lats",))
        b = d.createVariable("bathy", "i2", ("lats", "lons"), zlib=True, complevel=1)
        b.set_auto_maskandscale(False)
        x[:] = lon
        y[:] = lat
        b.add_offset = np.float32(offset)
        b.scale_factor = np.float32(scale)
        b.missing_value = np.int16(-32768)
        b.decompressed_missing_value = np.float32(nodata)
        b[:] = packed


def composite(lon, lat):
    """{band: median Rrs over the scenes}, and the scenes' mean tide (GOT4.10 at the box centre)."""
    from scipy.spatial import cKDTree
    LON, LAT = np.meshgrid(lon, lat)
    stacks = {wl: [] for wl in tp.L8_BANDS}
    tides = []
    tide_at = g.tide_for({"box": (g.LIMIT[1], g.LIMIT[3], g.LIMIT[0], g.LIMIT[2], CELL)})
    if g.REGION != "palm":
        import netCDF4
        for l2w in projected_l2w():
            with netCDF4.Dataset(l2w) as d:
                for wl in tp.L8_BANDS:
                    v = np.ma.filled(d["Rrs_%d" % wl][:].astype(np.float32), np.nan)[::-1]
                    v[~(v > 0)] = np.nan
                    stacks[wl].append(v)
                tides.append(tide_at(d.getncattr("isodate")))
        return {wl: np.nanmedian(np.array(v), axis=0) for wl, v in stacks.items()}, float(np.mean(tides)), tides, stacks
    for l2w in sorted(glob.glob(os.path.join(DATA, "derived", "tier4", g.ACOLITE, "*", "*L2W.nc"))):
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
    return {wl: np.nanmedian(np.array(v), axis=0) for wl, v in stacks.items()}, float(np.mean(tides)), tides, stacks


def water_mask():
    """Cells of the composite where every band has data (saved by prepare; else from the .asc
    files, whose nodata is -9999)."""
    if os.path.exists(os.path.join(OUT, "water.npy")):
        return np.load(os.path.join(OUT, "water.npy"))
    return np.all([np.loadtxt(os.path.join(OUT, "b%d.asc" % b), skiprows=6, ndmin=2)[::-1] != -9999.0
                   for b in tp.L8_BANDS], axis=0)


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
    bands, mean_tide, tides, stacks = composite(lon, lat)
    names = []
    for wl, b in bands.items():
        name = "b%d" % wl
        write_band(os.path.join(OUT, name + EXT), b, lon, lat)
        names.append(name)
    # Each scene's own bands too, for the joint model (all scenes' bands as predictors at once).
    for wl, v in stacks.items():
        for k, b in enumerate(v):
            write_band(os.path.join(OUT, "s%d_b%d%s" % (k + 1, wl, EXT)), b, lon, lat)
    water = np.all([np.isfinite(b) for b in bands.values()], axis=0)
    np.save(os.path.join(OUT, "water.npy"), water)
    del stacks
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
                n_scenes=len(tides),
                seed=SEED, passes_calibration=n_cal, passes_test=n_test,
                cells_calibration=int(np.isfinite(cal_depth).sum()), cells_test=int(np.isfinite(test_depth).sum()))
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
    print(json.dumps(meta, indent=1))


def write_xyz(path, depth, lon, lat, tide):
    i, j = np.nonzero(np.isfinite(depth))
    with open(path, "w") as f:
        for a, b in zip(i, j):
            f.write("%.7f %.7f %.3f\n" % (lon[b], lat[a], depth[a, b] + tide))


def run_photic(rundir, soundings, bands, lsm, radius, extra="", plot=None, name="H", inputs="composite"):
    """MODEL lyzenga_log_linear on the composite, or on every scene's bands jointly; returns the
    depth grid (positive down, nan where none) and photic's printed fit."""
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    if inputs == "joint":
        n_scenes = json.load(open(os.path.join(OUT, "meta.json")))["n_scenes"]
        names = ["s%d_b%d" % (k + 1, b) for k in range(n_scenes) for b in bands]
    else:
        names = ["b%d" % b for b in bands]
    lines = ["READ GRID %s %s" % (n, os.path.join(OUT, n + EXT)) for n in names]
    lines.append("READ SOUNDINGS %s" % soundings)
    words = "LSM " + " ".join("0" for _ in names) if lsm == "none" else "PERCENTILE 1"
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
    with NC_LOCK:
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
    water = water_mask()
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
        write_xyz(path.replace(".xyz", "_msl.xyz"), train, lon, lat, 0.0)
        folds.append((path, test))
    # The joint model has many predictors (up to all four bands of each of six scenes, 24): only
    # ln L (Lsm per scene and band is not estimated), radius 2 or 3; the composite keeps the full search.
    configs = [(b, l, r, pw, "composite") for b in BAND_SETS for l in LSM_CHOICES for r in RADII for pw in POWER_CHOICES]
    configs += [(b, "none", r, pw, "joint") for b in BAND_SETS for r in (2, 3) for pw in POWER_CHOICES]

    def one(job):
        (bands, lsm, radius, pw, inputs), f = job
        # The joint model's target is depth below MSL (the tides are in its intercept).
        tide = meta["mean_tide"] if inputs == "composite" else 0.0
        soundings = folds[f][0] if inputs == "composite" else folds[f][0].replace(".xyz", "_msl.xyz")
        H, _ = run_photic(os.path.join(WORK, "cv", "%s_%s_%s_r%d%s_f%d" % (inputs, "-".join(map(str, bands)), lsm,
                                                                            radius, "_power" if pw else "", f)),
                          soundings, bands, lsm, radius, " TIDE %.4f%s" % (tide, " POWER" if pw else ""),
                          inputs=inputs)
        test = folds[f][1]
        m = np.isfinite(test) & np.isfinite(H)
        e = H[m] - test[m]
        return (bands, lsm, radius, pw, inputs), dict(n=int(np.isfinite(test).sum()), got=int(m.sum()), sse=float(np.sum(e**2)),
                                          hits=int(np.sum(np.abs(e) <= 0.5 + 0.1*test[m])))
    jobs = [(c, f) for c in configs for f in range(FOLDS)]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(one, jobs))
    table = []
    for c in configs:
        r = [x for k, x in results if k == c]
        n, got = sum(x["n"] for x in r), sum(x["got"] for x in r)
        table.append(dict(bands=list(c[0]), lsm=c[1], radius=c[2], power=c[3], inputs=c[4], hits=sum(x["hits"] for x in r)/n,
                          with_depth=got/n, rmse=float(np.sqrt(sum(x["sse"] for x in r)/max(got, 1)))))
    table.sort(key=lambda t: -t["hits"])
    for t in table:
        print("%-9s %-22s %-10s r%d %-5s  hits %.3f  with depth %.3f  rmse %.2f" % (
            t["inputs"], t["bands"], t["lsm"], t["radius"], "power" if t["power"] else "", t["hits"], t["with_depth"],
            t["rmse"]))
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


def fit(override=None):
    """The chosen settings (cv.json), or override = (inputs, bands, radius, power) for a named
    alternative, written to runs/gulf_lyzenga/<inputs>_<bands>_r<radius>[_power]/ and
    score_<same>.json."""
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    cv_file = os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM
    c = json.load(open(cv_file))["chosen"]
    rundir, score_file = os.path.join(WORK, "final"), os.path.join(OUT, "score.json")
    if override:
        inputs, bands, radius, power = override
        c = dict(inputs=inputs, bands=bands, lsm="none", radius=radius, power=power, override=True)
        tag = "%s_%s_r%d%s" % (inputs, "-".join(map(str, bands)), radius, "_power" if power else "")
        rundir, score_file = os.path.join(WORK, tag), os.path.join(OUT, "score_%s.json" % tag)
    # As the 2019 product's figure (dubai_to_abu_dhabi_gmt.png): photic's own plot.
    plot = ["SET BACKGROUND WHITE", "SET PAGESIZE 100", "SET LINEWIDTH 20",
            "PLOT GRID Dubai_to_Abu_Dhabi PALETTE gmt MIN %g MAX %g" % COLOUR_RANGE]
    joint = c.get("inputs") == "joint"
    soundings = os.path.join(OUT, "calibration.xyz")
    if joint:
        soundings = os.path.join(OUT, "calibration_msl.xyz")
        write_xyz(soundings, ref["calibration"], lon, lat, 0.0)
    H, fit_log = run_photic(rundir, soundings, c["bands"], c["lsm"], c["radius"],
                            " TIDE %.4f%s ERROR H_error" % (0.0 if joint else meta["mean_tide"],
                                                           " POWER" if c.get("power") else ""), plot,
                            name="Dubai_to_Abu_Dhabi", inputs=c.get("inputs", "composite"))
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
    json.dump(result, open(score_file, "w"), indent=1)
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
    c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    fig.suptitle("Dubai to Abu Dhabi: Lyzenga SDB (Landsat 8%s) against ICESat-2"
                 % (", power-law rescaled beyond the crossover" if c.get("power") else ""), fontsize=12)
    out = os.path.join(rundir, "correlation.png")
    fig.savefig(out, dpi=150)
    print(out)


def points():
    """Every ICESat-2 reference photon (both halves), binned to POINTS_CELL cells (the median of at
    least 3 photons spreading at most 0.5 m), drawn by photic as the bathymetry map is."""
    rundir = os.path.join(WORK, "final")
    glon, glat = grid()
    s, w, n, e = glat[0], glon[0], glat[-1], glon[-1]
    lon = w + POINTS_CELL*np.arange(int(round((e - w)/POINTS_CELL)) + 1)
    lat = s + POINTS_CELL*np.arange(int(round((n - s)/POINTS_CELL)) + 1)
    R = np.load(g.REFERENCE)
    j = np.rint((R["lon"] - lon[0])/POINTS_CELL).astype(int)
    i = np.rint((R["lat"] - lat[0])/POINTS_CELL).astype(int)
    m = (i >= 0) & (i < len(lat)) & (j >= 0) & (j < len(lon))
    key, z = i[m]*len(lon) + j[m], R["depth"][m]
    order = np.argsort(key, kind="stable")
    key, z = key[order], z[order]
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    depth = np.full(len(lat)*len(lon), np.nan)
    for a, b in zip(starts, np.r_[starts[1:], len(key)]):
        sel = z[a:b]
        if len(sel) >= 3 and np.std(sel) <= 0.5 and np.median(sel) > 0:
            depth[key[a]] = -np.median(sel)  # negative down, as photic's bathymetry
    depth = depth.reshape(len(lat), len(lon))
    # Only over water in the imagery, as the reference used for scoring (land photons are misclassified).
    water = water_mask()
    wi = np.clip(np.rint((lat - lat[0])/CELL).astype(int), 0, water.shape[0] - 1)
    wj = np.clip(np.rint((lon - lon[0])/CELL).astype(int), 0, water.shape[1] - 1)
    depth[~water[np.ix_(wi, wj)]] = np.nan
    tp.write_asc(os.path.join(rundir, "icesat2.asc"), depth, lon, lat, POINTS_CELL)
    lines = ["READ GRID ICESat-2_ATL24 icesat2.asc", "SET BACKGROUND WHITE", "SET PAGESIZE 100", "SET LINEWIDTH 20",
             "PLOT GRID ICESat-2_ATL24 PALETTE gmt MIN %g MAX %g" % COLOUR_RANGE, "EXIT"]
    open(os.path.join(rundir, "points.bam"), "w").write("\n".join(lines) + "\n")
    with open(os.path.join(rundir, "points.log"), "w") as log:
        subprocess.run([PHOTIC, "points.bam"], cwd=rundir, input=b"icesat2.png/PNG\n", stdout=log,
                       stderr=subprocess.STDOUT)
    print(os.path.join(rundir, "icesat2.png"), int(np.isfinite(depth).sum()), "cells")


if __name__ == "__main__":
    commands = {"prepare": prepare, "cv": cv, "fit": fit, "correlation": correlation, "points": points}
    if len(sys.argv) == 6 and sys.argv[1] == "fit":
        # fit composite|joint 443,483,561,655 /radius/ power|nopower
        fit((sys.argv[2], [int(b) for b in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5] == "power"))
    elif len(sys.argv) == 2 and sys.argv[1] in commands:
        commands[sys.argv[1]]()
    else:
        raise SystemExit(__doc__)
