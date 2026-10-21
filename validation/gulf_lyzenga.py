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
  correlation [/run/]  the fitted depths against ICESat-2, test and calibration halves, as density
            plots (runs/gulf_lyzenga/final/correlation.png, or a fit's run, e.g.
            composite_443-483-561-655_r3_power)
  reference_maps  the ICESat-2 cells by depth: all, and the training cells only after the ratio
            filter (runs/gulf_lyzenga/final/icesat2_reference.png)
  points    the ICESat-2 reference itself as coloured points, in the same format and colour scale as
            the bathymetry map (runs/gulf_lyzenga/final/icesat2.png)
  ratio_filter  the calibration soundings filtered by their agreement with the blue/green log
            ratio (shallower than the ratio implies; or both ways, but not where it saturates), scored
            out of fold and on the test half, unfiltered and filtered (ratio_filter.json)
  errors    ERROR's calibration with Lyzenga 1978's noise (eq. 16, NOISE) and gain (eq. 17, GAIN),
            out of fold on cv's folds (errors.json)
  calibrate  Lyzenga et al. 2006's calibration condition (eq. 12) imposed on the joint fit
            (CONSTRAIN GAIN [BRIGHTNESS], NORMALISE), against the free fit, out of fold on cv's
            folds (calibrate.json)
  attenuation [/bands/]  MODEL lyzenga_log_linear's ATTENUATION with the chosen joint model: K of every scene
            and band over the box, and the depth from it, scored against Lyzenga's out of fold on
            cv's folds (the calibration half only; attenuation.json), then the maps from the whole
            calibration half (runs/gulf_lyzenga/attenuation/: attenuation.png, and K of every scene and band
            as maps, k_maps.png, and spectra against pure water, k_spectra.png; Lyzenga 1978's
            depth-invariant bottom indices, bottom_indices.png)
  consistent [/bands/]  ATTENUATION ... CONSISTENT (one bottom in every scene, the scenes' differences of
            attenuation from every cell) against stage 2's flat bottom, with Lyzenga's free fit and
            CONSTRAIN GAIN beside them, out of fold on cv's folds; then the whole calibration half
            scored on the test half, with the indices' depth and between-scene correlations
            (consistent.json; default 4 and 3 bands)
  consistent_plots  the indices of consistent's final runs (443/483/561) as maps, with stage 1's and
            CONSISTENT's attenuation and the one bottom's (bottom_indices_consistent.png)

Every run passes ATTENUATION (each scene's tide and ACOLITE's sun and view zeniths), so photic
also writes K and the attenuation's depth alongside Lyzenga's.

Inputs: the composite (one image at the scenes' mean tide), or every scene at once: one Lyzenga
regression with each scene's log-bands as predictors, depth = a0 + sum_s sum_b a_sb ln L_sb. The
joint model includes every weighted sum of per-scene fits as a special case, so one least-squares
solve finds the best combination; the scenes' tides (constant per scene) go into the intercept, so
the target is depth below MSL (TIDE 0). A cell needs every scene.

The split is by ICESat-2 pass (granule), alternating by date after a seeded shuffle, so the two
halves are separate tracks rather than neighbouring photons of one track. A cell is the median of
at least 3 photons spreading at most 0.5 m (tier4_gulf.binned); a test cell that also holds
calibration photons is dropped. Depths are below MSL, and so are the soundings files; the composite
is at the scenes' mean tide, which photic's TIDE adds to the soundings (the depth the light saw) and
removes again from the modelled depths.
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
# percentile of the image, or none: ln L, with no deep-water signal subtracted) and the window radius.
BAND_SETS = [(443, 483, 561), (483, 561), (443, 483, 561, 655), (483, 561, 655), (443, 561)]
LSM_CHOICES = ["percentile", "none"]
RADII = [1, 2, 3]
INPUTS = ["composite", "joint"]
POWER_CHOICES = [False, True]  # MODEL lyzenga_log_linear ... POWER: the inverted power law, deep side
COLOUR_RANGE = (-25.0, 0.0)  # m, the maps' colour scale (PLOT GRID ... MIN MAX), shared by both maps
NC_LOCK = threading.Lock()  # netCDF4 is not thread-safe
POINTS_CELL = 0.0006  # deg, the ICESat-2 points map's cells: each a visible dot
ATT_NODES = 0.05  # deg, the attenuation's node spacing (ATTENUATION ... NODES)
# The training soundings (calibration.xyz and every fold's) keep only the cells that agree with the
# blue/green log ratio (RatioFilter: mode, k; ratio_filter() for the evidence); None keeps all. The
# test half and the folds' held-out cells are never filtered.
RATIO_FILTER = ("both", 3.0)


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
    write_xyz(os.path.join(OUT, "calibration.xyz"), training(cal_depth), lon, lat)
    np.savez(os.path.join(OUT, "reference.npz"), lon=lon, lat=lat, calibration=cal_depth, test=test_depth,
             mean_tide=mean_tide)
    meta = dict(cell=CELL, shape=[len(lat), len(lon)], bands=names, scene_tides=tides, mean_tide=mean_tide,
                n_scenes=len(tides),
                seed=SEED, passes_calibration=n_cal, passes_test=n_test,
                cells_calibration=int(np.isfinite(cal_depth).sum()), cells_test=int(np.isfinite(test_depth).sum()))
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
    print(json.dumps(meta, indent=1))


def training(depth):
    """The training soundings: depth without the cells RATIO_FILTER drops (fitted on these cells)."""
    if RATIO_FILTER is None:
        return depth
    ratio = log_ratio()
    return np.where(RatioFilter(depth, ratio, *RATIO_FILTER).keep(depth, ratio), depth, np.nan)


def write_xyz(path, depth, lon, lat):
    """photic's soundings: lon lat depth below MSL (positive down). The tide at the imagery's time
    is photic's TIDE, which adds it; adding it here as well would count it twice."""
    i, j = np.nonzero(np.isfinite(depth))
    with open(path, "w") as f:
        for a, b in zip(i, j):
            f.write("%.7f %.7f %.3f\n" % (lon[b], lat[a], depth[a, b]))


def scene_geometry():
    """Each scene's sun and view zenith (deg, ACOLITE's L2W attributes), in the composite's order."""
    import netCDF4
    files = projected_l2w() if g.REGION != "palm" else \
        sorted(glob.glob(os.path.join(DATA, "derived", "tier4", g.ACOLITE, "*", "*L2W.nc")))
    out = []
    with NC_LOCK:
        for f in files:
            with netCDF4.Dataset(f) as d:
                out.append((float(d.getncattr("sza")), float(d.getncattr("vza"))))
    return out


def attenuation_words(inputs):
    """ATTENUATION K: each scene's tide above MSL (the joint model's soundings are below MSL) and
    geometry; the composite is one scene at the mean tide and mean geometry (TIDE)."""
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    geo = scene_geometry()
    if len(geo) != meta["n_scenes"]:
        raise SystemExit("%d L2W files for %d scenes" % (len(geo), meta["n_scenes"]))
    if inputs == "joint":
        return " ATTENUATION K SCENES %d TIDES %s SUNZENITH %s VIEWZENITH %s NODES %g" % (
            len(geo), " ".join("%.4f" % h for h in meta["scene_tides"]), " ".join("%.3f" % s for s, _ in geo),
            " ".join("%.3f" % v for _, v in geo), ATT_NODES)
    return " ATTENUATION K SCENES 1 SUNZENITH %.3f VIEWZENITH %.3f NODES %g" % (
        np.mean([s for s, _ in geo]), np.mean([v for _, v in geo]), ATT_NODES)


def run_photic(rundir, soundings, bands, lsm, radius, extra="", plot=None, name="H", inputs="composite", writes=()):
    """MODEL lyzenga_log_linear on the composite, or on every scene's bands jointly, with
    ATTENUATION; returns the depth grid (positive down, nan where none) and photic's printed fit.
    writes: further (grid, file) pairs to write, e.g. ("K_depth", "K_depth.nc")."""
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
    lines.append("MODEL lyzenga_log_linear GRIDS %s BATHYMETRY %s %s RADIUS %d DEPTHS 0 30%s%s"
                 % (" ".join(names), name, words, radius, extra, attenuation_words(inputs)))
    lines.append("WRITE GRID %s H.nc" % name)
    lines += ["WRITE GRID %s %s" % w for w in writes]
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


def cv_folds(lon, lat, meta):
    """The cross-validation folds by pass within the calibration half: each fold's training
    soundings (cv_train_<f>.xyz, below MSL, after RATIO_FILTER) and its held-out cells (unfiltered)."""
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
        write_xyz(path, training(train), lon, lat)
        folds.append((path, test))
    return folds


def cv():
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    folds = cv_folds(lon, lat, meta)
    # The joint model has many predictors (up to all four bands of each of six scenes, 24): only
    # ln L (Lsm per scene and band is not estimated), radius 2 or 3; the composite keeps the full search.
    configs = [(b, l, r, pw, "composite") for b in BAND_SETS for l in LSM_CHOICES for r in RADII for pw in POWER_CHOICES]
    configs += [(b, "none", r, pw, "joint") for b in BAND_SETS for r in (2, 3) for pw in POWER_CHOICES]

    def one(job):
        (bands, lsm, radius, pw, inputs), f = job
        # The joint model's target is depth below MSL (the tides are in its intercept).
        tide = meta["mean_tide"] if inputs == "composite" else 0.0
        soundings = folds[f][0]
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


# Lyzenga et al. 2006, eq. 12: the variants of calibrate(), (name, deep-water signal, options).
CALIBRATE_VARIANTS = [
    ("lnL_power", "none", " POWER"),
    ("lnL", "none", ""),
    ("lnL_gain", "none", " CONSTRAIN GAIN"),
    ("sub_power", "percentile", " POWER"),
    ("sub", "percentile", ""),
    ("sub_gain", "percentile", " CONSTRAIN GAIN"),
    ("sub_gain_bright", "percentile", " CONSTRAIN GAIN BRIGHTNESS"),
    ("norm", "percentile", " NORMALISE"),
    ("norm_gain", "percentile", " NORMALISE CONSTRAIN GAIN"),
    ("norm_gain_bright", "percentile", " NORMALISE CONSTRAIN GAIN BRIGHTNESS"),
    ("norm_gain_power", "percentile", " NORMALISE CONSTRAIN GAIN POWER"),
    ("norm_deep1", "percentile", " NORMALISE DEEP 1"),
    ("norm_deep0", "percentile", " NORMALISE DEEP 0"),
]


ERROR_BINS = [(0, 2), (2, 5), (5, 10), (10, 15), (15, 20), (20, 30)]


def error_scores(ref, H, S):
    """How well a 1-sigma error grid S describes the errors of H against ref: the share within
    1 and 2 sigma, the rms of e/sigma and the median sigma, overall and by depth."""
    m = np.isfinite(ref) & np.isfinite(H) & np.isfinite(S) & (S > 0)
    def one(k):
        e, s = H[k] - ref[k], S[k]
        return dict(n=int(k.sum()), within_1=float(np.mean(np.abs(e) <= s)), within_2=float(np.mean(np.abs(e) <= 2*s)),
                    rms_z=float(np.sqrt(np.mean((e/s)**2))), median_sigma=float(np.median(s)),
                    rmse=float(np.sqrt(np.mean(e**2))))
    out = one(m)
    out["by_depth"] = {"%d-%d" % b: one(m & (ref >= b[0]) & (ref < b[1])) for b in ERROR_BINS if (m & (ref >= b[0]) & (ref < b[1])).any()}
    return out


ERROR_VARIANTS = [("error", ""), ("noise", " NOISE"), ("gain", " GAIN"), ("noise_gain", " NOISE GAIN")]


def errors():
    """ERROR's calibration with Lyzenga 1978's eqs. 16 (NOISE) and 17 (GAIN): the chosen model, out
    of fold on cv's folds (the calibration half only), each 1-sigma grid scored on the held-out
    cells (errors.json)."""
    import tier3_scenes
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    folds = cv_folds(lon, lat, meta)
    jobs = [(v, f) for v in ERROR_VARIANTS for f in range(FOLDS)]

    def one(job):
        (name, words), f = job
        rundir = os.path.join(WORK, "errors", "%s_f%d" % (name, f))
        H, _ = run_photic(rundir, folds[f][0], c["bands"], c["lsm"], c["radius"],
                          " TIDE 0%s ERROR H_error%s" % (" POWER" if c.get("power") else "", words), inputs="joint",
                          writes=[("H_error", "H_error.nc")])
        with NC_LOCK:
            S = tier3_scenes.read_nc(os.path.join(rundir, "H_error.nc"))
        return job, H, S
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(one, jobs))
    table = {}
    for name, words in ERROR_VARIANTS:
        held = np.full(ref["calibration"].shape, np.nan)
        Hf, Sf = held.copy(), held.copy()
        for ((v, _), f), H, S in results:
            if v != name:
                continue
            m = np.isfinite(folds[f][1])
            held[m], Hf[m], Sf[m] = folds[f][1][m], H[m], S[m]
        r = table[name] = dict(options=("ERROR" + words), **error_scores(held, Hf, Sf))
        print("%-11s within 1 sigma %.0f%%, 2 sigma %.0f%%, rms z %.2f, median sigma %.2f m | %s" % (
            name, 100*r["within_1"], 100*r["within_2"], r["rms_z"], r["median_sigma"],
            "  ".join("%s %.0f%% (%.1f m)" % (k, 100*x["within_1"], x["median_sigma"]) for k, x in r["by_depth"].items())))
    json.dump(dict(settings=c, folds=FOLDS, variants=table), open(os.path.join(OUT, "errors.json"), "w"), indent=1)


def log_ratio(blue=483, green=561):
    """Stumpf et al.'s blue/green log ratio, ln(1000 R_blue)/ln(1000 R_green), of the median
    composite, over 3 x 3 cells as photic's window (nan where either band is missing)."""
    from scipy.ndimage import uniform_filter
    import rasterio
    out = []
    for wl in (blue, green):
        with rasterio.open(os.path.join(OUT, "b%d%s" % (wl, EXT))) as d:
            x = d.read(1).astype(float)[::-1]
            x[x == d.nodata] = np.nan
        m = np.isfinite(x) & (x > 0)
        xb = uniform_filter(np.where(m, x, 0), 3)/np.maximum(uniform_filter(m.astype(float), 3), 1e-9)
        xb[~m] = np.nan
        out.append(xb)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.log(1000*out[0])/np.log(1000*out[1])


class RatioFilter:
    """The soundings' depth against the log ratio: the running median and robust spread
    (1.4826 MAD, at least FLOOR m) in equal-count bins, fitted on one set of cells and applied to
    any (the spread held at its value where the ratio saturates). mode "shallow" drops cells more than k spreads shallower than the ratio implies (an
    ICESat-2 seabed where the image sees deep water); "both" also drops those deeper, except where
    the ratio is saturated (its top 10%), as the ratio cannot tell 15 m from 25 m."""
    FLOOR, BINS = 0.5, 40

    def __init__(self, depth, ratio, mode="shallow", k=3.0):
        m = np.isfinite(depth) & np.isfinite(ratio)
        z, r = depth[m], ratio[m]
        o = np.argsort(r)
        e = np.linspace(0, len(o), self.BINS + 1).astype(int)
        parts = [o[a:b] for a, b in zip(e[:-1], e[1:])]
        self.r = np.array([np.median(r[p]) for p in parts])
        self.z = np.array([np.median(z[p]) for p in parts])
        self.s = np.array([1.4826*np.median(np.abs(z[p] - np.median(z[p]))) for p in parts])
        self.saturated = np.percentile(r, 90)
        # Where the ratio saturates, the shallow outliers themselves would inflate the spread: it is
        # held at its value at the saturation point.
        sat = self.r >= self.saturated
        if sat.any() and (~sat).any():
            self.s[sat] = self.s[~sat][-1]
        self.mode, self.k = mode, k

    def residual(self, depth, ratio):
        """(depth - median)/spread, nan where either is missing."""
        med = np.interp(ratio, self.r, self.z)
        spr = np.maximum(np.interp(ratio, self.r, self.s), self.FLOOR)
        return (depth - med)/spr

    def keep(self, depth, ratio):
        with np.errstate(invalid="ignore"):
            res = self.residual(depth, ratio)
            drop = res < -self.k
            if self.mode == "both":
                drop |= (res > self.k) & (ratio < self.saturated)
        return np.isfinite(depth) & np.isfinite(ratio) & ~drop | (np.isfinite(depth) & ~np.isfinite(ratio))


RATIO_VARIANTS = [("none", None, None), ("shallow_k3", "shallow", 3.0), ("both_k3", "both", 3.0)]


def ratio_filter():
    """The calibration soundings filtered by their agreement with the blue/green log ratio
    (RatioFilter), for the chosen model: out of fold on cv's folds (each fold's filter fitted on its
    own training cells), scored on the held-out cells unfiltered (the honest score) and filtered
    (optimistic: the filter uses the image); then fitted on the whole calibration half and scored
    once on the test half, both ways (ratio_filter.json, runs/gulf_lyzenga/ratio_filter/)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    ratio = log_ratio()
    folds = cv_folds(lon, lat, meta)
    water = water_mask()
    cal, _, _, _ = split_photons()
    extra = " TIDE 0%s" % (" POWER" if c.get("power") else "")
    work = os.path.join(WORK, "ratio_filter")
    jobs = []
    for name, mode, k in RATIO_VARIANTS:
        for f in range(FOLDS):
            passes = np.random.default_rng(SEED + 1).permutation(np.unique(cal["granule"]))
            held = np.isin(cal["granule"], passes[f::FOLDS])
            train, _ = bin_cells({key: v[~held] for key, v in cal.items()}, lon, lat)
            train[~water] = np.nan
            filt = RatioFilter(train, ratio, mode, k) if mode else None
            if filt:
                train[~filt.keep(train, ratio)] = np.nan
            path = os.path.join(work, "train_%s_%d.xyz" % (name, f))
            os.makedirs(work, exist_ok=True)
            write_xyz(path, train, lon, lat)
            jobs.append((name, f, path, filt, int(np.isfinite(train).sum())))

    def one(job):
        name, f, path, filt, n = job
        H, _ = run_photic(os.path.join(work, "%s_f%d" % (name, f)), path, c["bands"], c["lsm"], c["radius"], extra,
                          inputs="joint")
        return job, H
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(one, jobs))
    table = {}
    for name, mode, k in RATIO_VARIANTS:
        held = np.full(ref["calibration"].shape, np.nan)
        Hf, kept = held.copy(), np.zeros(held.shape, bool)
        n_train = []
        for (nm, f, _, filt, n), H in results:
            if nm != name:
                continue
            m = np.isfinite(folds[f][1])
            held[m], Hf[m] = folds[f][1][m], H[m]
            kept[m] = filt.keep(folds[f][1], ratio)[m] if filt else True
            n_train.append(n)
        r_all, r_kept = scores(held, Hf), scores(np.where(kept, held, np.nan), Hf)
        table[name] = dict(mode=mode, k=k, train_cells=n_train, out_of_fold=r_all, out_of_fold_filtered=r_kept)
        print("%-11s train cells %s | held out, all %d: hits %.3f rmse %.2f bias %+.2f %s | filtered %d: hits %.3f" % (
            name, n_train, r_all["n_reference"], r_all["hits"], r_all["rmse"], r_all["bias"],
            " ".join("%s %2.0f" % (b, 100*x["hits"]) for b, x in r_all["by_depth"].items()), r_kept["n_reference"],
            r_kept["hits"]))
    # The whole calibration half, scored once on the test half.
    for name, mode, k in RATIO_VARIANTS:
        R = ref["calibration"].copy()
        filt = RatioFilter(R, ratio, mode, k) if mode else None
        if filt:
            R[~filt.keep(R, ratio)] = np.nan
        path = os.path.join(work, "calibration_%s.xyz" % name)
        write_xyz(path, R, lon, lat)
        H, _ = run_photic(os.path.join(work, "%s_final" % name), path, c["bands"], c["lsm"], c["radius"], extra,
                          inputs="joint")
        T = ref["test"]
        t_all = scores(T, H)
        t_kept = scores(np.where(filt.keep(T, ratio), T, np.nan), H) if filt else t_all
        table[name].update(calibration_cells=int(np.isfinite(R).sum()), test=t_all, test_filtered=t_kept)
        print("%-11s final: calibration cells %d | test, all: hits %.3f rmse %.2f bias %+.2f %s | filtered %d: hits %.3f rmse %.2f" % (
            name, np.isfinite(R).sum(), t_all["hits"], t_all["rmse"], t_all["bias"],
            " ".join("%s %2.0f" % (b, 100*x["hits"]) for b, x in t_all["by_depth"].items()), t_kept["n_reference"],
            t_kept["hits"], t_kept["rmse"]))
    json.dump(dict(settings=c, ratio="ln(1000 R483)/ln(1000 R561), median composite", variants=table),
              open(os.path.join(OUT, "ratio_filter.json"), "w"), indent=1)

    # The calibration half against the ratio, kept and dropped (both_k3), and where the dropped cells are.
    R = ref["calibration"]
    filt = RatioFilter(R, ratio, "both", 3.0)
    m = np.isfinite(R) & np.isfinite(ratio)
    keep = filt.keep(R, ratio)
    res = filt.residual(R, ratio)
    shallow = m & ~keep & (res < 0)
    deep = m & ~keep & (res > 0)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    ax = axes[0]
    ax.scatter(ratio[m & keep], R[m & keep], s=2, c="0.6", label="kept (%d)" % (m & keep).sum())
    ax.scatter(ratio[shallow], R[shallow], s=6, c="tab:red", label="dropped: shallower than the ratio (%d)" % shallow.sum())
    ax.scatter(ratio[deep], R[deep], s=6, c="tab:blue", label="dropped: deeper, ratio not saturated (%d)" % deep.sum())
    rr = np.linspace(filt.r[0], filt.r[-1], 200)
    med = np.interp(rr, filt.r, filt.z)
    spr = np.maximum(np.interp(rr, filt.r, filt.s), filt.FLOOR)
    ax.plot(rr, med, "k-", lw=1.5, label="running median")
    ax.plot(rr, med - 3*spr, "k--", lw=1, label="3 robust spreads")
    ax.plot(rr, med + 3*spr, "k--", lw=1)
    ax.axvline(filt.saturated, color="0.3", ls=":", label="saturated beyond (top 10% of the ratio)")
    ax.set_xlabel("ln(1000 R$_{483}$) / ln(1000 R$_{561}$), median composite")
    ax.set_ylabel("ICESat-2 ATL24 depth below MSL (m)")
    ax.set_ylim(0, 26)
    ax.set_xlim(np.nanpercentile(ratio[m], 0.2), np.nanpercentile(ratio[m], 99.8))
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Calibration half against the blue/green log ratio")
    ax = axes[1]
    ax.imshow(np.where(water, 1.0, np.nan), origin="lower", extent=(lon[0], lon[-1], lat[0], lat[-1]), cmap="Greys",
              vmin=0, vmax=4, interpolation="nearest")
    ii, jj = np.nonzero(m & keep)
    ax.scatter(lon[jj], lat[ii], s=1, c="0.5")
    for sel, col in ((shallow, "tab:red"), (deep, "tab:blue")):
        ii, jj = np.nonzero(sel)
        ax.scatter(lon[jj], lat[ii], s=8, c=col)
    ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
    ax.set_title("Where the dropped cells are (red: shallower than the ratio; blue: deeper)")
    out = os.path.join(work, "ratio_filter.png")
    fig.savefig(out, dpi=120)
    print(out)


def calibrate():
    """Lyzenga's calibration condition (eq. 12) imposed on the joint fit (CONSTRAIN GAIN
    [BRIGHTNESS], NORMALISE), against the free fit with and without POWER: out of fold on cv's
    folds (the calibration half only), for four bands and for 443/483/561 (calibrate.json)."""
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    folds = cv_folds(lon, lat, meta)
    jobs = [(bands, v, f) for bands in ([443, 483, 561, 655], [443, 483, 561]) for v in CALIBRATE_VARIANTS
            for f in range(FOLDS)]

    def one(job):
        bands, (name, lsm, words), f = job
        rundir = os.path.join(WORK, "calibrate", "%s_%s_f%d" % ("-".join(map(str, bands)), name, f))
        H, log = run_photic(rundir, folds[f][0], bands, lsm, c["radius"], " TIDE 0" + words, inputs="joint")
        gain = [l for l in log if l.startswith("  gain ")]
        return job, H, float(gain[0].split()[1]) if gain else None
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(one, jobs))
    table = []
    for bands in ([443, 483, 561, 655], [443, 483, 561]):
        maps = {}
        for name, lsm, words in CALIBRATE_VARIANTS:
            held = np.full(ref["calibration"].shape, np.nan)
            Hf = held.copy()
            gains = []
            for (b, v, f), H, gain in results:
                if b != bands or v[0] != name:
                    continue
                m = np.isfinite(folds[f][1])
                held[m], Hf[m] = folds[f][1][m], H[m]
                gains.append(gain)
            maps[name] = (held, Hf, gains)
        # ln(L - Lsm) loses the optically deep cells, which count as misses: also score every
        # variant on the cells where all of them have a depth.
        common = np.all([np.isfinite(Hf) for _, Hf, _ in maps.values()], axis=0)
        for name, lsm, words in CALIBRATE_VARIANTS:
            held, Hf, gains = maps[name]
            r, rc = scores(held, Hf), scores(np.where(common, held, np.nan), Hf)
            table.append(dict(bands=bands, variant=name, options=words.strip(), lsm=lsm, gain=gains, **r,
                              common=rc))
            print("%-16s %-17s hits %.3f  with depth %.3f  rmse %.2f | common %d: hits %.3f rmse %.2f bias %+.2f  %s" % (
                "-".join(map(str, bands)), name, r["hits"], r["with_depth"], r["rmse"], rc["n_reference"], rc["hits"],
                rc["rmse"], rc["bias"], " ".join("%s %2.0f" % (k, 100*x["hits"]) for k, x in rc["by_depth"].items())))
    json.dump(dict(radius=c["radius"], folds=FOLDS, variants=table), open(os.path.join(OUT, "calibrate.json"), "w"),
              indent=1)


def attenuation(bands=None):
    """ATTENUATION with the chosen joint model (or other bands, e.g. 443,483,561). Out of fold on
    cv's folds (the calibration half only), Lyzenga's depth and the attenuation's depth (stage 2)
    are scored on the same held-out cells; then the whole calibration half gives the K maps and
    the attenuation's depth map."""
    import tier3_scenes
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    if c.get("inputs") != "joint":
        raise SystemExit("attenuation needs the joint model (one set of bands per scene)")
    tag = ""
    if bands:
        c = dict(c, bands=bands, override=True)
        tag = "_" + "-".join(map(str, bands))
    work = os.path.join(WORK, "attenuation" + tag)
    extra = " TIDE 0%s" % (" POWER" if c.get("power") else "")
    folds = cv_folds(lon, lat, meta)
    held = np.full(ref["calibration"].shape, np.nan)
    H_lyz, H_att, S_att = held.copy(), held.copy(), held.copy()

    def one(f):
        rundir = os.path.join(work, "fold%d" % f)
        H, _ = run_photic(rundir, folds[f][0], c["bands"], c["lsm"], c["radius"], extra,
                          inputs="joint", writes=[("K_depth", "K_depth.nc"), ("K_depth_sigma", "K_depth_sigma.nc")])
        with NC_LOCK:
            A = -tier3_scenes.read_nc(os.path.join(rundir, "K_depth.nc"))
            S = tier3_scenes.read_nc(os.path.join(rundir, "K_depth_sigma.nc"))
        return f, H, A, S
    with ThreadPoolExecutor(max_workers=FOLDS) as pool:
        for f, H, A, S in pool.map(one, range(FOLDS)):
            m = np.isfinite(folds[f][1])
            held[m], H_lyz[m], H_att[m], S_att[m] = folds[f][1][m], H[m], A[m], S[m]
    both = np.isfinite(H_lyz) & np.isfinite(H_att)
    e = H_att[both] - held[both]
    within = float(np.mean(np.abs(e) <= S_att[both])) if both.any() else None
    result = dict(settings=c, nodes=ATT_NODES, out_of_fold=dict(lyzenga=scores(held, H_lyz), attenuation=scores(held, H_att),
                  attenuation_within_sigma=within))
    # The whole calibration half: every K grid, the depth, its sigma and the bottom term.
    rundir = os.path.join(work, "final")
    names = ["s%d_b%d" % (k + 1, b) for k in range(meta["n_scenes"]) for b in c["bands"]]
    nb = len(c["bands"])
    index_names = ["K_index%d_s%d" % (i + 1, k + 1) for k in range(meta["n_scenes"]) for i in range(nb - 1)]
    outs = ["K_" + n for n in names] + ["K_depth", "K_depth_sigma", "K_bottom"] + index_names
    _, fit_log = run_photic(rundir, os.path.join(OUT, "calibration.xyz"), c["bands"], c["lsm"], c["radius"],
                            extra + " INDICES", inputs="joint", writes=[(o, o + ".nc") for o in outs])
    result["fit"] = fit_log
    json.dump(result, open(os.path.join(OUT, "attenuation%s.json" % tag), "w"), indent=1)
    attenuation_maps(rundir, lon, lat, meta, c, names)
    attenuation_k_plots(rundir, lon, lat, meta, c, names)
    result["indices"] = attenuation_index_plots(rundir, lon, lat, meta, c, ref)
    json.dump(result, open(os.path.join(OUT, "attenuation%s.json" % tag), "w"), indent=1)
    print("\n".join(l for l in fit_log if "soundings, c" in l or "attenuation" in l or "Attenuation" in l))
    print(json.dumps(result["out_of_fold"], indent=1))


def index_checks(rundir, nb, ns, R, water):
    """The depth-invariant indices: each scene's (K_index<i>_s<k>) correlation with the depth R and
    their median correlation between scenes, and with CONSISTENT the one bottom's (K_index<i>)."""
    import tier3_scenes
    out = []
    for i in range(nb - 1):
        Y = []
        for k in range(ns):
            with NC_LOCK:
                y = tier3_scenes.read_nc(os.path.join(rundir, "K_index%d_s%d.nc" % (i + 1, k + 1)))
            y[~water] = np.nan
            Y.append(y)
        depth_r = [float(np.corrcoef(R[m], y[m])[0, 1]) for y in Y for m in [np.isfinite(R) & np.isfinite(y)]]
        pair_r = [float(np.corrcoef(Y[a][m], Y[b][m])[0, 1]) for a in range(ns) for b in range(a + 1, ns)
                  for m in [np.isfinite(Y[a]) & np.isfinite(Y[b])]]
        row = dict(index=i + 1, depth_correlation=depth_r, between_scenes_median=float(np.median(pair_r)))
        path = os.path.join(rundir, "K_index%d.nc" % (i + 1))
        if os.path.exists(path):
            with NC_LOCK:
                y = tier3_scenes.read_nc(path)
            y[~water] = np.nan
            m = np.isfinite(R) & np.isfinite(y)
            row["one_bottom_depth_correlation"] = float(np.corrcoef(R[m], y[m])[0, 1])
            # How far each scene's own index agrees with the one bottom's.
            row["scene_vs_one_bottom"] = [float(np.corrcoef(Y[k][mk], y[mk])[0, 1]) for k in range(ns)
                                          for mk in [np.isfinite(Y[k]) & np.isfinite(y)]]
        out.append(row)
    return out


def consistent_index_plots(bands=(443, 483, 561)):
    """The indices of consistent's final run (whole calibration half): each scene's, with the
    consistent K, and the one bottom's (K_index<i>), as anomalies from their medians, beside each
    scene's with stage 1's K (bottom_indices_consistent.png, in the consistent run's directory)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tier3_scenes
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    water = water_mask()
    nb, ns = len(bands), meta["n_scenes"]
    tag = "-".join(map(str, bands))
    dirs = {r: os.path.join(WORK, "consistent", "final_%s_%s" % (tag, r)) for r in ("free", "consistent")}
    dates = scene_dates()

    def read(path):
        with NC_LOCK:
            y = tier3_scenes.read_nc(path)
        y[~water] = np.nan
        return y - np.nanmedian(y)
    checks = {r: index_checks(d, nb, ns, ref["calibration"], water) for r, d in dirs.items()}
    # Columns: for each index, stage 1's K, then the consistent K; rows: the scenes, then the one bottom.
    fig, axes = plt.subplots(ns + 1, 2*(nb - 1), figsize=(3.6*2*(nb - 1) + 1, 2.3*(ns + 1) + 1.0),
                             constrained_layout=True, squeeze=False)
    extent = (lon[0], lon[-1], lat[0], lat[-1])
    for i in range(nb - 1):
        maps = {r: [read(os.path.join(d, "K_index%d_s%d.nc" % (i + 1, k + 1))) for k in range(ns)] for r, d in dirs.items()}
        one = read(os.path.join(dirs["consistent"], "K_index%d.nc" % (i + 1)))
        lim = np.nanpercentile(np.abs(np.array(maps["free"])), 98)*np.array([-1, 1])
        for col, r in enumerate(("free", "consistent")):
            ch = checks[r][i]
            ax_col = 2*i + col
            for k in range(ns + 1):
                ax = axes[k, ax_col]
                if k == ns and r == "free":
                    ax.axis("off")
                    continue
                z = one if k == ns else maps[r][k]
                im = ax.imshow(z, origin="lower", extent=extent, cmap="RdBu_r", vmin=lim[0], vmax=lim[1],
                               interpolation="nearest")
                ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
                ax.tick_params(labelsize=6)
                if k == 0:
                    ax.set_title("index %d, %s\n(r with depth %+.2f to %+.2f; between scenes %.2f)" % (
                        i + 1, "stage 1 K" if r == "free" else "consistent K", min(ch["depth_correlation"]),
                        max(ch["depth_correlation"]), ch["between_scenes_median"]), fontsize=8)
                if k == ns:
                    ax.set_title("index %d of the one bottom (r with depth %+.2f)" % (
                        i + 1, ch["one_bottom_depth_correlation"]), fontsize=8)
                if ax_col == 0:
                    ax.set_ylabel(dates[k], fontsize=8)
        fig.colorbar(im, ax=axes[:, 2*i:2*i + 2], location="bottom", shrink=0.8, aspect=30,
                     label="index %d, anomaly from the median (ln units)" % (i + 1))
    fig.suptitle("Dubai to Abu Dhabi: depth-invariant bottom indices (%s nm), with stage 1's attenuation and with "
                 "CONSISTENT's (one bottom in every scene)" % "/".join(map(str, bands)), fontsize=11)
    out = os.path.join(dirs["consistent"], "bottom_indices_consistent.png")
    fig.savefig(out, dpi=100, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(out)


def consistent(bands_list=None):
    """CONSISTENT (one bottom in every scene; the scenes' differences of attenuation from every cell)
    against stage 2's flat bottom, with Lyzenga's free fit (the chosen settings) and its
    CONSTRAIN GAIN variant beside them: out of fold on cv's folds (the calibration half only), all
    scored on the same held-out cells; then fitted on the whole calibration half and scored on the
    test half, with the depth-invariant indices (consistent.json)."""
    import tier3_scenes
    lon, lat = grid()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    ref = np.load(os.path.join(OUT, "reference.npz"))
    water = water_mask()
    c = json.load(open(os.path.join(OUT, "cv.json")))["chosen"]
    if c.get("inputs") != "joint":
        raise SystemExit("consistent needs the joint model (one set of bands per scene)")
    bands_list = bands_list or [[443, 483, 561, 655], [443, 483, 561]]
    folds = cv_folds(lon, lat, meta)
    # (name, lsm, extra words, which depth): the attenuation's depth (K_depth), or Lyzenga's (H).
    # PRIOR 1 pulls the fit towards stage 1 (the attenuation and the scenes' constants), with the
    # weight of an average parameter's data.
    variants = [("flat", "K_depth"), ("consistent", "K_depth"), ("consistent_prior1", "K_depth"),
                ("lyzenga", "H"), ("lyzenga_gain", "H")]
    runs = {"flat": "free", "lyzenga": "free", "consistent": "consistent", "consistent_prior1": "consistent_prior1",
            "lyzenga_gain": "gain"}
    words = {"free": (c["lsm"], ""), "consistent": (c["lsm"], " CONSISTENT"),
             "consistent_prior1": (c["lsm"], " CONSISTENT PRIOR 1"), "gain": ("percentile", " CONSTRAIN GAIN")}
    jobs = [(tuple(b), r, f) for b in bands_list for r in words for f in range(FOLDS)]

    def one(job):
        bands, run, f = job
        rundir = os.path.join(WORK, "consistent", "%s_%s_f%d" % ("-".join(map(str, bands)), run, f))
        lsm, w = words[run]
        H, log = run_photic(rundir, folds[f][0], list(bands), lsm, c["radius"], " TIDE 0" + w, inputs="joint",
                            writes=[("K_depth", "K_depth.nc")])
        with NC_LOCK:
            A = -tier3_scenes.read_nc(os.path.join(rundir, "K_depth.nc"))
        return job, {"H": H, "K_depth": A}, [l for l in log if "CONSISTENT" in l]
    # The CONSISTENT runs first, as they are the slow ones.
    jobs.sort(key=lambda j: not j[1].startswith("consistent"))
    with ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(one, jobs))
    table = []
    for bands in map(tuple, bands_list):
        maps = {}
        for name, which in variants:
            held = np.full(ref["calibration"].shape, np.nan)
            Hf = held.copy()
            for (b, run, f), out, _ in results:
                if b != bands or run != runs[name]:
                    continue
                m = np.isfinite(folds[f][1])
                held[m], Hf[m] = folds[f][1][m], out[which][m]
            maps[name] = (held, Hf)
        common = np.all([np.isfinite(Hf) for _, Hf in maps.values()], axis=0)
        for name, _ in variants:
            held, Hf = maps[name]
            r, rc = scores(held, Hf), scores(np.where(common, held, np.nan), Hf)
            table.append(dict(bands=list(bands), variant=name, out_of_fold=r, common=rc))
            print("%-16s %-13s hits %.3f  with depth %.3f  rmse %.2f | common %d: hits %.3f rmse %.2f bias %+.2f  %s" % (
                "-".join(map(str, bands)), name, r["hits"], r["with_depth"], r["rmse"], rc["n_reference"], rc["hits"],
                rc["rmse"], rc["bias"], " ".join("%s %2.0f" % (k, 100*x["hits"]) for k, x in rc["by_depth"].items())))
    logs = {"%s_%s_f%d" % ("-".join(map(str, b)), r, f): l for (b, r, f), _, l in results if l}
    # The whole calibration half, scored on the test half, with the indices.
    final = []
    for bands in bands_list:
        nb, ns = len(bands), meta["n_scenes"]
        outs = ["K_depth", "K_bottom"] + ["K_index%d_s%d" % (i + 1, k + 1) for k in range(ns) for i in range(nb - 1)]
        for run in ("free", "consistent"):
            rundir = os.path.join(WORK, "consistent", "final_%s_%s" % ("-".join(map(str, bands)), run))
            lsm, w = words[run]
            o = outs + (["K_index%d" % (i + 1) for i in range(nb - 1)] if run == "consistent" else [])
            H, log = run_photic(rundir, os.path.join(OUT, "calibration.xyz"), bands, lsm, c["radius"],
                                " TIDE 0" + w + " INDICES", inputs="joint", writes=[(x, x + ".nc") for x in o])
            with NC_LOCK:
                A = -tier3_scenes.read_nc(os.path.join(rundir, "K_depth.nc"))
            row = dict(bands=bands, run=run, test=dict(attenuation=scores(ref["test"], A), lyzenga=scores(ref["test"], H)),
                       indices=index_checks(rundir, nb, ns, ref["calibration"], water),
                       log=[l for l in log if "CONSISTENT" in l or "depth from the attenuation" in l])
            final.append(row)
            t = row["test"]["attenuation"]
            print("test half, %s %-10s attenuation depth: hits %.3f with depth %.3f rmse %.2f bias %+.2f  %s" % (
                "-".join(map(str, bands)), run, t["hits"], t["with_depth"], t["rmse"], t["bias"],
                " ".join("%s %2.0f" % (k, 100*x["hits"]) for k, x in t["by_depth"].items())))
            for ch in row["indices"]:
                print("   index %d: r with depth %s; between scenes %.2f%s" % (
                    ch["index"], " ".join("%+.2f" % x for x in ch["depth_correlation"]), ch["between_scenes_median"],
                    "; one bottom: r with depth %+.2f, with each scene's %s" % (
                        ch["one_bottom_depth_correlation"], " ".join("%.2f" % x for x in ch["scene_vs_one_bottom"]))
                    if "one_bottom_depth_correlation" in ch else ""))
    json.dump(dict(radius=c["radius"], folds=FOLDS, variants=table, final=final, logs=logs),
              open(os.path.join(OUT, "consistent.json"), "w"), indent=1)


def attenuation_maps(rundir, lon, lat, meta, c, names):
    """K (the median over the scenes) for each band, and the attenuation's depth and sigma
    (runs/gulf_lyzenga/attenuation/final/attenuation.png)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tier3_scenes
    read = lambda n: tier3_scenes.read_nc(os.path.join(rundir, n + ".nc"))
    nb = len(c["bands"])
    panels = [("K, %d nm (median of %d scenes), m$^{-1}$" % (b, meta["n_scenes"]),
               np.nanmedian([read("K_" + n) for n in names[k::nb]], axis=0), "viridis", None) for k, b in enumerate(c["bands"])]
    panels.append(("Depth from the attenuation, m", -read("K_depth"), "viridis_r", (0, 25)))
    panels.append(("Its 1-sigma error, m", read("K_depth_sigma"), "magma", (0, 3)))
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    extent = (lon[0], lon[-1], lat[0], lat[-1])
    for ax, (title, z, cmap, lim) in zip(axes.ravel(), panels):
        lim = lim or tuple(np.nanpercentile(z, [2, 98]))
        im = ax.imshow(z, origin="lower", extent=extent, cmap=cmap, vmin=lim[0], vmax=lim[1], interpolation="nearest")
        ax.set_title(title, fontsize=10)
        ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
        fig.colorbar(im, ax=ax, shrink=0.8)
    nodes = np.genfromtxt(os.path.join(rundir, "attenuation.csv"), delimiter=",", names=True)
    for ax in axes.ravel()[:nb]:
        ax.scatter(nodes["LON"], nodes["LAT"], s=4 + 2*np.sqrt(nodes["SOUNDINGS"]), facecolors="none", edgecolors="w", lw=0.5)
    fig.suptitle("Dubai to Abu Dhabi: attenuation from MODEL lyzenga_log_linear ... ATTENUATION (nodes %g deg; "
                 "circles: soundings near each node)" % ATT_NODES, fontsize=12)
    out = os.path.join(rundir, "attenuation.png")
    fig.savefig(out, dpi=110)
    print(out)


def scene_dates():
    """Each scene's date (ACOLITE's L2W isodate), in the composite's order."""
    import netCDF4
    files = projected_l2w() if g.REGION != "palm" else \
        sorted(glob.glob(os.path.join(DATA, "derived", "tier4", g.ACOLITE, "*", "*L2W.nc")))
    with NC_LOCK:
        out = []
        for f in files:
            with netCDF4.Dataset(f) as d:
                out.append(d.getncattr("isodate")[:10])
    return out


def attenuation_k_plots(rundir, lon, lat, meta, c, names):
    """K of every scene and band as maps (k_maps.png, one colour scale per band) and as spectra
    (k_spectra.png: each scene's median over the nodes with soundings, the bars their 10th to 90th
    percentiles), with pure water's absorption, below which K cannot physically fall."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tier3_scenes
    nb, ns = len(c["bands"]), meta["n_scenes"]
    dates, geo = scene_dates(), scene_geometry()
    K = np.array([tier3_scenes.read_nc(os.path.join(rundir, "K_%s.nc" % n)) for n in names])
    K[:, ~water_mask()] = np.nan
    fig, axes = plt.subplots(ns, nb, figsize=(3.6*nb + 1, 2.5*ns + 0.8), constrained_layout=True, squeeze=False)
    extent = (lon[0], lon[-1], lat[0], lat[-1])
    for k, b in enumerate(c["bands"]):
        lim = np.nanpercentile(K[k::nb], [2, 98])
        for s in range(ns):
            ax = axes[s, k]
            im = ax.imshow(K[s*nb + k], origin="lower", extent=extent, cmap="viridis", vmin=lim[0], vmax=lim[1],
                           interpolation="nearest")
            ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
            ax.tick_params(labelsize=7)
            if s == 0:
                ax.set_title("%d nm" % b, fontsize=11)
            if k == 0:
                ax.set_ylabel("%s\nsun %.0f\u00b0" % (dates[s], geo[s][0]), fontsize=9)
        fig.colorbar(im, ax=axes[:, k], location="bottom", shrink=0.8, aspect=30, label="K, m$^{-1}$")
    fig.suptitle("Dubai to Abu Dhabi: diffuse attenuation K of every scene and band (ATTENUATION, nodes %g deg)"
                 % ATT_NODES, fontsize=12)
    fig.savefig(os.path.join(rundir, "k_maps.png"), dpi=100)
    plt.close(fig)

    # The spectra at the nodes the soundings support (elsewhere K is the smoothing's extrapolation).
    nodes = np.genfromtxt(os.path.join(rundir, "attenuation.csv"), delimiter=",", names=True)
    held = nodes["SOUNDINGS"] > 0
    aw = np.loadtxt([l for l in open(os.path.join(os.path.dirname(PHOTIC), "data", "aw_pope_fry_1997.csv"))
                     if l[0].isdigit()], delimiter=",")
    wl = np.array(c["bands"], dtype=float)
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    for s in range(ns):
        v = np.array([nodes[n][held] for n in names[s*nb:(s + 1)*nb]])
        mid, lo, hi = np.median(v, axis=1), np.percentile(v, 10, axis=1), np.percentile(v, 90, axis=1)
        x = wl + 2.0*(s - (ns - 1)/2.0)
        ax.errorbar(x, mid, yerr=[mid - lo, hi - mid], marker="o", ms=5, capsize=3, lw=1.2,
                    label="%s (sun %.0f\u00b0)" % (dates[s], geo[s][0]))
    t = np.linspace(wl.min() - 15, wl.max() + 15, 200)
    ax.plot(t, np.interp(t, aw[:, 0], aw[:, 1]), color="k", ls="--", lw=1.2,
            label="pure water absorption $a_w$ (Pope and Fry 1997): K cannot be lower")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("K (m$^{-1}$)")
    ax.set_yscale("log")
    ax.set_xticks(wl)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("K at the %d nodes with soundings: median and 10th to 90th percentiles, by scene" % held.sum(),
                 fontsize=10)
    fig.savefig(os.path.join(rundir, "k_spectra.png"), dpi=130)
    plt.close(fig)
    print(os.path.join(rundir, "k_maps.png"), os.path.join(rundir, "k_spectra.png"))


def attenuation_index_plots(rundir, lon, lat, meta, c, ref):
    """Lyzenga 1978's depth-invariant bottom indices (INDICES) of every scene as maps
    (bottom_indices.png, one colour scale per index, and the median over the scenes), with two
    checks: the correlation of each index with the ICESat-2 depth (calibration half), against that
    of a band's own ln(L - Lsm), and the median correlation of each index between scenes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tier3_scenes
    nb, ns = len(c["bands"]), meta["n_scenes"]
    water = water_mask()
    dates, geo = scene_dates(), scene_geometry()
    Y = np.array([[tier3_scenes.read_nc(os.path.join(rundir, "K_index%d_s%d.nc" % (i + 1, k + 1)))
                   for i in range(nb - 1)] for k in range(ns)])
    Y[:, :, ~water] = np.nan
    R = ref["calibration"]
    checks = []
    for i in range(nb - 1):
        depth_r = []
        for k in range(ns):
            m = np.isfinite(R) & np.isfinite(Y[k, i])
            depth_r.append(float(np.corrcoef(R[m], Y[k, i][m])[0, 1]))
        pair_r = []
        for k in range(ns):
            for k2 in range(k + 1, ns):
                m = np.isfinite(Y[k, i]) & np.isfinite(Y[k2, i])
                pair_r.append(float(np.corrcoef(Y[k, i][m], Y[k2, i][m])[0, 1]))
        checks.append(dict(index=i + 1, depth_correlation=depth_r, between_scenes_median=float(np.median(pair_r))))
    # For comparison: a band's own signal, ln(L - Lsm) of the first band's 1st percentile, against depth.
    import rasterio
    band_r = []
    for k in range(ns):
        with rasterio.open(os.path.join(OUT, "s%d_b%d%s" % (k + 1, c["bands"][1], EXT))) as d:
            x = d.read(1).astype(float)[::-1]
            x[x == d.nodata] = np.nan
        v = x - np.nanpercentile(x[water], 1)
        X = np.log(np.where(v > 0, v, np.nan))
        m = np.isfinite(R) & np.isfinite(X)
        band_r.append(float(np.corrcoef(R[m], X[m])[0, 1]))
    # Each scene's index carries that scene's constants (sum_j A_ij c_j: the bands' bottom colour,
    # deep-water and atmospheric terms), so the maps are anomalies from each scene's median.
    A = Y - np.nanmedian(Y, axis=(2, 3), keepdims=True)
    fig, axes = plt.subplots(ns + 1, nb - 1, figsize=(4.2*(nb - 1) + 1, 2.5*(ns + 1) + 0.8), constrained_layout=True,
                             squeeze=False)
    extent = (lon[0], lon[-1], lat[0], lat[-1])
    for i in range(nb - 1):
        lim = np.nanpercentile(np.abs(A[:, i]), 98)*np.array([-1, 1])
        rows = [A[k, i] for k in range(ns)] + [np.nanmedian(A[:, i], axis=0)]
        for k, z in enumerate(rows):
            ax = axes[k, i]
            im = ax.imshow(z, origin="lower", extent=extent, cmap="RdBu_r", vmin=lim[0], vmax=lim[1],
                           interpolation="nearest")
            ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
            ax.tick_params(labelsize=7)
            if k == 0:
                ax.set_title("index %d (r with depth %.2f to %.2f;\nbetween scenes %.2f)" % (
                    i + 1, min(checks[i]["depth_correlation"]), max(checks[i]["depth_correlation"]),
                    checks[i]["between_scenes_median"]), fontsize=9)
            if i == 0:
                ax.set_ylabel("median of %d scenes" % ns if k == ns else "%s\nsun %.0f\u00b0" % (dates[k], geo[k][0]),
                              fontsize=9)
        fig.colorbar(im, ax=axes[:, i], location="bottom", shrink=0.8, aspect=30, label="index, anomaly from the scene's median (ln units)")
    fig.suptitle("Dubai to Abu Dhabi: depth-invariant bottom indices (Lyzenga 1978, with the local attenuation; "
                 "%s nm)\nfor comparison, ln(L - Lsm) at %d nm has r = %.2f to %.2f with depth" % (
                     "/".join(map(str, c["bands"])), c["bands"][1], min(band_r), max(band_r)), fontsize=11)
    out = os.path.join(rundir, "bottom_indices.png")
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print(out)
    for ch in checks:
        print("index %d: r with depth %s; between scenes (median) %.2f" % (
            ch["index"], " ".join("%.2f" % r for r in ch["depth_correlation"]), ch["between_scenes_median"]))
    print("ln(L - Lsm) at %d nm: r with depth %s" % (c["bands"][1], " ".join("%.2f" % r for r in band_r)))
    return dict(checks=checks, band_depth_correlation=band_r)


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
    grids = ["s%d_b%d" % (k + 1, b) for k in range(meta["n_scenes"]) for b in c["bands"]] if joint else \
        ["b%d" % b for b in c["bands"]]
    H, fit_log = run_photic(rundir, soundings, c["bands"], c["lsm"], c["radius"],
                            " TIDE %.4f%s ERROR H_error NOISE GAIN" % (0.0 if joint else meta["mean_tide"],
                                                           " POWER" if c.get("power") else ""), plot,
                            name="Dubai_to_Abu_Dhabi", inputs=c.get("inputs", "composite"),
                            writes=[("K_" + n, "K_%s.nc" % n) for n in grids] + [("H_error", "H_error.nc")])
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


def correlation(run="final"):
    """The run's depths (runs/gulf_lyzenga/<run>/, default the chosen model's) against ICESat-2."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    import tier3_scenes
    # White -> yellow -> orange -> red on a log scale from half a cell, so white is no cells and a
    # single cell is already pale yellow.
    cmap = LinearSegmentedColormap.from_list("white_yellow_orange_red", ["white", "yellow", "orange", "red"])
    rundir = os.path.join(WORK, run)
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
        hb = ax.hexbin(x, y, gridsize=70, extent=(0, top, 0, top), mincnt=1, cmap=cmap, linewidths=0)
        hb.set_norm(LogNorm(vmin=0.5, vmax=hb.get_array().max()))
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
        ticks = [n for n in (1, 3, 10, 30, 100, 300, 1000) if n <= hb.get_array().max()]
        cb = fig.colorbar(hb, ax=ax, shrink=0.8, label="cells")
        cb.set_ticks(ticks)
        cb.set_ticklabels([str(n) for n in ticks])
        cb.minorticks_off()
    if run == "final":
        c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    else:
        c = json.load(open(os.path.join(OUT, "score_%s.json" % run)))["settings"]
    fig.suptitle("Dubai to Abu Dhabi: Lyzenga SDB against ICESat-2\n(Landsat 8, %s, %s nm%s)"
                 % ("all scenes jointly" if c.get("inputs") == "joint" else "median composite",
                    "/".join(map(str, c["bands"])), ", power-law rescaled beyond the crossover" if c.get("power") else ""),
                 fontsize=12)
    out = os.path.join(rundir, "correlation.png")
    fig.savefig(out, dpi=150)
    print(out)


def reference_maps():
    """The ICESat-2 reference as cells coloured by depth: (a) every cell, calibration and test
    halves; (b) only the training cells (the calibration half after RATIO_FILTER), with the cells
    the filter drops (runs/gulf_lyzenga/final/icesat2_reference.png)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lon, lat = grid()
    ref = np.load(os.path.join(OUT, "reference.npz"))
    water = water_mask()
    cal, test = ref["calibration"], ref["test"]
    train = training(cal)
    dropped = np.isfinite(cal) & ~np.isfinite(train)
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.4), constrained_layout=True)
    extent = (lon[0], lon[-1], lat[0], lat[-1])
    panels = [("(a) ICESat-2 ATL24, every cell: calibration %d + test %d" % (np.isfinite(cal).sum(), np.isfinite(test).sum()),
               np.where(np.isfinite(cal), cal, test), None),
              ("(b) Training cells only: the calibration half after the ratio filter (%d; %d dropped, black)" % (
                  np.isfinite(train).sum(), dropped.sum()), train, dropped)]
    for ax, (title, z, drop) in zip(axes, panels):
        ax.imshow(np.where(water, 1.0, np.nan), origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=8,
                  interpolation="nearest")
        ii, jj = np.nonzero(np.isfinite(z))
        sc = ax.scatter(lon[jj], lat[ii], c=z[ii, jj], s=2.5, cmap="viridis_r", vmin=0, vmax=25, linewidths=0)
        if drop is not None:
            ii, jj = np.nonzero(drop)
            ax.scatter(lon[jj], lat[ii], c="k", s=5, marker="x", linewidths=0.6)
        ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
        ax.set_xlim(lon[0], lon[-1])
        ax.set_ylim(lat[0], lat[-1])
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
    fig.colorbar(sc, ax=axes, shrink=0.75, label="depth below MSL (m)")
    fig.suptitle("Dubai to Abu Dhabi: the ICESat-2 reference (cells of 0.0003\u00b0, the median of at least 3 photons "
                 "spreading at most 0.5 m)", fontsize=12)
    out = os.path.join(WORK, "final", "icesat2_reference.png")
    fig.savefig(out, dpi=130)
    print(out)


def sdb_maps():
    """The chosen model's depth over every water cell, (a) the depth and (b) its 1-sigma error
    estimate (ERROR ... NOISE GAIN), cropped to the water (runs/gulf_lyzenga/final/sdb_map.png)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tier3_scenes
    lon, lat = grid()
    rundir = os.path.join(WORK, "final")
    water = water_mask()
    H = -tier3_scenes.read_nc(os.path.join(rundir, "H.nc"))
    S = tier3_scenes.read_nc(os.path.join(rundir, "H_error.nc"))
    H[~water], S[~water] = np.nan, np.nan
    S[~np.isfinite(H)] = np.nan
    c = json.load(open(os.path.join(OUT, "cv.json") if os.path.exists(os.path.join(OUT, "cv.json")) else CV_PALM))["chosen"]
    fig, axes = plt.subplots(2, 1, figsize=(10, 11), constrained_layout=True)
    extent = (lon[0], lon[-1], lat[0], lat[-1])
    ii, jj = np.nonzero(np.isfinite(H))
    pad = 20*CELL
    panels = [("(a) Depth below MSL", H, "viridis_r", 0.0, 20.0, "depth below MSL (m)"),
              ("(b) Its 1$\\sigma$ error estimate (ERROR NOISE GAIN)", S, "magma_r", 0.0, 3.0, "1$\\sigma$ error (m)")]
    for ax, (title, z, cmap, lo, hi, label) in zip(axes, panels):
        ax.set_facecolor("0.85")
        im = ax.imshow(z, origin="lower", extent=extent, cmap=cmap, vmin=lo, vmax=hi, interpolation="nearest")
        ax.set_aspect(1.0/np.cos(np.radians(np.mean(lat))))
        ax.set_xlim(lon[jj.min()] - pad, lon[jj.max()] + pad)
        ax.set_ylim(lat[ii.min()] - pad, lat[ii.max()] + pad)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        fig.colorbar(im, ax=ax, shrink=0.85, label=label, extend="max")
    fig.suptitle("Dubai to Abu Dhabi: Lyzenga SDB (Landsat 8, %s, %s nm, radius %d; land grey)"
                 % ("all scenes jointly" if c.get("inputs") == "joint" else "median composite",
                    "/".join(map(str, c["bands"])), c["radius"]), fontsize=12)
    out = os.path.join(rundir, "sdb_map.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(out, int(np.isfinite(H).sum()), "cells; depth max %.2f m; error median, 90th, 99th %s m"
          % (np.nanmax(H), np.round(np.nanpercentile(S, [50, 90, 99]), 2)))


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
    commands = {"prepare": prepare, "cv": cv, "fit": fit, "correlation": correlation, "points": points, "calibrate": calibrate, "errors": errors,
                "ratio_filter": ratio_filter, "reference_maps": reference_maps, "sdb_maps": sdb_maps,
                "attenuation": attenuation, "consistent": consistent,
                "consistent_plots": consistent_index_plots}
    if len(sys.argv) == 3 and sys.argv[1] == "attenuation":
        attenuation([int(b) for b in sys.argv[2].split(",")])
    elif len(sys.argv) == 3 and sys.argv[1] == "consistent":
        consistent([[int(b) for b in sys.argv[2].split(",")]])
    elif len(sys.argv) == 6 and sys.argv[1] == "fit":
        # fit composite|joint 443,483,561,655 /radius/ power|nopower
        fit((sys.argv[2], [int(b) for b in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5] == "power"))
    elif len(sys.argv) == 3 and sys.argv[1] == "correlation":
        correlation(sys.argv[2])
    elif len(sys.argv) == 2 and sys.argv[1] in commands:
        commands[sys.argv[1]]()
    else:
        raise SystemExit(__doc__)
