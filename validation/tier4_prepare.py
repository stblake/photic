"""Tier 4: prepare Exmouth Gulf inputs for photic.

  acolite    ACOLITE (DSF, residual glint correction) on each Sentinel-2 scene, over the
             bounding box of the validation sites, at 20 m. Scene-mean geometry
             (output_geometry False) keeps ACOLITE on its fast standard LUT path; the
             sites span < 60 km, so sun and view angles barely vary.
  sites      for each site: resample the ACOLITE Rrs bands to a regular WGS84 grid, write
             photic .asc bands and a command file (one SCENE per image, with HTIDE from
             the Exmouth 02 gauge), and grid the DoT singlebeam soundings as the
             reference depth.

Vertical datums: DoT soundings are relative to AHD; the gauge is relative to LAT. HTIDE
is the observed level minus the gauge's long-term mean (MSL), so depths and tides are
both relative to MSL ~ AHD (the AHD-MSL offset, typically < 0.1 m, is not corrected).
"""

import glob
import json
import os
import sys

import numpy as np

DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))

SITES = {  # west, east, south, north, cell size (deg; 0.0002 ~ 21 m)
    "marina": (114.120, 114.170, -21.985, -21.925, 0.0002),
    "bundegi": (114.160, 114.205, -21.845, -21.800, 0.0002),
    "tantabiddi": (113.955, 113.995, -21.925, -21.880, 0.0002),
    "gulf_entrance": (114.340, 114.440, -21.730, -21.630, 0.0006),  # sparse 1997 lines, deep
}
LIMIT = [-21.99, 113.95, -21.62, 114.45]  # S, W, N, E for ACOLITE


def run_acolite(safe, out):
    sys.path.insert(0, os.path.join(DATA, "tools", "acolite"))
    import acolite as ac
    settings = dict(inputfile=safe, output=out, limit=LIMIT, s2_target_res=20,
                    output_geometry=False, ancillary_data=False, dsf_aot_estimate="fixed",
                    dsf_residual_glint_correction=True, l2w_parameters=["Rrs_*"],
                    rgb_rhot=False, rgb_rhos=False, map_l2w=False,
                    l1r_delete_netcdf=True, l2r_delete_netcdf=True)
    ac.acolite.acolite_run(settings=settings)


def acolite_all():
    for safe in sorted(glob.glob(os.path.join(DATA, "raw", "s2", "*.SAFE"))):
        name = os.path.basename(safe)[:-5]
        out = os.path.join(DATA, "derived", "tier4", "acolite", name)
        if glob.glob(os.path.join(out, "*L2W.nc")):
            print("done:", name)
            continue
        print("ACOLITE:", name, flush=True)
        run_acolite(safe, out)


def site_grid(site):
    """Node-registered WGS84 grid (photic convention: row 0 south) for a site."""
    w, e, s, n, cell = SITES[site]
    ncols = int(round((e - w)/cell)) + 1
    nrows = int(round((n - s)/cell)) + 1
    lon = w + cell*np.arange(ncols)
    lat = s + cell*np.arange(nrows)
    return lon, lat, cell


def reference_depth(site, min_count=3, max_std=0.5):
    """Median depth below AHD (m, positive down) of the DoT singlebeam soundings in each
    cell, using only the most recent survey that covers the cell; cells with fewer than
    min_count soundings or a spread above max_std (steep or mixed) are nodata."""
    P = np.load(os.path.join(DATA, "derived", "tier4", "soundings_lonlat_zAHD_year.npy"))
    lon, lat, cell = site_grid(site)
    j = np.rint((P[:, 0] - lon[0])/cell).astype(int)
    i = np.rint((P[:, 1] - lat[0])/cell).astype(int)
    m = (i >= 0) & (i < len(lat)) & (j >= 0) & (j < len(lon)) & (P[:, 2] < 0)
    i, j, z, yr = i[m], j[m], -P[m, 2], P[m, 3]
    depth = np.full((len(lat), len(lon)), np.nan)
    count = np.zeros(depth.shape, dtype=int)
    year = np.zeros(depth.shape, dtype=int)
    order = np.lexsort((-yr, j, i))
    i, j, z, yr = i[order], j[order], z[order], yr[order]
    key = i*len(lon) + j
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    ends = np.r_[starts[1:], len(key)]
    for a, b in zip(starts, ends):
        latest = yr[a]
        sel = z[a:b][yr[a:b] == latest]
        if len(sel) >= min_count and np.std(sel) <= max_std:
            depth[i[a], j[a]] = np.median(sel)
            count[i[a], j[a]] = len(sel)
            year[i[a], j[a]] = int(latest)
    return depth, count, year


PHOTIC_BANDS = [443.0, 492.0, 560.0, 665.0]  # nominal S2 B1-B4 centres
RED_EDGE_BANDS = PHOTIC_BANDS + [704.0]       # plus B5 (EXPERIMENTS.md B7)
N_WATER = 1.34


def gauge():
    """Exmouth 02 observed water levels: {YYYYMMDD.hhmm (AWST): cm above LAT}, and MSL."""
    levels = {}
    for f in glob.glob(os.path.join(DATA, "raw", "tides", "EXEXM02_0*", "*.txt")):
        for line in open(f, errors="replace"):
            if "," not in line or line.startswith("#") or line.startswith("H("):
                continue
            h, t = line.strip().rstrip("/").split(",")
            h = int(h)
            if h != -9999:
                levels[t] = h
    return levels, float(np.mean(list(levels.values())))


def scene_tide(isodate, levels, msl):
    """Water level above MSL (m) at the overpass, from the nearest 5-minute record."""
    import datetime
    t = datetime.datetime.fromisoformat(isodate.replace("Z", "+00:00")) + datetime.timedelta(hours=8)
    for dm in [0, 5, -5, 10, -10]:
        u = t + datetime.timedelta(minutes=dm)
        key = u.strftime("%Y%m%d.%H") + "%02d" % (u.minute - u.minute % 5)
        if key in levels:
            return (levels[key] - msl)/100.0
    return None


def read_l2w(path):
    import netCDF4
    with netCDF4.Dataset(path) as d:
        g = {k: d.getncattr(k) for k in d.ncattrs()}
        lon = d.variables["lon"][:].astype(float)
        lat = d.variables["lat"][:].astype(float)
        rrs = {}
        for v in d.variables:
            if v.startswith("Rrs_"):
                rrs[float(v[4:])] = np.ma.filled(d.variables[v][:].astype(float), np.nan)
    return g, np.ma.filled(lon, np.nan), np.ma.filled(lat, np.nan), rrs


def write_asc(path, grid, lon, lat, cell, nodata=-9999.0):
    g = np.where(np.isfinite(grid), grid, nodata)
    with open(path, "w") as f:
        f.write("ncols %d\nnrows %d\nxllcorner %.7f\nyllcorner %.7f\ncellsize %.7f\nNODATA_value %g\n"
                % (len(lon), len(lat), lon[0], lat[0], cell, nodata))
        for i in range(len(lat) - 1, -1, -1):
            f.write(" ".join("%.7g" % v for v in g[i]) + "\n")


def prepare_sites(wavelengths=PHOTIC_BANDS, sites_name="sites"):
    """photic inputs and reference grids for every site from all ACOLITE scenes, written to
    derived/tier4/<sites_name>/ (tier4_score.py reads $PHOTIC_TIER4_SITES, default "sites")."""
    from scipy.spatial import cKDTree
    levels, msl = gauge()
    scenes = []
    for l2w in sorted(glob.glob(os.path.join(DATA, "derived", "tier4", "acolite", "*", "*L2W.nc"))):
        g, lon, lat, rrs = read_l2w(l2w)
        ok = np.isfinite(lon) & np.isfinite(lat)
        tree = cKDTree(np.column_stack([lon[ok], lat[ok]]))
        flat = {w: v[ok] for w, v in rrs.items()}
        tide = scene_tide(g["isodate"], levels, msl)
        if tide is None:
            print("no tide for", l2w)
            continue
        scenes.append(dict(file=l2w, isodate=g["isodate"], sza=float(g["sza"]), vza=float(g["vza"]),
                           tide=tide, tree=tree, rrs=flat))
    meta = dict(msl_cm_above_lat=msl, scenes=[{k: s[k] for k in ("file", "isodate", "sza", "vza", "tide")}
                                              for s in scenes], sites={})
    for site in SITES:
        lon, lat, cell = site_grid(site)
        out = os.path.join(DATA, "derived", "tier4", sites_name, site)
        os.makedirs(out, exist_ok=True)
        LON, LAT = np.meshgrid(lon, lat)
        # Resample every scene; photic needs every scene valid at a pixel, so drop scenes
        # that cover less than half as much of the site as the best one (swath edges).
        resampled = []
        for s in scenes:
            dist, idx = s["tree"].query(np.column_stack([LON.ravel(), LAT.ravel()]))
            far = (dist > 1.5*0.0002).reshape(LON.shape)
            bands = {}
            for wl in wavelengths:
                src = min(s["rrs"], key=lambda w: abs(w - wl))
                band = s["rrs"][src][idx].reshape(LON.shape)
                band[far] = np.nan
                band[band <= 0] = np.nan
                bands[wl] = band
            cover = float(np.mean(np.all([np.isfinite(b) for b in bands.values()], axis=0)))
            resampled.append((s, bands, cover))
        best = max(c for _, _, c in resampled)
        kept = [(s, b) for s, b, c in resampled if c >= 0.5*best]
        for f in glob.glob(os.path.join(out, "s*_b*.asc")):
            os.remove(f)
        lines, scene_lines = [], []
        for k, (s, bands) in enumerate(kept):
            names = []
            for wl in wavelengths:
                name = "s%d_b%d" % (k + 1, wl)
                write_asc(os.path.join(out, name + ".asc"), bands[wl], lon, lat, cell)
                lines.append("READ GRID %s %s.asc" % (name, name))
                names.append(name)
            tw = np.degrees(np.arcsin(np.sin(np.radians(s["sza"]))/N_WATER))
            tv = np.degrees(np.arcsin(np.sin(np.radians(s["vza"]))/N_WATER))
            scene_lines.append("SCENE NAME s%d BANDS %s WAVELENGTHS %s THETAVIEW %.4f THETASUN %.4f HTIDE %.3f" % (
                k + 1, " ".join(names), " ".join("%d" % w for w in wavelengths), tv, tw, s["tide"]))
        model = "MODEL Lee_Semi_Analytical IN %s OUT H H_sigma H_error albedo sand seagrass coral kmin btype iod" % (
            " ".join("s%d" % (k + 1) for k in range(len(kept))))
        open(os.path.join(out, "run.bam"), "w").write("\n".join(lines + scene_lines + [model, "EXIT"]) + "\n")
        depth, count, year = reference_depth(site)
        np.savez(os.path.join(out, "reference.npz"), depth=depth, count=count, year=year, lon=lon, lat=lat)
        meta["sites"][site] = dict(shape=list(depth.shape), cell=cell, reference_cells=int(np.isfinite(depth).sum()),
                                   scenes=[s["isodate"] for s, _ in kept],
                                   coverage=[round(c, 3) for _, _, c in resampled])
        print("%-14s %dx%d, %d of %d scenes, %d reference cells" % (site, len(lon), len(lat), len(kept), len(scenes),
                                                             meta["sites"][site]["reference_cells"]))
    json.dump(meta, open(os.path.join(DATA, "derived", "tier4", sites_name, "meta.json"), "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1:] == ["acolite"]:
        acolite_all()
    elif sys.argv[1:] == ["sites"]:
        prepare_sites()
    elif sys.argv[1:] == ["sites5"]:
        prepare_sites(RED_EDGE_BANDS, "sites5")
    else:
        raise SystemExit("usage: tier4_prepare.py acolite|sites|sites5")
