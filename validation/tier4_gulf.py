"""Tier 4, second region: the Arabian Gulf between Abu Dhabi and Dubai, scored against ICESat-2.

A held-out test of sounding-free models: nothing here was used to choose any setting. Landsat 8
only: few scenes here are clear enough to model, and the Sentinel-2 scenes first chosen (by cloud,
sun and tide) were too hazy for ACOLITE's SWIR water mask.

  screen     SWIR over open water for every candidate scene (haze): gulf_screen_l8.json
  select     the scenes, by the rule below, as Collection 1 names for tier4_fetch_l8.GULF_SCENES
  reference  the ATL24 seabed photons (tier4_fetch_atl24.py) as depths below MSL:
             derived/tier4/atl24_gulf.npz (lon, lat, depth, year, granule: the pass index)
  acolite    ACOLITE (as tier4_prepare.py: DSF, fixed AOT, residual glint) on the gulf scenes
             (tier4_fetch_l8.py gulf) over the site's box
  sites      photic inputs on 0.1 deg tiles: sites_gulf_l8

Scene selection, in order: haze first (the six clearest by SWIR surface reflectance over open
water), among scenes with the sun at 38-56 deg (lower: little light reaches depth; higher: sun
glint) and a tide changing by at most 0.2 m/h over the hour around the overpass.

Photons: class 40 (seabed), not low_confidence_flag and not sensor_depth_exceeded, and consistent
with the other passes nearby (consistent(): some passes classify a near-surface layer as seabed). ATL24's
ortho_h is refraction-corrected and relative to the EGM2008 geoid. Mean sea level is found from
the passes themselves: each pass's sea surface (surface_h, median over its seabed photons) minus
the GOT4.10 tide at that time, averaged over all passes (GOT4.10 follows the passes' surface with
r = 0.86 and a non-tidal residual of 0.20 m). The depth below MSL is then
msl - ortho_h, without each pass's non-tidal residual.

Scenes' HTIDE is GOT4.10 (pyTMD, cutoff 50 km) at the tile's centre at the overpass. Tiles are
0.1 deg squares over the box, kept when they hold at least MIN_REFERENCE reference cells; a
reference cell is the median of at least 3 photons with a spread of at most 0.5 m (the rule for
the Exmouth soundings), over water in the imagery.
"""

import datetime
import glob
import json
import os
import sys

import numpy as np

import tier4_prepare as tp

DATA = tp.DATA
# The region, from $GULF_REGION: "palm", the box around the Palm Jebel Ali first studied (the
# default), or "scene", the whole Landsat 8 scene (path 160 row 43). Each has its own data.
REGIONS = {"palm": [24.75, 54.7, 25.1, 55.2], "scene": [23.47, 54.05, 25.63, 56.37]}  # S, W, N, E
# The scene's water (the coast from Abu Dhabi to Sharjah and the open Gulf), for its imagery: the
# whole scene box at 0.0003 deg would be 56 M cells per band.
WATER_LIMITS = {"scene": [24.0, 54.05, 25.63, 55.7]}
REGION = os.environ.get("GULF_REGION", "palm")
if REGION not in REGIONS:
    raise SystemExit("GULF_REGION must be one of %s" % ", ".join(REGIONS))
SUFFIX = "" if REGION == "palm" else "_" + REGION
ATL24 = os.path.join(DATA, "raw", "atl24" + SUFFIX)
REFERENCE = os.path.join(DATA, "derived", "tier4", "atl24_gulf%s.npz" % SUFFIX)
ACOLITE = "acolite_gulf_l8" + SUFFIX
GOT = os.path.join(DATA, "raw", "tides", "got")
LIMIT = REGIONS[REGION]  # S, W, N, E for ACOLITE and the ICESat-2 crop
TILE = 0.1
CELL = tp.L8_CELL  # ~31 m
MIN_REFERENCE = 200
EPOCH_2000 = datetime.datetime(2000, 1, 1)
ATLAS_EPOCH = datetime.datetime(2018, 1, 1)  # ATL24 delta_time (GPS seconds; 18 s from UTC ignored)


def got_tide(lon, lat, seconds_since_2000):
    """GOT4.10 tide (m above MSL) at points and times."""
    import pyTMD.compute
    lon, lat, t = (np.atleast_1d(np.asarray(x, dtype=float)) for x in (lon, lat, seconds_since_2000))
    return np.asarray(pyTMD.compute.tide_elevations(lon, lat, t, directory=GOT, model="GOT4.10", type="drift",
                                                    extrapolate=True, cutoff=50)).ravel()


CONSISTENCY_RADIUS = 0.003  # deg (~300 m): other passes' cells compared with each pass's cell
CONSISTENCY_CELL = 0.0003    # deg: cells of one pass (median depth) for the comparison
CONSISTENCY_MIN = 5          # fewest neighbouring cells of other passes for a comparison
CONSISTENCY_TOL = (3.0, 0.3)  # m and fraction: a cell disagreeing by more than max of these is dropped


def consistent(lon, lat, depth, granule):
    """Photons whose pass agrees with the other passes nearby. Some ATL24 passes classify a layer
    about 2 m below the surface as seabed across deep water (passes 6-7.5 m shallower than every
    pass crossing them, with normal confidence, by night and by day). Each pass's cells (median
    depth) are compared with the median of other passes' cells within CONSISTENCY_RADIUS; a cell
    disagreeing by more than CONSISTENCY_TOL is dropped. Cells without enough neighbours are kept."""
    from scipy.spatial import cKDTree
    key = (np.rint(lon/CONSISTENCY_CELL).astype(np.int64)*1000000 + np.rint(lat/CONSISTENCY_CELL).astype(np.int64))*1000 + granule
    u, inv = np.unique(key, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    starts = np.flatnonzero(np.r_[True, inv[order][1:] != inv[order][:-1]])
    ends = np.r_[starts[1:], len(order)]
    D = np.array([np.median(depth[order[a:b]]) for a, b in zip(starts, ends)])
    X = np.array([lon[order[a]] for a in starts])
    Y = np.array([lat[order[a]] for a in starts])
    G = np.array([granule[order[a]] for a in starts])
    tree = cKDTree(np.column_stack([X, Y]))
    bad = np.zeros(len(D), dtype=bool)
    for c, nb in enumerate(tree.query_ball_point(np.column_stack([X, Y]), CONSISTENCY_RADIUS)):
        nb = [n for n in nb if G[n] != G[c]]
        if len(nb) >= CONSISTENCY_MIN:
            ref = np.median(D[nb])
            bad[c] = abs(D[c] - ref) > max(CONSISTENCY_TOL[0], CONSISTENCY_TOL[1]*ref)
    keep = np.ones(len(depth), dtype=bool)
    for c in np.flatnonzero(bad):
        keep[order[starts[c]:ends[c]]] = False
    return keep, int(bad.sum()), len(D)


def reference():
    offset = (ATLAS_EPOCH - EPOCH_2000).total_seconds()
    passes, photons, granules = [], [], []
    for f in sorted(glob.glob(os.path.join(ATL24, "*.npz"))):
        d = np.load(f)
        if "class_ph" not in d.files:
            continue
        m = (d["class_ph"] == 40) & (d["low_confidence_flag"] == 0) & (d["sensor_depth_exceeded"] == 0)
        if m.sum() < 20:
            continue
        t = float(np.median(d["delta_time"][m])) + offset
        passes.append((np.median(d["lon_ph"][m]), np.median(d["lat_ph"][m]), t, np.median(d["surface_h"][m])))
        year = (ATLAS_EPOCH + datetime.timedelta(seconds=t - offset)).year
        photons.append(np.column_stack([d["lon_ph"][m], d["lat_ph"][m], d["ortho_h"][m], np.full(m.sum(), year)]))
        granules.append(np.full(m.sum(), len(passes) - 1))
    P = np.array(passes)
    tide = got_tide(P[:, 0], P[:, 1], P[:, 2])
    msl = float(np.mean(P[:, 3] - tide))
    r = float(np.corrcoef(P[:, 3], tide)[0, 1])
    X = np.concatenate(photons)
    G = np.concatenate(granules)
    keep, n_bad, n_cells = consistent(X[:, 0], X[:, 1], msl - X[:, 2], G)
    np.savez(REFERENCE, lon=X[keep, 0], lat=X[keep, 1], depth=msl - X[keep, 2], year=X[keep, 3].astype(int),
             granule=G[keep], msl_above_geoid=msl, n_passes=len(P), tide_r=r,
             nontidal_sd=float(np.std(P[:, 3] - tide)), dropped_inconsistent=int((~keep).sum()))
    print("%d passes, %d photons; MSL %.3f m above EGM2008; GOT4.10 r %.2f, non-tidal sd %.3f m"
          % (len(P), len(X), msl, r, np.std(P[:, 3] - tide)))
    print("consistency with other passes: %d of %d pass cells dropped (%d photons, %.1f%%)"
          % (n_bad, n_cells, (~keep).sum(), 100.0*(~keep).mean()))
    for gnum in np.unique(G):
        f = np.mean(~keep[G == gnum])
        if f > 0.2:
            print("  pass %d: %.0f%% of its photons dropped" % (gnum, 100*f))


def tiles():
    """{name: (W, E, S, N, cell)} for every 0.1 deg tile of the box."""
    s, w, n, e = LIMIT
    out = {}
    for x in np.arange(w, e - 1e-9, TILE):
        for y in np.arange(s, n - 1e-9, TILE):
            out["t%d_%d" % (round(x*100), round(y*100))] = (round(x, 4), round(min(x + TILE, e), 4),
                                                            round(y, 4), round(min(y + TILE, n), 4), CELL)
    return out


def binned(lon, lat, cell, min_count=3, max_std=0.5):
    """The median depth of the photons in each cell (nodata with too few or too spread)."""
    R = np.load(REFERENCE)
    j = np.rint((R["lon"] - lon[0])/cell).astype(int)
    i = np.rint((R["lat"] - lat[0])/cell).astype(int)
    m = (i >= 0) & (i < len(lat)) & (j >= 0) & (j < len(lon))
    i, j, z = i[m], j[m], R["depth"][m]
    depth = np.full((len(lat), len(lon)), np.nan)
    count = np.zeros(depth.shape, dtype=int)
    key = i*len(lon) + j
    order = np.argsort(key, kind="stable")
    key, z = key[order], z[order]
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    for a, b in zip(starts, np.r_[starts[1:], len(key)]):
        sel = z[a:b]
        if len(sel) >= min_count and np.std(sel) <= max_std and np.median(sel) > 0:
            depth.flat[key[a]] = np.median(sel)
            count.flat[key[a]] = len(sel)
    return depth, count


def atl24_reference(site, lon, lat, cell):
    depth, count = binned(lon, lat, cell)
    return dict(depth=depth, count=count)


def tide_for(sites):
    """tide(isodate) at the centre of the sites' box."""
    w = min(v[0] for v in sites.values()); e = max(v[1] for v in sites.values())
    s = min(v[2] for v in sites.values()); n = max(v[3] for v in sites.values())
    def fn(isodate):
        t = datetime.datetime.fromisoformat(isodate.replace("Z", "+00:00")).replace(tzinfo=None)
        return float(got_tide((w + e)/2, (s + n)/2, (t - EPOCH_2000).total_seconds())[0])
    return fn


WATER_BOX = (54.75, 25.0, 54.85, 25.08)  # W, S, N, E: open water 13-20 m deep, for the haze screen


SUN = (38.0, 56.0)
MAX_TIDE_RATE = 0.2  # m over the hour around the overpass
N_SCENES = 6


def screen_items():
    """STAC items over the water box with little cloud: Landsat 8 C2 L2 (Planetary Computer, path
    160 row 43, 2018-10 to 2021, the years of Collection 1 on Google Cloud, which ACOLITE reads)."""
    import urllib.request
    url = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
    body = dict(collections=["landsat-c2-l2"], datetime="2018-10-01T00:00:00Z/2021-12-31T23:59:59Z",
                query={"eo:cloud_cover": {"lt": 5}, "platform": {"eq": "landsat-8"},
                       "landsat:wrs_path": {"eq": "160"}, "landsat:wrs_row": {"eq": "043"}})
    w, s, e, n = WATER_BOX
    body.update(intersects={"type": "Point", "coordinates": [(w + e)/2, (s + n)/2]}, limit=200)
    items = []
    while True:
        req = urllib.request.Request(url, json.dumps(body).encode(), {"Content-Type": "application/json"})
        page = json.load(urllib.request.urlopen(req))
        items += page["features"]
        nxt = [l for l in page.get("links", []) if l.get("rel") == "next"]
        if not nxt or not page["features"]:
            return items
        body = nxt[0].get("body", body) if nxt[0].get("method", "GET") == "POST" else body
        if nxt[0].get("method", "GET") != "POST":
            url, body = nxt[0]["href"], None
            page_items = []
            while url:
                p = json.load(urllib.request.urlopen(url))
                items += p["features"]
                url = next((l["href"] for l in p.get("links", []) if l.get("rel") == "next"), None) if p["features"] else None
            return items


def window_median(href, asset):
    """Median reflectance of an asset over the water box (a ranged read of a cloud-optimised GeoTIFF)."""
    import rasterio
    from rasterio.warp import transform
    band = (asset.get("raster:bands") or [{}])[0]
    with rasterio.open(href) as src:
        xs, ys = transform("EPSG:4326", src.crs, [WATER_BOX[0], WATER_BOX[2]], [WATER_BOX[1], WATER_BOX[3]])
        win = rasterio.windows.from_bounds(min(xs), min(ys), max(xs), max(ys), src.transform)
        a = src.read(1, window=win, boundless=True, fill_value=band.get("nodata", 0)).astype(float)
    ok = a != band.get("nodata", 0)
    if ok.mean() < 0.9:
        return None
    return float(np.median(a[ok]))*band.get("scale", 1.0) + band.get("offset", 0.0)


def screen():
    """SWIR (1.6 um) surface reflectance over open water for every candidate scene: water is black
    there, so what remains is haze, glint or adjacency that the processor did not remove (ACOLITE
    masks water where TOA reflectance at 1.6 um exceeds 0.0215). Also the tide and its rate at
    the overpass and the sun elevation, for the selection."""
    from concurrent.futures import ThreadPoolExecutor
    import urllib.request
    items = screen_items()
    token = "?" + json.load(urllib.request.urlopen(
        "https://planetarycomputer.microsoft.com/api/sas/v1/token/landsateuwest/landsat-c2"))["token"]

    def one(it):
        a = it["assets"]["swir16"]
        try:
            v = window_median(a["href"] + token, a)
        except Exception as ex:  # noqa: BLE001 (a missing or unreadable asset: skip the scene)
            return dict(id=it["id"], error=str(ex)[:100])
        p = it["properties"]
        return dict(id=it["id"], datetime=p["datetime"],
                    cloud=p.get("eo:cloud_cover"), sun_elevation=p.get("view:sun_elevation"), swir=v)

    with ThreadPoolExecutor(max_workers=16) as pool:
        rows = [r for r in pool.map(one, items)]
    ok = [r for r in rows if r.get("swir") is not None]
    times = [datetime.datetime.fromisoformat(r["datetime"].replace("Z", "+00:00")).replace(tzinfo=None) for r in ok]
    t = np.array([(x - EPOCH_2000).total_seconds() for x in times])
    w, s, e, n = WATER_BOX
    lon, lat = np.full(t.size, (w + e)/2), np.full(t.size, (s + n)/2)
    h = got_tide(lon, lat, t)
    rate = got_tide(lon, lat, t + 1800) - got_tide(lon, lat, t - 1800)
    for r, a, b in zip(ok, h, rate):
        r.update(tide=float(a), tide_rate=float(b))
    out = os.path.join(DATA, "derived", "tier4", "gulf_screen_l8.json")
    json.dump(sorted(rows, key=lambda r: r["datetime"]), open(out, "w"), indent=1)
    v = np.array([r["swir"] for r in ok])
    print("%d items, %d read; SWIR over water: %s (5/25/50/75/95%%)"
          % (len(rows), len(ok), np.round(np.percentile(v, [5, 25, 50, 75, 95]), 4)))


def select():
    """The six clearest screened scenes with the sun and tide rate in range, as Collection 1
    product names on Google Cloud (tier4_fetch_l8.GULF_SCENES)."""
    import urllib.request
    rows = json.load(open(os.path.join(DATA, "derived", "tier4", "gulf_screen_l8.json")))
    ok = [r for r in rows if r.get("swir") is not None and SUN[0] <= r["sun_elevation"] <= SUN[1]
          and abs(r["tide_rate"]) <= MAX_TIDE_RATE]
    chosen = sorted(ok, key=lambda r: r["swir"])[:N_SCENES]
    listing = json.load(urllib.request.urlopen("https://storage.googleapis.com/storage/v1/b/gcp-public-data-landsat/o?"
                                               "prefix=LC08/01/160/043/&delimiter=/&maxResults=1000"))["prefixes"]
    for r in sorted(chosen, key=lambda r: r["tide"]):
        date = r["id"].split("_")[3]
        c1 = sorted(p.split("/")[-2] for p in listing if p.split("/")[-2].split("_")[3] == date
                    and p.split("/")[-2].endswith("_T1"))
        print('    "%s",  # SWIR %.4f, %.0f deg, %+.2f m' % (c1[-1] if c1 else "none for " + date, r["swir"],
                                                          r["sun_elevation"], r["tide"]))


def acolite():
    if REGION == "palm":
        tp.acolite_all("l8_gulf", ACOLITE, limit=LIMIT)
        return
    # The whole scene: ACOLITE reprojects L2W to one regular WGS84 grid (0.0003 deg, the same
    # limit for every scene), so the scenes stack without resampling; only the four bands photic
    # uses are kept, and the unprojected L2W is deleted (about 0.5 GB per scene remains).
    from concurrent.futures import ProcessPoolExecutor
    limit = WATER_LIMITS[REGION]
    jobs = [d for d in sorted(glob.glob(os.path.join(DATA, "raw", "l8_gulf", "*"))) if os.path.isdir(d)]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for out in pool.map(acolite_projected, jobs, [limit]*len(jobs)):
            print(out, flush=True)


def acolite_projected(scene, limit):
    import contextlib
    import io
    out = os.path.join(DATA, "derived", "tier4", ACOLITE, os.path.basename(scene))
    if glob.glob(os.path.join(out, "*L2W_projected.nc")):
        return out + " (done)"
    with contextlib.redirect_stdout(io.StringIO()):
        tp.run_acolite(scene, out, limit, l2w_parameters=["Rrs_443", "Rrs_483", "Rrs_561", "Rrs_655"],
                       output_projection=True, output_projection_epsg="4326",
                       output_projection_resolution=[tp.L8_CELL, tp.L8_CELL], output_projection_limit=limit,
                       reproject_outputs="L2W")
    for f in glob.glob(os.path.join(out, "*.nc")):
        if not f.endswith("L2W_projected.nc"):
            os.remove(f)
    return out


def sites():
    if REGION != "palm":
        raise SystemExit("sites (DROP's tiles) are built for the palm region only")
    wavelengths, name, acolite_dir, min_cell = tp.L8_BANDS, "sites_gulf_l8", ACOLITE, tp.L8_CELL
    # Tiles are chosen from the reference alone (before any model is run), one per tile with
    # its own centre's tide.
    chosen = {}
    for t, box in tiles().items():
        lon, lat, cell = tp.site_grid(t, {t: box}, min_cell)
        n = int(np.isfinite(binned(lon, lat, cell)[0]).sum())
        if n >= MIN_REFERENCE:
            chosen[t] = box
    root = os.path.join(DATA, "derived", "tier4", name)
    metas = {}
    for t, box in chosen.items():
        tp.prepare_sites(wavelengths, name, acolite=acolite_dir, sites={t: box}, min_cell=min_cell,
                         tide=tide_for({t: box}), reference=atl24_reference)
        metas[t] = json.load(open(os.path.join(root, "meta.json")))
        # Reference cells over land or masked water in every scene are not part of the test.
        out = os.path.join(root, t)
        ref = dict(np.load(os.path.join(out, "reference.npz")))
        water = np.zeros(ref["depth"].shape, dtype=bool)
        for f in glob.glob(os.path.join(out, "s*_b%d.asc" % wavelengths[-1])):
            g = np.loadtxt(f, skiprows=6)[::-1]
            water |= g != -9999.0
        ref["depth"] = np.where(water, ref["depth"], np.nan)
        np.savez(os.path.join(out, "reference.npz"), **ref)
        print("  %s: %d reference cells over water" % (t, np.isfinite(ref["depth"]).sum()))
    meta = dict(scenes={t: m["scenes"] for t, m in metas.items()},
                sites={t: m["sites"][t] for t, m in metas.items()})
    json.dump(meta, open(os.path.join(root, "meta.json"), "w"), indent=1)


if __name__ == "__main__":
    cmd = sys.argv[1:]
    commands = {"screen": screen, "select": select, "reference": reference, "acolite": acolite, "sites": sites}
    if len(cmd) == 1 and cmd[0] in commands:
        commands[cmd[0]]()
    else:
        raise SystemExit(__doc__)
