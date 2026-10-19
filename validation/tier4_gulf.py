"""Tier 4, second region: the Arabian Gulf between Abu Dhabi and Dubai, scored against ICESat-2.

A held-out test of sounding-free models: nothing here was used to choose any setting.

  reference      the ATL24 seabed photons (tier4_fetch_atl24.py) as depths below MSL:
                 derived/tier4/atl24_gulf.npz (lon, lat, depth, year)
  acolite s2|l8  ACOLITE (as tier4_prepare.py: DSF, fixed AOT, residual glint) on the gulf scenes
                 (tier4_fetch_s2.py/tier4_fetch_l8.py gulf) over the site's box
  sites s2|l8    photic inputs on 0.1 deg tiles: sites_gulf, sites_gulf_l8

Photons: class 40 (seabed), not low_confidence_flag and not sensor_depth_exceeded. ATL24's
ortho_h is refraction-corrected and relative to the EGM2008 geoid. Mean sea level is found from
the passes themselves: each pass's sea surface (surface_h, median over its seabed photons) minus
the GOT4.10 tide at that time, averaged over all passes (GOT4.10 follows the passes' surface with
r = 0.84, slope 0.85 and a non-tidal residual of 0.22 m). The depth below MSL is then
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
ATL24 = os.path.join(DATA, "raw", "atl24")
REFERENCE = os.path.join(DATA, "derived", "tier4", "atl24_gulf.npz")
GOT = os.path.join(DATA, "raw", "tides", "got")
LIMIT = [24.75, 54.7, 25.1, 55.2]  # S, W, N, E for ACOLITE (tier4_fetch_atl24.BOX)
TILE = 0.1
CELL = 0.0002  # as the Exmouth sites; Landsat 8 uses tp.L8_CELL
MIN_REFERENCE = 200
EPOCH_2000 = datetime.datetime(2000, 1, 1)
ATLAS_EPOCH = datetime.datetime(2018, 1, 1)  # ATL24 delta_time (GPS seconds; 18 s from UTC ignored)


def got_tide(lon, lat, seconds_since_2000):
    """GOT4.10 tide (m above MSL) at points and times."""
    import pyTMD.compute
    lon, lat, t = (np.atleast_1d(np.asarray(x, dtype=float)) for x in (lon, lat, seconds_since_2000))
    return np.asarray(pyTMD.compute.tide_elevations(lon, lat, t, directory=GOT, model="GOT4.10", type="drift",
                                                    extrapolate=True, cutoff=50)).ravel()


def reference():
    offset = (ATLAS_EPOCH - EPOCH_2000).total_seconds()
    passes, photons = [], []
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
    P = np.array(passes)
    tide = got_tide(P[:, 0], P[:, 1], P[:, 2])
    msl = float(np.mean(P[:, 3] - tide))
    r = float(np.corrcoef(P[:, 3], tide)[0, 1])
    X = np.concatenate(photons)
    np.savez(REFERENCE, lon=X[:, 0], lat=X[:, 1], depth=msl - X[:, 2], year=X[:, 3].astype(int),
             msl_above_geoid=msl, n_passes=len(P), tide_r=r, nontidal_sd=float(np.std(P[:, 3] - tide)))
    print("%d passes, %d photons; MSL %.3f m above EGM2008; GOT4.10 r %.2f, non-tidal sd %.3f m"
          % (len(P), len(X), msl, r, np.std(P[:, 3] - tide)))


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


def acolite(sensor):
    raw, derived = ("s2_gulf", "acolite_gulf") if sensor == "s2" else ("l8_gulf", "acolite_gulf_l8")
    tp.acolite_all(raw, derived, limit=LIMIT)


def sites(sensor):
    if sensor == "s2":
        wavelengths, name, acolite_dir, min_cell = tp.PHOTIC_BANDS, "sites_gulf", "acolite_gulf", 0.0
    else:
        wavelengths, name, acolite_dir, min_cell = tp.L8_BANDS, "sites_gulf_l8", "acolite_gulf_l8", tp.L8_CELL
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
    if cmd == ["reference"]:
        reference()
    elif len(cmd) == 2 and cmd[0] in ("acolite", "sites") and cmd[1] in ("s2", "l8"):
        {"acolite": acolite, "sites": sites}[cmd[0]](cmd[1])
    else:
        raise SystemExit(__doc__)
