"""Tier 4: quality indicators for each ACOLITE scene over each site.

Over water the survey says is deeper than DEEP m (where no bottom signal reaches the red or
near-infrared), for each scene and site:
  red     median Rrs(665): a turbid or offset red, which photic reads as a bright, shallow
          bottom (EXPERIMENTS.md, the tier 4 red excess);
  nir     median Rrs(865): residual glint or aerosol (water itself is nearly black here);
  swir    median Rrs(1614): as nir, with no water-leaving signal at all;
  texture standard deviation of Rrs(865): glint from waves, which varies from pixel to pixel.
Also the scene's sun elevation, tide and the tide's rate of change (cm/h over the hour around
the overpass, from the Exmouth gauge), a proxy for tidal currents.

  tier4_scene_quality.py [sites ...]    writes derived/tier4/scene_quality.json and prints a table
"""

import datetime
import glob
import json
import os
import sys

import numpy as np

import tier4_prepare as t4

DEEP = 10.0


def tide_rate(isodate, levels):
    t = datetime.datetime.fromisoformat(isodate.replace("Z", "+00:00")) + datetime.timedelta(hours=8)
    def level(u):
        return levels.get(u.strftime("%Y%m%d.%H") + "%02d" % (u.minute - u.minute % 5))
    a, b = level(t - datetime.timedelta(minutes=30)), level(t + datetime.timedelta(minutes=30))
    return None if a is None or b is None else float(b - a)


def main(sites):
    import netCDF4
    from scipy.spatial import cKDTree
    levels, msl = t4.gauge()
    refs = {}
    for site in sites:
        lon, lat, _ = t4.site_grid(site)
        depth, _, _ = t4.reference_depth(site)
        LON, LAT = np.meshgrid(lon, lat)
        deep = np.isfinite(depth) & (depth > DEEP)
        refs[site] = np.column_stack([LON[deep], LAT[deep]])
    out = []
    for l2w in sorted(glob.glob(os.path.join(t4.DATA, "derived", "tier4", "acolite", "*", "*L2W.nc"))):
        with netCDF4.Dataset(l2w) as d:
            g = {k: d.getncattr(k) for k in d.ncattrs()}
            lon = np.ma.filled(d["lon"][:].astype(float), np.nan)
            lat = np.ma.filled(d["lat"][:].astype(float), np.nan)
            have = {int(v[4:]): v for v in d.variables if v.startswith("Rrs_")}
            # Band centres differ slightly between S2A, S2B and S2C (864/865, 1610/1614).
            bands = {w: np.ma.filled(d[have[min(have, key=lambda x: abs(x - w))]][:].astype(float), np.nan)
                     for w in (665, 865, 1614)}
        ok = np.isfinite(lon) & np.isfinite(lat)
        tree = cKDTree(np.column_stack([lon[ok], lat[ok]]))
        flat = {w: b[ok] for w, b in bands.items()}
        tide = t4.scene_tide(g["isodate"], levels, msl)
        row = dict(scene=os.path.basename(os.path.dirname(l2w)), date=g["isodate"][:10],
                   sun_el=90.0 - float(g["sza"]), tide=tide, tide_rate=tide_rate(g["isodate"], levels), sites={})
        for site, pts in refs.items():
            if not len(pts):
                continue
            dist, idx = tree.query(pts)
            near = dist < 1.5*0.0002
            if near.sum() < 20:
                continue
            v = {w: flat[w][idx[near]] for w in flat}
            v = {w: x[np.isfinite(x)] for w, x in v.items()}
            row["sites"][site] = dict(n=int(near.sum()), red=float(np.median(v[665])), nir=float(np.median(v[865])),
                                      swir=float(np.median(v[1614])), texture=float(np.std(v[865])))
        out.append(row)
    json.dump(out, open(os.path.join(t4.DATA, "derived", "tier4", "scene_quality.json"), "w"), indent=1)
    print("%-10s %6s %6s %6s  %s" % ("date", "sun_el", "tide", "cm/h", "  ".join(
        "%-30s" % (s + ": red nir swir tex (1e-4)") for s in sites)))
    for r in out:
        cells = []
        for s in sites:
            q = r["sites"].get(s)
            cells.append("%-30s" % ("-" if q is None else "%5.1f %5.1f %5.1f %5.1f" % tuple(
                1e4*q[k] for k in ("red", "nir", "swir", "texture"))))
        print("%-10s %6.1f %+6.2f %+6.0f  %s" % (r["date"], r["sun_el"], r["tide"] if r["tide"] is not None else np.nan,
                                                r["tide_rate"] if r["tide_rate"] is not None else np.nan, "  ".join(cells)))


if __name__ == "__main__":
    main(sys.argv[1:] or list(t4.SITES))
