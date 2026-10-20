"""Tier 4, second region: download the ICESat-2 ATL24 bathymetry over the Arabian Gulf site.

ATL24 (version 002, NSIDC) gives classified ICESat-2 photons with refraction-corrected
orthometric heights. The granules crossing the site are found with CMR, downloaded with the
Earthdata login in ~/.netrc, cropped to the site's box and saved as one .npz per granule in
raw/atl24/ (seabed and sea-surface photons only, every field, with the beam name); the full
granule (about 100 MB) is then deleted. Downloads are idempotent (a granule with an .npz is skipped).

    [GULF_REGION=scene] tier4_fetch_atl24.py
"""

import json
import os
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))
import tier4_gulf  # noqa: E402 (the region: $GULF_REGION, palm or scene)
OUT = tier4_gulf.ATL24
BOX = (tier4_gulf.LIMIT[1], tier4_gulf.LIMIT[0], tier4_gulf.LIMIT[3], tier4_gulf.LIMIT[2])  # W, S, E, N
BEAMS = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
KEEP_CLASSES = (40, 41)  # ATL24 class_ph: 40 seabed, 41 sea surface (0 is unclassified)


def granules():
    q = dict(short_name="ATL24", version="002", bounding_box=",".join(map(str, BOX)), page_size=2000)
    with urllib.request.urlopen(CMR + "?" + urllib.parse.urlencode(q)) as r:
        entries = json.load(r)["feed"]["entry"]
    urls = []
    for e in entries:
        h5 = [l["href"] for l in e["links"] if l["href"].startswith("https://") and l["href"].endswith(".h5")]
        if h5:
            urls.append(h5[0])
    return sorted(set(urls))


def crop(h5, out):
    w, s, e, n = BOX
    keep = {}
    with h5py.File(h5) as f:
        for b in BEAMS:
            if b not in f or "lat_ph" not in f[b]:
                continue
            g = f[b]
            lat, lon = g["lat_ph"][:], g["lon_ph"][:]
            m = (lon >= w) & (lon <= e) & (lat >= s) & (lat <= n) & np.isin(g["class_ph"][:], KEEP_CLASSES)
            if not m.any():
                continue
            for k in g:
                if hasattr(g[k], "shape") and g[k].shape == lat.shape:
                    keep.setdefault(k, []).append(g[k][:][m])
            keep.setdefault("beam", []).append(np.full(int(m.sum()), b))
    np.savez_compressed(out, **{k: np.concatenate(v) for k, v in keep.items()})
    return sum(len(v) for v in keep.get("beam", []))


def shrink(path):
    """Drop the unclassified photons from an .npz written before KEEP_CLASSES."""
    d = dict(np.load(path))
    m = np.isin(d["class_ph"], KEEP_CLASSES) if "class_ph" in d else None
    if m is None or m.all():
        return
    np.savez_compressed(path[:-4] + ".tmp.npz", **{k: v[m] for k, v in d.items()})
    os.replace(path[:-4] + ".tmp.npz", path)


def fetch_one(url):
    name = os.path.basename(url)[:-3]
    out = os.path.join(OUT, name + ".npz")
    if os.path.exists(out):
        shrink(out)
        return name, "kept"
    with tempfile.TemporaryDirectory(dir=OUT) as tmp:
        h5, cookies = os.path.join(tmp, "g.h5"), os.path.join(tmp, "cookies")
        # curl -n reads the Earthdata login from ~/.netrc; the cookies carry it through NSIDC's redirect.
        subprocess.run(["curl", "-s", "-f", "-n", "-L", "-c", cookies, "-b", cookies, "-o", h5, url], check=True)
        n = crop(h5, out + ".tmp.npz")
        os.replace(out + ".tmp.npz", out)
    return name, "%d photons" % n


def main():
    os.makedirs(OUT, exist_ok=True)
    urls = granules()
    print("%d granules" % len(urls), flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for name, status in pool.map(fetch_one, urls):
            print(name, status, flush=True)


if __name__ == "__main__":
    if sys.argv[1:]:
        raise SystemExit(__doc__)
    main()
