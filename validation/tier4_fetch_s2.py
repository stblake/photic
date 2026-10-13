"""Tier 4: download the Sentinel-2 L1C scenes over Exmouth Gulf (tile 49KHR).

Scenes come from the public Google Cloud bucket gcp-public-data-sentinel-2 (no
credentials). Only the files ACOLITE needs are fetched: the product and tile metadata
and the band JP2s. Downloads are idempotent (existing files of the right size are kept).

The scenes were chosen from Earth Search STAC (cloud cover 0%, observed tide at the
Exmouth 02 gauge, sun elevation 38-56 deg to limit glint, a spread of tides from
-0.7 to +0.7 m about MSL); see derived/tier4/candidates.json.
"""

import json
import os
import shutil
import urllib.parse
import urllib.request

BUCKET = "gcp-public-data-sentinel-2"
TILE_PREFIX = "tiles/49/K/HR/"

SCENES = [
    "S2B_MSIL1C_20190428T024329_N0207_R003_T49KHR_20190428T054446.SAFE",
    "S2A_MSIL1C_20190612T024321_N0207_R003_T49KHR_20190612T064555.SAFE",
    "S2B_MSIL1C_20210723T023249_N0500_R103_T49KHR_20230127T060055.SAFE",
    "S2B_MSIL1C_20210914T024319_N0500_R003_T49KHR_20230108T130759.SAFE",
    "S2B_MSIL1C_20220708T023259_N0400_R103_T49KHR_20220708T045218.SAFE",
    "S2B_MSIL1C_20240620T024319_N0510_R003_T49KHR_20240620T054509.SAFE",
]


def list_objects(prefix):
    items, token = [], None
    while True:
        q = {"prefix": prefix, "fields": "items(name,size),nextPageToken"}
        if token:
            q["pageToken"] = token
        url = "https://storage.googleapis.com/storage/v1/b/%s/o?%s" % (BUCKET, urllib.parse.urlencode(q))
        d = json.load(urllib.request.urlopen(url))
        items += d.get("items", [])
        token = d.get("nextPageToken")
        if not token:
            return items


def wanted(name):
    base = os.path.basename(name)
    if base.endswith("_$folder$"):
        return False
    if "/IMG_DATA/" in name and base.endswith(".jp2"):
        return True
    return base in ("MTD_MSIL1C.xml", "MTD_TL.xml", "manifest.safe")


def fetch(dest_root):
    os.makedirs(dest_root, exist_ok=True)
    for safe in SCENES:
        items = [i for i in list_objects(TILE_PREFIX + safe + "/") if wanted(i["name"])]
        if not items:
            raise SystemExit("no objects for %s in gs://%s/%s" % (safe, BUCKET, TILE_PREFIX))
        total = sum(int(i["size"]) for i in items)
        print("%s: %d files, %.0f MB" % (safe, len(items), total/1e6), flush=True)
        for i in items:
            rel = i["name"][len(TILE_PREFIX):]
            path = os.path.join(dest_root, rel)
            if os.path.exists(path) and os.path.getsize(path) == int(i["size"]):
                continue
            os.makedirs(os.path.dirname(path), exist_ok=True)
            url = "https://storage.googleapis.com/%s/%s" % (BUCKET, urllib.parse.quote(i["name"]))
            tmp = path + ".part"
            with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(tmp, path)


if __name__ == "__main__":
    fetch(os.path.join(os.environ.get("PHOTIC_VALIDATION_DATA",
                                      os.path.expanduser("~/photic_dev/validation_data")), "raw", "s2"))
