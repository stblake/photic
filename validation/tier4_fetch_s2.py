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

# Ten more cloud-free scenes (2026-10-14) to test which conditions spoil a scene: they span
# sun elevation (glint) and the rate of change of the tide at the Exmouth gauge over the hour
# around the overpass (tidal currents), as well as the tide. Chosen at random within each
# group (seed 7) from derived/tier4/candidates.json, one per year where possible. Tier 4's
# default sites use the first six only; the sites16 set uses all.

EXTRA_SCENES = [
    "S2B_MSIL1C_20220412T024319_N0400_R003_T49KHR_20220412T055324.SAFE",  # low glint, slack; sun el. 51, tide +0.12 m, -3 cm/h
    "S2A_MSIL1C_20210906T023251_N0500_R103_T49KHR_20230112T232004.SAFE",  # low glint, slack; sun el. 51, tide +0.84 m, +2 cm/h
    "S2A_MSIL1C_20190503T024331_N0207_R003_T49KHR_20190503T054637.SAFE",  # low glint, slack; sun el. 46, tide +0.80 m, -5 cm/h
    "S2A_MSIL1C_20190523T024331_N0207_R003_T49KHR_20190523T064343.SAFE",  # low glint, fast tide; sun el. 41, tide +0.33 m, +40 cm/h
    "S2A_MSIL1C_20210402T024321_N0300_R003_T49KHR_20210402T051142.SAFE",  # low glint, fast tide; sun el. 53, tide +0.40 m, +51 cm/h
    "S2B_MSIL1C_20220519T023249_N0400_R103_T49KHR_20220519T055351.SAFE",  # low glint, fast tide; sun el. 41, tide +0.87 m, +40 cm/h
    "S2B_MSIL1C_20211031T023249_N0500_R103_T49KHR_20230105T051602.SAFE",  # high sun, slack; sun el. 67, tide -0.16 m, -5 cm/h
    "S2A_MSIL1C_20221004T024331_N0400_R003_T49KHR_20221004T055416.SAFE",  # high sun, slack; sun el. 62, tide -0.50 m, +2 cm/h
    "S2A_MSIL1C_20230112T024321_N0509_R003_T49KHR_20230112T054603.SAFE",  # high sun, fast tide; sun el. 65, tide -0.47 m, +36 cm/h
    "S2A_MSIL1C_20191116T023251_N0500_R103_T49KHR_20230526T101946.SAFE",  # high sun, fast tide; sun el. 68, tide -0.35 m, +47 cm/h
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


# The Arabian Gulf between Abu Dhabi and Dubai (tile 40RBN; tier4_gulf.py). Chosen from Earth Search
# STAC (2026-10-15): cloud cover 0%, relative orbit 63 (the whole tile), 2021-2025 (the years of the
# ICESat-2 reference), sun elevation 38-56 deg, a tide changing by at most 0.2 m/h and a spread of
# tides about MSL (GOT4.10 at 54.6 E, 25.0 N; there is no open gauge record).

GULF_TILE_PREFIX = "tiles/40/R/BN/"

GULF_SCENES = [
    "S2C_MSIL1C_20250204T070131_N0511_R063_T40RBN_20250204T090729.SAFE",  # sun el. 44, tide -0.38 m
    "S2B_MSIL1C_20220126T070139_N0400_R063_T40RBN_20220126T074954.SAFE",  # sun el. 42, tide -0.15 m
    "S2A_MSIL1C_20221018T065901_N0400_R063_T40RBN_20221018T074328.SAFE",  # sun el. 53, tide +0.12 m
    "S2B_MSIL1C_20221023T065929_N0400_R063_T40RBN_20221023T090520.SAFE",  # sun el. 52, tide +0.47 m
    "S2A_MSIL1C_20230205T070101_N0509_R063_T40RBN_20230205T074629.SAFE",  # sun el. 44, tide +0.66 m
    "S2A_MSIL1C_20220131T070121_N0400_R063_T40RBN_20220131T080317.SAFE",  # sun el. 43, tide +0.92 m
]


def fetch(dest_root, scenes=None, tile_prefix=TILE_PREFIX):
    os.makedirs(dest_root, exist_ok=True)
    for safe in scenes or SCENES + EXTRA_SCENES:
        items = [i for i in list_objects(tile_prefix + safe + "/") if wanted(i["name"])]
        if not items:
            raise SystemExit("no objects for %s in gs://%s/%s" % (safe, BUCKET, tile_prefix))
        total = sum(int(i["size"]) for i in items)
        print("%s: %d files, %.0f MB" % (safe, len(items), total/1e6), flush=True)
        for i in items:
            rel = i["name"][len(tile_prefix):]
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
    import sys
    raw = os.path.join(os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data")), "raw")
    if sys.argv[1:] == ["gulf"]:
        fetch(os.path.join(raw, "s2_gulf"), GULF_SCENES, GULF_TILE_PREFIX)
    else:
        fetch(os.path.join(raw, "s2"))
