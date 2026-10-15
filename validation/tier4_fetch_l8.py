"""Tier 4: download the Landsat 8 OLI L1TP scenes over Exmouth Gulf (WRS-2 path 115, row 75).

Scenes come from the public Google Cloud bucket gcp-public-data-landsat (Collection 1, no
credentials; it holds Landsat 8 to the end of 2021). Only the files ACOLITE needs are fetched:
the metadata, the angle file, the quality band and bands 1-7 and 9. Downloads are idempotent
(existing files of the right size are kept).

The scenes were chosen like the Sentinel-2 ones (tier4_fetch_s2.py): scene cloud cover <= 0.1%,
observed tide at the Exmouth 02 gauge, sun elevation 38-56 deg to limit glint and a spread of
tides from -0.7 to +0.9 m about MSL, in the same years as the Sentinel-2 scenes where possible.
"""

import json
import os
import shutil
import urllib.parse
import urllib.request

BUCKET = "gcp-public-data-landsat"
PREFIX = "LC08/01/115/075/"

SCENES = [  # sun elevation, tide about MSL at the overpass
    "LC08_L1TP_115075_20190809_20190820_01_T1",  # 41 deg, -0.67 m
    "LC08_L1TP_115075_20210830_20210909_01_T1",  # 47 deg, -0.62 m
    "LC08_L1TP_115075_20180822_20180829_01_T1",  # 45 deg, -0.08 m
    "LC08_L1TP_115075_20190910_20190917_01_T1",  # 51 deg, +0.01 m
    "LC08_L1TP_115075_20190403_20190421_01_T1",  # 49 deg, +0.59 m
    "LC08_L1TP_115075_20190505_20190520_01_T1",  # 42 deg, +0.87 m
]

# The Arabian Gulf between Abu Dhabi and Dubai (path 160, row 43; tier4_gulf.py): scene cloud cover
# <= 0.3%, sun elevation 38-56 deg, a tide changing by at most 0.2 m/h, a spread of tides (GOT4.10 at
# 54.6 E, 25.0 N), 2018-2021.

GULF_PREFIX = "LC08/01/160/043/"

GULF_SCENES = [
    "LC08_L1TP_160043_20180310_20180320_01_T1",  # 52 deg, -0.16 m
    "LC08_L1TP_160043_20211012_20211019_01_T1",  # 53 deg, -0.07 m
    "LC08_L1TP_160043_20191023_20191030_01_T1",  # 50 deg, +0.17 m
    "LC08_L1TP_160043_20201025_20201106_01_T1",  # 49 deg, +0.24 m
    "LC08_L1TP_160043_20191108_20191115_01_T1",  # 45 deg, +0.39 m
    "LC08_L1TP_160043_20181105_20181115_01_T1",  # 46 deg, +0.60 m
]

FILES = ["MTL.txt", "ANG.txt", "BQA.TIF"] + ["B%d.TIF" % b for b in (1, 2, 3, 4, 5, 6, 7, 9)]


def fetch(dest_root, scenes=SCENES, prefix=PREFIX):
    for scene in scenes:
        os.makedirs(os.path.join(dest_root, scene), exist_ok=True)
        for suffix in FILES:
            name = "%s_%s" % (scene, suffix)
            obj = prefix + scene + "/" + name
            meta = json.load(urllib.request.urlopen("https://storage.googleapis.com/storage/v1/b/%s/o/%s"
                                                    % (BUCKET, urllib.parse.quote(obj, safe=""))))
            path = os.path.join(dest_root, scene, name)
            if os.path.exists(path) and os.path.getsize(path) == int(meta["size"]):
                continue
            print(name, "%.0f MB" % (int(meta["size"])/1e6), flush=True)
            tmp = path + ".part"
            with urllib.request.urlopen("https://storage.googleapis.com/%s/%s" % (BUCKET, obj)) as r, \
                    open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(tmp, path)


if __name__ == "__main__":
    import sys
    raw = os.path.join(os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data")), "raw")
    if sys.argv[1:] == ["gulf"]:
        fetch(os.path.join(raw, "l8_gulf"), GULF_SCENES, GULF_PREFIX)
    else:
        fetch(os.path.join(raw, "l8"))
