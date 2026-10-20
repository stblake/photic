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

# The Arabian Gulf between Abu Dhabi and Dubai (path 160, row 43), chosen by tier4_gulf.py select:
# haze first (the six clearest by SWIR surface reflectance over open water, of 51 scenes with < 5%
# cloud, 2018-10 to 2021), among those with the sun at 38-56 deg and a tide changing by at most
# 0.2 m/h (GOT4.10).

GULF_PREFIX = "LC08/01/160/043/"

GULF_SCENES = [
    "LC08_L1TP_160043_20190313_20190325_01_T1",  # SWIR 0.0042, 53 deg, -0.49 m
    "LC08_L1TP_160043_20190209_20190221_01_T1",  # SWIR 0.0026, 43 deg, -0.07 m
    "LC08_L1TP_160043_20181004_20181010_01_T1",  # SWIR 0.0048, 55 deg, +0.10 m
    "LC08_L1TP_160043_20201110_20210316_01_T1",  # SWIR 0.0037, 44 deg, +0.19 m
    "LC08_L1TP_160043_20191108_20191115_01_T1",  # SWIR 0.0048, 45 deg, +0.38 m
    "LC08_L1TP_160043_20181207_20181211_01_T1",  # SWIR 0.0030, 39 deg, +0.64 m
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
