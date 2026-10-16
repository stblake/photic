"""Write photic's spectral response function files (model/data/srf/*.csv) from ACOLITE's RSR tables.

ACOLITE (https://github.com/acolite/acolite, data/RSR) tabulates the agencies' published
relative spectral responses:
  - Sentinel-2A/2B/2C MSI: ESA, COPE-GSEG-EOPG-TN-15-0007, Sentinel-2 Spectral Response
    Functions, version 4.0 (2024);
  - Landsat 8 OLI: NASA/USGS band-averaged RSR (Barsi et al. 2014, Remote Sens. 6, 10232);
  - Landsat 9 OLI-2: NASA, L9_OLI2_Ball_FPM_RSR.v1.0 (Barsi et al. 2019, Proc. SPIE 11127).

Only the bands photic's tables cover (400-800 nm) are written: MSI B1-B5, OLI B1-B4.

  python srf_files.py <directory with S2A_MSI_V4.txt S2B_MSI_V4.txt S2C_MSI_V4.txt L8_OLI.txt L9_OLI.txt>
"""

import os
import re
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "model", "data", "srf")

SENSORS = [  # ACOLITE file, output, bands, description
    ("S2A_MSI_V4.txt", "S2A_MSI.csv", ["1", "2", "3", "4", "5"],
     "Sentinel-2A MSI; ESA COPE-GSEG-EOPG-TN-15-0007, Sentinel-2 Spectral Response Functions 4.0 (2024)"),
    ("S2B_MSI_V4.txt", "S2B_MSI.csv", ["1", "2", "3", "4", "5"],
     "Sentinel-2B MSI; ESA COPE-GSEG-EOPG-TN-15-0007, Sentinel-2 Spectral Response Functions 4.0 (2024)"),
    ("S2C_MSI_V4.txt", "S2C_MSI.csv", ["1", "2", "3", "4", "5"],
     "Sentinel-2C MSI; ESA COPE-GSEG-EOPG-TN-15-0007, Sentinel-2 Spectral Response Functions 4.0 (2024)"),
    ("L8_OLI.txt", "L8_OLI.csv", ["1", "2", "3", "4"],
     "Landsat 8 OLI; NASA/USGS band-averaged relative spectral response (Barsi et al. 2014)"),
    ("L9_OLI.txt", "L9_OLI.csv", ["1", "2", "3", "4"],
     "Landsat 9 OLI-2; NASA L9_OLI2_Ball_FPM_RSR.v1.0 (Barsi et al. 2019)"),
]


def read_acolite(path):
    bands, cur = {}, None
    for line in open(path, encoding="latin-1"):
        s = line.strip()
        if s.startswith(("#", ";")):
            m = re.search(r"band\s+(\w+)", s, re.I)
            if m:
                cur = m.group(1)
                bands[cur] = []
            continue
        f = s.split()
        if cur and len(f) >= 2:
            bands[cur].append((float(f[0]), float(f[1])))
    out = {}
    for k, v in bands.items():
        a = np.array(v)
        wl = a[:, 0]*(1000.0 if a[:, 0].max() < 10 else 1.0)
        out[k] = (wl, np.clip(a[:, 1], 0.0, None))
    return out


def main(src):
    os.makedirs(OUT, exist_ok=True)
    for fname, out, keep, desc in SENSORS:
        rsr = read_acolite(os.path.join(src, fname))
        lo = int(np.floor(min(rsr[b][0][rsr[b][1] > 0].min() for b in keep)))
        hi = int(np.ceil(max(rsr[b][0][rsr[b][1] > 0].max() for b in keep)))
        lam = np.arange(lo, hi + 1, 1.0)
        cols = [np.interp(lam, *rsr[b], left=0.0, right=0.0) for b in keep]
        with open(os.path.join(OUT, out), "w") as f:
            f.write("# Relative spectral response, %s.\n" % desc)
            f.write("# Tabulated by ACOLITE (data/RSR/%s), resampled to 1 nm by validation/srf_files.py.\n" % fname)
            f.write("# Band centroids (nm): %s\n" % ", ".join(
                "B%s %.1f" % (b, np.sum(lam*c)/np.sum(c)) for b, c in zip(keep, cols)))
            f.write("wavelength," + ",".join("B" + b for b in keep) + "\n")
            for i, w in enumerate(lam):
                f.write("%g," % w + ",".join("%.6g" % c[i] for c in cols) + "\n")
        print(out, lo, hi)


if __name__ == "__main__":
    main(sys.argv[1])
