"""Phase G6, first step: fit each scene's aerosol optical thickness to the water.

ACOLITE's dark spectrum fit leaves a positive residual in the red and near infrared over clear,
deep water (EXPERIMENTS.md, Phase G: Rrs(865) 4-7e-4, Rrs(665) about 1e-3 too high). Here each
scene's AOT and aerosol model are chosen so that the deep water of a dark target is consistent:
  - at 665 nm, rrs equals what the water found by QAA at 560 nm predicts (the red offset of
    DROP's PROFILE_QAA 3 is zero), and
  - at 865 nm, where water is nearly black, rrs is zero.
The cost is delta665^2 + rrs865^2 over the target's deep cells (DROP's rule: lowest 20% of the
smoothed red, then lowest 15% of those by smoothed green, chosen once at the first AOT), for
each model on a grid of AOT, refined by a parabola through the three lowest points. The target
is clear ocean west of the Ningaloo reef, else north of the cape, else the gulf entrance.

    aot_fit.py fit s2|l8            # writes derived/tier4/aot_fit_<sensor>.json
    aot_fit.py acolite s2|l8        # full-crop ACOLITE with the fitted AOT: acolite_aot[_l8]
    aot_fit.py sites s2|l8          # photic inputs: sites_aot, sites_l8_aot
"""

import contextlib
import glob
import io
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.ndimage import uniform_filter

import tier4_prepare as tp

DATA = tp.DATA
WORK = os.path.join(DATA, "derived", "tier4", "aot_fit_work")
LUTS = ["ACOLITE-LUT-202110-MOD1", "ACOLITE-LUT-202110-MOD2"]
AOTS = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
TARGETS = [  # S, W, N, E; clear deep ocean first
    ("west", [-21.95, 113.95, -21.85, 113.962]),
    ("north", [-21.70, 113.95, -21.62, 114.10]),
    ("gulf_entrance", [-21.73, 114.34, -21.63, 114.44]),
]

# QAA at 560 nm (Lee et al. 2014, v6) with photic's g0, g1, and water from model/data (photic's
# default tables) at each sensor's band centres; bands are keyed by nominal 443, 490, 560, 665.
CENTRES = {"s2": {443: 443.0, 490: 492.0, 560: 560.0, 665: 665.0},
           "l8": {443: 443.0, 490: 483.0, 560: 561.0, 665: 655.0}}
G0, G1 = 0.084, 0.17
MODEL_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model", "data")


def table(name, column):
    rows = [l for l in open(os.path.join(MODEL_DATA, name)) if not l.startswith("#") and l.strip()]
    head = rows[0].strip().split(",")
    a = np.array([[float(x) for x in r.split(",")] for r in rows[1:] if r[0].isdigit()])
    return a[:, 0], a[:, head.index(column)]


def water(sensor):
    wa, aw = table("aw_pope_fry_1997.csv", "aw")
    wb, bbw = table("bbw_morel_1974.csv", "bbw")
    c = CENTRES[sensor]
    return ({k: float(np.interp(v, wa, aw)) for k, v in c.items()},
            {k: float(np.interp(v, wb, bbw)) for k, v in c.items()})


def rrs_of(R):
    return R/(0.52 + 1.7*R)


def red_offset(r, sensor):
    """The flat offset of rrs for which rrs(red) - offset is what the water found by QAA at 560 nm
    from rrs - offset predicts, as profile.c's red_offset (negative allowed here)."""
    AW, BBW = water(sensor)
    c = CENTRES[sensor]
    d = 0.0
    for _ in range(40):
        q = {k: max(v - d, 1e-6) for k, v in r.items()}
        chi = np.log10((q[443] + q[490])/(q[560] + 5*q[665]**2/q[490]))
        a560 = AW[560] + 10**(-1.146 - 1.366*chi - 0.469*chi*chi)
        u = (-G0 + np.sqrt(G0*G0 + 4*G1*q[560]))/(2*G1)
        bbp = max(u*a560/(1 - u) - BBW[560], 1e-6)
        Y = np.clip(2*(1 - 1.2*np.exp(-0.9*q[443]/q[560])), 0, 2.5)
        bb = BBW[665] + bbp*(c[560]/c[665])**Y
        u = bb/(AW[665] + bb)
        d = r[665] - (G0 + G1*u)*u
    return d


def bands(l2w, sensor):
    import netCDF4
    with netCDF4.Dataset(l2w) as d:
        names = [v for v in d.variables if v.startswith("Rrs_")]
        out = {}
        for k, wl in list(CENTRES[sensor].items()) + [(865, 865.0)]:
            v = min(names, key=lambda x: abs(float(x[4:]) - wl))
            out[k] = np.ma.filled(d[v][:].astype(float), np.nan)
    return out


def deep_cells(B):
    ok = np.all([np.isfinite(B[w]) & (B[w] > 0) for w in (443, 490, 560, 665)], axis=0)
    def sm(x):
        y = np.where(ok, x, 0.0)
        return uniform_filter(y, 7)/np.maximum(uniform_filter(ok.astype(float), 7), 1e-9)
    red, grn = sm(B[665]), sm(B[560])
    c1 = ok & (red <= np.percentile(red[ok], 20))
    return c1 & (grn <= np.percentile(grn[c1], 15))


def run_one(args):
    scene, limit, lut, aot, out = args
    if not glob.glob(os.path.join(out, "*L2W.nc")):
        with contextlib.redirect_stdout(io.StringIO()):
            tp.run_acolite(scene, out, limit, dsf_fixed_aot=aot, dsf_fixed_lut=lut, l2w_mask=False)
    f = glob.glob(os.path.join(out, "*L2W.nc"))
    return f[0] if f else None


def covered(scene_l2w, limit, n_min=5000):
    import netCDF4
    s, w, n, e = limit
    with netCDF4.Dataset(scene_l2w) as d:
        lon, lat = d["lon"][:], d["lat"][:]
        v = [x for x in d.variables if x.startswith("Rrs_6")][0]
        r = np.ma.filled(d[v][:].astype(float), np.nan)
    m = (lon >= w) & (lon <= e) & (lat >= s) & (lat <= n)
    return int(np.isfinite(r[m]).sum()) >= n_min


def scenes_for(sensor):
    if sensor == "s2":
        import tier4_fetch_s2
        raw = [os.path.join(DATA, "raw", "s2", x) for x in tier4_fetch_s2.SCENES]
        base = "acolite"
    else:
        raw = sorted(d for d in glob.glob(os.path.join(DATA, "raw", "l8", "*")) if os.path.isdir(d))
        base = "acolite_l8"
    out = []
    for r in raw:
        name = os.path.basename(r).replace(".SAFE", "")
        l2w = glob.glob(os.path.join(DATA, "derived", "tier4", base, name, "*L2W.nc"))
        out.append((r, name, l2w[0] if l2w else None))
    return out


def fit(sensor):
    jobs, plan = [], {}
    for raw, name, l2w in scenes_for(sensor):
        target = next((t for t in TARGETS if l2w is None or covered(l2w, t[1])), TARGETS[-1])
        plan[name] = (raw, target)
        for lut in LUTS:
            for aot in AOTS:
                out = os.path.join(WORK, sensor, "%s_%s_%s_%.3f" % (name, target[0], lut[-4:], aot))
                jobs.append((raw, target[1], lut, aot, out))
    with ProcessPoolExecutor(max_workers=6) as pool:
        files = dict(zip([(j[0], j[2], j[3]) for j in jobs], pool.map(run_one, jobs)))
    result = {}
    for name, (raw, target) in plan.items():
        rows, deep = {}, None
        for lut in LUTS:
            J = []
            for aot in AOTS:
                f = files[(raw, lut, aot)]
                B = bands(f, sensor)
                if deep is None:
                    deep = deep_cells(B)
                med = {w: float(np.nanmedian(B[w][deep])) for w in B}
                d665 = red_offset({w: rrs_of(med[w]) for w in (443, 490, 560, 665)}, sensor)
                r865 = med[865]/0.52
                J.append((d665**2 + r865**2, d665, r865))
            k = int(np.argmin([j[0] for j in J]))
            k = min(max(k, 1), len(AOTS) - 2)
            x, y = np.array(AOTS[k - 1:k + 2]), np.array([j[0] for j in J[k - 1:k + 2]])
            a, b, _ = np.polyfit(x, y, 2)
            best = float(np.clip(-b/(2*a), x[0], x[-1])) if a > 0 else float(x[np.argmin(y)])
            rows[lut] = dict(aot=best, cost=float(np.interp(best, AOTS, [j[0] for j in J])),
                             scan=[dict(aot=t, delta665=j[1], rrs865=j[2]) for t, j in zip(AOTS, J)])
        lut = min(rows, key=lambda m: rows[m]["cost"])
        result[name] = dict(target=target[0], lut=lut, aot=rows[lut]["aot"], n_deep=int(deep.sum()), models=rows)
        print("%-45s %-14s %s aot %.3f  cost %.2e (other %.2e)" % (name, target[0], lut[-4:], rows[lut]["aot"],
              rows[lut]["cost"], min(r["cost"] for m, r in rows.items() if m != lut)), flush=True)
    json.dump(result, open(os.path.join(DATA, "derived", "tier4", "aot_fit_%s.json" % sensor), "w"), indent=1)


def acolite(sensor):
    fits = json.load(open(os.path.join(DATA, "derived", "tier4", "aot_fit_%s.json" % sensor)))
    derived = "acolite_aot" if sensor == "s2" else "acolite_l8_aot"
    jobs = []
    for raw, name, _ in scenes_for(sensor):
        out = os.path.join(DATA, "derived", "tier4", derived, name)
        jobs.append((raw, tp.LIMIT, fits[name]["lut"], fits[name]["aot"], out))
    with ProcessPoolExecutor(max_workers=3) as pool:
        for f in pool.map(run_acolite_full, jobs):
            print(f, flush=True)


def run_acolite_full(args):
    scene, limit, lut, aot, out = args
    if not glob.glob(os.path.join(out, "*L2W.nc")):
        with contextlib.redirect_stdout(io.StringIO()):
            # With the AOT fitted so deep water is ~0 at 865 nm, about half the pixels are slightly
            # negative there (noise about zero); ACOLITE's negative-reflectance mask (400-900 nm)
            # would then drop them, although photic uses only 443-665 nm (prepare_sites rejects
            # non-positive values in those bands itself).
            tp.run_acolite(scene, out, limit, dsf_fixed_aot=aot, dsf_fixed_lut=lut, l2w_mask_negative_rhow=False)
    return out


def sites(sensor):
    if sensor == "s2":
        import tier4_fetch_s2
        tp.prepare_sites(sites_name="sites_aot", scenes_only=tier4_fetch_s2.SCENES, acolite="acolite_aot")
    else:
        tp.prepare_sites(tp.L8_BANDS, "sites_l8_aot", acolite="acolite_l8_aot", min_cell=tp.L8_CELL)


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("fit", "acolite", "sites") or sys.argv[2] not in ("s2", "l8"):
        raise SystemExit(__doc__)
    {"fit": fit, "acolite": acolite, "sites": sites}[sys.argv[1]](sys.argv[2])
