"""Tier 2: forward physics against exact radiative transfer (OSOAA).

For each case (wavelength, water type, depth, bottom, sun angle), OSOAA solves the full
radiative transfer problem for a water column whose absorption and backscattering are
exactly photic's, over a Lambertian bottom with photic's albedo. photic's semi-analytical
rrs/Rrs (photic_ref.forward, verified bit-exact against the C code in tier 1) is then
compared with OSOAA's.

OSOAA setup (scalar, no atmosphere, so the comparison is with Lee's scalar model):
  * tools/OSOAA_scalar: built with CTE_POLAR_SWITCHED_OFF 1 and
    EXPERT_MODE_FORCED_FATM_NULL .TRUE.; AOT 0 and near-zero pressure.
  * HYD.Model 3: user profile of hydrosol absorption/scattering (water removed) plus a
    Fournier-Forand phase function with a fixed backscattering ratio.
  * rrs(0-) = REFL/pi just below the surface at nadir (RESLUM_vsZ),
    Rrs(0+) = REFL/pi just above the surface at nadir (RESLUM_vsVZA, level 3).
"""

import csv
import hashlib
import json
import itertools
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import photic_ref
import tier1

N_WATER = 1.34
WIND = 1.0  # m/s
PRESSURE = 1013.0  # hPa. The black-sky build nulls the atmospheric source, so this only
# avoids the NaNs OSOAA produces for a near-zero Rayleigh optical thickness.
FF_N, FF_MU = 1.10, 3.5835  # Fournier-Forand parameters (backscattering ratio ~0.0183).

WATER_TYPES = {  # P, G, X, Y
    "clear": (0.02, 0.02, 0.004, 1.0),
    "moderate": (0.08, 0.15, 0.02, 1.0),
    "turbid": (0.30, 0.50, 0.08, 1.0),
}
DEPTHS = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, math.inf]
BOTTOMS = {"sand": 0, "seagrass": 1, "coral": 2}
B_ALBEDO = 0.6
THETA_AIR = [20.0, 40.0, 60.0]
BANDS = [443.0, 490.0, 560.0, 665.0, 705.0]

QUICK = dict(water=["clear", "moderate", "turbid"], depths=[1.0, 3.0, 8.0, math.inf],
             bottoms=["sand"], theta=[20.0, 40.0], bands=[490.0, 560.0])
FULL = dict(water=list(WATER_TYPES), depths=DEPTHS, bottoms=list(BOTTOMS),
            theta=THETA_AIR, bands=BANDS)


def fournier_forand(theta_deg):
    """Fournier-Forand phase function (Mobley et al. 2002), normalised to 1 over 4 pi."""
    th = np.radians(np.maximum(theta_deg, 1e-3))
    nu = (3.0 - FF_MU)/2.0
    d = 4.0/(3.0*(FF_N - 1.0)**2)*np.sin(th/2.0)**2
    d180 = 4.0/(3.0*(FF_N - 1.0)**2)
    p = (1.0/(4.0*np.pi*(1.0 - d)**2*d**nu))*(
        nu*(1.0 - d) - (1.0 - d**nu) + (d*(1.0 - d**nu) - nu*(1.0 - d))/np.sin(th/2.0)**2)
    p += (1.0 - d180**nu)/(16.0*np.pi*(d180 - 1.0)*d180**nu)*(3.0*np.cos(th)**2 - 1.0)
    return p


def phase_table():
    """Angles and F11 (normalised so that integral F11 sin(theta) dtheta = 2) with its
    backscattering ratio, on <= 200 angles (CTE_MAXNB_ANG_EXT)."""
    ang = np.unique(np.concatenate([[0.0], np.geomspace(0.01, 10.0, 80),
                                    np.linspace(10.5, 180.0, 119)]))
    f11 = 4.0*np.pi*fournier_forand(ang)
    th = np.radians(ang)
    norm = np.trapezoid(f11*np.sin(th), th)
    f11 *= 2.0/norm
    back = th >= np.pi/2
    btilde = np.trapezoid(f11[back]*np.sin(th[back]), th[back])/2.0
    return ang, f11, btilde


def water_coeffs(osoaa_root):
    """OSOAA's pure seawater absorption (table) and Morel (1974) scattering."""
    data = np.loadtxt(os.path.join(osoaa_root, "fic", "OSOAA_SEA_MOL_COEFFS_JUNE_2013.txt"),
                      comments=("/", "!"))
    def aw(lam):
        return float(np.interp(lam, data[:, 0], data[:, 1]))
    def bw(lam):
        return 0.00288*(lam/500.0)**(-4.32)
    return aw, bw


def photic_tables(driver, bands):
    """photic's interpolated a0, a1, a_w, b_bw and bottom tables at the bands."""
    cases = {k: [] for k in ["id", "lam", "P", "G", "X", "B", "bottom", "H", "tide",
                             "theta_sun", "theta_view", "R440", "R490", "R750"]}
    n = 0
    for lam in bands:
        for b in BOTTOMS.values():
            for k, v in zip(cases, [n, lam, 0.05, 0.1, 0.01, 1.0, b, 5.0, 0.0, 30.0, 0.0,
                                    0.01, 0.01, 0.0]):
                cases[k].append(v)
            n += 1
    cases = {k: np.array(v, dtype=float) for k, v in cases.items()}
    dump = tier1.run_driver(driver, cases)
    t = {}
    for i in range(n):
        lam, b = cases["lam"][i], int(cases["bottom"][i])
        t[(lam, b)] = dict(a0=dump["a_0"][i], a1=dump["a_1"][i], aw=dump["a_w"][i],
                           bbw=dump["b_bw"][i], rho=dump["rho_table"][i])
    return t


def make_cases(sel):
    cases = []
    for water, H, bottom, theta_air, lam in itertools.product(
            sel["water"], sel["depths"], sel["bottoms"], sel["theta"], sel["bands"]):
        if math.isinf(H) and bottom != sel["bottoms"][0]:
            continue  # the bottom does not matter in optically deep water
        cases.append(dict(water=water, H=H, bottom=bottom, theta_air=theta_air, lam=lam))
    return cases


def case_key(c):
    s = "%s|%r|%s|%r|%r" % (c["water"], c["H"], c["bottom"], c["theta_air"], c["lam"])
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def photic_side(c, tables):
    P, G, X, Y = WATER_TYPES[c["water"]]
    t = tables[(c["lam"], BOTTOMS[c["bottom"]])]
    theta_w = math.degrees(math.asin(math.sin(math.radians(c["theta_air"]))/N_WATER))
    H = 1e4 if math.isinf(c["H"]) else c["H"]
    r = photic_ref.forward(c["lam"], P, G, X, B_ALBEDO, t["rho"], H, 0.0, theta_w, 0.0,
                           0.01, 0.01, 0.0, t["a0"], t["a1"], t["aw"], t["bbw"], Y=Y)
    out = {k: float(np.asarray(v)) for k, v in r.items()}
    out["rho_table"] = float(t["rho"])
    return out, theta_w


def run_osoaa(c, ph, exe, osoaa_root, phase, workdir, dbdir, exe_black_ocean=None, albedo=None, aot=0.0):
    """OSOAA for one case: rrs(0-) and water-leaving Rrs(0+) at nadir, and the IOPs used.

    The radiance just above the surface includes sun glint (and any sky reflection), which
    is not water-leaving radiance. If exe_black_ocean (OSOAA built with
    EXPERT_MODE_FORCED_FSEA_NULL) is given, a second run with no marine scattering and a
    black bottom gives that surface-reflected part, which is subtracted. albedo, if given,
    overrides the bottom albedo (photic's B * rho_table). aot is the aerosol optical thickness
    at 550 nm (maritime, RH 80%); it only matters with a build that has atmospheric scattering
    (OSOAA_sky), where the sky then lights the water.
    """
    ang, f11, btilde = phase
    aw_o, bw_o = water_coeffs(osoaa_root)
    lam = c["lam"]
    a_hyd = ph["a"] - aw_o(lam)
    bb_hyd = ph["bb"] - 0.5*bw_o(lam)
    # Where photic's tabulated a_w is a little below OSOAA's (clear water, red bands) the
    # hydrosol absorption would be negative; clamp it and record the mismatch.
    a_mismatch = max(0.0, -a_hyd)/ph["a"]
    a_hyd = max(a_hyd, 1e-5)  # OSOAA returns NaN for exactly zero hydrosol absorption
    if bb_hyd < 0:
        raise ValueError("case %s: negative hydrosol backscattering (%g)" % (case_key(c), bb_hyd))
    b_hyd = bb_hyd/btilde
    c_tot = ph["a"] + b_hyd + bw_o(lam)
    if math.isinf(c["H"]):
        depth, alb = 25.0/c_tot, 0.0  # optical thickness 25 over a black bottom
    else:
        depth, alb = c["H"], B_ALBEDO*ph["rho_table"]
    if albedo is not None:
        alb = albedo

    os.makedirs(workdir, exist_ok=True)
    ext = os.path.join(workdir, "hyd_phase.txt")
    with open(ext, "w") as f:
        f.write("EXTINCTION_COEF : 1.0\nSCATTERING_COEF : 1.0\nNB_LINES : %d\n" % len(ang))
        f.write("Angle F11 -F12/F11 F22/F11 F33/F11\n")
        for a_, p_ in zip(ang[::-1], f11[::-1]):
            f.write("%.4f %.6E 0.0 1.0 1.0\n" % (a_, p_))
    prof = os.path.join(workdir, "hyd_profile.txt")
    with open(prof, "w") as f:
        f.write("photic tier 2\ndepth a_hyd b_hyd\n(m) (1/m) (1/m)\n-\n-\n")
        f.write("0.0 %.8E %.8E\n%.4f %.8E %.8E" % (a_hyd, b_hyd, depth + 1.0, a_hyd, b_hyd))

    rrs_0m, refl_0p = osoaa_once(c, exe, osoaa_root, workdir, dbdir, ext, prof, depth, alb, aot)
    surface = 0.0
    if exe_black_ocean is not None:
        _, surface = osoaa_once(c, exe_black_ocean, osoaa_root, os.path.join(workdir, "black_ocean"),
                                dbdir, ext, prof, depth, 0.0, aot)
    return dict(rrs_osoaa=rrs_0m, Rrs_osoaa=(refl_0p - surface)/math.pi,
                Rrs_surface_osoaa=surface/math.pi, a_hyd=a_hyd, b_hyd=b_hyd, depth_osoaa=depth,
                bottom_albedo=alb, a_rel_mismatch=a_mismatch)


def osoaa_once(c, exe, osoaa_root, workdir, dbdir, ext, prof, depth, alb, aot=0.0):
    """Run OSOAA; return rrs(0-) and REFL(0+) = pi L(0+)/Ed(0+) at nadir."""
    os.makedirs(workdir, exist_ok=True)
    res = os.path.join(workdir, "res")
    shutil.rmtree(res, ignore_errors=True)
    cmd = [exe, "-OSOAA.ResRoot", res, "-OSOAA.Log", "Main.Log", "-OSOAA.Wa", "%.4f" % (c["lam"]/1000.0),
           "-ANG.Thetas", "%.4f" % c["theta_air"],
           "-AP.Pressure", "%.3f" % PRESSURE, "-AP.HR", "8.0", "-AP.HA", "2.0",
           "-AER.Waref", "0.550", "-AER.AOTref", "%.4f" % aot, "-AER.Model", "2",
           "-AER.SF.Model", "3", "-AER.SF.RH", "98." if aot == 0.0 else "80.",
           "-AER.DirMie", os.path.join(dbdir, "MIE_AER"),
           "-PHYTO.Chl", "0.0", "-SED.Csed", "0.0", "-YS.Abs440", "0.0", "-DET.Abs440", "0.0",
           "-HYD.Model", "3", "-HYD.ExtData", ext, "-HYD.UserProfile", prof,
           "-HYD.DirMie", os.path.join(dbdir, "MIE_HYD"),
           "-SEA.Depth", "%.4f" % depth, "-SEA.Dir", os.path.join(dbdir, "SURF"),
           "-SEA.Ind", "%.3f" % N_WATER, "-SEA.Wind", "%.2f" % WIND, "-SEA.SurfAlb", "0.0",
           "-SEA.BotType", "1", "-SEA.BotAlb", "%.6f" % alb,
           "-OSOAA.View.Phi", "0.0", "-OSOAA.View.Level", "3", "-OSOAA.View.VZA", "0.0",
           "-OSOAA.ResFile.vsVZA", "vsVZA.txt", "-OSOAA.ResFile.vsZ", "vsZ.txt"]
    env = dict(os.environ, OSOAA_ROOT=osoaa_root)
    p = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    vsz = os.path.join(res, "Standard_outputs", "vsZ.txt")
    vsv = os.path.join(res, "Standard_outputs", "vsVZA.txt")
    if p.returncode != 0 or not os.path.exists(vsz):
        raise RuntimeError("OSOAA failed for %s:\n%s" % (case_key(c), (p.stdout + p.stderr)[-2000:]))
    rrs_0m = None
    for x, refl in read_table(vsz, "Z"):
        if abs(x) < 1e-9:
            rrs_0m = refl/math.pi
            break
    refl_0p = None
    for x, refl in read_table(vsv, "VZA"):
        if abs(x) < 1e-6:
            refl_0p = refl
    if rrs_0m is None or refl_0p is None or refl_0p < 0:
        log = ""
        for root, _, files in os.walk(res):
            for name in files:
                if name.endswith(".Log") or name.endswith(".log"):
                    log += open(os.path.join(root, name), errors="replace").read()[-1500:]
        raise RuntimeError("bad OSOAA output for %s (rrs %s, REFL(0+) %s):\n%s\n%s" % (
            case_key(c), rrs_0m, refl_0p, p.stdout[-1500:], log))
    shutil.rmtree(res, ignore_errors=True)
    return rrs_0m, refl_0p


def read_table(path, first_column):
    """(first column, REFL) pairs from an OSOAA standard output file, after its header."""
    rows, started = [], False
    for line in open(path):
        f = line.split()
        if not started:
            started = len(f) > 3 and f[0] == first_column and f[3] == "REFL"
            continue
        if len(f) >= 4:
            rows.append((float(f[0]), float(f[3])))
    return rows


FIELDS = ["key", "water", "H", "bottom", "theta_air", "theta_w", "lam", "a", "bb", "u", "K",
          "kappaH", "bottom_albedo", "rrs_photic", "Rrs_photic", "rrs_osoaa", "Rrs_osoaa",
          "rrs_rel_err", "Rrs_rel_err", "Rrs_surface_osoaa", "a_rel_mismatch", "a_hyd", "b_hyd", "depth_osoaa"]


def run(driver, osoaa_root, workroot, quick=False, jobs=12, black_ocean_root=None, aot=0.0):
    sel = QUICK if quick else FULL
    exe = os.path.join(osoaa_root, "exe", "OSOAA_MAIN.exe")
    exe_black_ocean = None
    if black_ocean_root is not None:
        exe_black_ocean = os.path.join(black_ocean_root, "exe", "OSOAA_MAIN.exe")
    dbdir = os.path.join(workroot, "db")
    for d in ["MIE_AER", "MIE_HYD", "SURF"]:
        os.makedirs(os.path.join(dbdir, d), exist_ok=True)
    phase = phase_table()
    tables = photic_tables(driver, sel["bands"])
    cases = make_cases(sel)

    out_csv = os.path.join(workroot, "tier2_cases.csv")
    done = {}
    if os.path.exists(out_csv):
        for row in csv.DictReader(open(out_csv)):
            done[row["key"]] = row

    def one(c):
        key = case_key(c)
        cache = os.path.join(workroot, "cases", key, "result.json")
        if os.path.exists(cache):
            return json.load(open(cache))
        ph, theta_w = photic_side(c, tables)
        r = run_osoaa(c, ph, exe, osoaa_root, phase, os.path.join(workroot, "cases", key), dbdir,
                      exe_black_ocean=exe_black_ocean, aot=aot)
        H = c["H"]
        row = dict(key=key, water=c["water"], H=H, bottom=c["bottom"], theta_air=c["theta_air"],
                    theta_w=theta_w, lam=c["lam"], a=ph["a"], bb=ph["bb"], u=ph["u"], K=ph["K"],
                    kappaH=(math.inf if math.isinf(H) else ph["K"]*H),
                    rrs_photic=ph["rrs"], Rrs_photic=ph["Rrs"],
                    rrs_rel_err=(ph["rrs"] - r["rrs_osoaa"])/r["rrs_osoaa"],
                    Rrs_rel_err=(ph["Rrs"] - r["Rrs_osoaa"])/r["Rrs_osoaa"], **r)
        row = {k: (str(v) if isinstance(v, float) and math.isinf(v) else v) for k, v in row.items()}
        json.dump(row, open(cache, "w"))
        return row

    todo = [c for c in cases if case_key(c) not in done]
    # Surface matrices (per sun angle, db/SURF) and aerosol Mie files (per wavelength,
    # db/MIE_AER) are computed once and cached; compute one case for each serially first so
    # that parallel runs do not race to write (or read half-written) files.
    for theta in sel["theta"]:
        for lam in sel["bands"]:
            first = [c for c in todo if c["theta_air"] == theta and c["lam"] == lam][:1]
            for c in first:
                row = one(c)
                done[row["key"]] = row
                todo.remove(c)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for row in pool.map(one, todo):
            done[row["key"]] = row

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for c in cases:
            w.writerow(done[case_key(c)])

    rows = [done[case_key(c)] for c in cases]
    return summarise(rows, phase[2])


def summarise(rows, btilde):
    # photic's bottom spectra are normalised at 550 nm and exceed 3 at the red edge, so
    # B * rho_table can be above 1 (physically impossible). Those cases are reported but
    # left out of the statistics.
    unphysical = [r for r in rows if float(r["bottom_albedo"]) > 1.0]
    rows = [r for r in rows if float(r["bottom_albedo"]) <= 1.0]
    def arr(name, subset=None):
        return np.array([float(r[name]) for r in (subset or rows)])

    def err_stats(subset):
        e = np.abs(arr("rrs_rel_err", subset))
        E = np.abs(arr("Rrs_rel_err", subset))
        return dict(n=len(subset), rrs_median=float(np.median(e)), rrs_p95=float(np.percentile(e, 95)),
                    rrs_max=float(e.max()), Rrs_median=float(np.median(E)),
                    Rrs_p95=float(np.percentile(E, 95)), Rrs_max=float(E.max()),
                    rrs_bias=float(np.mean(arr("rrs_rel_err", subset))))

    out = dict(n_cases=len(rows), ff_backscatter_ratio=btilde, overall=err_stats(rows))
    out["excluded_albedo_above_1"] = dict(
        n=len(unphysical), max_albedo=max([float(r["bottom_albedo"]) for r in unphysical] or [0.0]),
        Rrs_max_rel_err=max([abs(float(r["Rrs_rel_err"])) for r in unphysical] or [0.0]))
    mism = arr("a_rel_mismatch")
    out["absorption_clamped"] = dict(n=int(np.sum(mism > 0)), max_rel=float(mism.max()))
    out["by_water"] = {}
    for g in WATER_TYPES:
        sub = [r for r in rows if r["water"] == g]
        if sub:
            out["by_water"][g] = err_stats(sub)
    for key, field, groups in [("by_theta_air", "theta_air", THETA_AIR), ("by_band", "lam", BANDS)]:
        out[key] = {}
        for g in groups:
            sub = [r for r in rows if float(r[field]) == g]
            if sub:
                out[key]["%g" % g] = err_stats(sub)
    out["by_depth"] = {}
    for H in DEPTHS:
        sub = [r for r in rows if float(r["H"]) == H]
        if sub:
            out["by_depth"]["inf" if math.isinf(H) else "%g" % H] = err_stats(sub)
    out["by_kappaH"] = {}
    for lo, hi in [(0, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 1e9)]:
        sub = [r for r in rows if lo <= float(r["kappaH"]) < hi]
        if sub:
            out["by_kappaH"]["%g-%g" % (lo, hi)] = err_stats(sub)
    return out
