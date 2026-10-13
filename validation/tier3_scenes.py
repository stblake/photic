"""Tier 3: inversion on RT-simulated scenes. Can photic retrieve the truth?

Scenes are built from OSOAA, not from photic's own forward model, so photic is tested
against the real radiative transfer it approximates. For each scene (IOPs, sun angle,
tide) and band, OSOAA is run at depth nodes with bottom albedo 0, 0.5 and 1. For fixed
water and depth, Rrs(rho) = R0 + A rho/(1 - c rho): the bottom term saturates because
upwelling light is partly reflected back down at the surface and re-reflected by the
bottom. The three albedos determine R0, A and c; in depth, R0 and c are interpolated
linearly and A log-linearly. These lookup tables are cached in the data directory, so
after a change to photic's core only photic is re-run.

Each scenario (water type, number of scenes, noise, depth prior) writes .asc bands and a
command file, runs MODEL Lee_Semi_Analytical, and scores depth, albedo and K_min against
the truth. A control set uses photic's own forward model without noise, which separates
optimiser failure from model-structure error.
"""

import json
import math
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import photic_ref
import tier2_osoaa as t2

NROWS = NCOLS = 32
WLON, SLAT, CELL = 115.0, -21.0, 0.0003
BANDS = [443.0, 490.0, 560.0, 665.0]
DEPTH_NODES = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0,
               18.0, 22.0, 27.0, 32.0]
B_ALBEDO = 0.6
ALBEDO_NODES = (0.0, 0.5, 1.0)
BOTTOM_COLUMNS = [(0, 16, "sand"), (16, 24, "seagrass"), (24, 32, "coral")]
DEPTH_BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 30)]

# Up to four scenes per water type: P, G, X, Y, sun zenith in air (deg), tide (m).
SCENE_SETS = {
    "clear": [(0.02, 0.02, 0.004, 1.0, 30.0, 0.0), (0.03, 0.03, 0.005, 1.0, 45.0, 0.4),
              (0.015, 0.025, 0.0035, 1.0, 25.0, -0.3), (0.025, 0.02, 0.0045, 1.0, 50.0, 0.8)],
    "moderate": [(0.08, 0.15, 0.02, 1.0, 30.0, 0.0), (0.10, 0.18, 0.025, 1.0, 45.0, 0.4),
                 (0.07, 0.13, 0.018, 1.0, 25.0, -0.3), (0.09, 0.16, 0.022, 1.0, 50.0, 0.8)],
    "turbid": [(0.30, 0.50, 0.08, 1.0, 30.0, 0.0), (0.35, 0.60, 0.09, 1.0, 45.0, 0.4),
               (0.28, 0.45, 0.075, 1.0, 25.0, -0.3), (0.32, 0.55, 0.085, 1.0, 50.0, 0.8)],
}

FULL = dict(water=["clear", "moderate", "turbid"], n_scenes=[1, 2, 4], noise=[0.0, 5e-4],
            prior=[False, True])
QUICK = dict(water=["moderate"], n_scenes=[2], noise=[5e-4], prior=[False])


def truth():
    i = np.arange(NROWS)[:, None]
    j = np.arange(NCOLS)[None, :]
    H = 0.8 + 24.5*(i/(NROWS - 1.0))**1.3 + 0.3*np.sin(2*np.pi*j/16.0)
    bottom = np.empty((NROWS, NCOLS), dtype=object)
    for j0, j1, name in BOTTOM_COLUMNS:
        bottom[:, j0:j1] = name
    return H, bottom


def theta_w(theta_air):
    return math.degrees(math.asin(math.sin(math.radians(theta_air))/t2.N_WATER))


def scene_key(water, k):
    return "%s_s%d" % (water, k + 1)


def build_lut(water, k, tables, osoaa, workroot, jobs):
    """OSOAA Rrs(0+) at the depth nodes for albedo 0 and 1, per band, cached as JSON."""
    P, G, X, Y, theta_air, tide = SCENE_SETS[water][k]
    t2.WATER_TYPES[scene_key(water, k)] = (P, G, X, Y)
    path = os.path.join(workroot, "lut", scene_key(water, k) + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    jobs_list = []
    for lam in BANDS:
        for H in DEPTH_NODES:
            for alb in ALBEDO_NODES:
                jobs_list.append((lam, H, alb))

    def one(args):
        lam, H, alb = args
        c = dict(water=scene_key(water, k), H=H, bottom="sand", theta_air=theta_air, lam=lam)
        ph, _ = t2.photic_side(c, tables)
        wd = os.path.join(workroot, "osoaa", scene_key(water, k), "%g_%g_%g" % (lam, H, alb))
        black = osoaa["black_ocean"] if (H == DEPTH_NODES[0] and alb == 0.0) else None
        r = t2.run_osoaa(c, ph, osoaa["exe"], osoaa["root"], osoaa["phase"], wd, osoaa["db"],
                         exe_black_ocean=black, albedo=alb)
        return args, r

    # Surface matrices for this sun angle first (avoids parallel writes of the same file).
    first, rest = jobs_list[:1], jobs_list[1:]
    results = [one(a) for a in first]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results += list(pool.map(one, rest))

    surface = {}
    for (lam, H, alb), r in results:
        if H == DEPTH_NODES[0] and alb == 0.0:
            surface[lam] = r["Rrs_surface_osoaa"]
    lut = {}
    for lam in BANDS:
        R0, A, C = [], [], []
        for H in DEPTH_NODES:
            v = {}
            for (l_, H_, a_), r in results:
                if l_ == lam and H_ == H:
                    # Water-leaving Rrs: total REFL(0+) minus this band's surface reflection.
                    v[a_] = r["Rrs_osoaa"] + r["Rrs_surface_osoaa"] - surface[lam]
            a, b = v[0.5] - v[0.0], v[1.0] - v[0.0]
            c = (b - 2.0*a)/(b - a) if abs(b - a) > 1e-15 else 0.0
            R0.append(v[0.0])
            A.append(b*(1.0 - c))
            C.append(c)
        lut["%g" % lam] = dict(R0=R0, A=A, C=C)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(lut, open(path, "w"))
    return lut


def lut_rrs(lut, lam, Hs, albedo):
    """Rrs from the lookup table at scene depth Hs and bottom albedo."""
    nodes = np.array(DEPTH_NODES)
    t = lut["%g" % lam]
    Hs = np.clip(Hs, nodes[0], nodes[-1])
    r0 = np.interp(Hs, nodes, np.array(t["R0"]))
    A = np.exp(np.interp(Hs, nodes, np.log(np.maximum(np.array(t["A"]), 1e-12))))
    c = np.interp(Hs, nodes, np.array(t["C"]))
    return r0 + A*albedo/(1.0 - c*albedo)


def write_asc(path, grid, nodata=-9999.0):
    with open(path, "w") as f:
        f.write("ncols %d\nnrows %d\nxllcorner %.6f\nyllcorner %.6f\ncellsize %.6f\nNODATA_value %g\n"
                % (NCOLS, NROWS, WLON, SLAT, CELL, nodata))
        for i in range(NROWS - 1, -1, -1):  # north row first; photic's row 0 is south
            f.write(" ".join("%.8f" % v for v in grid[i]) + "\n")


def scene_bands(water, k, H, bottom, tables, lut, control):
    """Rrs grids for one scene: from the OSOAA lookup table, or (control) photic's own model."""
    P, G, X, Y, theta_air, tide = SCENE_SETS[water][k]
    out = {}
    if not control:
        for lam in BANDS:
            rho = np.zeros(H.shape)
            for name, b in t2.BOTTOMS.items():
                rho[bottom == name] = B_ALBEDO*tables[(lam, b)]["rho"]
            out[lam] = lut_rrs(lut, lam, H + tide, rho)
        return out
    # Control: photic's forward model with its own chi -> Y relation (fixed point, as in
    # model/tests/synth_gen.c).
    tw = theta_w(theta_air)
    rrs = {lam: np.zeros(H.shape) for lam in BANDS}
    for name, b in t2.BOTTOMS.items():
        m = bottom == name
        R440 = np.full(m.sum(), 0.01)
        R490 = np.full(m.sum(), 0.01)
        for _ in range(100):
            vals = {}
            for lam in BANDS:
                t = tables[(lam, b)]
                vals[lam] = photic_ref.forward(lam, P, G, X, B_ALBEDO, t["rho"], H[m], tide, tw, 0.0,
                                               R440, R490, 0.0, t["a0"], t["a1"], t["aw"], t["bbw"])["Rrs"]
            spec = np.array([vals[l] for l in BANDS])
            # photic's interp_1d extrapolates linearly from the first two bands (440 < 443).
            slope = (spec[1] - spec[0])/(BANDS[1] - BANDS[0])
            new440 = spec[0] + slope*(440.0 - BANDS[0])
            new490 = np.array([np.interp(490.0, BANDS, spec[:, n]) for n in range(spec.shape[1])])
            change = np.max(np.abs(new440 - R440) + np.abs(new490 - R490))
            R440, R490 = new440, new490
            if change < 1e-12:
                break
        for lam in BANDS:
            rrs[lam][m] = vals[lam]
    return rrs


def read_nc(path):
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        name = [v for v in ds.variables if v not in ("lons", "lats")][0]
        return np.ma.filled(ds.variables[name][:].astype(float), np.nan)


def score(H_true, H_model, albedo_model, kmin_model, kmin_true, elapsed):
    valid = np.isfinite(H_model) & (H_model > 0)
    err = H_model - H_true
    e = err[valid]
    out = dict(frac_modelled=float(valid.mean()), seconds=elapsed)
    if e.size == 0:
        return out
    tol = 0.5 + 0.1*H_true[valid]
    out.update(bias=float(e.mean()), rmse=float(np.sqrt(np.mean(e**2))), mae=float(np.mean(np.abs(e))),
               median_abs=float(np.median(np.abs(e))), frac_within_tol=float(np.mean(np.abs(e) <= tol)))
    out["by_depth"] = {}
    for lo, hi in DEPTH_BINS:
        m = valid & (H_true >= lo) & (H_true < hi)
        if m.any():
            d = err[m]
            out["by_depth"]["%g-%g" % (lo, hi)] = dict(
                n=int(m.sum()), bias=float(d.mean()), rmse=float(np.sqrt(np.mean(d**2))),
                frac_within_tol=float(np.mean(np.abs(d) <= 0.5 + 0.1*H_true[m])))
    a = albedo_model[valid]
    a = a[np.isfinite(a)]
    if a.size:
        out["albedo_bias"] = float(np.mean(a - B_ALBEDO))
        out["albedo_mae"] = float(np.mean(np.abs(a - B_ALBEDO)))
    k = kmin_model[valid]
    k = k[np.isfinite(k)]
    if k.size:
        out["kmin_rel_bias"] = float(np.mean(k - kmin_true)/kmin_true)
    return out


def run_photic(photic, rundir, water, n, H, bottom, tables, luts, noise, prior, control):
    os.makedirs(rundir, exist_ok=True)
    rng = np.random.default_rng(12345)
    lines = []
    scenes = []
    for k in range(n):
        P, G, X, Y, theta_air, tide = SCENE_SETS[water][k]
        bands = scene_bands(water, k, H, bottom, tables, None if control else luts[k], control)
        names = []
        for lam in BANDS:
            g = bands[lam] + (rng.normal(0.0, noise, H.shape) if noise > 0 else 0.0)
            name = "s%d_b%d" % (k + 1, lam)
            write_asc(os.path.join(rundir, name + ".asc"), g)
            lines.append("READ GRID %s %s.asc" % (name, name))
            names.append(name)
        scenes.append("SCENE NAME s%d BANDS %s WAVELENGTHS %s THETAVIEW 0 THETASUN %.4f HTIDE %.3f" % (
            k + 1, " ".join(names), " ".join("%d" % l for l in BANDS), theta_w(theta_air), tide))
    model = "MODEL Lee_Semi_Analytical IN %s OUT H H_sigma H_error albedo sand seagrass coral kmin btype iod" % (
        " ".join("s%d" % (k + 1) for k in range(n)))
    if prior:
        # A rough empirical depth estimate (+-20%), negative down as photic expects.
        est = -H*(1.0 + rng.normal(0.0, 0.2, H.shape))
        write_asc(os.path.join(rundir, "prior.asc"), est)
        lines.append("READ GRID prior prior.asc")
        model += " DEPTHS prior"
    lines += scenes + [model, "EXIT"]
    open(os.path.join(rundir, "run.bam"), "w").write("\n".join(lines) + "\n")
    t0 = time.time()
    with open(os.path.join(rundir, "photic.log"), "w") as log:
        subprocess.run([photic, "run.bam"], cwd=rundir, stdin=subprocess.DEVNULL, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    elapsed = time.time() - t0
    Hm = -read_nc(os.path.join(rundir, "modelled_H.nc"))
    alb = read_nc(os.path.join(rundir, "modelled_albedo.nc"))
    kmin = read_nc(os.path.join(rundir, "modelled_min_K.nc"))
    kmin_true = min(min(photic_side_K(water, k, lam, tables) for lam in BANDS) for k in range(n))
    out = score(H, Hm, alb, kmin, kmin_true, elapsed)
    if prior:
        # The prior on its own: photic only adds skill if it beats this.
        p = score(H, -est, np.full(H.shape, np.nan), np.full(H.shape, np.nan), kmin_true, 0.0)
        out["prior_only"] = {k: p[k] for k in ("bias", "rmse", "mae", "frac_within_tol")}
        if out.get("rmse"):
            out["rmse_skill_vs_prior"] = 1.0 - out["rmse"]/p["rmse"]
    return out


def photic_side_K(water, k, lam, tables):
    P, G, X, Y, theta_air, tide = SCENE_SETS[water][k]
    t = tables[(lam, 0)]
    r = photic_ref.forward(lam, P, G, X, B_ALBEDO, t["rho"], 5.0, 0.0, 20.0, 0.0, 0.01, 0.01, 0.0,
                           t["a0"], t["a1"], t["aw"], t["bbw"], Y=Y)
    return float(r["K"])


def run(driver, photic, data, workroot, quick=False, jobs=12):
    sel = QUICK if quick else FULL
    osoaa = dict(root=os.path.join(data, "tools", "OSOAA_scalar"),
                 exe=os.path.join(data, "tools", "OSOAA_scalar", "exe", "OSOAA_MAIN.exe"),
                 black_ocean=os.path.join(data, "tools", "OSOAA_blackocean", "exe", "OSOAA_MAIN.exe"),
                 db=os.path.join(workroot, "db"), phase=t2.phase_table())
    for d in ["MIE_AER", "MIE_HYD", "SURF"]:
        os.makedirs(os.path.join(osoaa["db"], d), exist_ok=True)
    tables = t2.photic_tables(driver, BANDS)
    H, bottom = truth()
    results = {}
    for water in sel["water"]:
        n_max = max(sel["n_scenes"])
        luts = [build_lut(water, k, tables, osoaa, workroot, jobs) for k in range(n_max)]
        for n in sel["n_scenes"]:
            for noise in sel["noise"]:
                for prior in sel["prior"]:
                    name = "%s_n%d_noise%g%s" % (water, n, noise, "_prior" if prior else "")
                    results[name] = run_photic(photic, os.path.join(workroot, "photic", name), water, n,
                                               H, bottom, tables, luts, noise, prior, control=False)
        name = "%s_n%d_control" % (water, 2)
        results[name] = run_photic(photic, os.path.join(workroot, "photic", name), water, 2, H, bottom,
                                   tables, None, 0.0, False, control=True)
    results["lut_check"] = lut_check(tables, osoaa, workroot, sel["water"][0])
    return dict(grid="%dx%d" % (NROWS, NCOLS), depth_range=[float(H.min()), float(H.max())],
                scenarios=results)


def lut_check(tables, osoaa, workroot, water):
    """Interpolation error of the lookup table at depths between the nodes."""
    path = os.path.join(workroot, "lut", scene_key(water, 0) + ".json")
    lut = json.load(open(path))
    P, G, X, Y, theta_air, tide = SCENE_SETS[water][0]
    errs = []
    for lam in [490.0, 560.0]:
        for H in [0.6, 2.5, 7.0, 13.5]:
            c = dict(water=scene_key(water, 0), H=H, bottom="sand", theta_air=theta_air, lam=lam)
            ph, _ = t2.photic_side(c, tables)
            alb = B_ALBEDO*tables[(lam, 0)]["rho"]
            wd = os.path.join(workroot, "osoaa_check", "%g_%g" % (lam, H))
            r = t2.run_osoaa(c, ph, osoaa["exe"], osoaa["root"], osoaa["phase"], wd, osoaa["db"],
                             exe_black_ocean=osoaa["black_ocean"], albedo=alb)
            errs.append(abs(lut_rrs(lut, lam, H, alb) - r["Rrs_osoaa"])/r["Rrs_osoaa"])
    return dict(max_rel_err=float(np.max(errs)), median_rel_err=float(np.median(errs)), n=len(errs))
