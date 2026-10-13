"""Tier 1: component physics. Does photic's C code implement the equations?

1. Equations: photic's samodel_Rrs (via drivers/forward_dump) against the independent
   reference (photic_ref.forward, mode "photic"), given the same table values.
2. Departures from the paper: the size of each deliberate difference, as the relative
   change in Rrs over the case grid.
3. Tables: photic's interpolated a_w, b_bw and a0/a1 against their published sources.
4. Cross-implementation: SAMBUCA (CSIRO sambuca_core) on the same IOPs.
"""

import itertools
import math
import re
import subprocess

import numpy as np

import photic_ref

BANDS = [443.0, 490.0, 560.0, 665.0, 705.0]  # Sentinel-2 B1-B5 centres (nm).

GRID = dict(
    P=[0.01, 0.05, 0.2], G=[0.01, 0.1, 0.5], X=[0.001, 0.01, 0.05], B=[0.2, 0.6],
    bottom=[0, 1, 2], H=[0.5, 2.0, 5.0, 15.0], tide=[0.0, 0.5], theta_sun=[20.0, 45.0],
    theta_view=[0.0, 10.0],
    chi_inputs=[(0.008, 0.010, 0.0005), (0.020, 0.018, 0.0020)],
    lam=BANDS)

TOL_EQUATIONS = 1e-9


def make_cases(grid=GRID):
    keys = ["P", "G", "X", "B", "bottom", "H", "tide", "theta_sun", "theta_view",
            "chi_inputs", "lam"]
    rows = []
    for n, combo in enumerate(itertools.product(*[grid[k] for k in keys])):
        c = dict(zip(keys, combo))
        R440, R490, R750 = c.pop("chi_inputs")
        c.update(id=n, R440=R440, R490=R490, R750=R750)
        rows.append(c)
    cols = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
    return cols


def run_driver(driver, cases):
    lines = []
    for i in range(len(cases["id"])):
        lines.append("%d %.17g %.17g %.17g %.17g %.17g %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g" % (
            cases["id"][i], cases["lam"][i], cases["P"][i], cases["G"][i], cases["X"][i],
            cases["B"][i], cases["bottom"][i], cases["H"][i], cases["tide"][i],
            cases["theta_sun"][i], cases["theta_view"][i], cases["R440"][i],
            cases["R490"][i], cases["R750"][i]))
    out = subprocess.run([driver], input="\n".join(lines) + "\n", capture_output=True,
                         text=True, check=True).stdout
    header, *body = out.strip().split("\n")
    names = header.split(",")
    data = np.array([[float(v) for v in l.split(",")] for l in body])
    return {k: data[:, i] for i, k in enumerate(names)}


def rel(a, b):
    return np.abs(a - b)/np.maximum(np.abs(b), 1e-300)


def stats(x):
    return dict(median=float(np.median(x)), p95=float(np.percentile(x, 95)),
                max=float(np.max(x)))


def reference(cases, dump, **kw):
    return photic_ref.forward(
        cases["lam"], cases["P"], cases["G"], cases["X"], cases["B"], dump["rho_table"],
        cases["H"], cases["tide"], cases["theta_sun"], cases["theta_view"],
        cases["R440"], cases["R490"], cases["R750"],
        dump["a_0"], dump["a_1"], dump["a_w"], dump["b_bw"], **kw)


def blake_a0_a1(paper_tex):
    """The CSIRO-derived a0/a1 table the paper says photic uses (Appendix A)."""
    text = open(paper_tex).read()
    start = text.index("440nm-normalised")
    end = text.index("\\end{tabular}", start)
    table = {}
    for line in text[start:end].split("\n"):
        nums = re.findall(r"-?\d+\.?\d*", line)
        if len(nums) % 3 == 0 and len(nums) >= 3 and "&" in line:
            for k in range(0, len(nums), 3):
                table[float(nums[k])] = (float(nums[k + 1]), float(nums[k + 2]))
    wl = np.array(sorted(table))
    return wl, np.array([table[w][0] for w in wl]), np.array([table[w][1] for w in wl])


def pope_fry(path):
    wl, a = [], []
    for line in open(path, encoding="latin-1"):
        parts = line.split()
        if len(parts) == 2:
            try:
                wl.append(float(parts[0]))
                a.append(100.0*float(parts[1]))  # 1/cm -> 1/m
            except ValueError:
                pass
    return np.array(wl), np.array(a)


def sambuca_check(cases, dump, ref):
    """SAMBUCA on the same IOPs (theta_view = 0, tide = 0 subset)."""
    import sambuca_core as sc
    n_w = 1.33784
    sel = np.where((cases["theta_view"] == 0.0) & (cases["tide"] == 0.0))[0]
    out_same_bbw, out_own_bbw = [], []
    for i in sel:
        lam = cases["lam"][i]
        Y = ref["Y"][i]
        theta_air = math.degrees(math.asin(n_w*math.sin(math.radians(cases["theta_sun"][i]))))
        r = sc.forward_model(
            chl=1.0, cdom=cases["G"][i], nap=1.0, depth=cases["H"][i],
            substrate1=np.array([cases["B"][i]*dump["rho_table"][i]]),
            wavelengths=np.array([lam]), a_water=np.array([dump["a_w"][i]]),
            a_ph_star=np.array([ref["a_phi"][i]]), num_bands=1,
            a_cdom_slope=photic_ref.S_CDOM, lambda0cdom=440.0, a_cdom_lambda0cdom=1.0,
            a_nap_lambda0nap=0.0, x_ph_lambda0x=0.0, lambda0x=546.0,
            x_nap_lambda0x=cases["X"][i]*(440.0/546.0)**Y, bb_nap_slope=Y,
            water_refractive_index=n_w, theta_air=theta_air, off_nadir=0.0)
        out_own_bbw.append(rel(r.rrs[0], dump["rrs"][i]))
        # Same comparison with photic's b_bw replaced by SAMBUCA's, to isolate the physics.
        bbw_s = r.bb_water[0]
        alt = photic_ref.forward(
            lam, cases["P"][i], cases["G"][i], cases["X"][i], cases["B"][i],
            dump["rho_table"][i], cases["H"][i], 0.0, cases["theta_sun"][i], 0.0,
            cases["R440"][i], cases["R490"][i], cases["R750"][i],
            dump["a_0"][i], dump["a_1"][i], dump["a_w"][i], bbw_s, clamp_kappa=False)
        out_same_bbw.append(rel(r.rrs[0], alt["rrs"]))
    return dict(n=int(len(sel)), rrs_rel_diff=stats(np.array(out_own_bbw)),
                rrs_rel_diff_same_bbw=stats(np.array(out_same_bbw)))


def run(driver, paper_tex, pope_fry_path):
    cases = make_cases()
    dump = run_driver(driver, cases)
    assert np.array_equal(dump["id"], cases["id"])

    # 1. Equations.
    ref = reference(cases, dump)
    eq = {k: float(np.max(rel(dump[k], ref[k]))) for k in ["K", "rrs_dp", "rrs", "Rrs"]}
    eq["rrs_B"] = float(np.max(rel(dump["rrs_B"], ref["rrs_B"])))
    equations_pass = max(eq.values()) <= TOL_EQUATIONS

    # 2. Departures from the paper (each switched on alone, then all together).
    departures = {}
    for name, kw in [("chi_750nm_correction", dict(chi_750=True)),
                     ("no_Y_clamp", dict(clamp_Y=False)),
                     ("no_kappa_clamp", dict(clamp_kappa=False)),
                     ("paper_all", dict(chi_750=True, clamp_Y=False, clamp_kappa=False))]:
        alt = reference(cases, dump, **kw)
        departures[name] = stats(rel(alt["Rrs"], ref["Rrs"]))
    wl_b, a0_b, a1_b = blake_a0_a1(paper_tex)
    a0_bl = np.interp(cases["lam"], wl_b, a0_b)
    a1_bl = np.interp(cases["lam"], wl_b, a1_b)
    alt = photic_ref.forward(
        cases["lam"], cases["P"], cases["G"], cases["X"], cases["B"], dump["rho_table"],
        cases["H"], cases["tide"], cases["theta_sun"], cases["theta_view"],
        cases["R440"], cases["R490"], cases["R750"], a0_bl, a1_bl, dump["a_w"], dump["b_bw"])
    departures["a0_a1_paper_appendix_A"] = stats(rel(alt["Rrs"], ref["Rrs"]))

    # 3. Tables against their sources, at the S2 band centres.
    band = {lam: np.where(cases["lam"] == lam)[0][0] for lam in BANDS}
    wl_pf, a_pf = pope_fry(pope_fry_path)
    tables = {}
    for lam in BANDS:
        i = band[lam]
        entry = dict(a_w=float(dump["a_w"][i]), b_bw=float(dump["b_bw"][i]),
                     a0=float(dump["a_0"][i]), a1=float(dump["a_1"][i]),
                     b_bw_lee_formula=float(photic_ref.bbw_lee(lam)),
                     a0_appendix_A=float(np.interp(lam, wl_b, a0_b)),
                     a1_appendix_A=float(np.interp(lam, wl_b, a1_b)))
        if wl_pf[0] <= lam <= wl_pf[-1]:
            entry["a_w_pope_fry"] = float(np.interp(lam, wl_pf, a_pf))
            entry["a_w_rel_diff"] = float(rel(entry["a_w"], entry["a_w_pope_fry"]))
        entry["b_bw_rel_diff"] = float(rel(entry["b_bw"], entry["b_bw_lee_formula"]))
        tables["%d" % lam] = entry

    # 4. SAMBUCA.
    sambuca = sambuca_check(cases, dump, ref)

    return dict(
        n_cases=int(len(cases["id"])),
        equations=dict(max_rel_diff=eq, tolerance=TOL_EQUATIONS, passed=bool(equations_pass)),
        departures_rel_Rrs=departures,
        tables=tables,
        sambuca=sambuca)
