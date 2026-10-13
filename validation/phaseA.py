"""Phase A of EXPERIMENTS.md: diagnose photic's inversion one pixel at a time.

Uses drivers/invert_lab (photic's own model_data, samodel_error, nelmin and
samodel_optimise) on pixels of the tier 3 scenes, where the truth is known:

  A0  photic's per-pixel optimiser, with no LUT and no carry-over between pixels;
      A0x is the same search without the early exit, keeping the best of all depth starts.
  A1  error against depth: conditional (everything else true) and profiled (everything
      else re-fitted at each fixed depth).
  A2  Nelder-Mead started at the truth, with photic's steps.
  A3  Nelder-Mead from photic's start with the water (PGX), the bottom (B, q) or the
      depth fixed at the truth.

  phaseA.py [--water moderate] [--scenes 2] [--osoaa] [--noise 0] [--nspatial 2]

Writes a JSON result per scenario to $PHOTIC_VALIDATION_DATA/runs/phaseA/.
"""

import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tier2_osoaa as t2
import tier3_scenes as t3

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))
LAB = os.environ.get("PHOTIC_INVERT_LAB", os.path.join(HERE, "drivers", "invert_lab"))
DRIVER = os.path.join(HERE, "drivers", "forward_dump")
N_BOTTOMS = 3  # photic's default NBOTTOMS: sand, seagrass, coral (t2.BOTTOMS indexes)
H_STARTS = [40.0, 30.0, 20.0, 15.0, 10.0, 7.5, 2.5, 1.0]  # samodel_optimise_one_bottom_combination
NM_PHOTIC = (1.0e-2, 100, 5000)  # reqmin, konvge, kcount
PROFILE_H = [0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 26, 30, 35, 40]


def scene_grids(water, n, osoaa, noise, seed=12345):
    """Rrs grids [scene][band] for the tier 3 truth, as tier3_scenes.run_photic makes them."""
    tables = t2.photic_tables(DRIVER, t3.BANDS)
    H, bottom = t3.truth()
    luts = None
    if osoaa:
        work = os.path.join(DATA, "runs", "tier3")
        luts = [t3.build_lut(water, k, tables, None, work, 1) for k in range(n)]
    rng = np.random.default_rng(seed)
    grids = []
    for k in range(n):
        bands = t3.scene_bands(water, k, H, bottom, tables, luts[k] if osoaa else None, not osoaa)
        grids.append([bands[lam] + (rng.normal(0.0, noise, H.shape) if noise > 0 else 0.0)
                      for lam in t3.BANDS])
    return H, bottom, grids


def window(i, j, nspatial):
    """Region cells in extract_Rrs_data order (edges clamped), and the origin index."""
    cells, origin = [], 0
    s = max(nspatial, 1)
    for ii in range(1 - s, s):
        for jj in range(1 - s, s):
            c = (min(max(i + ii, 0), t3.NROWS - 1), min(max(j + jj, 0), t3.NCOLS - 1))
            if c == (i, j):
                origin = len(cells)
            cells.append(c)
    return cells, origin


def pixel_setup(water, n, i, j, H, bottom, grids, nspatial, noise=0.0):
    cells, origin = window(i, j, nspatial)
    nr = len(cells)
    lines = ["SETUP %d %d %d %d %d" % (n, len(t3.BANDS), nr, origin, N_BOTTOMS)]
    for k in range(n):
        P, G, X, Y, theta_air, tide = t3.SCENE_SETS[water][k]
        lines.append("SCENE %d %.10g 0 %.10g %s" % (k, t3.theta_w(theta_air), tide,
                                                    " ".join("%g" % l for l in t3.BANDS)))
    for r, (a, b) in enumerate(cells):
        for k in range(n):
            lines.append("RRS %d %d %s" % (r, k, " ".join("%.12g" % grids[k][m][a, b]
                                                          for m in range(len(t3.BANDS)))))
    if noise > 0:
        for k in range(n):
            lines.append("NOISE %d %s" % (k, " ".join(["%g" % noise]*len(t3.BANDS))))
    lines.append("END")
    # Truth as physical values; packed into photic's layout by the driver (PACK).
    q = []
    for a, b in cells:
        row = [0.0]*N_BOTTOMS
        row[t2.BOTTOMS[bottom[a, b]]] = 1.0
        q += row
    pgx = []
    for k in range(n):
        pgx += list(t3.SCENE_SETS[water][k][:3])
    phys = [H[a, b] for a, b in cells] + [t3.albedo_for(bottom[a, b]) for a, b in cells] + q + pgx
    setup = "\n".join(lines) + "\n"
    truth = lab(setup, ["PACK " + vec(phys)])[0]["PACKED"]
    return setup, np.array(truth), nr, origin, np.array(phys)


def shifted(setup, phys, nr, origin, depths):
    """Truth vectors with all regional depths shifted so the origin is at each depth."""
    cmds = []
    for h in depths:
        p = phys.copy()
        p[:nr] = np.maximum(phys[:nr] + (h - phys[origin]), 0.01)
        cmds.append("PACK " + vec(p))
    return [np.array(r["PACKED"]) for r in lab(setup, cmds)]



def lab(setup, commands):
    """Run invert_lab; returns one parsed dict per command."""
    p = subprocess.run([LAB], input=setup + "\n".join(commands) + "\n", capture_output=True, text=True)
    out, cur = [], {}
    for line in p.stdout.splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] == "DONE":
            out.append(cur)
            cur = {}
        elif w[0] in ("START", "STEP", "PACKED", "params", "depths", "PGX", "bottom_percent"):
            cur[w[0]] = [float(v) for v in w[1:]]
        elif w[0] == "Rrs_modelled":
            cur.setdefault("Rrs_modelled", []).append([float(v) for v in w[3:]])
        elif w[0] in ("model_error", "NM", "PHOTIC", "LM"):
            kv = w[1:] if w[0] != "model_error" else w
            for a, b in zip(kv[0::2], kv[1::2]):
                cur[a] = float(b)
        elif w[0] == "ERROR:":
            raise RuntimeError(line)
    if len(out) != len(commands):
        raise RuntimeError("invert_lab: %d results for %d commands\n%s" % (len(out), len(commands), p.stdout[-2000:]))
    return out


def vec(v):
    return " ".join("%.12g" % x for x in v)


def nm(start, step, budget=NM_PHOTIC):
    return "NM %g %d %d %s %s" % (budget[0], budget[1], budget[2], vec(start), vec(step))


def depth_of(fit, nr, origin):
    """photic's reported depth from an NM result: origin weighted sqrt(n_regions) in the
    regional mean of the depths the driver unpacked."""
    H = np.asarray(fit["depths"])
    w = np.ones(nr)
    w[origin] = np.sqrt(nr)
    return float(np.sum(w*H)/np.sum(w))


def analyse_pixel(args):
    water, n, i, j, H, bottom, grids, nspatial = args
    setup, truth, nr, origin, phys = pixel_setup(water, n, i, j, H, bottom, grids, nspatial)
    shifts = shifted(setup, phys, nr, origin, PROFILE_H)
    nb = N_BOTTOMS
    iH = slice(0, nr)
    iB = slice(nr, nr*(1 + nb))
    iW = slice(nr*(1 + nb), len(truth))
    res = dict(i=i, j=j, H_true=float(H[i, j]), bottom=str(bottom[i, j]))

    # Starting points (photic's) for every depth start, and the truth's own error.
    starts = lab(setup, ["START %g" % h for h in H_STARTS] + ["START %g" % H[i, j], "EVAL " + vec(truth)])
    res["truth_error"] = starts[-1]["model_error"]
    res["truth_components"] = {k: starts[-1][k] for k in ("Rrs_error", "depth_error", "bottom_error", "K_error")}
    truth_step = np.array(starts[-2]["STEP"])

    cmds = ["PHOTIC -1"]                                                     # A0
    cmds += [nm(s["START"], s["STEP"]) for s in starts[:len(H_STARTS)]]     # A0x
    # A1 conditional: shift all regional depths together, everything else true.
    for v in shifts:
        cmds.append("EVAL " + vec(v))
    # A1 profiled: depths fixed at the shifted values; the rest re-fitted from the truth
    # and from photic's start (for that depth).
    prof_idx = []
    for h, v in zip(PROFILE_H, shifts):
        st = truth_step.copy()
        st[iH] = 0.0
        cmds.append(nm(v, st, (1e-6, 100, 20000)))
        p_start = lab(setup, ["START %g" % h])[0]
        w = np.array(p_start["START"])
        w[iH] = v[iH]
        st2 = np.array(p_start["STEP"])
        st2[iH] = 0.0
        cmds.append(nm(w, st2, (1e-6, 100, 20000)))
        prof_idx.append(h)
    # A2: from the truth with photic's steps.
    cmds.append(nm(truth, truth_step))
    # A3: one block fixed at the truth, the rest from photic's best-guess start (10 m).
    s10 = lab(setup, ["START 10"])[0]
    for name, blk in (("water_fixed", iW), ("bottom_fixed", iB), ("depth_fixed", iH)):
        v = np.array(s10["START"])
        st = np.array(s10["STEP"])
        v[blk] = truth[blk]
        st[blk] = 0.0
        cmds.append(nm(v, st, (1e-6, 100, 20000)))
    r = lab(setup, cmds)

    k = 0
    a0 = r[k]; k += 1
    res["A0"] = dict(depth=a0["depth"], albedo=a0["albedo"], K_min=a0["K_min"], error=a0["model_error"],
                     Rrs_error=a0["Rrs_error"], K_error=a0["K_error"])
    fits = r[k:k + len(H_STARTS)]; k += len(H_STARTS)
    res["A0x_starts"] = [dict(h_start=h, error=f["error"], depth=depth_of(f, nr, origin))
                         for h, f in zip(H_STARTS, fits)]
    best = min(fits, key=lambda f: f["error"])
    res["A0x"] = dict(depth=depth_of(best, nr, origin), error=best["error"])
    first_ok = next((f for f in fits if f["error"] < 2.5*n), None)
    res["A0_early_exit_depth"] = depth_of(first_ok, nr, origin) if first_ok else None
    res["A1_conditional"] = [dict(h=h, error=r[k + m]["model_error"], Rrs_error=r[k + m]["Rrs_error"],
                                  K_error=r[k + m]["K_error"], depth_error=r[k + m]["depth_error"])
                             for m, h in enumerate(PROFILE_H)]
    k += len(PROFILE_H)
    prof = []
    for h in prof_idx:
        a, b = r[k], r[k + 1]
        k += 2
        f = a if a["error"] <= b["error"] else b
        prof.append(dict(h=h, error=f["error"], Rrs_error=f["Rrs_error"], K_error=f["K_error"],
                         depth_error=f["depth_error"], bottom_error=f["bottom_error"]))
    res["A1_profiled"] = prof
    a2 = r[k]; k += 1
    res["A2"] = dict(depth=depth_of(a2, nr, origin), error=a2["error"],
                     moved=float(np.max(np.abs(np.asarray(a2["depths"]) - np.array([H[c] for c in window(i, j, nspatial)[0]])))))
    for name in ("water_fixed", "bottom_fixed", "depth_fixed"):
        f = r[k]; k += 1
        res["A3_" + name] = dict(depth=depth_of(f, nr, origin), error=f["error"],
                                 Rrs_error=f["Rrs_error"], K_error=f["K_error"])
    return res


def summarise(rows):
    def stats(key):
        sel = [d for d in flat if d[key] is not None]
        e = np.array([d[key] - d["H_true"] for d in sel])
        t = np.array([d["H_true"] for d in sel])
        return dict(n=len(sel), bias=float(e.mean()), rmse=float(np.sqrt(np.mean(e**2))),
                    within_tol=float(np.mean(np.abs(e) <= 0.5 + 0.1*t)))
    flat = []
    for r in rows:
        d = dict(H_true=r["H_true"], A0=r["A0"]["depth"], A0x=r["A0x"]["depth"], A2=r["A2"]["depth"],
                 early=r["A0_early_exit_depth"])
        for name in ("water_fixed", "bottom_fixed"):
            d[name] = r["A3_" + name]["depth"]
        cond = min(r["A1_conditional"], key=lambda x: x["error"])
        prof = min(r["A1_profiled"], key=lambda x: x["error"])
        d["argmin_conditional"] = cond["h"]
        d["argmin_profiled"] = prof["h"]
        d["truth_beats_profile_min"] = r["truth_error"] <= prof["error"] + 1e-9
        flat.append(d)
    out = {k: stats(k) for k in ("A0", "early", "A0x", "A2", "water_fixed", "bottom_fixed", "argmin_conditional",
                                 "argmin_profiled")}
    out["frac_truth_is_global_min"] = float(np.mean([d["truth_beats_profile_min"] for d in flat]))
    return out, flat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--water", default="moderate")
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("--osoaa", action="store_true")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--nspatial", type=int, default=2)
    ap.add_argument("--rows", default="0,3,6,9,12,15,18,21,24,27,31")
    ap.add_argument("--cols", default="8,20,28")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--tag", default="", help="appended to the output name, e.g. the branch")
    a = ap.parse_args()
    H, bottom, grids = scene_grids(a.water, a.scenes, a.osoaa, a.noise)
    pix = [(a.water, a.scenes, i, j, H, bottom, grids, a.nspatial)
           for i in map(int, a.rows.split(",")) for j in map(int, a.cols.split(","))]
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        rows = list(pool.map(analyse_pixel, pix))
    summary, flat = summarise(rows)
    name = "%s_n%d_%s_noise%g_ns%d" % (a.water, a.scenes, "osoaa" if a.osoaa else "control", a.noise, a.nspatial)
    if a.tag:
        name += "_" + a.tag
    out = os.path.join(DATA, "runs", "phaseA")
    os.makedirs(out, exist_ok=True)
    json.dump(dict(summary=summary, pixels=rows), open(os.path.join(out, name + ".json"), "w"), indent=1)
    print(name)
    print(json.dumps(summary, indent=1))
    print("%6s %8s %7s %7s %7s %7s %7s %7s %7s %5s" % ("H", "bottom", "A0", "A0x", "A2", "Wfix", "Bfix",
                                                     "argC", "argP", "truth"))
    for r, d in zip(rows, flat):
        print("%6.2f %8s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %5s" % (
            d["H_true"], r["bottom"], d["A0"], d["A0x"], d["A2"], d["water_fixed"], d["bottom_fixed"],
            d["argmin_conditional"], d["argmin_profiled"], "min" if d["truth_beats_profile_min"] else "-"))


if __name__ == "__main__":
    main()
