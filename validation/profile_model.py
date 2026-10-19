"""Phase F prototype of DROP, Deep-water Referenced Optical Profiling (MODEL DROP; EXPERIMENTS.md, Phase F).

Lee's shallow-water model with the water observed rather than fitted. Each scene's deep-water
reflectance rrs_inf comes from its optically deep pixels, and kappa and u from rrs_inf by QAA.
Then, per pixel and over all scenes at once,

    rrs_s = rrs_inf_s [1 - exp(-kappa_s M^C_s Z_s)] + exp(-kappa_s M^B_s Z_s) b rho_c / pi,

with Z_s = H + tide_s and the bottom one class c with brightness b in the class range (B8),
shared by the scenes. b (and, with --glint, a flat offset per scene) is linear, so for every
depth on a grid and every class it is solved exactly: the profile chi2(H) is exact and its
global minimum is found by construction. A depth is reported only when the optically deep
alternative is worse by more than --t-deep and no depth outside the tolerance scores within
--t-amb of the best.

Runs on ab.py's points and compares with a photic variant on the same points, with ab.py's
paired scoring:

    profile_model.py [--split train] [--n 300] [--seed 1] [--against HEAD] [--ref red]
                     [--eps 0.03] [--t-deep 9] [--t-amb 0] [--glint] [--sites sites5]

--ref: where rrs_inf comes from. "red" (default) the pixels with the lowest red (and 705 nm)
reflectance, no soundings; "oracle" tier 3's deepest rows (truth >= 22 m) or tier 4's sounded
cells in the deepest decile, a diagnostic of identifiability only.
"""

import argparse
import datetime
import json
import os
import re
import time

import numpy as np

import ab
import tier2_osoaa as t2
from tier4_prior import read_asc

# Bottom classes: photic's table index and brightness range (samodel.c bright_lo/hi, B8).
CLASSES = {"sand": (0, 0.15, 0.75), "seagrass": (1, 0.03, 0.15), "coral": (2, 0.03, 0.25)}
H_GRID = np.concatenate([np.linspace(0.05, 2.0, 40), np.linspace(2.05, 10.0, 80), np.linspace(10.2, 30.0, 100)])
G0, G1 = 0.084, 0.170  # Lee 1999, photic's
GLINT_SIGMA = 5e-4     # prior on the per-scene offset in Rrs, as GLINT_OFFSET


def rrs_of(R):
    return R/(0.5 + 1.5*R)


# Inputs.

def read_case(case_dir):
    """Scenes of a case: list of dict(lam, theta_w, theta_v, tide, R (L, ny, nx))."""
    files, scenes = {}, []
    for line in open(os.path.join(case_dir, "run.bam")):
        tok = line.split()
        if line.startswith("READ GRID"):
            files[tok[2]] = tok[3]
        elif line.startswith("SCENE NAME"):
            i = tok.index("BANDS")
            j = tok.index("WAVELENGTHS")
            bands = tok[i + 1:j]
            lam = [float(x) for x in tok[j + 1:j + 1 + len(bands)]]
            get = lambda key, d=0.0: float(tok[tok.index(key) + 1]) if key in tok else d
            R = np.array([read_asc(os.path.join(case_dir, files[b])) for b in bands])
            scenes.append(dict(name=tok[2], lam=lam, theta_w=get("THETASUN"), theta_v=get("THETAVIEW"),
                               tide=get("HTIDE"), R=R))
    return scenes


def smoothed(scene, size):
    """Valid mask and the scene's red (+705 nm) and green reflectance, smoothed over size x size
    pixels (a nan-aware mean), so that ranking pixels by them does not select noise."""
    from scipy import ndimage
    R = scene["R"]
    valid = np.all(np.isfinite(R) & (R > 0), axis=0)
    lam = np.array(scene["lam"])
    red = np.where(valid, R[lam >= 660.0].sum(axis=0), 0.0)
    green = np.where(valid, R[int(np.argmin(np.abs(lam - 560.0)))], 0.0)
    if size > 1:
        n = ndimage.uniform_filter(valid.astype(float), size, mode="constant")
        red = ndimage.uniform_filter(red, size, mode="constant")/np.maximum(n, 1e-9)
        green = ndimage.uniform_filter(green, size, mode="constant")/np.maximum(n, 1e-9)
    return valid, red, green


def deep_reference(scene, truth, a, fields, region=None):
    """rrs_inf (L), the deep pixels' noise in rrs (L) and their number, for one scene, from the
    pixels in region (a mask; default the whole site)."""
    R = scene["R"]
    valid, red, green = fields
    if region is not None:
        valid = valid & region
    if a.ref == "oracle":
        t = np.where(valid, truth, np.nan)
        deep = valid & (t >= (22.0 if np.nanmax(t) > 22.0 else np.nanpercentile(t, 90)))
    else:
        # Every bottom is far brighter than deep water in the red and near infrared, so the
        # pixels with the lowest (smoothed) red are deeper than the red's reach. Among them, the
        # lowest green are the deepest, where the water is clear enough for the bottom to show
        # in the green beyond the red's reach.
        cand = valid & (red <= np.percentile(red[valid], 100*a.red_fraction))
        deep = cand & (green <= np.percentile(green[cand], 100*a.green_fraction))
    r = rrs_of(R[:, deep])
    med = np.median(r, axis=1)
    noise = 1.4826*np.median(np.abs(r - med[:, None]), axis=1)
    return med, np.maximum(noise, 1e-5), int(deep.sum())


def qaa(rrs_inf, lam, tables):
    """kappa and u at each band from deep-water rrs (QAA-style, as samodel_h0): u from
    rrs = (g0 + g1 u) u; at 665 nm a ~ a_w fixes b_b; Y from QAA's 443/560 ratio."""
    u = (-G0 + np.sqrt(G0*G0 + 4*G1*rrs_inf))/(2*G1)
    lam = np.asarray(lam)
    i_ref = int(np.argmin(np.abs(lam - 665.0)))
    i443, i560 = int(np.argmin(np.abs(lam - 443.0))), int(np.argmin(np.abs(lam - 560.0)))
    t = [tables[(l, 0)] for l in lam]
    bbp_ref = max(u[i_ref]*t[i_ref]["aw"]/(1 - u[i_ref]) - t[i_ref]["bbw"], 1e-6)
    Y = np.clip(2.0*(1 - 1.2*np.exp(-0.9*rrs_inf[i443]/rrs_inf[i560])), 0.0, 2.5)
    bb = np.array([tk["bbw"] for tk in t]) + bbp_ref*(lam[i_ref]/lam)**Y
    a = (1 - u)*bb/u
    return a + bb, u


# The model.

def profile(rrs, w, ref, glint):
    """Profile chi2 over H_GRID for one pixel.

    rrs, w: (S, L) observed rrs and weights 1/sigma^2 (0 where missing). ref: per scene
    dict(rrs_inf, kappa, MC, MB, tide), and shapes {class: rho (S, L)}. Returns chi2 (nH,),
    the class (nH,) and brightness (nH,) at each depth, and chi2 of the optically deep model."""
    S, L = rrs.shape
    rinf = np.array([r["rrs_inf"] for r in ref])
    delta = np.array([r.get("delta", 0.0) for r in ref])[:, None]
    kap = np.array([r["kappa"] for r in ref])
    MC = np.array([r["MC"] for r in ref])
    MB = np.array([r["MB"] for r in ref])
    tide = np.array([r["tide"] for r in ref])
    Z = np.maximum(H_GRID[:, None] + tide[None, :], 0.0)[:, :, None]   # (nH, S, 1)
    y = rrs[None] - delta[None] - (rinf - delta)[None]*(1 - np.exp(-kap[None]*MC[None]*Z))  # (nH, S, L)
    T = np.exp(-kap[None]*MB[None]*Z)
    best = np.full(H_GRID.size, np.inf)
    cls = np.zeros(H_GRID.size, dtype=int)
    bri = np.zeros(H_GRID.size)
    # Offsets: rrs offset per scene d_s (flat), with a Gaussian prior; in rrs units ~ 2 Rrs.
    gw = 1.0/(2*GLINT_SIGMA)**2
    for ci, (name, (_, lo, hi)) in enumerate(CLASSES.items()):
        lo, hi = ref[0].get("bounds", {}).get(name, (lo, hi))
        x = T*ref[0]["shapes"][name][None]/np.pi                         # (nH, S, L)
        if not glint:
            sxx = (w*x*x).sum(axis=(1, 2))
            sxy = (w*x*y).sum(axis=(1, 2))
            b = np.clip(sxy/np.maximum(sxx, 1e-30), lo, hi)
            chi = (w*(y - b[:, None, None]*x)**2).sum(axis=(1, 2))
        else:
            # Minimise sum w (y - b x - d_s)^2 + gw sum d_s^2: for fixed b, d_s is explicit; b
            # is then a 1-D quadratic, clipped to its range (exact for one bounded variable).
            sw = w.sum(axis=1) + gw                                        # (S,)
            xs = (w[None]*x).sum(axis=2)/sw[None]                          # (nH, S)
            ys = (w[None]*y).sum(axis=2)/sw[None]
            xr = x - xs[:, :, None]
            yr = y - ys[:, :, None]
            # chi(b) = sum w (yr - b xr)^2 + gw sum (ys - b xs)^2
            sxx = (w*xr*xr).sum(axis=(1, 2)) + gw*(xs*xs).sum(axis=1)
            sxy = (w*xr*yr).sum(axis=(1, 2)) + gw*(xs*ys).sum(axis=1)
            b = np.clip(sxy/np.maximum(sxx, 1e-30), lo, hi)
            chi = (w*(yr - b[:, None, None]*xr)**2).sum(axis=(1, 2)) + \
                gw*((ys - b[:, None]*xs)**2).sum(axis=1)
        better = chi < best
        best[better], cls[better], bri[better] = chi[better], ci, b[better]
    yd = rrs - rinf
    if not glint:
        chi_deep = float((w*yd*yd).sum())
    else:
        sw = w.sum(axis=1) + gw
        d = (w*yd).sum(axis=1)/sw
        chi_deep = float((w*(yd - d[:, None])**2).sum() + gw*(d*d).sum())
    return best, cls, bri, chi_deep


def decide(chi, chi_deep, t_deep, t_amb, n):
    """Depth, sigma and why not, from the profile."""
    k = int(np.argmin(chi))
    H = H_GRID[k]
    # Parabolic refinement on the grid neighbours.
    if 0 < k < H_GRID.size - 1:
        x0, x1, x2 = H_GRID[k - 1:k + 2]
        y0, y1, y2 = chi[k - 1:k + 2]
        den = (x0 - x1)*(x0 - x2)*(x1 - x2)
        A = (x2*(y1 - y0) + x1*(y0 - y2) + x0*(y2 - y1))/den
        B = (x2*x2*(y0 - y1) + x1*x1*(y2 - y0) + x0*x0*(y1 - y2))/den
        if A > 0:
            H = float(np.clip(-B/(2*A), x0, x2))
    scale = max(1.0, chi[k]/max(n - 2, 1))
    inside = H_GRID[chi <= chi[k] + scale]
    sigma = 0.5*(inside.max() - inside.min()) if inside.size else np.nan
    sigma = max(sigma, 0.5*(H_GRID[min(k + 1, H_GRID.size - 1)] - H_GRID[max(k - 1, 0)]))
    if chi_deep - chi[k] <= t_deep*scale:
        return np.nan, sigma, "deep"
    if t_amb > 0:
        far = np.abs(H_GRID - H) > 0.5 + 0.1*H
        if far.any() and chi[far].min() - chi[k] <= t_amb*scale:
            return np.nan, sigma, "ambiguous"
    return H, sigma, "ok"


# Image-calibrated attenuation (EXPERIMENTS.md, "The depth ceiling is the attenuation"; after 4SM,
# Morel & Favoretto 2017). QAA from the absolute deep-water rrs turns any offset left by the
# atmospheric correction into backscatter, and so into attenuation at every band. The bottom
# contrast c = rrs - rrs_inf is free of a flat offset: for one bright bottom at varying depth,
# ln c_i against ln c_j is a line of slope K_i/K_j (Lyzenga; 4SM's brightest pixels line).

# Jerlov's Kd (m^-1; Austin & Petzold's update for the oceanic types), as model/jerlov.c.
JERLOV_WL = np.array([400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700.0])
JERLOV = np.array([
    [0.0217, 0.0185, 0.0176, 0.0184, 0.0280, 0.0504, 0.0640, 0.0931, 0.2408, 0.3174, 0.3559, 0.4372, 0.6513],
    [0.0316, 0.0280, 0.0257, 0.0250, 0.0332, 0.0545, 0.0674, 0.0960, 0.2437, 0.3206, 0.3601, 0.4410, 0.6530],
    [0.0438, 0.0395, 0.0355, 0.0330, 0.0396, 0.0596, 0.0715, 0.0995, 0.2471, 0.3245, 0.3652, 0.4457, 0.6550],
    [0.0878, 0.0814, 0.0714, 0.0620, 0.0627, 0.0779, 0.0863, 0.1122, 0.2595, 0.3389, 0.3837, 0.4626, 0.6623],
    [0.1697, 0.1594, 0.1381, 0.1160, 0.1056, 0.1120, 0.1139, 0.1359, 0.2826, 0.3655, 0.4181, 0.4942, 0.6760],
    [0.2516, 0.2374, 0.2048, 0.1700, 0.1486, 0.1461, 0.1415, 0.1596, 0.3057, 0.3922, 0.4525, 0.5257, 0.6896],
    [0.78, 0.54, 0.39, 0.29, 0.22, 0.2, 0.19, 0.21, 0.33, 0.4, 0.46, 0.56, 0.71],
    [1.1, 0.78, 0.56, 0.43, 0.36, 0.31, 0.3, 0.33, 0.4, 0.48, 0.54, 0.65, 0.8],
    [1.6, 1.2, 0.89, 0.71, 0.58, 0.49, 0.46, 0.46, 0.48, 0.54, 0.63, 0.78, 0.92],
    [2.4, 1.9, 1.6, 1.23, 0.99, 0.78, 0.63, 0.58, 0.6, 0.65, 0.76, 0.92, 1.1]])  # I IA IB II III 1C 3C 5C 7C 9C


def bpl_ratio(rrs, rinf, i, j, q=95.0, nbins=14, min_count=50):
    """K_i/K_j from the brightest pixels line: the upper envelope (q-th percentile of ln c_i in
    bins of ln c_j) of the bottom contrast, fitted by a line. None if there are too few bins."""
    ci, cj = rrs[i] - rinf[i], rrs[j] - rinf[j]
    m = np.isfinite(ci) & np.isfinite(cj) & (ci > 0) & (cj > 0)
    if m.sum() < nbins*min_count:
        return None
    Xi, Xj = np.log(ci[m]), np.log(cj[m])
    edges = np.percentile(Xj, np.linspace(5, 95, nbins + 1))
    env = [(0.5*(e0 + e1), np.percentile(Xi[(Xj >= e0) & (Xj < e1)], q))
           for e0, e1 in zip(edges[:-1], edges[1:]) if ((Xj >= e0) & (Xj < e1)).sum() >= min_count]
    if len(env) < 4:
        return None
    env = np.array(env)
    return float(np.polyfit(env[:, 0], env[:, 1], 1)[0])


def jerlov_kd(lam, t):
    """Kd at the bands for a continuous Jerlov type index t (0 = I ... 9 = 9C)."""
    k = np.array([np.interp(lam, JERLOV_WL, row) for row in JERLOV])  # types x bands
    t = float(np.clip(t, 0, len(JERLOV) - 1))
    i = min(int(t), len(JERLOV) - 2)
    f = t - i
    return np.exp((1 - f)*np.log(k[i]) + f*np.log(k[i + 1]))


def calibrate(s, rinf, lam, tables, a, region=None):
    """Attenuation calibration from the scene's bottom contrast, over the cells in region (default
    the whole site), with that region's rrs_inf. Returns dict(mode, ratios, delta | t), or None."""
    if a.kappa == "qaa":
        return None
    lam = np.asarray(lam)
    ig = int(np.argmin(np.abs(lam - 560.0)))
    blue = [int(np.argmin(np.abs(lam - l))) for l in (443.0, 492.0)]
    rrs = rrs_of(s["R"])
    if region is not None:
        rrs = np.where(region[None], rrs, np.nan)
    ratios = [bpl_ratio(rrs, rinf, i, ig) for i in blue]
    if any(r is None or not (0.05 < r < 3.0) for r in ratios):
        return None
    target = np.log(ratios)
    # The flat offset delta that makes QAA's kappa ratios match the image's (used for kappa in
    # mode qaa_offset, and by the offset-aware forward model in both modes).
    best = (np.inf, 0.0)
    for d in np.linspace(0.0, 0.999*rinf.min(), 200):
        k, _ = qaa(rinf - d, lam, tables)
        c = ((np.log(k[blue]/k[ig]) - target)**2).sum()
        if c < best[0]:
            best = (c, d)
    if a.kappa == "jerlov":
        ts = np.linspace(0, len(JERLOV) - 1, 901)
        cost = [((np.log(jerlov_kd(lam[blue], t)/jerlov_kd(lam[[ig]], t)) - target)**2).sum() for t in ts]
        return dict(mode="jerlov", ratios=ratios, delta=best[1], t=float(ts[int(np.argmin(cost))]))
    return dict(mode="qaa_offset", ratios=ratios, delta=best[1])


def bright_bounds(scenes, ref, lam, tables, a):
    """Bottom brightness ranges from the image (after 4SM's brightest bottom at null depth): near
    the surface the bottom term is b rho/pi, so the brightest bottoms have b ~ pi (rrs - rrs_inf)/rho
    at 560 nm. The cap is 1.2 x the largest scene's 99.5th percentile of that. "clip": each class's
    upper bound (and lower, if above) at most the cap; "scale": every range scaled by cap/0.75
    (if below 1). Default: the class ranges (CLASSES)."""
    mode = getattr(a, "bright", "classes")
    if mode == "classes":
        return {}
    lam = np.asarray(lam)
    ig = int(np.argmin(np.abs(lam - 560.0)))
    rho = tables[(lam[ig], CLASSES["sand"][0])]["rho"]
    caps = []
    for s, r in zip(scenes, ref):
        b = np.pi*(rrs_of(s["R"][ig]) - r["rrs_inf"][ig])/rho
        b = b[np.isfinite(b)]
        if b.size:
            caps.append(np.percentile(b, 99.5))
    if not caps:
        return {}
    cap = 1.2*max(caps)
    if mode == "clip":
        return {n: (min(lo, cap), min(hi, cap)) for n, (_, lo, hi) in CLASSES.items()}
    f = min(cap/CLASSES["sand"][2], 1.0)
    return {n: (lo*f, hi*f) for n, (_, lo, hi) in CLASSES.items()}


def scene_reference(s, lam, tables, truth, a, fields, region=None):
    rinf, noise, n_deep = deep_reference(s, truth, a, fields, region)
    cal = s.get("cal")
    if region is not None and getattr(a, "kappa_local", False):
        cal = calibrate(s, rinf, lam, tables, a, region) or cal
    if cal and cal["mode"] == "qaa_offset":
        kap, u = qaa(np.maximum(rinf - cal["delta"], 1e-6), lam, tables)
    else:
        kap, u = qaa(rinf, lam, tables)
        if cal and cal["mode"] == "jerlov":
            kap = jerlov_kd(np.asarray(lam), cal["t"])
    sec_w = 1.0/np.cos(np.radians(s["theta_w"]))
    sec_v = 1.0/np.cos(np.radians(s["theta_v"]))
    MC = sec_w + 1.03*np.sqrt(1 + 2.4*u)*sec_v
    MB = sec_w + 1.04*np.sqrt(1 + 5.4*u)*sec_v
    # Offset-aware model: the offset delta in rrs_inf is not water-column light, so it is not
    # attenuated: rrs = delta + (rrs_inf - delta)(1 - exp(-kappa M_C Z)) + bottom.
    delta = cal["delta"] if cal and getattr(a, "offset_model", False) else 0.0
    return dict(rrs_inf=rinf, noise=noise, kappa=kap, MC=MC, MB=MB, tide=s["tide"], n_deep=n_deep, delta=delta)


def run_case(case_dir, truth, pts, tables_for, a):
    scenes = read_case(case_dir)
    lam = scenes[0]["lam"]
    tables = tables_for(tuple(lam))
    shapes = {name: np.repeat(np.array([tables[(l, idx)]["rho"] for l in lam])[None], len(scenes), axis=0)
              for name, (idx, _, _) in CLASSES.items()}
    fields = [smoothed(s, a.smooth) for s in scenes]
    for s, f in zip(scenes, fields):
        s["cal"] = None
        s["cal"] = calibrate(s, scene_reference(s, lam, tables, truth, a, f)["rrs_inf"], lam, tables, a)
    ref = [scene_reference(s, lam, tables, truth, a, f) for s, f in zip(scenes, fields)]
    bounds = bright_bounds(scenes, ref, lam, tables, a)
    info = [dict(scene=s["name"], n_deep=r["n_deep"], rrs_inf=r["rrs_inf"].round(5).tolist(), cal=s["cal"],
                 noise=r["noise"].round(6).tolist(), kappa=r["kappa"].round(3).tolist())
            for s, r in zip(scenes, ref)]
    ny, nx = truth.shape
    cell_km = float(open(os.path.join(case_dir, "s1_b443.asc")).readlines()[4].split()[1])*111.0
    rad = int(round(a.radius/cell_km)) if a.radius > 0 else 0
    depth, sigma, why = [], [], []
    for i, j in pts:
        here = ref
        if rad:
            # The water within the radius: rrs_inf from the deep pixels near the point.
            region = np.zeros((ny, nx), dtype=bool)
            region[max(i - rad, 0):i + rad + 1, max(j - rad, 0):j + rad + 1] = True
            here = []
            for s, f, r in zip(scenes, fields, ref):
                if (f[0] & region).sum() >= a.min_local:
                    here.append(scene_reference(s, lam, tables, truth, a, f, region))
                else:
                    here.append(r)
        for r in here:
            r["shapes"] = shapes
            r["bounds"] = bounds
        rrs = np.array([rrs_of(s["R"][:, i, j]) for s in scenes])
        ok = np.isfinite(rrs) & (rrs > 0)
        sig2 = np.array([r["noise"] for r in here])**2 + (a.eps*np.where(ok, rrs, 0.0))**2
        w = np.where(ok, 1.0/sig2, 0.0)
        rrs = np.where(ok, rrs, 0.0)
        if ok.sum() < 3:
            depth.append(np.nan), sigma.append(np.nan), why.append("no data")
            continue
        chi, cls, bri, chi_deep = profile(rrs, w, here, a.glint)
        H, sg, reason = decide(chi, chi_deep, a.t_deep, a.t_amb, int(ok.sum()))
        depth.append(H), sigma.append(sg), why.append(reason)
    return np.array(depth), np.array(sigma), why, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", choices=sorted(ab.SPLITS), default="train")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--against", default="HEAD", help="the photic variant to compare with (ab.py syntax)")
    ap.add_argument("--sites", default="sites", help="tier 4 site set for the profile model (sites5: 705 nm)")
    ap.add_argument("--ref", choices=["red", "oracle"], default="red")
    ap.add_argument("--eps", type=float, default=0.03, help="relative model error")
    ap.add_argument("--t-deep", type=float, default=9.0)
    ap.add_argument("--t-amb", type=float, default=0.0)
    ap.add_argument("--glint", action="store_true")
    ap.add_argument("--kappa", choices=["qaa", "qaa_offset", "jerlov"], default="qaa",
                    help="attenuation: QAA from rrs_inf, QAA after a per-scene offset matched to the image's "
                         "brightest pixels line, or Jerlov's Kd for the type matching it (4SM)")
    ap.add_argument("--bright", choices=["classes", "clip", "scale"], default="classes",
                    help="bottom brightness ranges: the in situ classes, or capped from the image's brightest bottoms")
    ap.add_argument("--offset-model", action="store_true",
                    help="keep the calibration's offset outside the attenuated water-column term")
    ap.add_argument("--kappa-local", action="store_true",
                    help="calibrate the attenuation within each local window (fallback: the site-wide calibration)")
    ap.add_argument("--smooth", type=int, default=7, help="smoothing window (pixels) for ranking deep pixels")
    ap.add_argument("--red-fraction", type=float, default=0.2, help="share of pixels kept by the red ranking")
    ap.add_argument("--green-fraction", type=float, default=0.15, help="share of those kept by the green ranking")
    ap.add_argument("--radius", type=float, default=0.0, help="local water within this radius (km); 0: whole site")
    ap.add_argument("--min-local", type=int, default=400, help="fewest valid pixels for a local reference")
    ap.add_argument("--min-effect", type=float, default=0.02)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--tier", type=int, choices=[3, 4], help="only this tier")
    a = ap.parse_args()

    split = ab.SPLITS[a.split]
    cases = ab.tier3_cases(split["water"]) + ab.tier4_cases(split["sites"])
    points = ab.sample(cases, a.n, a.seed)
    if a.tier:
        cases = [c for c in cases if c["tier"] == a.tier]
    tag = "%s_n%d_s%d" % (a.split, a.n, a.seed)

    # photic, as ab.py runs it (cached).
    bspec, sites_set, scenes, start, model = ab.parse_variant(a.against)
    key, photic = ab.build(bspec)
    key += ("%" + sites_set if sites_set else "") + ("#" + scenes.replace(",", "+") if scenes else "") + \
        ("@" + start if start else "") + ("~" + model if model else "")
    alt = {c["site"]: c for c in ab.tier4_cases(split["sites"], sites_set)} if sites_set else {}
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(a.jobs) as ex:
        base = list(ex.map(lambda c: ab.run_case(photic, alt.get(c.get("site"), c) if c["tier"] == 4 else c,
                                                  points[c["name"]], os.path.join(ab.AB, "runs", key, tag, c["name"]),
                                                  False, start, scenes, model), cases))
    base = {c["name"]: o for c, o in zip(cases, base)}

    # The profile model.
    pm_alt = {c["site"]: c for c in ab.tier4_cases(split["sites"], a.sites)} if a.sites != "sites" else {}
    cache = {}

    def tables_for(lam):
        if lam not in cache:
            cache[lam] = t2.photic_tables(ab.DRIVER, list(lam))
        return cache[lam]

    t0 = time.time()
    pm, infos, reasons = {}, {}, {}
    for c in cases:
        cc = pm_alt.get(c.get("site"), c) if c["tier"] == 4 else c
        truth = np.load(os.path.join(cc["dir"], "truth.npy"))
        d, s, why, info = run_case(cc["dir"], truth, points[c["name"]], tables_for, a)
        pm[c["name"]] = dict(depth=d, sigma=s, misfit=np.full(d.size, np.nan))
        infos[c["name"]], reasons[c["name"]] = info, {k: why.count(k) for k in set(why)}
    seconds = time.time() - t0

    rng = np.random.default_rng(0)
    groups = {}
    sel = {c["name"]: [c["name"]] for c in cases}
    sel.update({"tier 3": [c["name"] for c in cases if c["tier"] == 3],
                "tier 4": [c["name"] for c in cases if c["tier"] == 4],
                "all": [c["name"] for c in cases]})
    sel = {k: v for k, v in sel.items() if v}
    for g, names in sel.items():
        cat = lambda d, k: np.concatenate([d[nm][k] for nm in names])
        groups[g] = ab.paired(cat(base, "truth"), cat(base, "depth"), cat(pm, "depth"),
                              cat(base, "sigma"), cat(pm, "sigma"), rng)
    label = "profile(ref %s, eps %g, t_deep %g, t_amb %g%s%s%s%s)" % (
        a.ref, a.eps, a.t_deep, a.t_amb, "" if a.kappa == "qaa" else ", " + a.kappa + (" local" if a.kappa_local else "") + (" offset model" if a.offset_model else "") + ("" if a.bright == "classes" else ", bright " + a.bright), ", glint" if a.glint else "",
        "" if a.sites == "sites" else ", " + a.sites,
        ", r %g km" % a.radius if a.radius else "")
    res = dict(A=a.against, B=label, split=a.split, seed=a.seed, min_effect=a.min_effect,
               n_points=sum(len(points[c["name"]]) for c in cases),
               date=datetime.datetime.now().isoformat(timespec="seconds"),
               seconds_A=sum(base[c["name"]]["seconds"] for c in cases), seconds_B=seconds,
               groups=groups, reasons=reasons, references=infos, args=vars(a))
    os.makedirs(os.path.join(ab.AB, "results"), exist_ok=True)
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", "profile_%s_vs_%s_%s" % (key, label, tag))
    path = os.path.join(ab.AB, "results", name + ".json")
    json.dump(res, open(path, "w"), indent=1, default=float)
    print(ab.report(res))
    print("\nabstentions (profile model): " + "; ".join("%s %s" % (k, v) for k, v in reasons.items()))
    print("results: %s" % path)


if __name__ == "__main__":
    main()
