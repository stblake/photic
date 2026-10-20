"""Iterative (alternating) inversion prototype: depth, water and bottom improved in turn.

A stand-in for the photic executable: `iterative.py run.bam` reads the command file's grids and
SCENE lines in the working directory and writes modelled_H.nc, modelled_H_sigma.nc,
modelled_albedo.nc and modelled_min_K.nc, so bench.py and ab.py run and score it unchanged
(bench.py run --photic validation/iterative_photic). The forward model is Lee's, as photic's:

    rrs_spl = rinf_sl [1 - exp(-kappa_sl MC_sl Z_ps)] + b_p rho_c(p),l/pi exp(-kappa_sl MB_sl Z_ps),

Z_ps = H_p + tide_s. Each pass:

  1. Pixels (water and bottom colours fixed): for every pixel, the exact profile chi2 over a depth
     grid, with the bottom class and brightness b (bounded, shared by the scenes) solved in closed
     form. One H and one bottom per pixel for every band and scene; the global minimum by
     construction. The optically deep test and sigma as the profile model (profile_model.py).
  2. Water (pixels fixed), per scene, pooled over the confident pixels: kappa from photic's IOPs,
     a = a_w + a_phi(P) + a_g(G), b_b = b_bw + X (440/lambda)^Y, with P, G, X fitted and Y fixed;
     rinf_sl free per band (it absorbs an atmospheric offset; profiled out in closed form). The
     known pure-water absorption anchors the depth scale.
  3. Bottom colours (pixels and water fixed; PHOTIC_ITER_BOTTOM=learn): each class's spectrum from
     its pixels, normalised to the library's mean brightness and damped towards it.

Pass 0 with QAA water is the profile model (Lee_Profile, site-wide reference, no offset).

Environment (bench.py passes it through):
  PHOTIC_ITER_N        passes after the first (default 4; 0: the profile model)
  PHOTIC_ITER_WATER    pooled (default) | fixed
  PHOTIC_ITER_BOTTOM   library (default) | learn
  PHOTIC_ITER_EPS      relative model error (default 0.03)
  PHOTIC_ITER_TDEEP    optically deep threshold (default 9)
  PHOTIC_ITER_BMAX     upper bound of sand's brightness, if above the class's 0.75
  PHOTIC_ITER_POOL     all (default: every pixel, as the objective) | confident: the water from those only
  PHOTIC_ITER_CONF     interior (default: pixels with b inside its bounds) | all: pixels pooled for the water
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ab  # noqa: E402
import photic_ref  # noqa: E402
import tier2_osoaa as t2  # noqa: E402

C = photic_ref.LEE_1999
CLASSES = {"sand": (0, 0.15, 0.75), "seagrass": (1, 0.03, 0.15), "coral": (2, 0.03, 0.25)}
H_GRID = np.concatenate([np.linspace(0.05, 2.0, 40), np.linspace(2.05, 10.0, 80), np.linspace(10.2, 30.0, 100)])
CHUNK = 8000


def env(key, default, cast=str):
    return cast(os.environ.get(key, default))


def log(*a):
    print(*a, flush=True)


# Inputs and outputs.

def read_asc(path):
    with open(path) as f:
        hdr = {}
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
        a = np.loadtxt(f)
    a = a[::-1].astype(float)  # photic's row 0 is south
    a[a == hdr["nodata_value"]] = np.nan
    return a, hdr


def read_bam(path):
    files, scenes, hdr = {}, [], None
    for line in open(path):
        tok = line.split()
        if line.startswith("READ GRID"):
            files[tok[2]] = tok[3]
        elif line.startswith("SCENE NAME"):
            bands = tok[tok.index("BANDS") + 1:tok.index("WAVELENGTHS")]
            j = tok.index("WAVELENGTHS")
            lam = [float(x) for x in tok[j + 1:j + 1 + len(bands)]]
            get = lambda key, d=0.0: float(tok[tok.index(key) + 1]) if key in tok else d
            R = []
            for b in bands:
                g, hdr = read_asc(files[b])
                R.append(g)
            scenes.append(dict(name=tok[2], lam=lam, theta_w=get("THETASUN"), theta_v=get("THETAVIEW"),
                               tide=get("HTIDE"), R=np.array(R)))
    return scenes, hdr


def write_nc(path, grid, hdr):
    import netCDF4
    ny, nx = grid.shape
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lons", nx)
        ds.createDimension("lats", ny)
        lo = ds.createVariable("lons", "f8", ("lons",))
        la = ds.createVariable("lats", "f8", ("lats",))
        v = ds.createVariable("bathy", "f4", ("lats", "lons"), fill_value=np.float32(-99999.0))
        lo[:] = hdr["xllcorner"] + hdr["cellsize"]*np.arange(nx)
        la[:] = hdr["yllcorner"] + hdr["cellsize"]*np.arange(ny)
        v[:] = np.ma.masked_invalid(grid)


# Water.

def rrs_of(R):
    return R/(C["ra"] + C["rb"]*R)


def deep_cells(scene, red_fraction=0.2, green_fraction=0.15, size=7):
    """The optically deep cells (rinf.py's rule): darkest smoothed red, then darkest green."""
    from scipy import ndimage
    R = scene["R"]
    lam = np.array(scene["lam"])
    valid = np.all(np.isfinite(R) & (R > 0), axis=0)
    n = ndimage.uniform_filter(valid.astype(float), size, mode="constant")
    sm = lambda f: ndimage.uniform_filter(np.where(valid, f, 0.0), size, mode="constant")/np.maximum(n, 1e-9)
    red = sm(R[lam >= 660.0].sum(axis=0))
    green = sm(R[int(np.argmin(np.abs(lam - 560.0)))])
    cand = valid & (red <= np.percentile(red[valid], 100*red_fraction))
    return cand & (green <= np.percentile(green[cand], 100*green_fraction))


class Water:
    """One scene's water: IOPs (P, G, X; Y fixed) give kappa, u and the path factors; rinf is free."""

    def __init__(self, scene, tables, rinf, noise):
        self.lam = np.array(scene["lam"])
        self.t = [tables[(l, 0)] for l in self.lam]
        self.aw = np.array([t["aw"] for t in self.t])
        self.bbw = np.array([t["bbw"] for t in self.t])
        self.a0 = np.array([t["a0"] for t in self.t])
        self.a1 = np.array([t["a1"] for t in self.t])
        self.sec_w = 1.0/np.cos(np.radians(scene["theta_w"]))
        self.sec_v = 1.0/np.cos(np.radians(scene["theta_v"]))
        self.tide = scene["tide"]
        self.rinf = rinf.copy()
        self.noise = noise
        # Y from photic's ratio on the deep water (R440/R490); then P, G, X matched to QAA's kappa.
        i443, i490 = int(np.argmin(np.abs(self.lam - 443))), int(np.argmin(np.abs(self.lam - 490)))
        chi = rinf[i443]/rinf[i490]
        self.Y = float(np.clip(3.44*(1.0 - 3.17*np.exp(-2.01*chi)), 0.0, 2.5))
        self.theta = np.log([0.02, 0.02, 0.005])
        k_qaa, _ = qaa(rinf, self.lam, tables)
        from scipy.optimize import least_squares
        r = least_squares(lambda th: np.log(self.kappa_u(th)[0]) - np.log(k_qaa), self.theta,
                          bounds=(np.log([1e-4, 1e-4, 1e-5]), np.log([1.0, 2.0, 0.5])))
        self.theta = r.x
        self.update()

    def kappa_u(self, theta):
        P, G, X = np.exp(theta)
        a = self.aw + (self.a0 + self.a1*np.log(P))*P + G*np.exp(-C["S"]*(self.lam - 440.0))
        bb = self.bbw + X*(440.0/self.lam)**self.Y
        return a + bb, bb/(a + bb)

    def paths(self, u):
        MC = self.sec_w + C["c0"]*np.sqrt(1.0 + C["c1"]*u)*self.sec_v
        MB = self.sec_w + C["b0"]*np.sqrt(1.0 + C["b1"]*u)*self.sec_v
        return MC, MB

    def update(self):
        self.kappa, self.u = self.kappa_u(self.theta)
        self.MC, self.MB = self.paths(self.u)


def qaa(rrs_inf, lam, tables):
    """As profile_model.qaa: kappa and u at each band from deep-water rrs."""
    G0, G1 = C["g0"], C["g1"]
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


# Step 1: pixels.

def pixel_step(rrs, w, waters, shapes, bounds, t_deep):
    """rrs, w: (P, S, L). Returns H (nan where not given), sigma, class, brightness, chi2 at the
    optimum, chi2 of the optically deep model and the reason."""
    P, S, L = rrs.shape
    kap = np.array([wt.kappa for wt in waters])
    MC = np.array([wt.MC for wt in waters])
    MB = np.array([wt.MB for wt in waters])
    rinf = np.array([wt.rinf for wt in waters])
    tide = np.array([wt.tide for wt in waters])
    Z = np.maximum(H_GRID[:, None] + tide[None, :], 0.0)[:, :, None]          # (nH, S, 1)
    Wc = (rinf[None]*(1 - np.exp(-kap[None]*MC[None]*Z))).reshape(H_GRID.size, S*L)
    T = np.exp(-kap[None]*MB[None]*Z)
    nH = H_GRID.size
    best = np.full((P, nH), np.inf)
    cls = np.zeros((P, nH), dtype=np.int8)
    bri = np.zeros((P, nH))
    rf = rrs.reshape(P, S*L)
    wf = w.reshape(P, S*L)
    A = wf*rf
    syy0 = (A*rf).sum(1)[:, None] - 2*A @ Wc.T + wf @ (Wc*Wc).T                # sum w (rrs - Wc)^2
    for ci, name in enumerate(CLASSES):
        lo, hi = bounds[name]
        x = (T*shapes[name][None, None, :]/np.pi).reshape(nH, S*L)
        sxy = A @ x.T - wf @ (x*Wc).T
        sxx = wf @ (x*x).T
        b = np.clip(sxy/np.maximum(sxx, 1e-30), lo, hi)
        chi = syy0 - 2*b*sxy + b*b*sxx
        m = chi < best
        best[m], cls[m], bri[m] = chi[m], ci, b[m]
    chi_deep = (wf*(rf - rinf.reshape(1, -1))**2).sum(1)
    n = (wf > 0).sum(1)
    k = np.argmin(best, axis=1)
    idx = np.arange(P)
    cmin = best[idx, k]
    H = H_GRID[k].copy()
    inner = (k > 0) & (k < nH - 1)
    kk = np.clip(k, 1, nH - 2)
    x0, x1, x2 = H_GRID[kk - 1], H_GRID[kk], H_GRID[kk + 1]
    y0, y1, y2 = best[idx, kk - 1], best[idx, kk], best[idx, kk + 1]
    den = (x0 - x1)*(x0 - x2)*(x1 - x2)
    Aq = (x2*(y1 - y0) + x1*(y0 - y2) + x0*(y2 - y1))/den
    Bq = (x2*x2*(y0 - y1) + x1*x1*(y2 - y0) + x0*x0*(y1 - y2))/den
    ok_q = inner & (Aq > 0)
    H[ok_q] = np.clip(-Bq[ok_q]/(2*Aq[ok_q]), x0[ok_q], x2[ok_q])
    scale = np.maximum(1.0, cmin/np.maximum(n - 2, 1))
    inside = best <= (cmin + scale)[:, None]
    Hin = np.where(inside, H_GRID[None, :], np.nan)
    sigma = 0.5*(np.nanmax(Hin, 1) - np.nanmin(Hin, 1))
    step = 0.5*(H_GRID[np.minimum(k + 1, nH - 1)] - H_GRID[np.maximum(k - 1, 0)])
    sigma = np.maximum(sigma, step)
    deep = chi_deep - cmin <= t_deep*scale
    nodata = n < 3
    reason = np.where(nodata, "nodata", np.where(deep, "deep", "ok"))
    Hout = np.where(reason == "ok", H, np.nan)
    return dict(H=Hout, Hraw=H, sigma=sigma, cls=cls[idx, k], b=bri[idx, k], chi=cmin, chi_deep=chi_deep,
                n=n, reason=reason)


# Step 2: water.

def water_step(wt, s, rrs, w, px, shapes, deep_rrs, deep_w, max_px=6000, rng=None, pool_mode="all"):
    """Refit one scene's P, G, X and rinf to the confident pixels (and the deep cells, which see
    rinf only)."""
    from scipy.optimize import least_squares
    names = list(CLASSES)
    # The pixels pooled: by default every pixel with data (optically deep ones see rinf only), so
    # the step lowers the same objective as the pixel step; PHOTIC_ITER_CONF=confident: only those.
    pool = px["conf"] if pool_mode == "confident" else px["reason"] != "nodata"
    sel = np.flatnonzero(pool)
    if sel.size < 50:
        return False
    if sel.size > max_px:
        sel = rng.choice(sel, max_px, replace=False)
    y = rrs[sel, s]                                   # (Q, L)
    ww = w[sel, s]
    deep_px = px["reason"][sel] != "ok"
    Z = np.where(deep_px, 1e3, np.maximum(px["Hraw"][sel] + wt.tide, 0.0))[:, None]
    bshape = np.array([shapes[names[c]] for c in px["cls"][sel]])*np.where(deep_px, 0.0, px["b"][sel])[:, None]/np.pi

    def resid(theta, want_rinf=False):
        kap, u = wt.kappa_u(theta)
        MC, MB = wt.paths(u)
        g = 1 - np.exp(-kap*MC*Z)                     # (Q, L)
        bot = bshape*np.exp(-kap*MB*Z)
        # rinf per band in closed form: y - bot = rinf g (pixels) and deep = rinf (deep cells).
        num = (ww*g*(y - bot)).sum(0) + (deep_w*deep_rrs).sum(0)
        den = (ww*g*g).sum(0) + deep_w.sum(0)
        rinf = num/np.maximum(den, 1e-30)
        if want_rinf:
            return rinf
        r = np.sqrt(ww)*(y - rinf*g - bot)
        rd = np.sqrt(deep_w)*(deep_rrs - rinf)
        return np.concatenate([r.ravel(), rd.ravel()])

    r = least_squares(resid, wt.theta, bounds=(np.log([1e-4, 1e-4, 1e-5]), np.log([1.0, 2.0, 0.5])),
                      x_scale=1.0, max_nfev=200)
    wt.theta = r.x
    wt.rinf = resid(r.x, True)
    wt.update()
    return True


# Step 3: bottom colours.

def bottom_step(rrs, w, px, waters, shapes, library, damp=0.5):
    names = list(CLASSES)
    conf = px["conf"]
    out = {}
    for ci, name in enumerate(names):
        sel = np.flatnonzero(conf & (px["cls"] == ci))
        if sel.size < 50:
            out[name] = shapes[name]
            continue
        num = np.zeros(len(library[name]))
        den = np.zeros_like(num)
        for s, wt in enumerate(waters):
            Z = np.maximum(px["Hraw"][sel] + wt.tide, 0.0)[:, None]
            y = rrs[sel, s] - wt.rinf*(1 - np.exp(-wt.kappa*wt.MC*Z))
            x = px["b"][sel, None]/np.pi*np.exp(-wt.kappa*wt.MB*Z)
            num += (w[sel, s]*x*y).sum(0)
            den += (w[sel, s]*x*x).sum(0)
        fit = num/np.maximum(den, 1e-30)
        lib = library[name]
        ratio = np.clip(fit/lib, 0.2, 5.0)
        ratio /= np.exp(np.mean(np.log(ratio)))       # brightness stays with b
        out[name] = lib*ratio**damp
    return out


# The loop.

def run(bam="run.bam"):
    t0 = time.time()
    n_pass = env("PHOTIC_ITER_N", 4, int)
    water_mode = env("PHOTIC_ITER_WATER", "pooled")
    bottom_mode = env("PHOTIC_ITER_BOTTOM", "library")
    eps = env("PHOTIC_ITER_EPS", 0.03, float)
    t_deep = env("PHOTIC_ITER_TDEEP", 9.0, float)
    conf_mode = env("PHOTIC_ITER_CONF", "interior")
    pool_mode = env("PHOTIC_ITER_POOL", "all")
    scenes, hdr = read_bam(bam)
    lam = scenes[0]["lam"]
    tables = t2.photic_tables(ab.DRIVER, list(lam))
    library = {name: np.array([tables[(l, idx)]["rho"] for l in lam]) for name, (idx, _, _) in CLASSES.items()}
    shapes = dict(library)
    bmax = env("PHOTIC_ITER_BMAX", 0.0, float)
    bounds = {name: (lo, max(hi, bmax) if name == "sand" else hi) for name, (_, lo, hi) in CLASSES.items()}
    ny, nx = scenes[0]["R"].shape[1:]
    S, L = len(scenes), len(lam)

    # Observations and the deep reference of each scene.
    rrs_g = np.array([rrs_of(s["R"]) for s in scenes])                          # (S, L, ny, nx)
    ok = np.isfinite(rrs_g) & (rrs_g > 0)
    waters, deep_obs = [], []
    for s in scenes:
        deep = deep_cells(s)
        r = rrs_of(s["R"][:, deep])
        med = np.median(r, axis=1)
        noise = np.maximum(1.4826*np.median(np.abs(r - med[:, None]), axis=1), 1e-5)
        waters.append(Water(s, tables, med, noise))
        deep_obs.append(r.T)
    valid = ok.sum(axis=(0, 1)) >= 3
    pix = np.flatnonzero(valid.ravel())
    rrs = np.where(ok, rrs_g, 0.0).reshape(S, L, -1)[:, :, pix].transpose(2, 0, 1)   # (P, S, L)
    okp = ok.reshape(S, L, -1)[:, :, pix].transpose(2, 0, 1)
    noise = np.array([wt.noise for wt in waters])
    w = np.where(okp, 1.0/(noise[None]**2 + (eps*rrs)**2), 0.0)
    rng = np.random.default_rng(0)

    def pixels():
        parts = [pixel_step(rrs[i:i + CHUNK], w[i:i + CHUNK], waters, shapes, bounds, t_deep)
                 for i in range(0, rrs.shape[0], CHUNK)]
        px = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        lo = np.array([bounds[n][0] for n in CLASSES])[px["cls"]]
        hi = np.array([bounds[n][1] for n in CLASSES])[px["cls"]]
        red = px["chi"]/np.maximum(px["n"] - 2, 1)
        okr = px["reason"] == "ok"
        lim = 4*np.median(red[okr]) if okr.any() else np.inf
        px["conf"] = okr & (px["Hraw"] > 0.3) & (px["Hraw"] < 20) & (red <= lim)
        if conf_mode == "interior":
            px["conf"] &= (px["b"] > lo*1.01) & (px["b"] < hi*0.99)
        return px

    px = pixels()
    total = lambda px: float(np.where(px["reason"] == "ok", px["chi"], px["chi_deep"]).sum())
    log("pass 0: objective %.4g, depths %d of %d, confident %d" % (total(px), np.isfinite(px["H"]).sum(),
                                                                   px["H"].size, px["conf"].sum()))
    save = lambda: [(wt.theta.copy(), wt.rinf.copy()) for wt in waters]

    def restore(st, f=1.0, new=None):
        # f = 1: st; else the point a fraction f of the way from st to new.
        for wt, (th, ri), nw in zip(waters, st, new or st):
            wt.theta = th + f*(nw[0] - th) if new else th.copy()
            wt.rinf = ri + f*(nw[1] - ri) if new else ri.copy()
            wt.update()

    for it in range(1, n_pass + 1):
        before, old_shapes, old_total = save(), shapes, total(px)
        if water_mode == "pooled":
            for s, wt in enumerate(waters):
                dw = 1.0/(wt.noise**2 + (eps*deep_obs[s])**2)
                water_step(wt, s, rrs, w, px, shapes, deep_obs[s], dw, rng=rng, pool_mode=pool_mode)
        if bottom_mode == "learn":
            shapes = bottom_step(rrs, w, px, waters, shapes, library)
        after, new_shapes = save(), shapes
        # Accept only a step that lowers the total objective; otherwise backtrack (halving).
        f, accepted = 1.0, False
        for _ in range(4):
            restore(before, f, after)
            shapes = {n: old_shapes[n]*(new_shapes[n]/old_shapes[n])**f for n in CLASSES}
            cand = pixels()
            if total(cand) < old_total:
                px, accepted = cand, True
                break
            f *= 0.5
        if not accepted:
            restore(before)
            shapes = old_shapes
            log("pass %d: no step lowers the objective (%.4g); stopping" % (it, old_total))
            break
        log("pass %d (step %.3g): objective %.4g, depths %d, confident %d; water P,G,X,Y = %s" % (
            it, f, total(px), np.isfinite(px["H"]).sum(), px["conf"].sum(),
            "; ".join("%.3g,%.3g,%.3g,%.2f" % (*np.exp(wt.theta), wt.Y) for wt in waters)))
        for wt in waters:
            log("   kappa %s rinf %s" % (np.round(wt.kappa, 3), np.round(wt.rinf, 5)))
    if bottom_mode == "learn":
        log("bottom shapes / library: %s" % "; ".join("%s %s" % (n, np.round(shapes[n]/library[n], 2)) for n in CLASSES))

    def grid(v):
        g = np.full(ny*nx, np.nan)
        g[pix] = v
        return g.reshape(ny, nx)

    write_nc("modelled_H.nc", -grid(px["H"]), hdr)
    write_nc("modelled_H_sigma.nc", grid(np.where(np.isfinite(px["H"]), px["sigma"], np.nan)), hdr)
    write_nc("modelled_albedo.nc", grid(np.where(np.isfinite(px["H"]), px["b"], np.nan)), hdr)
    kmin = min(float(np.min(wt.kappa)) for wt in waters)
    write_nc("modelled_min_K.nc", grid(np.where(np.isfinite(px["H"]), kmin, np.nan)), hdr)
    log("done in %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "run.bam")
