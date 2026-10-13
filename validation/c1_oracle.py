"""Phase C1 of EXPERIMENTS.md: can a better optimiser recover depth from photic's own
reflectance misfit?

The per-pixel problem is photic's, at the same tier 3 pixels as phaseA.py: a 3x3 window
(NSPATIAL 2) with three bottom types, fitted jointly over the scenes. The forward model is
photic_ref.forward, which tier 1 shows matches samodel_Rrs exactly. The misfit is
photic's Rrs error: the RMS difference over regions, scenes and bands, relative to the
mean measured Rrs. The K prior and the coherence terms are left out. The optimiser is SciPy's bounded
least squares (trust-region reflective), on the same search as B5: a profile over a
grid of common depths, then a refinement of all parameters from the two best.

  c1_oracle.py [--water moderate] [--osoaa] [--fix-depth]

--fix-depth fits everything but depth, at the true depths: a test of the optimiser alone.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

import phaseA as A
import photic_ref
import tier2_osoaa as t2
import tier3_scenes as t3

NB = 3
H_GRID = [0.5, 1.0, 2.0, 3.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0]
BOUNDS = dict(H=(0.05, 60.0), P=(0.0005, 1.0), G=(0.0005, 2.0), X=(0.0001, 0.5))
TABLES = None


def interp_photic(wl, spec, x):
    """photic's interp_1d: linear, extrapolating from the end segments."""
    wl = np.asarray(wl)
    k = np.clip(np.searchsorted(wl, x) - 1, 0, len(wl) - 2)
    return spec[..., k] + (spec[..., k + 1] - spec[..., k])*(x - wl[k])/(wl[k + 1] - wl[k])


class Pixel:
    def __init__(self, water, n, i, j, H, bottom, grids):
        cells, self.origin = A.window(i, j, 2)
        self.nr, self.n = len(cells), n
        self.lam = np.array(t3.BANDS)
        # Measured Rrs [region, scene, band].
        self.R = np.array([[[grids[k][m][a, b] for m in range(len(t3.BANDS))] for k in range(n)]
                           for a, b in cells])
        self.R440 = interp_photic(self.lam, self.R, 440.0)
        self.R490 = interp_photic(self.lam, self.R, 490.0)
        self.R550 = interp_photic(self.lam, self.R, 550.0)
        self.R640 = interp_photic(self.lam, self.R, 640.0)
        self.R440 = np.where(self.R440 < 0, 0.0001, self.R440)
        sets = t3.SCENE_SETS[water]
        self.tide = np.array([sets[k][5] for k in range(n)])
        self.theta_sun = np.array([t3.theta_w(sets[k][4]) for k in range(n)])
        t = TABLES
        self.rho = np.array([[t[(l, b)]["rho"] for l in t3.BANDS] for b in range(NB)])  # [bottom, band]
        self.a0 = np.array([t[(l, 0)]["a0"] for l in t3.BANDS])
        self.a1 = np.array([t[(l, 0)]["a1"] for l in t3.BANDS])
        self.aw = np.array([t[(l, 0)]["aw"] for l in t3.BANDS])
        self.bbw = np.array([t[(l, 0)]["bbw"] for l in t3.BANDS])
        self.H_true = np.array([H[c] for c in cells])
        self.mean = self.R.mean()
        self.scale = 100.0/self.mean/np.sqrt(self.R.size)

    def unpack(self, x):
        nr, n = self.nr, self.n
        H = x[:nr]
        alb = x[nr:2*nr]
        z = np.concatenate([np.zeros((nr, 1)), x[2*nr:2*nr + nr*(NB - 1)].reshape(nr, NB - 1)], axis=1)
        q = np.exp(z - z.max(axis=1, keepdims=True))
        q /= q.sum(axis=1, keepdims=True)
        w = np.exp(x[2*nr + nr*(NB - 1):].reshape(n, 3))
        return H, alb, q, w[:, 0], w[:, 1], w[:, 2]

    def residuals(self, x):
        H, alb, q, P, G, X = self.unpack(x)
        rho_tab = q @ self.rho  # [region, band]
        out = photic_ref.forward(self.lam[None, None, :], P[None, :, None], G[None, :, None], X[None, :, None],
                                 alb[:, None, None], rho_tab[:, None, :], H[:, None, None],
                                 self.tide[None, :, None], self.theta_sun[None, :, None], 0.0,
                                 self.R440[:, :, None], self.R490[:, :, None], 0.0,
                                 self.a0, self.a1, self.aw, self.bbw)["Rrs"]
        return ((out - self.R)*self.scale).ravel()

    def start(self, h):
        nr, n = self.nr, self.n
        m490, m550, m640 = self.R490.mean(0), self.R550.mean(0), self.R640.mean(0)
        P = np.clip(0.072*(m490/m550)**-1.7, *BOUNDS["P"])
        G = np.clip(1.5*P, *BOUNDS["G"])
        aw640 = np.interp(640.0, [630.0, 640.0, 650.0], [0.319, 0.329, 0.349])  # Pope & Fry
        X = np.clip(30.0*aw640*m640, *BOUNDS["X"])
        pgx = np.log(np.column_stack([P, G, X])).ravel()
        return np.concatenate([np.full(nr, h), np.full(nr, 0.3), np.zeros(nr*(NB - 1)), pgx])

    def bounds(self):
        nr, n = self.nr, self.n
        lo = [BOUNDS["H"][0]]*nr + [0.0]*nr + [-20.0]*(nr*(NB - 1))
        hi = [BOUNDS["H"][1]]*nr + [1.0]*nr + [20.0]*(nr*(NB - 1))
        for _ in range(n):
            lo += [np.log(BOUNDS[k][0]) for k in "PGX"]
            hi += [np.log(BOUNDS[k][1]) for k in "PGX"]
        return np.array(lo), np.array(hi)

    def fit(self, x0, free):
        lo, hi = self.bounds()
        x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

        def f(y):
            x = x0.copy()
            x[free] = y
            return self.residuals(x)
        r = least_squares(f, x0[free], bounds=(lo[free], hi[free]), method="trf", x_scale="jac",
                          max_nfev=2000)
        x = x0.copy()
        x[free] = r.x
        return x, float(np.sqrt(np.sum(r.fun**2))), r.nfev

    def depth(self, x):
        w = np.ones(self.nr)
        w[self.origin] = np.sqrt(self.nr)
        return float(np.sum(w*x[:self.nr])/w.sum())


def solve(args):
    water, n, i, j, H, bottom, grids, fix_depth = args
    p = Pixel(water, n, i, j, H, bottom, grids)
    allfree = np.ones(len(p.start(1.0)), dtype=bool)
    notH = allfree.copy()
    notH[:p.nr] = False
    nfev = 0
    if fix_depth:
        x0 = p.start(1.0)
        x0[:p.nr] = p.H_true
        x, err, k = p.fit(x0, notH)
        return dict(H_true=float(H[i, j]), depth=p.depth(x), error=err, nfev=k)
    prof = []
    for h in H_GRID:
        x, err, k = p.fit(p.start(h), notH)
        nfev += k
        prof.append((err, h, x))
    prof.sort(key=lambda t: t[0])
    best = None
    for err0, h, x0 in prof[:2]:
        x, err, k = p.fit(x0, allfree)
        nfev += k
        if best is None or err < best[0]:
            best = (err, x)
    return dict(H_true=float(H[i, j]), depth=p.depth(best[1]), error=best[0], nfev=nfev,
                albedo=float(best[1][p.nr + p.origin]))


def init(tables):
    global TABLES
    TABLES = tables


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--water", default="moderate")
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("--osoaa", action="store_true")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--fix-depth", action="store_true")
    ap.add_argument("--rows", default="0,3,6,9,12,15,18,21,24,27,31")
    ap.add_argument("--cols", default="8,20,28")
    a = ap.parse_args()
    tables = t2.photic_tables(A.DRIVER, t3.BANDS)
    H, bottom, grids = A.scene_grids(a.water, a.scenes, a.osoaa, a.noise)
    pix = [(a.water, a.scenes, i, j, H, bottom, grids, a.fix_depth)
           for i in map(int, a.rows.split(",")) for j in map(int, a.cols.split(","))]
    with ProcessPoolExecutor(max_workers=12, initializer=init, initargs=(tables,)) as ex:
        rows = list(ex.map(solve, pix))
    e = np.array([r["depth"] - r["H_true"] for r in rows])
    t = np.array([r["H_true"] for r in rows])
    err = np.array([r["error"] for r in rows])
    name = "%s_n%d_%s%s" % (a.water, a.scenes, "osoaa" if a.osoaa else "control", "_fixdepth" if a.fix_depth else "")
    print("%-28s bias %6.2f rmse %6.2f tol %.2f | Rrs error median %.3f%% | evaluations/pixel %d" % (
        name, e.mean(), np.sqrt(np.mean(e**2)), np.mean(np.abs(e) <= 0.5 + 0.1*t), np.median(err),
        np.mean([r["nfev"] for r in rows])))
    for r in rows:
        print("  H %5.2f  fit %6.2f  err %.3f" % (r["H_true"], r["depth"], r["error"]))


if __name__ == "__main__":
    main()
