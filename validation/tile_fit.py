"""Decoupling, step 3: water shared over a tile, fitted with its pixels (EXPERIMENTS.md).

Within one image the water varies slowly in space, while depth and bottom vary pixel to pixel;
many pixels of different depths under one water separate the water from the depth. For each
tile of TILE x TILE pixels and each scene, one water (ln P, ln G, ln X, Y) is shared by the
tile's pixels, and the fit alternates:

  pixel step   water fixed: each pixel fits its depth and bottom alone (photic's profile over
               common depths, then a refinement; LM through invert_lab, one region per pixel)
  water step   depths and bottoms fixed: the tile's water is fitted to all its pixels at once
               (invert_lab with the tile's pixels as regions, only the water free)

starting from photic's starting water (no extra information) or from the water fitted to each
scene's deep-water Rrs (as decouple.py). Scored on the Phase A pixels, by depth.

Needs an invert_lab built with Y_MODEL=3 and BOTTOM_BOUNDS=1 ($PHOTIC_INVERT_LAB).

  tile_fit.py [--realistic] [--noise 0.0005] [--init start|deep] [--rounds 3] [water ...]
"""

import argparse
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import decouple as DC
import error_budget as EB
import phaseA as A
import tier3_scenes as t3

TILE = 8
N_WATER = 4
NB = 3
H_GRID = [0.5, 1.0, 2.0, 3.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0]


def scene_header(water, n, n_regions, noise):
    lines = ["SETUP %d %d %d 0 %d" % (n, len(t3.BANDS), n_regions, NB)]
    for k in range(n):
        P, G, X, Y, theta_air, tide = t3.SCENE_SETS[water][k]
        lines.append("SCENE %d %.10g 0 %.10g %s" % (k, t3.theta_w(theta_air), tide, " ".join("%g" % l for l in t3.BANDS)))
    if noise > 0:
        for k in range(n):
            lines.append("NOISE %d %s" % (k, " ".join(["%g" % noise]*len(t3.BANDS))))
    return lines


def rrs_lines(grids, cells, n):
    return ["RRS %d %d %s" % (r, k, " ".join("%.12g" % grids[k][m][a, b] for m in range(len(t3.BANDS))))
            for r, (a, b) in enumerate(cells) for k in range(n)]


def pixel_fit(water, n, grids, cell, wvec, noise):
    """Depth and bottom of one pixel, fitted in its 3x3 window (as photic) with the water fixed.
    Returns (the pixel's own parameters [H, bottom block, water], its reported depth, Rrs error)."""
    cells, origin = A.window(cell[0], cell[1], 2)
    nr = len(cells)
    lines = scene_header(water, n, nr, noise)
    lines[0] = "SETUP %d %d %d %d %d" % (n, len(t3.BANDS), nr, origin, NB)
    setup = "\n".join(lines + rrs_lines(grids, cells, n) + ["END"]) + "\n"
    starts = A.lab(setup, ["START %g 0.3" % h for h in H_GRID])
    cmds = []
    for s_ in starts:
        v = np.array(s_["START"])
        v[-N_WATER*n:] = wvec
        f = np.ones(len(v))
        f[-N_WATER*n:] = 0.0
        f[:nr] = 0.0                       # depths held on the grid
        cmds.append("LM %s %s" % (A.vec(v), A.vec(f)))
    prof = A.lab(setup, cmds)
    order = np.argsort([p["rrs_error"] for p in prof])[:2]
    cmds = []
    for k in order:
        v = np.array(prof[k]["params"])
        f = np.ones(len(v))
        f[-N_WATER*n:] = 0.0
        cmds.append("LM %s %s" % (A.vec(v), A.vec(f)))
    ref = A.lab(setup, cmds)
    best = min(ref, key=lambda r: r["rrs_error"])
    p = np.array(best["params"])
    block = 2*NB
    own = np.concatenate([[p[origin]], p[nr + origin*block:nr + (origin + 1)*block], p[-N_WATER*n:]])
    return own, A.depth_of(best, nr, origin), best["rrs_error"]


def water_fit(water, n, grids, cells, pix, wvec, noise):
    """The tile's water with every pixel's depth and bottom fixed."""
    nr = len(cells)
    setup = "\n".join(scene_header(water, n, nr, noise) + rrs_lines(grids, cells, n) + ["END"]) + "\n"
    block = 2*NB
    v = np.concatenate([[p[0] for p in pix]] + [p[1:1 + block] for p in pix] + [wvec])
    f = np.zeros(len(v))
    f[-N_WATER*n:] = 1.0
    r = A.lab(setup, ["LM %s %s" % (A.vec(v), A.vec(f))])[0]
    return np.array(r["params"][-N_WATER*n:])


def start_water(water, n, grids, cells, noise):
    """photic's starting water for the tile (samodel_initial over the tile's pixels)."""
    setup = "\n".join(scene_header(water, n, len(cells), noise) + rrs_lines(grids, cells, n) + ["END"]) + "\n"
    return np.array(A.lab(setup, ["START 5 0.3"])[0]["START"][-N_WATER*n:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--realistic", action="store_true")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--scenes", type=int, default=2)
    ap.add_argument("--init", choices=["start", "deep"], default="start")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--shallow-only", type=float, default=None,
                    help="water step uses only pixels whose fitted depth is below this (m)")
    ap.add_argument("water", nargs="*", default=["clear", "moderate", "turbid"])
    a = ap.parse_args()
    if a.realistic:
        t3.ALBEDO_BY_BOTTOM = EB.REALISTIC
    n = a.scenes
    sample = {(i, j) for i in range(0, 32, 2) for j in (8, 20, 28)}
    for w in a.water:
        H, b, g = A.scene_grids(w, n, True, a.noise)
        deep = (np.concatenate([DC.deep_water(w, k, R) for k, R in enumerate(DC.deep_rrs(w, n))])
                if a.init == "deep" else None)
        # Tiles holding the sample pixels.
        tiles = sorted({(i//TILE, j//TILE) for i, j in sample})
        depth = {}
        history = []

        def run_tile(tile):
            ti, tj = tile
            cells = [(i, j) for i in range(ti*TILE, (ti + 1)*TILE) for j in range(tj*TILE, (tj + 1)*TILE)]
            wvec = deep.copy() if deep is not None else start_water(w, n, g, cells, a.noise)
            log = []
            for rnd in range(a.rounds):
                fits = [pixel_fit(w, n, g, c, wvec, a.noise) for c in cells]
                pix = [f[0] for f in fits]
                reported = [f[1] for f in fits]
                log.append(float(np.mean([f[2] for f in fits])))
                if rnd < a.rounds - 1:
                    use = [k for k, p in enumerate(pix) if a.shallow_only is None or p[0] < a.shallow_only]
                    if len(use) >= 8:
                        wvec = water_fit(w, n, g, [cells[k] for k in use], [pix[k] for k in use], wvec, a.noise)
            return {c: d for c, d in zip(cells, reported)}, wvec, log

        with ThreadPoolExecutor(12) as ex:
            for d, wv, log in ex.map(run_tile, tiles):
                depth.update(d)
                history.append(log)
        rows = np.array([(H[c], depth[c]) for c in sorted(sample)])
        print("%s, %d scenes, OSOAA%s, noise %g, tile %dx%d, water from %s, %d rounds%s; mean Rrs error by round %s" % (
            w, n, ", realistic bottoms" if a.realistic else "", a.noise, TILE, TILE, a.init, a.rounds,
            ", water from pixels < %g m" % a.shallow_only if a.shallow_only else "",
            " -> ".join("%.2f" % v for v in np.mean(history, axis=0))))
        cells = []
        for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 15), (15, 20), (20, 30)):
            m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
            e = rows[m, 1] - rows[m, 0]
            cells.append("%+6.2f (%5.2f)" % (e.mean(), np.sqrt((e**2).mean())))
        print("   tile    " + "".join("%17s" % c for c in cells), flush=True)


if __name__ == "__main__":
    main()
