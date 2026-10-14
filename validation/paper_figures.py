"""Figures and numbers for the profile-model paper (run as paper_profile.py figures)."""

import glob
import json
import os
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import Rectangle

import ab
import profile_model as pm
import tier2_osoaa as t2
import tier3_scenes
import tier4_prepare
import tier4_score
from paper_profile import DATA, RUNS, OUT, SITES_DIR, FIGS, SITES, SITE_LABEL, ABLATIONS

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "savefig.dpi": 300,
                     "font.family": "serif"})
DEPTH_MAX = 20.0
C_SA, C_PM, C_PUB = "#888888", "#1f5fa8", "#d07a1c"


def hits(m):
    return 100*(m.get("frac_within_tol") or 0)*m["n_compared"]/m["n_reference"]


def site_outputs(run, site, name="H"):
    path = os.path.join(RUNS, "tier4", run, site, "modelled_%s.nc" % name)
    if not os.path.exists(path):
        return None
    return tier3_scenes.read_nc(path)


def depth_of(run, site):
    H = site_outputs(run, site)
    return None if H is None else -H


def published(site, ref):
    pub = glob.glob(os.path.join(DATA, "raw", "published", "exmouth_gulf", "*.tif"))
    return tier4_score.published_on_grid(pub[0], ref["lon"], ref["lat"])


def sa_run(a, site):
    """HEAD's tier 4 run if it has this site, else candidate 2's (the same default model)."""
    for run in (a.sa, "head-d9d14d6", "a54973d"):
        if site_outputs(run, site) is not None:
            return run
    raise SystemExit("no semi-analytical run for " + site)


# Numbers.

def tier4_table(a):
    rows = {}
    for site in SITES:
        ref = np.load(os.path.join(SITES_DIR, site, "reference.npz"))
        D = ref["depth"]
        sa = sa_run(a, site)
        out = dict(sa_run=sa)
        for key, H, sig in (("sa", depth_of(sa, site), site_outputs(sa, site, "H_sigma")),
                            ("profile", depth_of("profile-wt", site), site_outputs("profile-wt", site, "H_sigma")),
                            ("published", published(site, ref), None)):
            m = tier4_score.metrics(D, H)
            m["hits"] = hits(m)
            m["coverage"] = 100*m["n_compared"]/m["n_reference"]
            v = np.isfinite(D)
            ok = v & np.isfinite(H) & (H > 0)
            hit = ok & (np.abs(np.where(ok, H, 0) - np.where(v, D, 0)) <= 0.5 + 0.1*np.where(v, D, 0))
            m["bins"] = {"%g-%g" % (lo, hi): dict(n=int((v & (D >= lo) & (D < hi)).sum()),
                         coverage=float(100*(ok & (D >= lo) & (D < hi)).sum()/max((v & (D >= lo) & (D < hi)).sum(), 1)),
                         hits=float(100*(hit & (D >= lo) & (D < hi)).sum()/max((v & (D >= lo) & (D < hi)).sum(), 1)))
                         for lo, hi in ((0, 5), (5, 10), (10, 30))}
            m["frac_H_over_5"] = float(np.mean(H[ok] > 5.0)) if ok.any() else None
            if sig is not None:
                v = np.isfinite(D) & np.isfinite(H) & (H > 0) & (np.abs(sig) > 0)
                e = np.abs(H[v] - D[v])
                s = np.abs(sig[v])
                m["sigma_cal"] = dict(within_1=float(np.mean(e <= s)), within_2=float(np.mean(e <= 2*s)),
                                      median_sigma=float(np.median(s)), median_err=float(np.median(e)))
            out[key] = m
        iod = site_outputs("profile-wt", site, "index_optical_depth")
        v = np.isfinite(D)
        out["profile_reasons"] = {k: float(np.mean(np.rint(iod[v]) == c)) for k, c in
                                  (("depth", 1), ("deep", 2), ("ambiguous", 3), ("no_data", 0))}
        out["reference"] = dict(n=int(v.sum()), min=float(np.nanmin(D)), max=float(np.nanmax(D)),
                                median=float(np.nanmedian(D)),
                                frac_by_bin={"%g-%g" % b: float(np.mean((D[v] >= b[0]) & (D[v] < b[1])))
                                             for b in tier4_score.DEPTH_BINS})
        out["seconds"] = dict(sa=json.load(open(os.path.join(RUNS, sa, "tier4.json")))["sites"][site]["seconds"]
                              if os.path.exists(os.path.join(RUNS, sa, "tier4.json")) else None,
                              profile=json.load(open(os.path.join(RUNS, "profile-wt", "tier4.json")))["sites"][site]["seconds"])
        rows[site] = out
    return rows


def tier3_table(a):
    A = json.load(open(os.path.join(RUNS, a.sa3, "tier3.json")))["scenarios"]
    B = json.load(open(os.path.join(RUNS, "profile-wt", "tier3.json")))["scenarios"]
    rows = {}
    for n in sorted(A):
        if n not in B or "rmse" not in A[n]:
            continue
        row = {}
        for key, v in (("sa", A[n]), ("profile", B[n])):
            row[key] = dict(coverage=100*v["frac_modelled"], tol=100*v["frac_within_tol"],
                            hits=100*v["frac_within_tol"]*v["frac_modelled"], rmse=v["rmse"], bias=v.get("bias"),
                            by_depth=v.get("by_depth"))
        rows[n] = row
    return rows


def ablation_table():
    rows = {}
    for name in ABLATIONS:
        p4 = os.path.join(RUNS, "paper-" + name, "tier4.json")
        p3 = os.path.join(RUNS, "paper-" + name, "tier3.json")
        if not os.path.exists(p4):
            continue
        t4 = json.load(open(p4))["sites"]
        row = {s: dict(hits=hits(t4[s]["photic"]), tol=100*(t4[s]["photic"].get("frac_within_tol") or 0),
                       coverage=100*t4[s]["photic"]["n_compared"]/t4[s]["photic"]["n_reference"],
                       rmse=t4[s]["photic"].get("rmse")) for s in SITES}
        if os.path.exists(p3):
            t3 = json.load(open(p3))["scenarios"]
            for w in ("clear", "moderate", "turbid"):
                hs = [100*v["frac_within_tol"]*v["frac_modelled"] for k, v in t3.items()
                      if k.startswith(w) and "prior" not in k and "rmse" in v]
                rm = [v["rmse"] for k, v in t3.items() if k.startswith(w) and "prior" not in k and v.get("rmse")]
                row["t3_" + w] = dict(hits=float(np.mean(hs)), rmse=float(np.median(rm)))
        rows[name] = row
    return rows


# Figures.

def fig_study_area():
    f = glob.glob(os.path.join(DATA, "derived", "tier4", "acolite", "*", "*2019_04_28*L2W.nc"))[0]
    import netCDF4
    d = netCDF4.Dataset(f)
    lon, lat = d["lon"][:], d["lat"][:]
    rgb = np.dstack([d[v][:].filled(np.nan) for v in ("Rrs_665", "Rrs_559", "Rrs_492")])
    land = ~np.isfinite(rgb).all(axis=2)
    img = np.clip((rgb/np.array([0.025, 0.035, 0.035]))**(1/2.2), 0, 1)
    img[land] = 0.75
    step = 3
    fig, ax = plt.subplots(figsize=(6.3, 4.6))
    ax.pcolormesh(lon[::step, ::step], lat[::step, ::step], img[::step, ::step], shading="nearest",
                  rasterized=True)
    for site, (w, e, s, n, _) in tier4_prepare.SITES.items():
        ax.add_patch(Rectangle((w, s), e - w, n - s, fill=False, ec="#ffd400", lw=1.2))
        ax.text(w, n + 0.004, SITE_LABEL[site], color="k", fontsize=7, fontweight="bold")
    ax.set_xlim(113.93, 114.47)
    ax.set_ylim(-22.0, -21.6)
    ax.set_aspect(1/np.cos(np.radians(21.8)))
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°)")
    ax.text(113.985, -21.75, "Indian Ocean\n(Ningaloo coast)", fontsize=7, ha="center", color="w")
    ax.text(114.06, -21.93, "North West\nCape", fontsize=7, ha="center")
    ax.text(114.30, -21.85, "Exmouth Gulf", fontsize=8, ha="center", style="italic")
    save(fig, "study_area")


def proto_args():
    return types.SimpleNamespace(ref="red", red_fraction=0.2, green_fraction=0.15, smooth=7, eps=0.05,
                                 glint=True, radius=1.5, min_local=400, t_deep=25.0, t_amb=1.0)


TABLES = {}


def tables_for(lam):
    if lam not in TABLES:
        TABLES[lam] = t2.photic_tables(ab.DRIVER, list(lam))
    return TABLES[lam]


def proto_cell(site, i, j, a=None):
    """The prototype's profile at one cell: chi2(H), chi2 of the deep model, the decision."""
    a = a or proto_args()
    case = os.path.join(SITES_DIR, site)
    scenes = pm.read_case(case)
    lam = scenes[0]["lam"]
    tables = tables_for(tuple(lam))
    shapes = {name: np.repeat(np.array([tables[(l, idx)]["rho"] for l in lam])[None], len(scenes), axis=0)
              for name, (idx, _, _) in pm.CLASSES.items()}
    fields = [pm.smoothed(s, a.smooth) for s in scenes]
    ny, nx = scenes[0]["R"].shape[1:]
    cell_km = float(open(os.path.join(case, "s1_b443.asc")).readlines()[4].split()[1])*111.0
    rad = int(round(a.radius/cell_km))
    region = np.zeros((ny, nx), dtype=bool)
    region[max(i - rad, 0):i + rad + 1, max(j - rad, 0):j + rad + 1] = True
    here = [pm.scene_reference(s, lam, tables, None, a, f, region) for s, f in zip(scenes, fields)]
    for r in here:
        r["shapes"] = shapes
    rrs = np.array([pm.rrs_of(s["R"][:, i, j]) for s in scenes])
    ok = np.isfinite(rrs) & (rrs > 0)
    sig2 = np.array([r["noise"] for r in here])**2 + (a.eps*np.where(ok, rrs, 0.0))**2
    w = np.where(ok, 1.0/sig2, 0.0)
    rrs = np.where(ok, rrs, 0.0)
    chi, cls, bri, chi_deep = pm.profile(rrs, w, here, a.glint)
    H, sg, why = pm.decide(chi, chi_deep, a.t_deep, a.t_amb, int(ok.sum()))
    return dict(chi=chi, chi_deep=chi_deep, cls=cls, bri=bri, H=H, sigma=sg, why=why,
                scale=max(1.0, chi.min()/max(int(ok.sum()) - 2, 1)), rrs=rrs, ok=ok, here=here, lam=lam,
                scenes=scenes, fields=fields, region=region)


def fig_profiles():
    """Three cells at Bundegi: a depth, optically deep and ambiguous."""
    ref = np.load(os.path.join(SITES_DIR, "bundegi", "reference.npz"))
    D = ref["depth"]
    iod = np.rint(site_outputs("profile-wt", "bundegi", "index_optical_depth"))
    H = depth_of("profile-wt", "bundegi")
    rng = np.random.default_rng(3)
    picks = {}
    good = np.argwhere(np.isfinite(D) & (iod == 1) & (np.abs(H - D) <= 0.5 + 0.1*D) & (D > 3) & (D < 6))
    picks["depth"] = tuple(good[rng.integers(len(good))])
    deep = np.argwhere(np.isfinite(D) & (iod == 2) & (D > 9))
    picks["deep"] = tuple(deep[rng.integers(len(deep))])
    amb = np.argwhere(np.isfinite(D) & (iod == 3) & (D > 2) & (D < 8))
    picks["ambiguous"] = tuple(amb[rng.integers(len(amb))])
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.2), sharey=False)
    titles = {"depth": "(a) depth reported", "deep": "(b) optically deep", "ambiguous": "(c) ambiguous"}
    out = {}
    for ax, (k, (i, j)) in zip(axs, picks.items()):
        c = proto_cell("bundegi", int(i), int(j))
        dchi = (c["chi"] - c["chi"].min())/c["scale"]
        ax.plot(pm.H_GRID, dchi, color=C_PM, lw=1)
        ax.axhline((c["chi_deep"] - c["chi"].min())/c["scale"], color="k", ls="--", lw=0.8)
        ax.axhline(25, color="#b03030", ls=":", lw=0.8)
        ax.axvline(D[i, j], color=C_PUB, lw=1)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(0, max(200, 2*(c["chi_deep"] - c["chi"].min())/c["scale"]))
        ax.set_xlim(0, 30)
        ax.set_xlabel("Depth H (m)")
        ax.set_title(titles[k])
        out[k] = dict(cell=[int(i), int(j)], sounding=float(D[i, j]), H=float(c["H"]) if np.isfinite(c["H"]) else None,
                      why=c["why"], dchi_deep=float((c["chi_deep"] - c["chi"].min())/c["scale"]))
    axs[0].set_ylabel(r"$\Delta\chi^2(H)$ / scale")
    fig.legend(["profile", "optically deep model", r"$T_\mathrm{deep}$", "sounding"], loc="lower center", ncol=4,
               frameon=False, bbox_to_anchor=(0.5, -0.16))
    save(fig, "profiles")
    return out


def fig_reference():
    """Deep cells chosen at Marina (scene 1) and every scene's deep-water spectrum at each site."""
    a = proto_args()
    case = os.path.join(SITES_DIR, "marina")
    scenes = pm.read_case(case)
    fields = pm.smoothed(scenes[0], a.smooth)
    valid, red, green = fields
    cand = valid & (red <= np.percentile(red[valid], 100*a.red_fraction))
    deep = cand & (green <= np.percentile(green[cand], 100*a.green_fraction))
    ref = np.load(os.path.join(SITES_DIR, "marina", "reference.npz"))
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8), gridspec_kw=dict(width_ratios=[1, 1.25]))
    R = scenes[0]["R"]
    img = np.where(valid, R[3], np.nan)
    ext = [ref["lon"][0], ref["lon"][-1], ref["lat"][0], ref["lat"][-1]]
    axs[0].imshow(img, origin="lower", extent=ext, cmap="Greys_r", vmin=0, vmax=0.01)
    m = np.ma.masked_where(~cand, np.ones_like(img))
    axs[0].imshow(m, origin="lower", extent=ext, cmap="autumn", alpha=0.35, vmin=0, vmax=1)
    m = np.ma.masked_where(~deep, np.ones_like(img))
    axs[0].imshow(m, origin="lower", extent=ext, cmap="cool", alpha=0.9, vmin=0, vmax=1)
    axs[0].set_aspect(1/np.cos(np.radians(21.9)))
    axs[0].set_title("(a) Marina, 2019-06-12")
    axs[0].ticklabel_format(useOffset=False)
    axs[0].set_xlabel("Longitude (°E)")
    axs[0].set_ylabel("Latitude (°)")
    axs[0].tick_params(labelsize=6)
    axs[0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    for site, col in zip(SITES, ("#1b9e77", "#d95f02", "#7570b3", "#e7298a")):
        sc = pm.read_case(os.path.join(SITES_DIR, site))
        for s in sc:
            f = pm.smoothed(s, a.smooth)
            rinf, noise, n = pm.deep_reference(s, None, a, f)
            axs[1].plot(s["lam"], rinf*1e3, color=col, lw=0.8, marker="o", ms=2,
                        label=SITE_LABEL[site] if s is sc[0] else None)
    axs[1].set_xlabel("Wavelength (nm)")
    axs[1].set_ylabel(r"$r_\mathrm{rs}^\infty$ ($10^{-3}$ sr$^{-1}$)")
    axs[1].set_title("(b) deep-water reference per scene")
    axs[1].legend(frameon=False)
    save(fig, "reference")


def fig_maps(a):
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("#e6e6e6")
    fig, axs = plt.subplots(4, 4, figsize=(6.8, 7.0))
    for r, site in enumerate(SITES):
        ref = np.load(os.path.join(SITES_DIR, site, "reference.npz"))
        ext = [ref["lon"][0], ref["lon"][-1], ref["lat"][0], ref["lat"][-1]]
        H_sa = depth_of(sa_run(a, site), site)
        H_pm = depth_of("profile-wt", site)
        iod = np.rint(site_outputs("profile-wt", site, "index_optical_depth"))
        cols = [("Soundings", ref["depth"]), ("Semi-analytical", np.where(H_sa > 0, H_sa, np.nan)),
                ("Profile model", np.where(H_pm > 0, H_pm, np.nan)), ("Profile: why not", iod)]
        for c, (title, X) in enumerate(cols):
            ax = axs[r, c]
            if c < 3:
                im = ax.imshow(X, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=DEPTH_MAX,
                               interpolation="nearest")
            else:
                from matplotlib.colors import ListedColormap
                cm = ListedColormap(["#ffffff", C_PM, "#203040", "#e0a030"])
                im2 = ax.imshow(np.where(np.isfinite(X), X, 0), origin="lower", extent=ext, cmap=cm, vmin=-0.5,
                                vmax=3.5, interpolation="nearest")
            ax.set_aspect(1/np.cos(np.radians(21.8)))
            ax.set_xticks([])
            ax.set_yticks([])
            if r == 0:
                ax.set_title(title)
            if c == 0:
                ax.set_ylabel(SITE_LABEL[site])
    fig.subplots_adjust(right=0.86, wspace=0.05, hspace=0.12)
    cax = fig.add_axes([0.88, 0.55, 0.015, 0.3])
    fig.colorbar(im, cax=cax, label="Depth below MSL (m)")
    cax2 = fig.add_axes([0.88, 0.12, 0.015, 0.3])
    cb = fig.colorbar(im2, cax=cax2, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(["no data", "depth", "optically\ndeep", "ambiguous"], fontsize=6)
    save(fig, "maps")


def fig_scatter(a):
    fig, axs = plt.subplots(2, 4, figsize=(6.8, 3.6), sharex=True, sharey=True)
    for c, site in enumerate(SITES):
        ref = np.load(os.path.join(SITES_DIR, site, "reference.npz"))
        D = ref["depth"]
        for r, (label, H) in enumerate((("Semi-analytical", depth_of(sa_run(a, site), site)),
                                        ("Profile model", depth_of("profile-wt", site)))):
            ax = axs[r, c]
            v = np.isfinite(D) & np.isfinite(H) & (H > 0)
            x = np.linspace(0, 25, 50)
            ax.fill_between(x, x - 0.5 - 0.1*x, x + 0.5 + 0.1*x, color="#cfe3cf", lw=0)
            ax.plot(x, x, color="k", lw=0.5)
            ax.scatter(D[v], H[v], s=1, alpha=0.25, color=C_SA if r == 0 else C_PM, rasterized=True, lw=0)
            ax.set_xlim(0, 25)
            ax.set_ylim(0, 25)
            m = tier4_score.metrics(D, H)
            ax.text(0.04, 0.96, "n=%d\nRMSE %.2f m\nhits %.0f%%" % (m["n_compared"], m.get("rmse") or np.nan, hits(m)),
                    transform=ax.transAxes, va="top", fontsize=6)
            if r == 0:
                ax.set_title(SITE_LABEL[site])
            if c == 0:
                ax.set_ylabel("%s depth (m)" % label.split()[0])
            if r == 1:
                ax.set_xlabel("Sounding depth (m)")
    fig.tight_layout()
    save(fig, "scatter")


def fig_depth_bins(rows):
    bins = ["0-2", "2-5", "5-10", "10-20"]
    fig, axs = plt.subplots(1, 4, figsize=(6.8, 2.1), sharey=True)
    for ax, site in zip(axs, SITES):
        x = np.arange(len(bins))
        for k, (key, col) in enumerate((("sa", C_SA), ("profile", C_PM), ("published", C_PUB))):
            bd = rows[site][key].get("by_depth", {})
            y = [bd.get(b, {}).get("bias", np.nan) for b in bins]
            e = [bd.get(b, {}).get("rmse", np.nan) for b in bins]
            ax.errorbar(x + (k - 1)*0.22, y, yerr=e, fmt="o", ms=3, color=col, lw=0.8, capsize=1.5,
                        label=["semi-analytical", "profile", "published"][k])
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(bins)
        ax.set_title(SITE_LABEL[site])
        ax.set_xlabel("Sounding depth (m)")
    axs[0].set_ylabel("Bias ± RMSE (m)")
    axs[0].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    save(fig, "depth_bins")


def fig_tier3(rows):
    groups = ["clear", "moderate", "turbid"]
    kinds = [("noise0", "uniform\nno noise"), ("noise0.0005", "uniform\nnoise"),
             ("noise0_realistic", "realistic\nno noise"), ("noise0.0005_realistic", "realistic\nnoise")]
    fig, axs = plt.subplots(2, 1, figsize=(6.8, 3.6), sharex=True)
    labels, sa_h, pm_h, sa_r, pm_r = [], [], [], [], []
    for w in groups:
        for suffix, name in kinds:
            ns = [n for n in rows if n.startswith(w + "_n") and n.endswith("_" + suffix)]
            if not ns:
                continue
            labels.append("%s\n%s" % (w, name))
            sa_h.append(np.mean([rows[n]["sa"]["hits"] for n in ns]))
            pm_h.append(np.mean([rows[n]["profile"]["hits"] for n in ns]))
            sa_r.append(np.median([rows[n]["sa"]["rmse"] for n in ns]))
            pm_r.append(np.median([rows[n]["profile"]["rmse"] for n in ns]))
    x = np.arange(len(labels))
    for ax, (A, B, yl) in zip(axs, ((sa_h, pm_h, "Hits, all cells (%)"), (sa_r, pm_r, "RMSE, covered (m)"))):
        ax.bar(x - 0.2, A, 0.4, color=C_SA, label="semi-analytical")
        ax.bar(x + 0.2, B, 0.4, color=C_PM, label="profile")
        ax.set_ylabel(yl)
        for k in (4, 8):
            ax.axvline(k - 0.5, color="k", lw=0.4)
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(labels, fontsize=5.5)
    axs[0].legend(frameon=False)
    fig.tight_layout()
    save(fig, "tier3")


def fig_scenes():
    path = os.path.join(OUT, "scenes.json")
    if not os.path.exists(path):
        return None
    res = json.load(open(path))
    fig, ax = plt.subplots(figsize=(3.3, 2.4))
    summary = {}
    for site, col in zip(SITES, ("#1b9e77", "#d95f02", "#7570b3", "#e7298a")):
        ks = sorted({r["k"] for r in res if r["site"] == site})
        mean, lo, hi = [], [], []
        for k in ks:
            h = [100*(r["frac_within_tol"] or 0)*r["n_compared"]/r["n_reference"] for r in res
                 if r["site"] == site and r["k"] == k]
            mean.append(np.mean(h)), lo.append(min(h)), hi.append(max(h))
        ax.plot(ks, mean, color=col, marker="o", ms=3, lw=1, label=SITE_LABEL[site])
        ax.fill_between(ks, lo, hi, color=col, alpha=0.15, lw=0)
        summary[site] = dict(k=ks, mean=mean, lo=lo, hi=hi)
    ax.set_xlabel("Number of scenes")
    ax.set_ylabel("Hits over all cells (%)")
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    save(fig, "scenes")
    return summary


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, name + ".png"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def main(a):
    t4 = tier4_table(a)
    t3 = tier3_table(a)
    numbers = dict(tier4=t4, tier3=t3, ablation=ablation_table())
    numbers["profiles"] = fig_profiles()
    numbers["scenes"] = fig_scenes()
    fig_study_area()
    fig_reference()
    fig_maps(a)
    fig_scatter(a)
    fig_depth_bins(t4)
    fig_tier3(t3)
    os.makedirs(OUT, exist_ok=True)
    json.dump(numbers, open(os.path.join(OUT, "numbers.json"), "w"), indent=1, default=float)
    print("numbers:", os.path.join(OUT, "numbers.json"))
