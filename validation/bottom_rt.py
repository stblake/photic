"""Bottom reflectance from first principles (EXPERIMENTS.md B8, step 1).

Physical models of benthic irradiance reflectance (submerged, 400-800 nm, 1 nm), each a
family over parameter ranges rather than one spectrum:

  sand      Kubelka-Munk semi-infinite particulate layer, R = 1 + q - sqrt(q^2 + 2q), with
            q = k/s from a mineral background, iron oxides and organic matter (exponential
            blue absorption), and microphytobenthos pigments (chlorophyll a, chlorophyll c2,
            fucoxanthin). Submersion lowers the grain-water refractive-index contrast and so
            the scattering s, which raises q.
  seagrass  PROSPECT-D leaves (Feret et al. 2017) with the leaf-water relative refractive
            index (submerged leaves), in a Kubelka-Munk canopy of leaf area index L over a
            sand background.
  coral     zooxanthellate tissue (chlorophyll a, chlorophyll c2, peridinin,
            diadinoxanthin) over an aragonite skeleton: a weakly scattering absorbing layer,
            R = R_t + T^2 R_sk/(1 - R_t R_sk), T = exp(-2 sum c_i a_i(lambda)).

Absorption spectra:
  - PROSPECT-D specific absorption and leaf refractive index (from the prosail package data,
    raw/bottom/prosail_pkg);
  - concentration-specific absorption of extracted phytoplankton pigments (Clementson &
    Wojtasiewicz 2019, Data in Brief 24:103875, CC BY 4.0; raw/bottom/pigments), in m^2 mg^-1,
    shifted to in vivo positions (chlorophyll a red band +12 nm, other bands +10 nm), an
    approximation adequate for Sentinel-2's broad bands.

The parameter ranges are assumptions to be checked against measured libraries (step 2).

  bottom_rt.py    prints each family's reflectance range at Sentinel-2 bands and compares
                  its shapes (normalised at 550 nm) with photic's tabulated bottoms.
"""

import os
import re

import numpy as np

DATA = os.environ.get("PHOTIC_VALIDATION_DATA", os.path.expanduser("~/photic_dev/validation_data"))
LAM = np.arange(400.0, 801.0)
N_WATER = 1.333
S2_BANDS = [443.0, 490.0, 560.0, 665.0, 705.0, 740.0]
RNG = np.random.default_rng(2026)


# Absorption spectra.

def prospect_d():
    d = np.loadtxt(os.path.join(DATA, "raw", "bottom", "prosail_pkg", "prosail", "prospect_d_spectra.txt"))
    m = (d[:, 0] >= LAM[0]) & (d[:, 0] <= LAM[-1])
    wl, nr, kab, kcar, kant, kbrown, kw, km = d[m].T
    return dict(nr=nr, kab=kab, kcar=kcar, kant=kant, kbrown=kbrown, kw=kw, km=km)


def pigments():
    """Specific absorption (m^2 mg^-1) of the pigments used, at in vivo positions, on LAM."""
    out = {}
    base = os.path.join(DATA, "raw", "bottom", "pigments", "supp")
    for f in ("mmc3.txt", "mmc2.txt"):
        names = [n.strip().strip('"') for n in open(os.path.join(base, f)).readline().rstrip("\n").split("\t")[1:]]
        d = np.genfromtxt(os.path.join(base, f), skip_header=1, delimiter="\t")
        wl = d[:, 0]
        order = np.argsort(wl)
        for k, name in enumerate(names):
            out[name.lower()] = (wl[order], np.clip(d[order, k + 1], 0.0, None))
    def shifted(name, red_shift=10.0, blue_shift=10.0, split=600.0):
        wl, a = out[name]
        # In vivo: shift bands to longer wavelengths (protein binding, packaging ignored).
        w = np.where(wl < split, wl + blue_shift, wl + red_shift)
        return np.interp(LAM, w, a, left=0.0, right=0.0)
    return {
        "chl_a": shifted("chlorophyll-a", red_shift=12.0, blue_shift=6.0),
        "chl_c2": shifted("chlorophyll c2"),
        "fucoxanthin": shifted("fucoxanthin"),
        "peridinin": shifted("peridinin"),
        "diadinoxanthin": shifted("diadinoxanthin"),
    }


# Interface optics (Stern 1964; Allen 1973), as in PROSPECT.

def tav(alpha, nr):
    n2 = nr*nr
    npx = n2 + 1
    nm = n2 - 1
    a = (nr + 1)*(nr + 1)/2.0
    k = -(n2 - 1)*(n2 - 1)/4.0
    sa = np.sin(np.radians(alpha))
    b1 = np.sqrt((sa*sa - npx/2)**2 + k) if alpha != 90 else 0.0
    b2 = sa*sa - npx/2
    b = b1 - b2
    ts = (k**2/(6*b**3) + k/b - b/2) - (k**2/(6*a**3) + k/a - a/2)
    tp1 = -2*n2*(b - a)/npx**2
    tp2 = -2*n2*npx*np.log(b/a)/nm**2
    tp3 = n2*(1/b - 1/a)/2
    tp4 = 16*n2**2*(n2**2 + 1)*np.log((2*npx*b - nm**2)/(2*npx*a - nm**2))/(npx**3*nm**2)
    tp5 = 16*n2**3*(1/(2*npx*b - nm**2) - 1/(2*npx*a - nm**2))/npx**3
    return (ts + tp1 + tp2 + tp3 + tp4 + tp5)/(2*sa**2)


def prospect_leaf(P, N, cab, car, cbrown, cw=0.0, cm=0.003, n_medium=N_WATER, alpha=40.0):
    """PROSPECT-D leaf reflectance and transmittance in a medium of index n_medium."""
    from scipy.special import expi
    nr = P["nr"]/n_medium
    kall = (cab*P["kab"] + car*P["kcar"] + cbrown*P["kbrown"] + cw*P["kw"] + cm*P["km"])/N
    tau = np.ones_like(kall)
    j = kall > 0
    tau[j] = (1 - kall[j])*np.exp(-kall[j]) + kall[j]**2*(-expi(-kall[j]))
    talf = tav(alpha, nr)
    t12 = tav(90, nr)
    ralf, r12 = 1 - talf, 1 - t12
    t21 = t12/(nr*nr)
    r21 = 1 - t21
    denom = 1 - r21*r21*tau*tau
    Ta = talf*tau*t21/denom
    Ra = ralf + r21*tau*Ta
    t = t12*tau*t21/denom
    r = r12 + r21*tau*t
    D = np.sqrt((1 + r + t)*(1 + r - t)*(1 - r + t)*(1 - r - t))
    a = (1 + r*r - t*t + D)/(2*r)
    b = (1 - r*r + t*t + D)/(2*t)
    bN1 = b**(N - 1)
    den = a*a*bN1*bN1 - 1
    Rs = a*(bN1*bN1 - 1)/den
    Ts = bN1*(a*a - 1)/den
    den2 = 1 - Rs*r
    return Ra + Ta*Rs*t/den2, Ta*Ts/den2


def km_layer(k, s, thickness, R_background):
    """Kubelka-Munk reflectance of a layer (absorption k, scattering s, per unit thickness)
    over a background."""
    a = 1 + k/s
    b = np.sqrt(a*a - 1)
    coth = 1/np.tanh(b*s*thickness)
    return (1 - R_background*(a - b*coth))/(a - R_background + b*coth)


# Bottom families.

def sand(n=400, P=None, pig=None):
    """Wet sediment, R = 1 + q - sqrt(q^2 + 2q), q = k/s."""
    pig = pig or pigments()
    out = []
    for _ in range(n):
        q0 = RNG.uniform(0.02, 0.15)        # mineral background k/s (carbonate/quartz, wet)
        fe = RNG.uniform(0.0, 0.25)         # iron-oxide k/s at 440 nm
        org = RNG.uniform(0.0, 0.15)        # organic (humic-like) k/s at 440 nm
        mpb = RNG.uniform(0.0, 30.0)        # microphytobenthos chl a, mg m^-2 (per unit s)
        q = (q0 + fe*np.exp(-0.012*(LAM - 440.0)) + org*np.exp(-0.011*(LAM - 440.0))
             + 0.02*mpb*(pig["chl_a"] + 0.3*pig["chl_c2"] + 0.6*pig["fucoxanthin"]))
        out.append(1 + q - np.sqrt(q*q + 2*q))
    return np.array(out)


def seagrass(n=400, P=None, background=None):
    P = P or prospect_d()
    bg = background if background is not None else np.median(sand(100), axis=0)
    out = []
    for _ in range(n):
        N = RNG.uniform(1.2, 2.0)
        cab = RNG.uniform(15.0, 60.0)       # ug cm^-2
        car = cab*RNG.uniform(0.15, 0.3)
        cbrown = RNG.uniform(0.0, 0.3)      # epiphytes and senescence
        L = RNG.uniform(1.0, 6.0)           # leaf area index
        rL, tL = prospect_leaf(P, N, cab, car, cbrown)
        k = np.clip(1 - rL - tL, 1e-6, None)
        s = np.clip(rL, 1e-6, None)
        out.append(km_layer(k, s, L, bg))
    return np.array(out)


def coral(n=400, pig=None):
    pig = pig or pigments()
    out = []
    for _ in range(n):
        chl = RNG.uniform(5.0, 60.0)        # zooxanthellae chl a, mg m^-2
        tau = chl*(pig["chl_a"] + 0.4*pig["chl_c2"] + 0.6*pig["peridinin"] + 0.3*pig["diadinoxanthin"])
        tau *= RNG.uniform(0.4, 0.8)        # pigment packaging in symbiont cells
        R_sk = RNG.uniform(0.6, 0.9)        # wet aragonite skeleton
        R_t = RNG.uniform(0.01, 0.05)       # tissue scattering
        T = np.exp(-2.0*tau)
        out.append(R_t + T*T*R_sk/(1 - R_t*R_sk))
    return np.array(out)


def photic_shapes():
    """photic's tabulated bottoms (model/samodel.c), normalised shapes on LAM."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "model", "samodel.c")).read()
    def arr(name):
        m = re.search(r"double %s\[\] = \{([^}]*)\}" % name, src)
        return np.array([float(x) for x in m.group(1).replace("\n", " ").split(",") if x.strip()])
    wl = arr("ref_wlens")
    return {k: np.interp(LAM, wl, arr(v)) for k, v in (("sand", "bottom_type_sand"),
            ("seagrass", "bottom_type_seagrass"), ("coral", "bottom_type_coral"))}


def at_bands(R):
    return np.array([[np.interp(b, LAM, r) for b in S2_BANDS] for r in np.atleast_2d(R)])


def main():
    P = prospect_d()
    pig = pigments()
    fams = {"sand": sand(400, pig=pig), "seagrass": seagrass(400, P=P), "coral": coral(400, pig=pig)}
    shapes = photic_shapes()
    print("Absolute reflectance (submerged) at Sentinel-2 bands: 5th / 50th / 95th percentile")
    print("  %-9s" % "" + "".join("%18s" % ("%g nm" % b) for b in S2_BANDS))
    for name, R in fams.items():
        B = at_bands(R)
        print("  %-9s" % name + "".join("%18s" % ("%.3f/%.3f/%.3f" % tuple(np.percentile(B[:, k], [5, 50, 95])))
                                         for k in range(len(S2_BANDS))))
    print("\nShape against photic's table (both normalised at 550 nm), at the bands:")
    for name, R in fams.items():
        i550 = int(np.argmin(np.abs(LAM - 550.0)))
        model = np.median(R/R[:, [i550]], axis=0)
        tab = shapes[name]/np.interp(550.0, LAM, shapes[name])
        print("  %-9s model " % name + " ".join("%.2f" % np.interp(b, LAM, model) for b in S2_BANDS)
              + "   photic " + " ".join("%.2f" % np.interp(b, LAM, tab) for b in S2_BANDS))


if __name__ == "__main__":
    main()


# Step 2: measured spectra and calibration of the families.

MEASURED = os.path.join(DATA, "raw", "bottom", "measured")
FIT_LAM = np.arange(410.0, 681.0, 5.0)  # Heron Reef reads 0 at 400 nm; signal is lost beyond ~680 nm


def heron_reef():
    """Roelfsema & Phinn (2012), PANGAEA.804589: in situ underwater reflectance, Heron Reef."""
    import collections
    txt = open(os.path.join(MEASURED, "pangaea_804589.tab"), encoding="utf-8").read()
    ev = collections.OrderedDict()
    for r in txt[txt.index("*/") + 2:].lstrip("\n").splitlines()[1:]:
        f = r.split("\t")
        e = ev.setdefault(f[0], dict(comment=f[8].split(" - ")[0], wl=[], R=[]))
        e["wl"].append(float(f[6]))
        e["R"].append(float(f[7]) if f[7] else np.nan)
    out = []
    for name, e in ev.items():
        c = e["comment"].lower()
        if "sand" in c and "algae" not in c:
            cls = "sediment"
        elif c.startswith("rock") or c.startswith("rubble") and "turf" not in c:
            cls = "sediment"
        elif c.startswith("live") and "coral" in c or c.startswith(("live montipora", "live acropora", "live porites", "live pocillopora")):
            cls = "coral"
        elif "(algae)" in c or c.startswith(("green algae", "algae", "turf")):
            cls = "algae"
        else:
            continue
        wl, R = np.array(e["wl"]), np.array(e["R"])
        ok = np.isfinite(R)
        out.append(dict(source="Heron Reef", name=name, comment=e["comment"], cls=cls,
                        R=np.interp(FIT_LAM, wl[ok], R[ok])))
    return out


def moreton_bay():
    """Roelfsema, Phinn & Joyce (2016), PANGAEA.864310: underwater (Diablo Passage, Deception
    Bay) seagrass and sediment reflectance, Moreton Bay. Samples measured in air are excluded."""
    import openpyxl
    ws = openpyxl.load_workbook(os.path.join(MEASURED, "pangaea_864310.xlsx"), read_only=True, data_only=True).worksheets[0]
    rows = [r for r in ws.iter_rows(values_only=True)]
    wl = np.array([r[0] for r in rows[17:]], dtype=float)
    out = []
    for c in range(2, len(rows[16])):
        name, place, fov = rows[16][c], str(rows[0][c]), str(rows[5][c])
        if name is None or "uw" not in fov:
            continue
        low = str(name).lower()
        if any(g in low for g in ("halophila", "syringodium", "zostera", "halodu", "cymad")):
            cls = "seagrass"
        elif "sand" in low:
            cls = "sediment"
        else:
            continue
        if "Roof" in place:
            continue
        R = np.array([r[c] if isinstance(r[c], (int, float)) else np.nan for r in rows[17:]], dtype=float)
        ok = np.isfinite(R)
        out.append(dict(source="Moreton Bay", name=str(name), comment=place, cls=cls,
                        R=np.interp(FIT_LAM, wl[ok], R[ok])))
    return out


class Families:
    """Model families as functions of a parameter vector, on FIT_LAM, with bounds."""

    def __init__(self):
        self.P = prospect_d()
        self.pig = pigments()
        self.Pf = {k: np.interp(FIT_LAM, LAM, v) for k, v in self.P.items()}
        self.pigf = {k: np.interp(FIT_LAM, LAM, v) for k, v in self.pig.items()}
        self.brown = self.Pf["kbrown"]/self.Pf["kbrown"].max()

    def sediment(self, x):
        q0, fe, org, mpb = x
        L = FIT_LAM
        q = (q0 + fe*np.exp(-0.012*(L - 440.0)) + org*np.exp(-0.011*(L - 440.0))
             + 0.02*mpb*(self.pigf["chl_a"] + 0.3*self.pigf["chl_c2"] + 0.6*self.pigf["fucoxanthin"]))
        return 1 + q - np.sqrt(q*q + 2*q)

    sediment_bounds = ([0.0, 0.0, 0.0, 0.0], [5.0, 3.0, 3.0, 200.0])
    sediment_names = ["q_mineral", "q_ironoxide_440", "q_organic_440", "mpb_chl"]

    def seagrass(self, x):
        N, cab, car_ratio, cbrown, L, cover, q_bg = x
        P = {k: v for k, v in self.Pf.items()}
        rL, tL = prospect_leaf(P, N, cab, cab*car_ratio, cbrown)
        k = np.clip(1 - rL - tL, 1e-6, None)
        s = np.clip(rL, 1e-6, None)
        bg = 1 + q_bg - np.sqrt(q_bg*q_bg + 2*q_bg)
        return cover*km_layer(k, s, L, bg) + (1 - cover)*bg

    seagrass_bounds = ([1.0, 2.0, 0.1, 0.0, 0.2, 0.0, 0.0], [2.5, 80.0, 0.4, 2.0, 8.0, 1.0, 5.0])
    seagrass_names = ["N", "Cab_ug_cm2", "Car/Cab", "Cbrown", "LAI", "cover", "q_background"]

    def coral(self, x):
        chl, R_sk, R_t, brown = x
        tau = chl*(self.pigf["chl_a"] + 0.4*self.pigf["chl_c2"] + 0.6*self.pigf["peridinin"]
                   + 0.3*self.pigf["diadinoxanthin"]) + brown*self.brown
        T = np.exp(-2.0*tau)
        return R_t + T*T*R_sk/(1 - R_t*R_sk)

    coral_bounds = ([0.0, 0.05, 0.0, 0.0], [150.0, 0.95, 0.1, 3.0])
    coral_names = ["chl_mg_m2 (x packaging)", "R_skeleton", "R_tissue", "brown_tau"]


def fit_spectrum(fam, kind, R):
    from scipy.optimize import least_squares
    f = getattr(fam, kind)
    lo, hi = getattr(fam, kind + "_bounds")
    w = 1.0/np.maximum(R, 0.01)
    best = None
    rng = np.random.default_rng(0)
    for _ in range(8):
        x0 = rng.uniform(lo, hi)
        r = least_squares(lambda x: (f(x) - R)*w, x0, bounds=(lo, hi))
        if best is None or r.cost < best.cost:
            best = r
    model = f(best.x)
    rel = np.sqrt(np.mean(((model - R)/np.maximum(R, 0.01))**2))
    return best.x, model, rel


def calibrate():
    fam = Families()
    meas, seen = [], []
    for m in heron_reef() + moreton_bay():
        if any(np.allclose(m["R"], r) for r in seen):
            continue  # duplicate spectrum (Moreton Bay lists one sand average twice)
        seen.append(m["R"])
        meas.append(m)
    kind_of = {"sediment": "sediment", "seagrass": "seagrass", "coral": "coral"}
    summary = {}
    for cls, kind in kind_of.items():
        sel = [m for m in meas if m["cls"] == cls]
        fits = [fit_spectrum(fam, kind, m["R"]) for m in sel]
        X = np.array([f[0] for f in fits])
        rel = np.array([f[2] for f in fits])
        B = np.array([[np.interp(b, FIT_LAM, m["R"]) for b in (443.0, 490.0, 560.0, 665.0)] for m in sel])
        summary[cls] = dict(n=len(sel), rel=rel, X=X, names=getattr(fam, kind + "_names"), B=B,
                            sources=sorted({m["source"] for m in sel}))
    return summary


def report_calibration():
    s = calibrate()
    for cls, d in s.items():
        print("%s: %d measured spectra (%s)" % (cls, d["n"], ", ".join(d["sources"])))
        print("  fit error (RMS relative, 400-680 nm): median %.1f%%, 90th percentile %.1f%%" % (
            100*np.median(d["rel"]), 100*np.percentile(d["rel"], 90)))
        for k, nm in enumerate(d["names"]):
            print("  %-26s 10th %8.3f  median %8.3f  90th %8.3f" % (nm, *np.percentile(d["X"][:, k], [10, 50, 90])))
        print("  measured R at 443/490/560/665 nm, 10th-90th percentile: " + "  ".join(
            "%.3f-%.3f" % tuple(np.percentile(d["B"][:, k], [10, 90])) for k in range(4)))
