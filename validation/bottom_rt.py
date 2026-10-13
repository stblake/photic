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
