"""Independent reference implementation of photic's semi-analytical forward model.

Written from the paper (paper/photic.tex, section 2) and Lee et al. (1998, 1999, 2002),
not from photic's C code, so that tier 1 can check the C implementation against the
equations. Two modes:

  mode="photic"  the choices photic's code makes on purpose (chi = Rrs440/Rrs490,
                 Y and kappa clamped to [0, 2.5]); the C code should match this to
                 rounding.
  mode="paper"   the equations exactly as the paper states them (chi with the 750 nm
                 correction, no clamps).

Each departure can also be switched individually to measure its effect.
"""

import numpy as np

S_CDOM = 0.015  # nm^-1, photic's fixed gelbstoff slope (paper: 0.011-0.021).


def forward(lam, P, G, X, B, rho_table, H, tide, theta_sun, theta_view,
            R440, R490, R750, a0, a1, aw, bbw,
            chi_750=False, clamp_Y=True, clamp_kappa=True, Y=None):
    """Subsurface rrs, Rrs and intermediates for arrays of cases.

    lam in nm; angles are subsurface, in degrees; B * rho_table is the bottom
    albedo rho(lambda); H is depth below datum and tide the scene's offset. Y, if
    given, overrides the backscattering slope derived from R440/R490 (tier 2 sets it
    directly so that photic and the RT model see the same b_bp).
    """
    P, G = np.abs(P), np.abs(G)

    # Absorption: a = a_w + a_phi + a_g (paper eqs. for a_phi(P), a_g(G)).
    a_phi = (a0 + a1*np.log(P))*P
    a_g = G*np.exp(-S_CDOM*(lam - 440.0))
    a = aw + a_phi + a_g

    # Backscattering: b_b = b_bw + X (440/lambda)^Y.
    if chi_750:
        chi = (R440 - R750)/(R490 - R750)
    else:
        chi = R440/R490
    if Y is None:
        Y = 3.44*(1.0 - 3.17*np.exp(-2.01*chi))
        if clamp_Y:
            Y = np.clip(Y, 0.0, 2.5)
    b_bp = X*(440.0/lam)**Y
    bb = bbw + b_bp

    # Lee et al. 1998/1999.
    u = bb/(a + bb)
    kappa = a + bb
    if clamp_kappa:
        kappa = np.clip(kappa, 0.0, 2.5)
    rrs_dp = (0.084 + 0.170*u)*u
    D_C = 1.03*np.sqrt(1.0 + 2.4*u)
    D_B = 1.04*np.sqrt(1.0 + 5.4*u)

    # Water column at the time of the scene (paper: H + H_tide).
    Hs = np.maximum(H + tide, 0.0)
    sec_sun = 1.0/np.cos(np.radians(theta_sun))
    sec_view = 1.0/np.cos(np.radians(theta_view))
    rho = B*rho_table
    rrs_C = rrs_dp*(1.0 - np.exp(-(sec_sun + D_C*sec_view)*kappa*Hs))
    rrs_B = rho/np.pi*np.exp(-(sec_sun + D_B*sec_view)*kappa*Hs)
    rrs = rrs_C + rrs_B

    # Lee et al. 2002: above-surface Rrs from subsurface rrs.
    Rrs = 0.5*rrs/(1.0 - 1.5*rrs)

    return dict(a_phi=a_phi, a_g=a_g, a=a, chi=chi, Y=Y, b_bp=b_bp, bb=bb, u=u,
                K=kappa, rrs_dp=rrs_dp, D_C=D_C, D_B=D_B, rrs_C=rrs_C, rrs_B=rrs_B,
                rrs=rrs, Rrs=Rrs)


def bbw_lee(lam):
    """Pure seawater backscattering, Morel (1974) as used by Lee: 0.0038 (400/lambda)^4.32."""
    return 0.0038*(400.0/lam)**4.32
