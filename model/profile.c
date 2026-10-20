//
//  MODEL Lee_Profile
//

// A profile model for multispectral imagery (validation/EXPERIMENTS.md, Phase F; the
// prototype is validation/profile_model.py). It is Lee's shallow-water model with the water
// observed rather than fitted, and a global search over the one non-linear unknown, depth.
//
// 1. Water, per scene. The deep-water reflectance rrs_inf is the median spectrum of the
//    optically deep pixels within PROFILE_RADIUS km. Every bottom is far brighter than deep
//    water in the red and near infrared, so the pixels with the lowest (smoothed) red are
//    deeper than the red's reach; among them, the lowest green are the deepest. kappa and u
//    come from rrs_inf by QAA (as samodel_h0), and the deep pixels' spread gives the noise.
//    With PROFILE_KAPPA, the attenuation is calibrated from the scene's bottom contrast instead
//    (after 4SM, Morel & Favoretto 2017): see calibrate() below.
// 2. Per pixel, all scenes at once, with Z_s = H + tide_s:
//      rrs_s = rrs_inf_s [1 - exp(-kappa_s M^C_s Z_s)] + exp(-kappa_s M^B_s Z_s) b rho_c/pi + d_s,
//    with one bottom class c of brightness b in the class range (B8) shared by the scenes, and
//    a flat offset d_s per scene for residual glint (PROFILE_GLINT). b and d_s are linear, so at
//    each depth of a fixed grid they are solved exactly: the profile chi2(H) is exact and its
//    global minimum is found by construction.
// 3. A depth is reported only if the optically deep model is worse by more than
//    PROFILE_T_DEEP, and no depth outside the tolerance (0.5 m + 10%) scores within
//    PROFILE_T_AMB of the best (both in units of the reduced chi2).
//
// Outputs as MODEL Lee_Semi_Analytical. Cells without a depth get 0; the index of optical depth
// says why: 1 a depth, 2 optically deep, 3 ambiguous, 0 no data. The bottom percentages are 100
// for the chosen class, the albedo is its brightness, and the model error is the rms relative
// misfit (%).

#include "profile.h"

// photic's tables (samodel.c).

extern int n_ref_wlens;
extern double ref_wlens[], aw_Pope_Fry1997[], bbw_Morel_1974[], bottom_type_sand[],
  bottom_type_seagrass[], bottom_type_coral[];

// Settings, tuned on the training split (EXPERIMENTS.md, Phase F).

#ifndef PROFILE_RADIUS
#define PROFILE_RADIUS 1.5         // km; 0 takes the deep water from the whole grid
#endif
#ifndef PROFILE_TILE
#define PROFILE_TILE 8             // cells: one local reference per tile
#endif
#ifndef PROFILE_MIN_LOCAL
#define PROFILE_MIN_LOCAL 400      // fewest valid cells for a local reference
#endif
#ifndef PROFILE_SMOOTH
#define PROFILE_SMOOTH 7           // cells: window for ranking cells by red and green
#endif
#ifndef PROFILE_RED_FRACTION
#define PROFILE_RED_FRACTION 0.2
#endif
#ifndef PROFILE_GREEN_FRACTION
#define PROFILE_GREEN_FRACTION 0.15
#endif
#ifndef PROFILE_EPS
#define PROFILE_EPS 0.05           // relative model error
#endif
#ifndef PROFILE_T_DEEP
#define PROFILE_T_DEEP 25.0
#endif
#ifndef PROFILE_T_AMB
#define PROFILE_T_AMB 1.0
#endif
#ifndef PROFILE_GLINT
#define PROFILE_GLINT 1
#endif
#ifndef PROFILE_GLINT_SIGMA
#define PROFILE_GLINT_SIGMA 5.0e-4 // prior on the offset in Rrs, as GLINT_OFFSET
#endif

#ifndef PROFILE_BRIGHT
#define PROFILE_BRIGHT 0           // bottom brightness ranges: 0 the in situ classes; capped from the
#endif                             // image's brightest bottoms: 1 clip each class, 2 scale every range
#ifndef PROFILE_OFFSET_MODEL
#define PROFILE_OFFSET_MODEL 0     // 1: keep the calibration's offset outside the attenuated water-column term
#endif
#ifndef PROFILE_KAPPA_LOCAL
#define PROFILE_KAPPA_LOCAL 0      // 1: calibrate the attenuation in each local window (fallback: site-wide)
#endif
#ifndef PROFILE_KAPPA
#define PROFILE_KAPPA 0            // attenuation: 0 QAA, 1 QAA after an offset matched to the image's
#endif                             // brightest pixels line, 2 Jerlov's Kd for the type matching it

#ifndef PROFILE_RATIO
#define PROFILE_RATIO 0            // 1: add the log band-ratio consistency rows to chi2 (below)
#endif
#ifndef PROFILE_RATIO_WEIGHT
#define PROFILE_RATIO_WEIGHT 1.0
#endif

#define N_CLASSES 3
#define MAX_BANDS 16
#define N_H 220
#define QAA_G0 0.084               // Lee 1999, photic's
#define QAA_G1 0.170

static const double class_lo[N_CLASSES] = {0.15, 0.03, 0.03};
static const double class_hi[N_CLASSES] = {0.75, 0.15, 0.25};

typedef struct {
  double rrs_inf[MAX_BANDS], noise[MAX_BANDS], kappa[MAX_BANDS], MC[MAX_BANDS], MB[MAX_BANDS];
  int n_deep, cal_mode;             // calibration: 0 none (QAA), 1 offset, 2 Jerlov
  double cal_delta, cal_t;
  bool ok;
} reference;

typedef struct {
  int n_scenes, nrows, ncols, n_bands[MAX_SCENES];
  double lam[MAX_SCENES][MAX_BANDS], tide[MAX_SCENES], sec_w[MAX_SCENES], sec_v[MAX_SCENES];
  double aw[MAX_SCENES][MAX_BANDS], bbw[MAX_SCENES][MAX_BANDS], rho[N_CLASSES][MAX_SCENES][MAX_BANDS];
  float **R[MAX_SCENES][MAX_BANDS], nodata[MAX_SCENES][MAX_BANDS];
  double **red[MAX_SCENES], **green[MAX_SCENES];
  char **valid[MAX_SCENES];
  double H[N_H], lo[N_CLASSES], hi[N_CLASSES];
} profile_data;

static double rrs_of(double R) {
  return R/(0.5 + 1.5*R);
}

static int compare_doubles(const void *a, const void *b) {
  double x = *(const double*) a, y = *(const double*) b;
  return (x > y) - (x < y);
}

// Percentile q (0-1) of v[0..n-1], interpolated linearly as numpy's default; v is sorted.

static double percentile(double *v, int n, double q) {
  double h, f;
  int k;

  qsort(v, n, sizeof(double), compare_doubles);
  h = q*(n - 1);
  k = (int) floor(h);
  f = h - k;
  return k + 1 < n ? v[k] + f*(v[k + 1] - v[k]) : v[n - 1];
}

static int nearest_band(profile_data *pd, int s, double lambda) {
  int b, best = 0;

  for (b = 1; b < pd->n_bands[s]; b++) {
    if (fabs(pd->lam[s][b] - lambda) < fabs(pd->lam[s][best] - lambda)) best = b;
  }
  return best;
}

// Valid cells of each scene (every band present and positive), and their red (+ near-infrared)
// and green reflectance smoothed over PROFILE_SMOOTH cells, from summed-area tables.

static void smooth_fields(profile_data *pd, int s) {
  int i, j, b, h = PROFILE_SMOOTH/2, nr = pd->nrows, nc = pd->ncols, g = nearest_band(pd, s, 560.0), i0, i1, j0, j1;
  double **sn, **sr, **sg, v, red, n;

  pd->valid[s] = (char**) malloc(nr*sizeof(char*));
  allocate_double_array_2d(&pd->red[s], nr, nc);
  allocate_double_array_2d(&pd->green[s], nr, nc);
  allocate_double_array_2d(&sn, nr + 1, nc + 1);
  allocate_double_array_2d(&sr, nr + 1, nc + 1);
  allocate_double_array_2d(&sg, nr + 1, nc + 1);
  for (i = 0; i <= nr; i++) {
    for (j = 0; j <= nc; j++) sn[i][j] = sr[i][j] = sg[i][j] = 0.0;
  }
  for (i = 0; i < nr; i++) {
    pd->valid[s][i] = (char*) malloc(nc);
    for (j = 0; j < nc; j++) {
      pd->valid[s][i][j] = 1;
      red = 0.0;
      for (b = 0; b < pd->n_bands[s]; b++) {
        v = pd->R[s][b][i][j];
        if (! isfinite(v) || approx_equal(v, pd->nodata[s][b], 1.0e-6) || v <= 0.0) {
          pd->valid[s][i][j] = 0;
          break ;
        }
        if (pd->lam[s][b] >= 660.0) red += v;
      }
      n = pd->valid[s][i][j];
      sn[i + 1][j + 1] = sn[i][j + 1] + sn[i + 1][j] - sn[i][j] + n;
      sr[i + 1][j + 1] = sr[i][j + 1] + sr[i + 1][j] - sr[i][j] + (n ? red : 0.0);
      sg[i + 1][j + 1] = sg[i][j + 1] + sg[i + 1][j] - sg[i][j] + (n ? pd->R[s][g][i][j] : 0.0);
    }
  }
  for (i = 0; i < nr; i++) {
    i0 = MAX(i - h, 0);
    i1 = MIN(i + h, nr - 1) + 1;
    for (j = 0; j < nc; j++) {
      j0 = MAX(j - h, 0);
      j1 = MIN(j + h, nc - 1) + 1;
      n = sn[i1][j1] - sn[i0][j1] - sn[i1][j0] + sn[i0][j0];
      pd->red[s][i][j] = (sr[i1][j1] - sr[i0][j1] - sr[i1][j0] + sr[i0][j0])/MAX(n, 1.0e-9);
      pd->green[s][i][j] = (sg[i1][j1] - sg[i0][j1] - sg[i1][j0] + sg[i0][j0])/MAX(n, 1.0e-9);
    }
  }
  free_double_array_2d(sn, nr + 1);
  free_double_array_2d(sr, nr + 1);
  free_double_array_2d(sg, nr + 1);
}

// QAA: u from rrs = (g0 + g1 u) u; at 665 nm a ~ a_w fixes b_b; Y from QAA's 443/560 ratio;
// then kappa = a + b_b with a = (1 - u) b_b/u.

static void qaa_spectrum(profile_data *pd, int s, const double *rrs_inf, double *kappa, double *u) {
  int b, i_ref = nearest_band(pd, s, 665.0), i443 = nearest_band(pd, s, 443.0), i560 = nearest_band(pd, s, 560.0);
  double bbp_ref, Y, bb, a;

  for (b = 0; b < pd->n_bands[s]; b++) {
    u[b] = (-QAA_G0 + sqrt(QAA_G0*QAA_G0 + 4.0*QAA_G1*rrs_inf[b]))/(2.0*QAA_G1);
  }
  bbp_ref = MAX(u[i_ref]*pd->aw[s][i_ref]/(1.0 - u[i_ref]) - pd->bbw[s][i_ref], 1.0e-6);
  Y = 2.0*(1.0 - 1.2*exp(-0.9*rrs_inf[i443]/rrs_inf[i560]));
  Y = MIN(MAX(Y, 0.0), 2.5);
  for (b = 0; b < pd->n_bands[s]; b++) {
    bb = pd->bbw[s][b] + bbp_ref*pow(pd->lam[s][i_ref]/pd->lam[s][b], Y);
    a = (1.0 - u[b])*bb/u[b];
    kappa[b] = a + bb;
  }
}

// Jerlov's Kd (m^-1; Austin & Petzold's update for the oceanic types), as jerlov.c: types I, IA,
// IB, II, III, 1C, 3C, 5C, 7C, 9C at 400-700 nm in steps of 25 nm.

#define N_JERLOV 10
#define N_JERLOV_WL 13

static const double jerlov_kd_table[N_JERLOV][N_JERLOV_WL] = {
  {0.0217, 0.0185, 0.0176, 0.0184, 0.0280, 0.0504, 0.0640, 0.0931, 0.2408, 0.3174, 0.3559, 0.4372, 0.6513},
  {0.0316, 0.0280, 0.0257, 0.0250, 0.0332, 0.0545, 0.0674, 0.0960, 0.2437, 0.3206, 0.3601, 0.4410, 0.6530},
  {0.0438, 0.0395, 0.0355, 0.0330, 0.0396, 0.0596, 0.0715, 0.0995, 0.2471, 0.3245, 0.3652, 0.4457, 0.6550},
  {0.0878, 0.0814, 0.0714, 0.0620, 0.0627, 0.0779, 0.0863, 0.1122, 0.2595, 0.3389, 0.3837, 0.4626, 0.6623},
  {0.1697, 0.1594, 0.1381, 0.1160, 0.1056, 0.1120, 0.1139, 0.1359, 0.2826, 0.3655, 0.4181, 0.4942, 0.6760},
  {0.2516, 0.2374, 0.2048, 0.1700, 0.1486, 0.1461, 0.1415, 0.1596, 0.3057, 0.3922, 0.4525, 0.5257, 0.6896},
  {0.78, 0.54, 0.39, 0.29, 0.22, 0.2, 0.19, 0.21, 0.33, 0.4, 0.46, 0.56, 0.71},
  {1.1, 0.78, 0.56, 0.43, 0.36, 0.31, 0.3, 0.33, 0.4, 0.48, 0.54, 0.65, 0.8},
  {1.6, 1.2, 0.89, 0.71, 0.58, 0.49, 0.46, 0.46, 0.48, 0.54, 0.63, 0.78, 0.92},
  {2.4, 1.9, 1.6, 1.23, 0.99, 0.78, 0.63, 0.58, 0.6, 0.65, 0.76, 0.92, 1.1}};

// Kd at wavelength lam for a continuous type index t (0 = I ... 9 = 9C): linear in wavelength
// (clamped at the ends, as numpy's interp), log-linear between types.

static double jerlov_kd(double lam, double t) {
  int i, w;
  double f, g, k0, k1;

  t = MIN(MAX(t, 0.0), N_JERLOV - 1.0);
  i = MIN((int) t, N_JERLOV - 2);
  f = t - i;
  g = (lam - 400.0)/25.0;
  w = (int) floor(g);
  if (w < 0) {
    w = 0;
    g = 0.0;
  } else if (w >= N_JERLOV_WL - 1) {
    w = N_JERLOV_WL - 2;
    g = 1.0;
  } else g -= w;
  k0 = jerlov_kd_table[i][w] + g*(jerlov_kd_table[i][w + 1] - jerlov_kd_table[i][w]);
  k1 = jerlov_kd_table[i + 1][w] + g*(jerlov_kd_table[i + 1][w + 1] - jerlov_kd_table[i + 1][w]);
  return exp((1.0 - f)*log(k0) + f*log(k1));
}

// kappa and u of a reference, with the scene's attenuation calibration if any.

static void qaa(profile_data *pd, int s, reference *r) {
  int b;
  double u[MAX_BANDS] = {0.0}, shifted[MAX_BANDS] = {0.0};

  if (r->cal_mode == 1) {
    for (b = 0; b < pd->n_bands[s]; b++) shifted[b] = MAX(r->rrs_inf[b] - r->cal_delta, 1.0e-6);
    qaa_spectrum(pd, s, shifted, r->kappa, u);
  } else {
    qaa_spectrum(pd, s, r->rrs_inf, r->kappa, u);
    if (r->cal_mode == 2) {
      for (b = 0; b < pd->n_bands[s]; b++) r->kappa[b] = jerlov_kd(pd->lam[s][b], r->cal_t);
    }
  }
  for (b = 0; b < pd->n_bands[s]; b++) {
    r->MC[b] = pd->sec_w[s] + 1.03*sqrt(1.0 + 2.4*u[b])*pd->sec_v[s];
    r->MB[b] = pd->sec_w[s] + 1.04*sqrt(1.0 + 5.4*u[b])*pd->sec_v[s];
  }
}

// K_i/K_j from the brightest pixels line (Lyzenga; 4SM): for one bright bottom at varying depth,
// ln c_i against ln c_j, with c = rrs - rrs_inf the bottom contrast, is a line of slope K_i/K_j.
// The line is fitted to the upper envelope: the 95th percentile of ln c_i in 14 bins of ln c_j
// between its 5th and 95th percentiles, using bins of at least 50 cells. Returns NAN if there are
// too few cells or bins. A flat offset in rrs cancels in c, unlike in QAA's absolute rrs_inf.

#define BPL_BINS 14
#define BPL_MIN 50

static double bpl_ratio(profile_data *pd, int s, int i0, int i1, int j0, int j1, const double *rrs_inf,
    int bi, int bj) {
  int i, j, k, n = 0, nb = 0, m;
  double *xi, *xj, *v, edges[BPL_BINS + 1], ex[BPL_BINS], ey[BPL_BINS], ci, cj, vi, vj, mx, my, sxx, sxy;

  xi = (double*) malloc((size_t) pd->nrows*pd->ncols*sizeof(double));
  xj = (double*) malloc((size_t) pd->nrows*pd->ncols*sizeof(double));
  v = (double*) malloc((size_t) pd->nrows*pd->ncols*sizeof(double));
  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {
      // Only the two bands need data here (red is often below the noise in clear water).
      vi = pd->R[s][bi][i][j];
      vj = pd->R[s][bj][i][j];
      if (! isfinite(vi) || ! isfinite(vj) || approx_equal(vi, pd->nodata[s][bi], 1.0e-6) ||
          approx_equal(vj, pd->nodata[s][bj], 1.0e-6)) continue ;
      ci = rrs_of(vi) - rrs_inf[bi];
      cj = rrs_of(vj) - rrs_inf[bj];
      if (ci > 0.0 && cj > 0.0) {
        xi[n] = log(ci);
        xj[n] = log(cj);
        n++;
      }
    }
  }
  if (n < BPL_BINS*BPL_MIN) {
    free(xi);
    free(xj);
    free(v);
    return NAN;
  }
  for (k = 0; k <= BPL_BINS; k++) {
    for (i = 0; i < n; i++) v[i] = xj[i];
    edges[k] = percentile(v, n, 0.05 + 0.9*k/BPL_BINS);
  }
  for (k = 0; k < BPL_BINS; k++) {
    m = 0;
    for (i = 0; i < n; i++) {
      if (xj[i] >= edges[k] && xj[i] < edges[k + 1]) v[m++] = xi[i];
    }
    if (m >= BPL_MIN) {
      ex[nb] = 0.5*(edges[k] + edges[k + 1]);
      ey[nb] = percentile(v, m, 0.95);
      nb++;
    }
  }
  free(xi);
  free(xj);
  free(v);
  if (nb < 4) return NAN;
  mx = my = sxx = sxy = 0.0;
  for (k = 0; k < nb; k++) {
    mx += ex[k]/nb;
    my += ey[k]/nb;
  }
  for (k = 0; k < nb; k++) {
    sxx += (ex[k] - mx)*(ex[k] - mx);
    sxy += (ex[k] - mx)*(ey[k] - my);
  }
  return sxy/sxx;
}

// The attenuation calibration (PROFILE_KAPPA 1 or 2) of reference r, from the cells in rows
// i0..i1 and columns j0..j1 and r's rrs_inf: the brightest pixels line gives K_443/K_560 and
// K_492/K_560; mode 1 finds the flat offset of rrs_inf for which QAA reproduces them, mode 2 the
// Jerlov type. Returns false, leaving r uncalibrated, if either ratio is missing or outside
// (0.05, 3).

static bool calibrate(profile_data *pd, int s, int i0, int i1, int j0, int j1, reference *r, bool verbose) {
  int b, k, best_k = 0, ig = nearest_band(pd, s, 560.0), blue[2];
  double ratio[2], cost, best = INFINITY, kappa[MAX_BANDS], u[MAX_BANDS], shifted[MAX_BANDS] = {0.0}, rmin, d, t;

  r->cal_mode = 0;
  if (PROFILE_KAPPA == 0) return false;
  blue[0] = nearest_band(pd, s, 443.0);
  blue[1] = nearest_band(pd, s, 492.0);
  for (k = 0; k < 2; k++) {
    ratio[k] = bpl_ratio(pd, s, i0, i1, j0, j1, r->rrs_inf, blue[k], ig);
    if (! isfinite(ratio[k]) || ratio[k] <= 0.05 || ratio[k] >= 3.0) return false;
  }
  // The flat offset of rrs_inf for which QAA reproduces the ratios: kappa in mode 1, and the
  // offset-aware forward model (PROFILE_OFFSET_MODEL) in both modes.
  rmin = INFINITY;
  for (b = 0; b < pd->n_bands[s]; b++) rmin = MIN(rmin, r->rrs_inf[b]);
  for (k = 0; k < 200; k++) {
    d = k*0.999*rmin/199.0;
    for (b = 0; b < pd->n_bands[s]; b++) shifted[b] = r->rrs_inf[b] - d;
    qaa_spectrum(pd, s, shifted, kappa, u);
    cost = 0.0;
    for (b = 0; b < 2; b++) {
      t = log(kappa[blue[b]]/kappa[ig]) - log(ratio[b]);
      cost += t*t;
    }
    if (cost < best) {
      best = cost;
      best_k = k;
    }
  }
  r->cal_mode = 1;
  r->cal_delta = best_k*0.999*rmin/199.0;
  if (PROFILE_KAPPA == 2) {
    best = INFINITY;
    best_k = 0;
    for (k = 0; k <= 900; k++) {
      t = k*(N_JERLOV - 1.0)/900.0;
      cost = 0.0;
      for (b = 0; b < 2; b++) {
        d = log(jerlov_kd(pd->lam[s][blue[b]], t)/jerlov_kd(pd->lam[s][ig], t)) - log(ratio[b]);
        cost += d*d;
      }
      if (cost < best) {
        best = cost;
        best_k = k;
      }
    }
    r->cal_mode = 2;
    r->cal_t = best_k*(N_JERLOV - 1.0)/900.0;
  }
  if (verbose) {
    printf("  calibration: K443/K560 %.3f, K492/K560 %.3f -> offset %.4f", ratio[0], ratio[1], r->cal_delta);
    if (r->cal_mode == 2) printf(", Jerlov type %.4f", r->cal_t);
    printf("\n");
  }
  return true;
}

// A scene's deep-water reference from the valid cells in rows i0..i1, columns j0..j1.

static void deep_reference(profile_data *pd, int s, int i0, int i1, int j0, int j1, reference *r) {
  int i, j, b, n = 0, n_cand = 0, n_deep = 0;
  double *v, *w, p_red, p_green, med;

  r->ok = false;
  r->cal_mode = 0;
  r->cal_delta = r->cal_t = 0.0;
  v = (double*) malloc((size_t) (i1 - i0 + 1)*(j1 - j0 + 1)*sizeof(double));
  w = (double*) malloc((size_t) (i1 - i0 + 1)*(j1 - j0 + 1)*sizeof(double));
  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {
      if (pd->valid[s][i][j]) v[n++] = pd->red[s][i][j];
    }
  }
  if (n == 0) {
    free(v);
    free(w);
    return ;
  }
  p_red = percentile(v, n, PROFILE_RED_FRACTION);
  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {
      if (pd->valid[s][i][j] && pd->red[s][i][j] <= p_red) v[n_cand++] = pd->green[s][i][j];
    }
  }
  p_green = percentile(v, n_cand, PROFILE_GREEN_FRACTION);
  for (b = 0; b < pd->n_bands[s]; b++) {
    n_deep = 0;
    for (i = i0; i <= i1; i++) {
      for (j = j0; j <= j1; j++) {
        if (pd->valid[s][i][j] && pd->red[s][i][j] <= p_red && pd->green[s][i][j] <= p_green) {
          v[n_deep++] = rrs_of(pd->R[s][b][i][j]);
        }
      }
    }
    med = percentile(v, n_deep, 0.5);
    for (i = 0; i < n_deep; i++) w[i] = fabs(v[i] - med);
    r->rrs_inf[b] = med;
    r->noise[b] = MAX(1.4826*percentile(w, n_deep, 0.5), 1.0e-5);
  }
  r->n_deep = n_deep;
  r->ok = n_deep > 0;
  free(v);
  free(w);
  if (r->ok) qaa(pd, s, r);
}

// The profile chi2 over the depth grid for one cell, and its decision. Returns 1 (a depth), 2
// (optically deep), 3 (ambiguous) or 0 (no data).

// Log band-ratio consistency (PROFILE_RATIO; EXPERIMENTS.md, "Ratio consistency in the fit"):
// for every pair of bands i < j of a scene, the chi2 of the log ratio of the bottom contrasts
// D_b = rrs_b - rrs_inf,b, modelled against measured,
//   [ln(D_i/D_j) modelled - ln(D_i/D_j) measured]^2/sigma_ln^2,
// sigma_ln^2 = (s_i/D_i)^2 + (s_j/D_j)^2, with s_b the band's error (1/sqrt(w)). res is the
// measured minus modelled rrs, so the modelled contrast is D - res. Contrasts are kept positive
// and smooth as 0.5 (D + sqrt(D^2 + 4 s^2)), so weak contrasts carry little weight. Every ratio
// is then matched by the one depth, bottom and water of the cell.

static double positive_contrast(double D, double s) {
  return 0.5*(D + sqrt(D*D + 4.0*s*s));
}

static double ratio_chi2(profile_data *pd, reference **ref, double rrs[][MAX_BANDS], double w[][MAX_BANDS],
    double res[][MAX_BANDS]) {
  int s, i, j;
  double lm[MAX_BANDS], lo[MAX_BANDS], rel[MAX_BANDS], sd, D, d, chi2 = 0.0;
  for (s = 0; s < pd->n_scenes; s++) {
    for (i = 0; i < pd->n_bands[s]; i++) {
      if (w[s][i] <= 0.0) continue;
      sd = 1.0/sqrt(w[s][i]);
      D = rrs[s][i] - ref[s]->rrs_inf[i];
      lo[i] = log(positive_contrast(D, sd));
      lm[i] = log(positive_contrast(D - res[s][i], sd));
      rel[i] = sd/positive_contrast(D, sd);
    }
    for (i = 0; i < pd->n_bands[s]; i++) {
      if (w[s][i] <= 0.0) continue;
      for (j = i + 1; j < pd->n_bands[s]; j++) {
        if (w[s][j] <= 0.0) continue;
        d = (lm[i] - lm[j]) - (lo[i] - lo[j]);
        chi2 += PROFILE_RATIO_WEIGHT*d*d/(rel[i]*rel[i] + rel[j]*rel[j]);
      }
    }
  }
  return chi2;
}

static int profile_cell(profile_data *pd, reference **ref, double rrs[][MAX_BANDS], char ok[][MAX_BANDS],
    double *H_out, double *sigma_out, double *misfit, int *class_out, double *b_out) {

  int s, b, h, c, k, n = 0, S = pd->n_scenes, best_c[N_H];
  double res[MAX_SCENES][MAX_BANDS], w[MAX_SCENES][MAX_BANDS], y[MAX_SCENES][MAX_BANDS], T[MAX_SCENES][MAX_BANDS], x[MAX_SCENES][MAX_BANDS],
    chi[N_H], best_b[N_H], sw[MAX_SCENES], xs[MAX_SCENES], ys[MAX_SCENES], gw, sig2, Z, sxx, sxy, bb, e,
    chi_deep, d, dl, scale, H, x0, x1, x2, y0, y1, y2, den, A, B, lo, hi, far_min, sum2;
  int n_rel;

  gw = 1.0/(4.0*PROFILE_GLINT_SIGMA*PROFILE_GLINT_SIGMA);
  for (s = 0; s < S; s++) {
    for (b = 0; b < pd->n_bands[s]; b++) {
      w[s][b] = 0.0;
      if (ok[s][b] && ref[s]->ok) {
        sig2 = ref[s]->noise[b]*ref[s]->noise[b] + PROFILE_EPS*PROFILE_EPS*rrs[s][b]*rrs[s][b];
        w[s][b] = 1.0/sig2;
        n++;
      }
    }
  }
  if (n < 3) return 0;

  for (h = 0; h < N_H; h++) {
    for (s = 0; s < S; s++) {
      Z = MAX(pd->H[h] + pd->tide[s], 0.0);
      for (b = 0; b < pd->n_bands[s]; b++) {
        // With PROFILE_OFFSET_MODEL the calibration's offset is not water-column light and is not
        // attenuated: rrs = delta + (rrs_inf - delta)(1 - exp(-kappa M_C Z)) + bottom.
        dl = PROFILE_OFFSET_MODEL && ref[s]->cal_mode ? ref[s]->cal_delta : 0.0;
        y[s][b] = rrs[s][b] - dl - (ref[s]->rrs_inf[b] - dl)*(1.0 - exp(-ref[s]->kappa[b]*ref[s]->MC[b]*Z));
        T[s][b] = exp(-ref[s]->kappa[b]*ref[s]->MB[b]*Z);
      }
    }
    chi[h] = INFINITY;
    for (c = 0; c < N_CLASSES; c++) {
      sxx = sxy = 0.0;
      for (s = 0; s < S; s++) {
        sw[s] = PROFILE_GLINT ? gw : 0.0;
        xs[s] = ys[s] = 0.0;
        for (b = 0; b < pd->n_bands[s]; b++) {
          x[s][b] = T[s][b]*pd->rho[c][s][b]/PI;
          sw[s] += w[s][b];
          xs[s] += w[s][b]*x[s][b];
          ys[s] += w[s][b]*y[s][b];
        }
        if (PROFILE_GLINT) {
          // Minimise sum w (y - b x - d_s)^2 + gw sum d_s^2: for fixed b, d_s is explicit; b
          // is then a 1-D quadratic, clipped to its range (exact for one bounded variable).
          xs[s] /= sw[s];
          ys[s] /= sw[s];
          for (b = 0; b < pd->n_bands[s]; b++) {
            sxx += w[s][b]*(x[s][b] - xs[s])*(x[s][b] - xs[s]);
            sxy += w[s][b]*(x[s][b] - xs[s])*(y[s][b] - ys[s]);
          }
          sxx += gw*xs[s]*xs[s];
          sxy += gw*xs[s]*ys[s];
        } else {
          for (b = 0; b < pd->n_bands[s]; b++) {
            sxx += w[s][b]*x[s][b]*x[s][b];
            sxy += w[s][b]*x[s][b]*y[s][b];
          }
        }
      }
      bb = MIN(MAX(sxy/MAX(sxx, 1.0e-30), pd->lo[c]), pd->hi[c]);
      e = 0.0;
      for (s = 0; s < S; s++) {
        d = PROFILE_GLINT ? ys[s] - bb*xs[s] : 0.0;
        for (b = 0; b < pd->n_bands[s]; b++) {
          res[s][b] = y[s][b] - bb*x[s][b] - d;
          e += w[s][b]*res[s][b]*res[s][b];
        }
        if (PROFILE_GLINT) e += gw*d*d;
      }
      if (PROFILE_RATIO) e += ratio_chi2(pd, ref, rrs, w, res);
      if (e < chi[h]) {
        chi[h] = e;
        best_c[h] = c;
        best_b[h] = bb;
      }
    }
  }

  // The optically deep model (no bottom).

  chi_deep = 0.0;
  for (s = 0; s < S; s++) {
    sw[s] = PROFILE_GLINT ? gw : 0.0;
    d = 0.0;
    for (b = 0; b < pd->n_bands[s]; b++) {
      sw[s] += w[s][b];
      d += w[s][b]*(rrs[s][b] - ref[s]->rrs_inf[b]);
    }
    d = PROFILE_GLINT ? d/sw[s] : 0.0;
    for (b = 0; b < pd->n_bands[s]; b++) {
      res[s][b] = rrs[s][b] - ref[s]->rrs_inf[b] - d;
      chi_deep += w[s][b]*res[s][b]*res[s][b];
    }
    if (PROFILE_GLINT) chi_deep += gw*d*d;
  }
  if (PROFILE_RATIO) chi_deep += ratio_chi2(pd, ref, rrs, w, res);

  // The best depth, refined by a parabola through the grid neighbours.

  k = 0;
  for (h = 1; h < N_H; h++) {
    if (chi[h] < chi[k]) k = h;
  }
  H = pd->H[k];
  if (k > 0 && k < N_H - 1) {
    x0 = pd->H[k - 1]; x1 = pd->H[k]; x2 = pd->H[k + 1];
    y0 = chi[k - 1]; y1 = chi[k]; y2 = chi[k + 1];
    den = (x0 - x1)*(x0 - x2)*(x1 - x2);
    A = (x2*(y1 - y0) + x1*(y0 - y2) + x0*(y2 - y1))/den;
    B = (x2*x2*(y0 - y1) + x1*x1*(y2 - y0) + x0*x0*(y1 - y2))/den;
    if (A > 0.0) H = MIN(MAX(-B/(2.0*A), x0), x2);
  }
  scale = MAX(1.0, chi[k]/MAX(n - 2, 1));
  lo = INFINITY;
  hi = -INFINITY;
  for (h = 0; h < N_H; h++) {
    if (chi[h] <= chi[k] + scale) {
      lo = MIN(lo, pd->H[h]);
      hi = MAX(hi, pd->H[h]);
    }
  }
  *sigma_out = MAX(0.5*(hi - lo), 0.5*(pd->H[MIN(k + 1, N_H - 1)] - pd->H[MAX(k - 1, 0)]));
  *H_out = H;
  *class_out = best_c[k];
  *b_out = best_b[k];

  // The rms relative misfit at the best grid depth (%).

  sum2 = 0.0;
  n_rel = 0;
  for (s = 0; s < S; s++) {
    Z = MAX(pd->H[k] + pd->tide[s], 0.0);
    sw[s] = PROFILE_GLINT ? gw : 0.0;
    d = 0.0;
    for (b = 0; b < pd->n_bands[s]; b++) {
      y[s][b] = rrs[s][b] - ref[s]->rrs_inf[b]*(1.0 - exp(-ref[s]->kappa[b]*ref[s]->MC[b]*Z))
        - best_b[k]*exp(-ref[s]->kappa[b]*ref[s]->MB[b]*Z)*pd->rho[best_c[k]][s][b]/PI;
      sw[s] += w[s][b];
      d += w[s][b]*y[s][b];
    }
    d = PROFILE_GLINT ? d/sw[s] : 0.0;
    for (b = 0; b < pd->n_bands[s]; b++) {
      if (w[s][b] > 0.0) {
        sum2 += (y[s][b] - d)*(y[s][b] - d)/(rrs[s][b]*rrs[s][b]);
        n_rel++;
      }
    }
  }
  *misfit = 100.0*sqrt(sum2/MAX(n_rel, 1));

  if (chi_deep - chi[k] <= PROFILE_T_DEEP*scale) return 2;
  if (PROFILE_T_AMB > 0.0) {
    far_min = INFINITY;
    for (h = 0; h < N_H; h++) {
      if (fabs(pd->H[h] - H) > 0.5 + 0.1*H) far_min = MIN(far_min, chi[h]);
    }
    if (far_min - chi[k] <= PROFILE_T_AMB*scale) return 3;
  }
  return 1;
}

// Bottom brightness ranges (PROFILE_BRIGHT) from the image, after 4SM's brightest bottom at null
// depth: near the surface the bottom term is b rho/pi, so the brightest bottoms have
// b ~ pi (rrs - rrs_inf)/rho at 560 nm. The cap is 1.2 x the largest scene's 99.5th percentile of
// that. Mode 1 caps each class's bounds; mode 2 scales every range by cap/0.75 (if below 1).

static void bright_bounds(profile_data *pd, const reference *global) {
  int s, i, j, c, n, ig;
  double *v, rho, cap = -INFINITY, f, p;

  for (c = 0; c < N_CLASSES; c++) {
    pd->lo[c] = class_lo[c];
    pd->hi[c] = class_hi[c];
  }
  if (PROFILE_BRIGHT == 0) return ;
  v = (double*) malloc((size_t) pd->nrows*pd->ncols*sizeof(double));
  for (s = 0; s < pd->n_scenes; s++) {
    if (! global[s].ok) continue ;
    ig = nearest_band(pd, s, 560.0);
    rho = pd->rho[0][s][ig];
    n = 0;
    for (i = 0; i < pd->nrows; i++) {
      for (j = 0; j < pd->ncols; j++) {
        double r = pd->R[s][ig][i][j];
        if (! isfinite(r) || approx_equal(r, pd->nodata[s][ig], 1.0e-6)) continue ;
        v[n] = PI*(rrs_of(r) - global[s].rrs_inf[ig])/rho;
        if (isfinite(v[n])) n++;
      }
    }
    if (n > 0) {
      p = percentile(v, n, 0.995);
      cap = MAX(cap, p);
    }
  }
  free(v);
  if (! isfinite(cap)) return ;
  cap *= 1.2;
  f = MIN(cap/class_hi[0], 1.0);
  for (c = 0; c < N_CLASSES; c++) {
    if (PROFILE_BRIGHT == 1) {
      pd->lo[c] = MIN(class_lo[c], cap);
      pd->hi[c] = MIN(class_hi[c], cap);
    } else {
      pd->lo[c] = class_lo[c]*f;
      pd->hi[c] = class_hi[c]*f;
    }
  }
  printf("  brightness cap %.3f: sand %.3f-%.3f, seagrass %.3f-%.3f, coral %.3f-%.3f\n", cap, pd->lo[0], pd->hi[0],
    pd->lo[1], pd->hi[1], pd->lo[2], pd->hi[2]);
}

void profile_model(
  // inputs
  scene scene_data[], geogrid gridded_data[], int *scene_indexes, int nscenes,
  bool points_present, geogrid points,
  // outputs
  float **depth, float **depth_sigma, float **model_error, float **bottom_albedo,
  float **bottom_sand, float **bottom_seagrass, float **bottom_coral,
  float **K_min, float **bottom_type, float **index_optical_depth) {

  int i, j, s, b, c, h, nr, nc, rad, n_ti, n_tj, ti, tj, ci, cj, n_valid, k_index, n, cls,
    counts[4] = {0, 0, 0, 0};
  double cell_km, H, sigma, misfit, bright, kmin, rrs[MAX_SCENES][MAX_BANDS], v;
  char ok[MAX_SCENES][MAX_BANDS], file_name[64];
  float *lons, *lats;
  reference global[MAX_SCENES], **tiles, *ref[MAX_SCENES];
  profile_data *pd;
  double *table[N_CLASSES] = {bottom_type_sand, bottom_type_seagrass, bottom_type_coral};
  clock_t t0 = clock();

  pd = (profile_data*) calloc(1, sizeof(profile_data));
  nr = pd->nrows = scene_data[scene_indexes[0]].nrows;
  nc = pd->ncols = scene_data[scene_indexes[0]].ncols;
  pd->n_scenes = nscenes;
  for (h = 0; h < N_H; h++) {
    // As validation/profile_model.py's H_GRID: 0.05-2 m (40), 2.05-10 m (80), 10.2-30 m (100).
    pd->H[h] = h < 40 ? 0.05 + h*(2.0 - 0.05)/39.0 :
      h < 120 ? 2.05 + (h - 40)*(10.0 - 2.05)/79.0 : 10.2 + (h - 120)*(30.0 - 10.2)/99.0;
  }

  printf("\nLee_Profile: radius %g km, tile %d, eps %g, t_deep %g, t_amb %g, glint %d, kappa %d (local %d)\n",
    (double) PROFILE_RADIUS, PROFILE_TILE, (double) PROFILE_EPS, (double) PROFILE_T_DEEP,
    (double) PROFILE_T_AMB, PROFILE_GLINT, PROFILE_KAPPA, PROFILE_KAPPA_LOCAL);

  for (s = 0; s < nscenes; s++) {
    scene *sc = &scene_data[scene_indexes[s]];
    pd->n_bands[s] = MIN(sc->n_bands, MAX_BANDS);
    pd->tide[s] = sc->H_tide;
    pd->sec_w[s] = 1.0/cos(sc->theta_w*PI/180.0);
    pd->sec_v[s] = 1.0/cos(sc->theta_v*PI/180.0);
    for (b = 0; b < pd->n_bands[s]; b++) {
      pd->lam[s][b] = sc->wavelengths[b];
      pd->R[s][b] = gridded_data[sc->band_indexes[b]].array;
      pd->nodata[s][b] = gridded_data[sc->band_indexes[b]].nodata_value;
      pd->aw[s][b] = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, pd->lam[s][b]);
      pd->bbw[s][b] = interp_1d(ref_wlens, bbw_Morel_1974, n_ref_wlens, pd->lam[s][b]);
      for (c = 0; c < N_CLASSES; c++) {
        pd->rho[c][s][b] = interp_1d(ref_wlens, table[c], n_ref_wlens, pd->lam[s][b]);
      }
    }
    smooth_fields(pd, s);
    deep_reference(pd, s, 0, nr - 1, 0, nc - 1, &global[s]);
    printf("  %s: %d deep cells, rrs_inf", sc->scene_name, global[s].n_deep);
    for (b = 0; b < pd->n_bands[s]; b++) printf(" %.5f", global[s].rrs_inf[b]);
    printf("\n");
    if (global[s].ok) {
      calibrate(pd, s, 0, nr - 1, 0, nc - 1, &global[s], true);
      qaa(pd, s, &global[s]);
    }
  }

  bright_bounds(pd, global);

  // Local references, one per tile, computed when first needed.

  k_index = scene_data[scene_indexes[0]].band_indexes[0];
  cell_km = gridded_data[k_index].cellsize*111.0;
  rad = PROFILE_RADIUS > 0.0 ? (int) lround(PROFILE_RADIUS/cell_km) : 0;
  n_ti = (nr + PROFILE_TILE - 1)/PROFILE_TILE;
  n_tj = (nc + PROFILE_TILE - 1)/PROFILE_TILE;
  tiles = (reference**) calloc((size_t) n_ti*n_tj, sizeof(reference*));

  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      depth[i][j] = depth_sigma[i][j] = model_error[i][j] = bottom_albedo[i][j] = 0.0;
      bottom_sand[i][j] = bottom_seagrass[i][j] = bottom_coral[i][j] = 0.0;
      K_min[i][j] = bottom_type[i][j] = index_optical_depth[i][j] = 0.0;

      if (points_present && (points.array[i][j] == 0.0 ||
          approx_equal(points.array[i][j], points.nodata_value, 1.0e-6))) continue ;

      n = 0;
      for (s = 0; s < nscenes; s++) {
        for (b = 0; b < pd->n_bands[s]; b++) {
          v = pd->R[s][b][i][j];
          ok[s][b] = isfinite(v) && ! approx_equal(v, pd->nodata[s][b], 1.0e-6) && v > 0.0;
          rrs[s][b] = ok[s][b] ? rrs_of(v) : 0.0;
          n += ok[s][b];
        }
      }
      if (n < 3) {
        counts[0]++;
        continue ;
      }

      for (s = 0; s < nscenes; s++) ref[s] = &global[s];
      if (rad > 0) {
        ti = i/PROFILE_TILE;
        tj = j/PROFILE_TILE;
        if (tiles[ti*n_tj + tj] == NULL) {
          tiles[ti*n_tj + tj] = (reference*) calloc(nscenes, sizeof(reference));
          ci = MIN(ti*PROFILE_TILE + PROFILE_TILE/2, nr - 1);
          cj = MIN(tj*PROFILE_TILE + PROFILE_TILE/2, nc - 1);
          for (s = 0; s < nscenes; s++) {
            n_valid = 0;
            for (int ii = MAX(ci - rad, 0); ii <= MIN(ci + rad, nr - 1); ii++) {
              for (int jj = MAX(cj - rad, 0); jj <= MIN(cj + rad, nc - 1); jj++) n_valid += pd->valid[s][ii][jj];
            }
            if (n_valid >= PROFILE_MIN_LOCAL) {
              reference *t = &tiles[ti*n_tj + tj][s];
              deep_reference(pd, s, MAX(ci - rad, 0), MIN(ci + rad, nr - 1), MAX(cj - rad, 0),
                MIN(cj + rad, nc - 1), t);
              // The scene's site-wide calibration, or (PROFILE_KAPPA_LOCAL) the window's own.
              if (t->ok && (! PROFILE_KAPPA_LOCAL || ! calibrate(pd, s, MAX(ci - rad, 0), MIN(ci + rad, nr - 1),
                  MAX(cj - rad, 0), MIN(cj + rad, nc - 1), t, false))) {
                t->cal_mode = global[s].cal_mode;
                t->cal_delta = global[s].cal_delta;
                t->cal_t = global[s].cal_t;
              }
              if (t->ok) qaa(pd, s, t);
            }
            if (! tiles[ti*n_tj + tj][s].ok) tiles[ti*n_tj + tj][s] = global[s];
          }
        }
        for (s = 0; s < nscenes; s++) ref[s] = &tiles[ti*n_tj + tj][s];
      }

      cls = profile_cell(pd, ref, rrs, ok, &H, &sigma, &misfit, &c, &bright);
      counts[cls]++;
      index_optical_depth[i][j] = cls;
      if (cls == 0) continue ;
      model_error[i][j] = misfit;
      kmin = INFINITY;
      for (s = 0; s < nscenes; s++) {
        for (b = 0; b < pd->n_bands[s]; b++) kmin = MIN(kmin, ref[s]->kappa[b]);
      }
      K_min[i][j] = kmin;
      if (cls == 1) {
        depth[i][j] = -H;  // negative down
        depth_sigma[i][j] = sigma;
        bottom_albedo[i][j] = bright;
        bottom_type[i][j] = c;
        bottom_sand[i][j] = c == 0 ? 100.0 : 0.0;
        bottom_seagrass[i][j] = c == 1 ? 100.0 : 0.0;
        bottom_coral[i][j] = c == 2 ? 100.0 : 0.0;
      }
    }
  }

  printf("Lee_Profile: %d cells with a depth, %d optically deep, %d ambiguous, %d without data (%.1f s)\n",
    counts[1], counts[2], counts[3], counts[0], (double) (clock() - t0)/CLOCKS_PER_SEC);

  // Write the outputs as MODEL Lee_Semi_Analytical does.

  lons = (float*) malloc(nc*sizeof(float));
  lats = (float*) malloc(nr*sizeof(float));
  for (j = 0; j < nc; j++) lons[j] = gridded_data[k_index].wlon + ((float) j)*gridded_data[k_index].cellsize;
  for (i = 0; i < nr; i++) lats[i] = gridded_data[k_index].slat + ((float) i)*gridded_data[k_index].cellsize;
  {
    const char *names[] = {"H", "H_sigma", "error", "albedo", "bottom_sand", "bottom_seagrass",
      "bottom_coral", "min_K", "bottom_type", "index_optical_depth"};
    float **grids[] = {depth, depth_sigma, model_error, bottom_albedo, bottom_sand, bottom_seagrass,
      bottom_coral, K_min, bottom_type, index_optical_depth};
    for (n = 0; n < 10; n++) {
      sprintf(file_name, "modelled_%s.nc", names[n]);
      write_nc(file_name, grids[n], nc, nr, lons, lats, gridded_data[k_index].nodata_value);
    }
  }
  free(lons);
  free(lats);

  for (n = 0; n < n_ti*n_tj; n++) free(tiles[n]);
  free(tiles);
  for (s = 0; s < nscenes; s++) {
    for (i = 0; i < nr; i++) free(pd->valid[s][i]);
    free(pd->valid[s]);
    free_double_array_2d(pd->red[s], nr);
    free_double_array_2d(pd->green[s], nr);
  }
  free(pd);
}
