//
//  MODEL lyzenga_log_linear
//

// Lyzenga's log-linear depth model (Lyzenga 1978, 1985; Lyzenga et al. 2006), fitted to soundings:
//
//   H = a_0 + sum_m a_m X_m,   X_m = ln(L_m - Lsm_m),
//
// where L_m is band m (any units; Rrs from ACOLITE, or radiance) averaged over a (2r-1)-wide
// window of valid cells, and Lsm_m is the band's deep-water signal. The model is linear in its
// coefficients, so the fit is a weighted least-squares solve (Householder QR): exact, fast and
// deterministic. It is made robust to outlying soundings by iteratively reweighted least squares
// with Huber weights (scale from the median absolute deviation).
//
// Cells are modelled where every band is valid (and, with LAND/SHALLOW, off land and in the
// shallow mask). A cell is optically deep, and gets no depth, where L_m - Lsm_m is below
// DEEP x sigma_m in any band; sigma_m is the robust spread of the darkest cells (those below
// the PERCENTILE that sets Lsm by default). Soundings go to the nearest cell (node-registered
// grids: lat = slat + i*cellsize). Only those in modelled, shallow cells with a depth in the
// fitting range enter the fit, and the counts of those dropped are printed.
//
// A least-squares fit of depth on the bands shrinks its depths towards the mean (regression
// dilution: the bands are noisy), so shallow water comes out too deep and deep water too shallow.
// POWER corrects the scale by classical calibration: the fit's depths h at the soundings are
// modelled as a function of the soundings' depths z, h = c z^d (robustly, in logs), and inverted,
// z = (h/c)^(1/d). The mapped depths are then unbiased against depth at every depth (slope 1),
// at the cost of more scatter; where the fit has saturated it cannot recover depth. The law
// crosses the identity at h_x = c^(1/(1-d)); by default it is applied only beyond the crossover
// (where the fit is compressed into too little depth), so the shallow depths keep the fit's
// smaller scatter (Arabian Gulf, Landsat 8 against ICESat-2: hits 64 -> 69% in cross-validation,
// against 58% when applied everywhere). POWER ALL applies it everywhere.
//
// Depth conventions: soundings are depths below the datum, positive down (a file whose median is
// negative is taken as negative down and negated). TIDE is the water level above the datum at
// the imagery's time: the fit uses depth + TIDE (the water the light saw) and the output is
// depth below the datum, negative down, as photic's other bathymetry.
//
// MODEL lyzenga_log_linear GRIDS /band 1/ ... /band n/ BATHYMETRY /output/ [optional, any order]
//   LAND /grid/                  cells where the grid has data are land
//   SHALLOW /grid/               only cells where the grid has data are modelled
//   LSM /Lsm 1/ ... /Lsm n/      deep-water signal per band (0 gives ln L; Lyzenga et al. 2006 subtract one)
//   PERCENTILE /p/               else Lsm is the p-th percentile of the window means (default 1)
//   DEEP /k/                     optically deep below k sigma above Lsm (default 2; 0: off)
//   RADIUS /r/                   window radius (default SET DEPTHRADIUS)
//   DEPTHS /min/ /max/           soundings used for the fit (default SET DEPTHMIN, DEPTHMAX)
//   TIDE /h/                     water level above the datum at the imagery's time (m)
//   HUBER /k/                    Huber threshold in robust sigmas (default 1.345; 0: least squares)
//   BALANCE /bin/                weight soundings by 1/count in depth bins of this width (m)
//   POWER [ALL]                  rescale the depths by the inverted power law above
//   ERROR /grid/ [NOISE] [GAIN]  output: the 1-sigma prediction error (m); see below
//
// ERROR is the prediction error of a new depth, s^2 + x'Cx: the fit's (robust) residual variance
// and the coefficients' covariance. NOISE adds the cell's own radiometric noise (Lyzenga 1978,
// eq. 16), sum a_m^2 sigma_m^2/(L_m - Lsm_m)^2; the soundings' mean of it is already in s^2, so it
// is taken out of s^2 first. GAIN (with ATTENUATION) divides by the local gain of the depth on the
// true depth, as Lyzenga 1978's eq. 17 divides by C: the bands' attenuation, -sum a_m g_m(x) f_m,
// where f_m = (L_m - L_inf,m)/(L_m - Lsm_m) is 1 for ln(L - L_inf) and falls towards deep water
// for ln L. Where the fit saturates the gain falls, and the error grows. The gain is capped at 1,
// so the shallow water keeps the fit's own scatter, and floored at 0.05; with POWER the larger of
// its rescaling and 1/gain is used, as both undo the same compression.
//
// Outputs model.bam (the fitted model as applied) and output.csv (every sounding, its cell, the
// model depth and whether it was used).
//
// INDICES (with ATTENUATION) adds Lyzenga 1978's N - 1 depth-invariant bottom indices of each
// scene of N bands, /prefix/_index/i/_s/scene/, computed with the local attenuation (see below).
//
// ATTENUATION /prefix/ also estimates the attenuation over the whole grid, with a physical depth
// from it. Lyzenga's model is the linearised form of
//
//   X_m = ln(L_m - Lsm_m) = B(x) + c_m - g_m(x) Z_s,   Z_s = depth + h_s,
//
// for grid m of scene s (the grids are scene-major: SCENES /n/ groups them, n bands each, in the
// same band order), with B the bottom's brightness (spectrally flat in logs, Lyzenga's assumption),
// c_m the band's bottom colour and deep-water constant, and g_m = K_m M_s the two-way attenuation:
// K_m the effective diffuse attenuation and M_s = sec(theta_w,sun) + sec(theta_w,view) the path
// in water (from SUNZENITH and VIEWZENITH in air, refracted with n = 1.34; default 0, so M_s = 2).
//   Stage 1, at the soundings (depth known): for each grid, X_m = c_m - g_m(x) Z_s with g_m(x)
// bilinear on a grid of nodes NODES degrees apart (default 0.05), smoothed by a penalty on
// neighbouring nodes' differences (SMOOTHING /w/, default 1, scaled to an average node's data;
// 0 needs data at every node), fitted by weighted least squares with Huber weights. Weights are
// 1/(sigma_m^2/(L_m - Lsm_m)^2 + 0.05^2), so a band counts only while it sees the bottom. B is
// left in the residuals, so K is unbiased only where the bottom's brightness does not change with
// depth across the soundings.
//   Stage 2, at every cell: with c_m and g_m(x) fixed, X_m - c_m + g_m h_s = B - g_m H over the
// grids above DEEP sigma is a straight line in g_m, and the weighted fit gives the depth H (from
// how the attenuation differs between bands and scenes, which is what makes it independent of
// B), its 1-sigma error and B. Cells with fewer than 3 such grids get none.
// Lsm is the LSM value where it is > 0, else the PERCENTILE (so LSM 0, the plain ln L, still
// subtracts a deep-water signal here). TIDES /h_1/ ... /h_n/ are the scenes' water levels above
// the datum (default TIDE for each). Outputs: /prefix/_/grid/, K of each grid (m^-1);
// /prefix/_depth (below the datum, negative down), /prefix/_depth_sigma and /prefix/_bottom (B);
// and attenuation.csv, K at every node with the number of soundings near it.
//
//   ATTENUATION /prefix/ SCENES /n/ [TIDES /h_1/ ... /h_n/] [SUNZENITH /theta_1/ ... /theta_n/]
//     [VIEWZENITH /theta_1/ ... /theta_n/] [NODES /deg/] [SMOOTHING /w/]   (SCENES before the
//     per-scene lists)
//
// CONSISTENT (with ATTENUATION, two or more scenes) asks for one bottom in every scene: at every
// cell, each band's bottom beta_b (any colour) is shared by the scenes, so only how the attenuation
// differs between the scenes carries the depth. It refits those differences, gamma_sb(x) = g_sb -
// gbar_b, and the scenes' constants over the image's cells, not only the soundings (see "The
// scenes' differences" below); gbar_b, the band's mean over the scenes, stays from stage 1, as the
// images cannot separate it from a bottom that changes with depth. Stage 2 then fits one bottom per
// band at every cell; /prefix/_bottom is their mean, and with INDICES /prefix/_index/i/ are the one
// bottom's N - 1 indices (Lyzenga 1978's rotation of beta_b with gbar_b). PRIOR /w/ pulls the fit
// towards stage 1's values, with w times the weight of an average parameter's data (default 0).
// At the Arabian Gulf it did not beat stage 2's flat bottom: depth explains only about half of
// the scenes' differences there, the rest being atmosphere and adjacency (validation/EXPERIMENTS.md).
//   CONSISTENT [PRIOR /w/]
//
// With ATTENUATION, stage 1 runs before Lyzenga's fit, and two options use it (see "The
// attenuation at the soundings" below): CONSTRAIN GAIN [BRIGHTNESS] imposes Lyzenga et al. 2006's
// calibration condition (eq. 12) on the fit, and NORMALISE divides each predictor by the local
// two-way attenuation. Both want a deep-water signal (PERCENTILE or LSM > 0), not ln L. At the
// Arabian Gulf neither beat the free fit on ln L with POWER (validation/EXPERIMENTS.md).

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lyzenga.h"

#define LYZ_MAX_BANDS 32  // e.g. 3 bands of each of 6 scenes, fitted jointly
#define LYZ_MAX_ITER 100
#define ATT_EPS 0.05      // the attenuation fit's model error in ln(L - Lsm)
#define ATT_MIN_ROWS 10   // fewest soundings for a grid's attenuation
#define ATT_MAX_PASSES 10 // reselections by the fitted signal
#define N_WATER 1.34

static int compare_doubles(const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

// The p-th percentile (0-100) of v[0..n-1], which is sorted in place.
static float percentile(float *v, int n, float p) {
  float r;
  int k;
  qsort(v, n, sizeof(float), compare_floats);
  r = p/100.0*(n - 1);
  k = (int) floor(r);
  if (k >= n - 1) return v[n - 1];
  return v[k] + (r - k)*(v[k + 1] - v[k]);
}

static double median_double(double *v, int n) {
  qsort(v, n, sizeof(double), compare_doubles);
  return n%2 ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
}

int lyzenga_wls(const double *A, const double *y, const double *w, int n, int p, double *x, double *R) {

  int i, j, k;
  double *Q, *b, norm, alpha, s, v0, scale = 0.0;

  Q = malloc((size_t) n*p*sizeof(double));
  b = malloc((size_t) n*sizeof(double));
  for (i = 0; i < n; i++) {
    s = sqrt(w[i]);
    for (j = 0; j < p; j++) Q[i*p + j] = s*A[i*p + j];
    b[i] = s*y[i];
  }

  // Householder reflections, applied to the columns and to b as they are formed.
  for (k = 0; k < p; k++) {
    norm = 0.0;
    for (i = k; i < n; i++) norm += Q[i*p + k]*Q[i*p + k];
    norm = sqrt(norm);
    scale = MAX(scale, norm);
    if (norm <= 1.0e-12*scale || norm == 0.0) {
      free(Q);
      free(b);
      return 1;
    }
    alpha = Q[k*p + k] > 0.0 ? -norm : norm;
    v0 = Q[k*p + k] - alpha;
    Q[k*p + k] = v0;  // the reflector v = (v0, Q[k+1..n-1][k]), with v'v = 2 norm |v0|
    for (j = k + 1; j < p; j++) {
      s = 0.0;
      for (i = k; i < n; i++) s += Q[i*p + k]*Q[i*p + j];
      s /= -alpha*v0;
      for (i = k; i < n; i++) Q[i*p + j] -= s*Q[i*p + k];
    }
    s = 0.0;
    for (i = k; i < n; i++) s += Q[i*p + k]*b[i];
    s /= -alpha*v0;
    for (i = k; i < n; i++) b[i] -= s*Q[i*p + k];
    Q[k*p + k] = alpha;
    for (i = k + 1; i < n; i++) Q[i*p + k] = 0.0;
  }

  // Back substitution.
  for (k = p - 1; k >= 0; k--) {
    s = b[k];
    for (j = k + 1; j < p; j++) s -= Q[k*p + j]*x[j];
    x[k] = s/Q[k*p + k];
  }
  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) R[i*p + j] = Q[i*p + j];
  }
  free(Q);
  free(b);
  return 0;
}

// The covariance factor (R'R)^-1 from the upper-triangular R (p x p).
static void covariance_factor(const double *R, int p, double *C) {

  int i, j, k;
  double Ri[LYZ_MAX_BANDS + 1][LYZ_MAX_BANDS + 1], s;

  // Ri = R^-1 (upper triangular).
  for (j = 0; j < p; j++) {
    for (i = 0; i < p; i++) Ri[i][j] = 0.0;
    Ri[j][j] = 1.0/R[j*p + j];
    for (i = j - 1; i >= 0; i--) {
      s = 0.0;
      for (k = i + 1; k <= j; k++) s += R[i*p + k]*Ri[k][j];
      Ri[i][j] = -s/R[i*p + i];
    }
  }
  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) {
      s = 0.0;
      for (k = MAX(i, j); k < p; k++) s += Ri[i][k]*Ri[j][k];
      C[i*p + j] = s;
    }
  }
}

// Iteratively reweighted least squares with Huber weights (scale from the median absolute
// deviation), from the least-squares solution: coef, R, the final weights w and residuals res.
// Returns the number of iterations, or -1 if the design is rank deficient.
static int huber_irls(const double *A, const double *y, const double *w0, int n, int p, double huber,
    double *coef, double *R, double *w, double *res) {

  int i, m, it;
  double prev[LYZ_MAX_BANDS + 1], s, change, *absr = malloc((size_t) n*sizeof(double));

  for (m = 0; m < p; m++) coef[m] = 0.0;
  for (i = 0; i < n; i++) w[i] = w0[i];
  for (it = 0; it < LYZ_MAX_ITER; it++) {
    for (m = 0; m < p; m++) prev[m] = coef[m];
    if (lyzenga_wls(A, y, w, n, p, coef, R) != 0) {
      free(absr);
      return -1;
    }
    for (i = 0; i < n; i++) {
      s = coef[0];
      for (m = 1; m < p; m++) s += coef[m]*A[i*p + m];
      res[i] = y[i] - s;
    }
    if (huber == 0.0) break;
    for (i = 0; i < n; i++) absr[i] = fabs(res[i]);
    s = 1.4826*median_double(absr, n);
    if (s <= 0.0) break;
    for (i = 0; i < n; i++) w[i] = w0[i]*(fabs(res[i]) <= huber*s ? 1.0 : huber*s/fabs(res[i]));
    if (it > 0) {
      change = 0.0;
      for (m = 0; m < p; m++) change = MAX(change, fabs(coef[m] - prev[m])/(1.0 + fabs(coef[m])));
      if (change < 1.0e-9) break;
    }
  }
  free(absr);
  return it + 1;
}

// The power-law rescaling of a fitted depth h (POWER): (h/c)^(1/d), with ln c = a; without ALL,
// only on the deep side of the crossover h_x = c^(1/(1-d)), and h shallower.
static double power_map(double h, double a, double d, bool all) {
  if (h <= 0.0) return h;
  if (! all && d != 1.0 && log(h) < a/(1.0 - d)) return h;
  return exp((log(h) - a)/d);
}

// The attenuation's nodes: nx by ny, spanning the grid, dx and dy cells apart.
typedef struct {
  int nx, ny;
  double dx, dy;
} att_nodes;

// The bilinear weights of cell (i, j) on its four nodes.
static void node_weights(const att_nodes *a, int i, int j, int *q, double *f) {
  double u = j/a->dx, v = i/a->dy, fu, fv;
  int b = MIN((int) floor(u), a->nx - 2), c = MIN((int) floor(v), a->ny - 2);
  fu = u - b;
  fv = v - c;
  q[0] = c*a->nx + b;
  q[1] = q[0] + 1;
  q[2] = q[0] + a->nx;
  q[3] = q[2] + 1;
  f[0] = (1.0 - fu)*(1.0 - fv);
  f[1] = fu*(1.0 - fv);
  f[2] = (1.0 - fu)*fv;
  f[3] = fu*fv;
}

// Solves N x = b for symmetric positive-definite N (n x n, row-major) by Cholesky; N is
// overwritten by its factor and b by x. Returns 1 if N is not positive definite.
static int cholesky_solve(double *N, double *b, int n) {

  int i, j, k;
  double s, scale = 0.0;

  for (j = 0; j < n; j++) scale = MAX(scale, N[j*n + j]);
  for (j = 0; j < n; j++) {
    s = N[j*n + j];
    for (k = 0; k < j; k++) s -= N[j*n + k]*N[j*n + k];
    if (s <= 1.0e-12*scale) return 1;
    N[j*n + j] = sqrt(s);
    for (i = j + 1; i < n; i++) {
      s = N[i*n + j];
      for (k = 0; k < j; k++) s -= N[i*n + k]*N[j*n + k];
      N[i*n + j] = s/N[j*n + j];
    }
  }
  for (i = 0; i < n; i++) {
    s = b[i];
    for (k = 0; k < i; k++) s -= N[i*n + k]*b[k];
    b[i] = s/N[i*n + i];
  }
  for (i = n - 1; i >= 0; i--) {
    s = b[i];
    for (k = i + 1; k < n; k++) s -= N[k*n + i]*b[k];
    b[i] = s/N[i*n + i];
  }
  return 0;
}

// The name of an attenuation output, /prefix/_/suffix/; false if it is too long or in use.
static bool att_name(char *name, const char *prefix, const char *suffix) {
  int n;
  if (strlen(prefix) + strlen(suffix) + 2 > MAX_ARG_SIZE ||
      snprintf(name, MAX_ARG_SIZE, "%s_%s", prefix, suffix) >= MAX_ARG_SIZE) {
    printf("\nERROR: ATTENUATION: the name '%s_%s' is too long.\n", prefix, suffix);
    return false;
  }
  for (n = 0; n < MAX_GRIDS; n++) {
    if (grid_name_taken(n, name)) {
      printf("\nERROR: grid '%s' already in use.\n", name);
      return false;
    }
  }
  return true;
}

// Stores an output grid under a name checked by att_name.
static bool store_grid(const char *name, float **array, const geogrid *g0) {
  int n;
  for (n = 0; n < MAX_GRIDS && allocated_grids[n]; n++) ;
  if (n == MAX_GRIDS) {
    printf("\nERROR: unable to store grid '%s'.\n", name);
    return false;
  }
  strcpy(grid_names[n], name);
  allocated_grids[n] = true;
  n_grids_in_use++;
  gridded_data[n] = *g0;
  gridded_data[n].array = array;
  meminuse += 1.0e-6*((float) g0->nrows*g0->ncols)*((float) sizeof(float));
  return true;
}

// ATTENUATION stage 1, from the soundings: each grid's c_m and two-way attenuation g_m(x) as a
// field on the nodes (see the header).
typedef struct {
  att_nodes nodes;
  int nq;
  double *theta, *c, *support;  // theta[m*nq + q]: g_m at node q
  bool *fitted;
  bool consistent;              // refitted by CONSISTENT: one bottom per band at every cell
} att_field;

static void att_free(att_field *af) {
  free(af->theta);
  free(af->c);
  free(af->support);
  free(af->fitted);
  af->theta = af->c = af->support = NULL;
  af->fitted = NULL;
}

// g_m at cell (i, j).
static double att_g(const att_field *af, int m, int i, int j) {
  int q[4], l;
  double f[4], g = 0.0;
  node_weights(&af->nodes, i, j, q, f);
  for (l = 0; l < 4; l++) g += f[l]*af->theta[m*af->nq + q[l]];
  return g;
}

static void att_stage1(att_field *af, const char *prefix, int nb, const int *bands, float ***Lbar,
    const geogrid *g0, const double *lsm, const double *sigma, double deep_k, int n_scenes, const double *tides,
    const double *path, double spacing, double smooth, double huber, double dmin, double dmax, double sign,
    const int *cell_i, const int *cell_j, const int *status) {

  int nr = g0->nrows, nc = g0->ncols, nbs = nb/n_scenes, m, s, n, i, j, l, q[4], nq, P, it, rows,
    *row_cell = NULL, *iters = NULL;
  double f[4], *N = NULL, *rhs = NULL, *theta = NULL, *c = NULL, *X = NULL, *Z = NULL, *w0 = NULL, *w = NULL,
    *res = NULL, *absr = NULL, *support = NULL, *scatter = NULL, lambda, sc, v, wzz, change, g, *Xc = NULL,
    *Zc = NULL, *vc = NULL, *vhat = NULL;
  bool *fitted = NULL, *active = NULL;
  int *cand_cell = NULL;
  att_nodes nodes;

  nodes.nx = MAX(2, (int) ceil((nc - 1)*g0->cellsize/spacing - 1.0e-9) + 1);
  nodes.ny = MAX(2, (int) ceil((nr - 1)*g0->cellsize/spacing - 1.0e-9) + 1);
  nodes.dx = MAX(nc - 1, 1)/((double) nodes.nx - 1);
  nodes.dy = MAX(nr - 1, 1)/((double) nodes.ny - 1);
  nq = nodes.nx*nodes.ny;
  P = nq + 1;

  N = malloc((size_t) P*P*sizeof(double));
  rhs = malloc((size_t) P*sizeof(double));
  theta = calloc((size_t) nb*nq, sizeof(double));
  c = calloc(nb, sizeof(double));
  scatter = calloc(nb, sizeof(double));
  iters = calloc(nb, sizeof(int));
  fitted = calloc(nb, sizeof(bool));
  support = calloc(nq, sizeof(double));
  X = malloc((size_t) nsoundings*sizeof(double));
  Z = malloc((size_t) nsoundings*sizeof(double));
  w0 = malloc((size_t) nsoundings*sizeof(double));
  w = malloc((size_t) nsoundings*sizeof(double));
  res = malloc((size_t) nsoundings*sizeof(double));
  absr = malloc((size_t) nsoundings*sizeof(double));
  row_cell = malloc((size_t) nsoundings*sizeof(int));
  Xc = malloc((size_t) nsoundings*sizeof(double));
  Zc = malloc((size_t) nsoundings*sizeof(double));
  vc = malloc((size_t) nsoundings*sizeof(double));
  vhat = malloc((size_t) nsoundings*sizeof(double));
  active = malloc((size_t) nsoundings*sizeof(bool));
  cand_cell = malloc((size_t) nsoundings*sizeof(int));

  printf("\nAttenuation (%s): %d scenes of %d bands, %d x %d nodes %.4f x %.4f deg apart, smoothing %g\n", prefix,
    n_scenes, nbs, nodes.nx, nodes.ny, nodes.dx*g0->cellsize, nodes.dy*g0->cellsize, smooth);

// Stage 1: each grid's c_m and g_m(x) at the soundings, X = c - g(x) Z.

  for (n = 0; n < nsoundings; n++) {
    if (status[n] == 1 || status[n] == 2) continue ;
    v = sign*soundings[n][2];
    if (v < dmin || v > dmax) continue ;
    node_weights(&nodes, cell_i[n], cell_j[n], q, f);
    for (l = 0; l < 4; l++) support[q[l]] += f[l];
  }
  // A band's signal sinks into the noise with depth. Selecting (and weighting) the soundings by
  // the observed signal keeps only the deep cells that happen to be bright, which flattens the
  // slope (K too small, most in the red); so after the first pass they are reselected and
  // reweighted by the fitted signal exp(c - g Z), until the selection settles.
  for (m = 0; m < nb; m++) {
    int n_cand = 0, pass, changed;
    double *th = &theta[m*nq], prev_c;
    s = m/nbs;
    for (n = 0; n < nsoundings; n++) {
      if (status[n] == 1 || status[n] == 2) continue ;
      v = sign*soundings[n][2];
      if (v < dmin || v > dmax || v + tides[s] <= 0.0) continue ;
      Zc[n_cand] = v + tides[s];
      v = Lbar[m][cell_i[n]][cell_j[n]] - lsm[m];
      if (v <= 0.0) continue ;
      Xc[n_cand] = log(v);
      vc[n_cand] = v;
      cand_cell[n_cand] = n;
      active[n_cand] = v >= deep_k*sigma[m];
      n_cand++;
    }
    rows = 0;
    for (pass = 0; pass < ATT_MAX_PASSES; pass++) {
      rows = 0;
      for (l = 0; l < n_cand; l++) {
        if (! active[l]) continue ;
        X[rows] = Xc[l];
        Z[rows] = Zc[l];
        v = pass == 0 ? vc[l] : vhat[l];
        w[rows] = w0[rows] = 1.0/(sigma[m]*sigma[m]/(v*v) + ATT_EPS*ATT_EPS);
        row_cell[rows] = cand_cell[l];
        rows++;
      }
      if (rows < ATT_MIN_ROWS) break;
      wzz = 0.0;
      for (n = 0; n < rows; n++) wzz += w0[n]*Z[n]*Z[n];
      lambda = smooth*wzz/nq;
      prev_c = 0.0;
    for (it = 0; it < LYZ_MAX_ITER; it++) {
      for (l = 0; l < P*P; l++) N[l] = 0.0;
      for (l = 0; l < P; l++) rhs[l] = 0.0;
      for (n = 0; n < rows; n++) {
        int k2, idx[5];
        double a[5];
        node_weights(&nodes, cell_i[row_cell[n]], cell_j[row_cell[n]], q, f);
        idx[0] = 0;
        a[0] = 1.0;
        for (l = 0; l < 4; l++) {
          idx[l + 1] = q[l] + 1;
          a[l + 1] = -Z[n]*f[l];
        }
        for (l = 0; l < 5; l++) {
          rhs[idx[l]] += w[n]*a[l]*X[n];
          for (k2 = 0; k2 < 5; k2++) N[idx[l]*P + idx[k2]] += w[n]*a[l]*a[k2];
        }
      }
      // First differences between neighbouring nodes.
      for (i = 0; i < nodes.ny && lambda > 0.0; i++) {
        for (j = 0; j < nodes.nx; j++) {
          int qa = i*nodes.nx + j + 1, qb;
          for (l = 0; l < 2; l++) {
            if (l == 0 && j == nodes.nx - 1) continue ;
            if (l == 1 && i == nodes.ny - 1) continue ;
            qb = l == 0 ? qa + 1 : qa + nodes.nx;
            N[qa*P + qa] += lambda;
            N[qb*P + qb] += lambda;
            N[qa*P + qb] -= lambda;
            N[qb*P + qa] -= lambda;
          }
        }
      }
      if (cholesky_solve(N, rhs, P) != 0) break;
      change = fabs(rhs[0] - prev_c)/(1.0 + fabs(rhs[0]));
      for (l = 0; l < nq; l++) change = MAX(change, fabs(rhs[l + 1] - th[l])/(1.0e-3 + fabs(rhs[l + 1])));
      c[m] = prev_c = rhs[0];
      for (l = 0; l < nq; l++) th[l] = rhs[l + 1];
      fitted[m] = true;
      for (n = 0; n < rows; n++) {
        node_weights(&nodes, cell_i[row_cell[n]], cell_j[row_cell[n]], q, f);
        g = 0.0;
        for (l = 0; l < 4; l++) g += f[l]*th[q[l]];
        res[n] = X[n] - (c[m] - g*Z[n]);
        absr[n] = fabs(res[n]);
      }
      sc = 1.4826*median_double(absr, rows);
      scatter[m] = sc;
      if (huber == 0.0 || sc <= 0.0 || (it > 0 && change < 1.0e-9)) break;
      for (n = 0; n < rows; n++) w[n] = w0[n]*(fabs(res[n]) <= huber*sc ? 1.0 : huber*sc/fabs(res[n]));
    }
      iters[m] = it + 1;
      if (! fitted[m]) break;
      changed = 0;
      for (l = 0; l < n_cand; l++) {
        bool now;
        node_weights(&nodes, cell_i[cand_cell[l]], cell_j[cand_cell[l]], q, f);
        g = 0.0;
        for (i = 0; i < 4; i++) g += f[i]*th[q[i]];
        vhat[l] = exp(c[m] - g*Zc[l]);
        now = vhat[l] >= deep_k*sigma[m];
        if (now != active[l]) changed++;
        active[l] = now;
      }
      if (changed == 0) break;
    }
    if (! fitted[m]) {
      printf("  %-20s %d soundings: %s\n", grid_names[bands[m]], rows, rows < ATT_MIN_ROWS ? "too few, no attenuation" :
        "the attenuation is not determined (nodes without soundings need SMOOTHING > 0)");
      continue ;
    }
    double kmin = INFINITY, kmax = -INFINITY, *kv = malloc((size_t) nq*sizeof(double));
    int nk = 0;
    for (l = 0; l < nq; l++) {
      kmin = MIN(kmin, th[l]/path[s]);
      kmax = MAX(kmax, th[l]/path[s]);
      if (support[l] >= 1.0) kv[nk++] = th[l]/path[s];
    }
    printf("  %-20s %5d soundings, c %8.4f, K %.4f (median of supported nodes; %.4f to %.4f) m^-1, "
      "robust sigma %.3f, %d passes\n", grid_names[bands[m]], rows, c[m], nk ? median_double(kv, nk) : NAN,
      kmin, kmax, scatter[m], pass + 1);
    free(kv);
  }

  af->nodes = nodes;
  af->nq = nq;
  af->theta = theta;
  af->c = c;
  af->support = support;
  af->fitted = fitted;
  free(N);
  free(rhs);
  free(scatter);
  free(iters);
  free(X);
  free(Z);
  free(w0);
  free(w);
  free(res);
  free(absr);
  free(row_cell);
  free(Xc);
  free(Zc);
  free(vc);
  free(vhat);
  free(active);
  free(cand_cell);
}

// The scenes' differences (CONSISTENT): the attenuation's differences between scenes from every
// cell, with one bottom for all scenes.
// The bottom is the same in every scene, so at a cell
//   X_sb = beta_b + c_sb - g_sb(x) (H + h_s)
// for band b of scene s, with beta_b free at every cell (any colour, not only a flat brightness).
// Removing each band's (weighted) mean over the scenes removes beta_b and gbar_b H, where gbar_b
// is the band's mean attenuation over the scenes, and leaves
//   r_sb = C_s[X_sb - c_sb + g_sb (H + h_s)],   g_sb = gbar_b(x) + gamma_sb(x),   sum_s gamma_sb = 0,
// so only the scenes' differences gamma_sb and the depth remain. gbar_b is not identified by the
// images: a bottom that darkens uniformly with depth is the same data as every scene's attenuation
// shifted by the same amount (the tides alone distinguish them, weakly), so it is kept from stage 1.
// gamma and c are fitted over the cells (every cell, or a regular subset of them, plus the
// soundings' cells, where H is fixed to the sounding), with H free at every other cell: a
// bilinear problem, solved by Gauss-Newton with each cell's H eliminated (its Schur complement)
// and a backtracking line search, H then re-solved exactly at each cell. There is no smoothing
// penalty: every node has cells, and a penalty would shrink gamma along the near-gauge gamma ->
// lambda gamma, H -> H/lambda that only the soundings and tides resist. gamma_sb = sum_k U_sk
// eta_bk(x) and c_sb = cbar_b + sum_k U_sk d_bk, with U the orthonormal (Helmert) basis of vectors
// over the scenes that sum to zero; eta_bk bilinear on stage 1's nodes. Residuals are weighted as
// stage 1, by the observed signal, with Huber weights; only grids above DEEP sigma enter, and a
// band only where at least two scenes see the bottom.

#define CON_MAX_CELLS 100000  // a regular subset of the cells beyond this
#define CON_MAX_ITER 50
#define CON_RIDGE 1.0e-8      // relative to the largest diagonal element, for nodes without cells

static void att_consistent(att_field *af, int nb, float ***Lbar, const geogrid *g0, const double *lsm,
    const double *sigma, double deep_k, int n_scenes, const double *tides, double huber, double prior, double dmin,
    double dmax, double sign, const int *cell_i, const int *cell_j, const int *status) {

  int nr = g0->nrows, nc = g0->ncols, S = n_scenes, nbs = nb/n_scenes, nq = af->nq, K1 = n_scenes - 1,
    nb1 = (n_scenes - 1)*(af->nq + 1), P = nbs*(n_scenes - 1)*(af->nq + 1), nloc = nbs*(n_scenes - 1)*5,
    n_valid = 0, stride, ncell = 0, n_fixed = 0, i, j, k, l, m, s, b, n, x, it, *cq = NULL, *loc = NULL;
  double U[LYZ_MAX_BANDS][LYZ_MAX_BANDS], *gbar = NULL, *cbar = NULL, *th = NULL, *trial = NULL, *N = NULL,
    *gr = NULL, *dth = NULL, *X = NULL, *w0 = NULL, *w = NULL, *Hfix = NULL, *Hcur = NULL, *cf = NULL,
    *Jl = NULL, *u = NULL, *JH = NULL, *ra = NULL, *th1 = NULL, *lam = NULL, J0 = 0.0, J1 = 0.0, alpha = 0.0, step = 0.0;
  char *mask = NULL;

  for (m = 0; m < nb; m++) {
    if (! af->fitted[m]) {
      printf("  CONSISTENT: grid %d has no attenuation from stage 1; not applied\n", m + 1);
      return ;
    }
  }
  if (S < 2) {
    printf("  CONSISTENT needs two or more scenes; not applied\n");
    return ;
  }

  // U: column k is (1, ..., 1, -(k + 1), 0, ...)/sqrt((k + 1)(k + 2)), with k + 1 ones.
  for (k = 0; k < K1; k++) {
    for (s = 0; s < S; s++) U[s][k] = s <= k ? 1.0 : (s == k + 1 ? -(k + 1.0) : 0.0);
    for (s = 0; s < S; s++) U[s][k] /= sqrt((k + 1.0)*(k + 2.0));
  }

// The cells: a regular subset of the valid ones, and the soundings' cells with H fixed.

  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      if (isfinite(Lbar[0][i][j])) n_valid++;
    }
  }
  stride = MAX(1, (int) ceil(sqrt((double) n_valid/CON_MAX_CELLS)));
  mask = calloc((size_t) nr*nc, 1);
  Hfix = malloc((size_t) nr*nc*sizeof(double));
  for (n = 0; n < nr*nc; n++) Hfix[n] = 0.0;
  for (n = 0; n < nsoundings; n++) {
    double v;
    if (status[n] == 1 || status[n] == 2) continue ;
    v = sign*soundings[n][2];
    if (v < dmin || v > dmax) continue ;
    x = cell_i[n]*nc + cell_j[n];
    mask[x] = mask[x] < 2 ? 2 : mask[x] + 1;  // 2 + soundings - 1
    Hfix[x] += v;
  }
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      x = i*nc + j;
      if (! isfinite(Lbar[0][i][j])) mask[x] = 0;
      else if (mask[x] == 0 && i % stride == 0 && j % stride == 0) mask[x] = 1;
      if (mask[x]) ncell++;
    }
  }
  cq = malloc((size_t) ncell*4*sizeof(int));
  cf = malloc((size_t) ncell*4*sizeof(double));
  X = malloc((size_t) ncell*nb*sizeof(double));
  w0 = malloc((size_t) ncell*nb*sizeof(double));
  w = malloc((size_t) ncell*nb*sizeof(double));
  Hcur = malloc((size_t) ncell*sizeof(double));
  double *Hgiven = malloc((size_t) ncell*sizeof(double));
  n = 0;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      x = i*nc + j;
      if (! mask[x]) continue ;
      node_weights(&af->nodes, i, j, &cq[4*n], &cf[4*n]);
      Hgiven[n] = mask[x] >= 2 ? Hfix[x]/(mask[x] - 1) : NAN;
      if (mask[x] >= 2) n_fixed++;
      for (m = 0; m < nb; m++) {
        double v = Lbar[m][i][j] - lsm[m];
        X[n*nb + m] = v > 0.0 ? log(v) : 0.0;
        w0[n*nb + m] = v > 0.0 && v >= deep_k*sigma[m] ? 1.0/(sigma[m]*sigma[m]/(v*v) + ATT_EPS*ATT_EPS) : 0.0;
        w[n*nb + m] = w0[n*nb + m];
      }
      n++;
    }
  }
  free(Hfix);
  Hfix = NULL;

// The parameters: per band b, eta_bk at the nodes (b*nb1 + k*nq + q), then d_bk (b*nb1 + K1*nq + k).

  gbar = calloc((size_t) nbs*nq, sizeof(double));
  cbar = calloc(nbs, sizeof(double));
  th = calloc(P, sizeof(double));
  trial = malloc((size_t) P*sizeof(double));
  dth = malloc((size_t) P*sizeof(double));
  gr = malloc((size_t) P*sizeof(double));
  N = malloc((size_t) P*P*sizeof(double));
  Jl = malloc((size_t) nb*nloc*sizeof(double));
  u = malloc((size_t) nloc*sizeof(double));
  loc = malloc((size_t) nloc*sizeof(int));
  JH = malloc((size_t) nb*sizeof(double));
  ra = malloc((size_t) nb*sizeof(double));
  for (b = 0; b < nbs; b++) {
    for (s = 0; s < S; s++) {
      cbar[b] += af->c[s*nbs + b]/S;
      for (l = 0; l < nq; l++) gbar[b*nq + l] += af->theta[(s*nbs + b)*nq + l]/S;
    }
    for (k = 0; k < K1; k++) {
      for (s = 0; s < S; s++) {
        for (l = 0; l < nq; l++) th[b*nb1 + k*nq + l] += U[s][k]*(af->theta[(s*nbs + b)*nq + l] - gbar[b*nq + l]);
        th[b*nb1 + K1*nq + k] += U[s][k]*(af->c[s*nbs + b] - cbar[b]);
      }
    }
  }

  th1 = malloc((size_t) P*sizeof(double));
  lam = calloc(P, sizeof(double));
  for (l = 0; l < P; l++) th1[l] = th[l];
  printf("  CONSISTENT: %d cells (every %d%s, and %d soundings' cells), %d parameters, prior %g\n", ncell, stride,
    stride == 1 ? "" : " in each direction", n_fixed, P, prior);

// One cell's centred residuals at depth H: ra (H = 0) and JH (the derivative), each band's weighted
// mean over its scenes removed and scaled by sqrt(w); grids of a band seen in fewer than two
// scenes get zero. g and Z use the parameters t.

#define CON_CELL(t, n) do { \
    for (b = 0; b < nbs; b++) { \
      double sw_ = 0.0, sa_ = 0.0, sg_ = 0.0, gv_[LYZ_MAX_BANDS], av_[LYZ_MAX_BANDS]; \
      int ns_ = 0; \
      for (s = 0; s < S; s++) { \
        m = s*nbs + b; \
        double g_ = 0.0, ga_ = 0.0, cc_ = cbar[b]; \
        for (l = 0; l < 4; l++) { \
          ga_ += cf[4*(n) + l]*gbar[b*nq + cq[4*(n) + l]]; \
          for (k = 0; k < K1; k++) g_ += U[s][k]*cf[4*(n) + l]*(t)[b*nb1 + k*nq + cq[4*(n) + l]]; \
        } \
        for (k = 0; k < K1; k++) cc_ += U[s][k]*(t)[b*nb1 + K1*nq + k]; \
        gv_[s] = ga_ + g_; \
        av_[s] = X[(n)*nb + m] - cc_ + gv_[s]*tides[s]; \
        if (w[(n)*nb + m] > 0.0) { \
          sw_ += w[(n)*nb + m]; \
          sa_ += w[(n)*nb + m]*av_[s]; \
          sg_ += w[(n)*nb + m]*gv_[s]; \
          ns_++; \
        } \
      } \
      for (s = 0; s < S; s++) { \
        m = s*nbs + b; \
        if (ns_ < 2 || w[(n)*nb + m] <= 0.0) { \
          ra[m] = JH[m] = 0.0; \
          continue ; \
        } \
        ra[m] = sqrt(w[(n)*nb + m])*(av_[s] - sa_/sw_); \
        JH[m] = sqrt(w[(n)*nb + m])*(gv_[s] - sg_/sw_); \
      } \
    } \
  } while (0)

  for (it = 0; it < CON_MAX_ITER; it++) {
    double hh, rH, scale, *absr = NULL;
    int nres = 0;

    // The current depths, the objective, then Huber weights from the residuals.
    J0 = 0.0;
    for (n = 0; n < ncell; n++) {
      CON_CELL(th, n);
      hh = rH = 0.0;
      for (m = 0; m < nb; m++) {
        hh += JH[m]*JH[m];
        rH += JH[m]*ra[m];
      }
      Hcur[n] = isfinite(Hgiven[n]) ? Hgiven[n] : (hh > 1.0e-12 ? -rH/hh : NAN);
      if (! isfinite(Hcur[n])) continue ;
      for (m = 0; m < nb; m++) J0 += (ra[m] + JH[m]*Hcur[n])*(ra[m] + JH[m]*Hcur[n]);
    }
    if (huber > 0.0) {
      absr = malloc((size_t) ncell*nb*sizeof(double));
      for (n = 0; n < ncell; n++) {
        if (! isfinite(Hcur[n])) continue ;
        CON_CELL(th, n);
        for (m = 0; m < nb; m++) {
          if (w[n*nb + m] > 0.0) absr[nres++] = fabs(ra[m] + JH[m]*Hcur[n])/sqrt(w[n*nb + m]/w0[n*nb + m]);
        }
      }
      scale = nres ? 1.4826*median_double(absr, nres) : 0.0;
      // A residual's Huber weight is from its size at the base weight; the next pass's centring uses it.
      for (n = 0; n < ncell && scale > 0.0; n++) {
        if (! isfinite(Hcur[n])) continue ;
        CON_CELL(th, n);
        for (m = 0; m < nb; m++) {
          double e;
          if (w[n*nb + m] <= 0.0) continue ;
          e = fabs(ra[m] + JH[m]*Hcur[n])/sqrt(w[n*nb + m]/w0[n*nb + m]);
          w[n*nb + m] = w0[n*nb + m]*(e <= huber*scale ? 1.0 : huber*scale/e);
        }
      }
      free(absr);
      // The objective again, at the new weights.
      J0 = 0.0;
      for (n = 0; n < ncell; n++) {
        CON_CELL(th, n);
        hh = rH = 0.0;
        for (m = 0; m < nb; m++) {
          hh += JH[m]*JH[m];
          rH += JH[m]*ra[m];
        }
        Hcur[n] = isfinite(Hgiven[n]) ? Hgiven[n] : (hh > 1.0e-12 ? -rH/hh : NAN);
        if (! isfinite(Hcur[n])) continue ;
        for (m = 0; m < nb; m++) J0 += (ra[m] + JH[m]*Hcur[n])*(ra[m] + JH[m]*Hcur[n]);
      }
    }

    // PRIOR: sum lambda_l (theta_l - theta_1,l)^2, towards stage 1's values, with lambda the prior
    // weight times the mean diagonal, in the first iteration's system, of the parameter's kind (the
    // attenuation's nodes, or the scenes' constants): an average parameter's data.
    for (l = 0; l < P; l++) J0 += lam[l]*(th[l] - th1[l])*(th[l] - th1[l]);

    // The Gauss-Newton system, each free cell's H eliminated.
    for (l = 0; l < P*P; l++) N[l] = 0.0;
    for (l = 0; l < P; l++) gr[l] = 0.0;
    for (n = 0; n < ncell; n++) {
      int a1, a2;
      double res[LYZ_MAX_BANDS];
      if (! isfinite(Hcur[n])) continue ;
      CON_CELL(th, n);
      // The local Jacobian: grid m (band b, scene s) depends on band b's eta_bk at the cell's four
      // nodes and its d_bk; columns ordered (b, k, then the four nodes and d), the band's mean over
      // its scenes removed as for the residual.
      for (l = 0; l < nb*nloc; l++) Jl[l] = 0.0;
      for (b = 0; b < nbs; b++) {
        double sw_ = 0.0, mj[5*LYZ_MAX_BANDS];
        int ns_ = 0, c0 = b*K1*5;
        for (l = 0; l < 5*K1; l++) mj[l] = 0.0;
        for (s = 0; s < S; s++) {
          m = s*nbs + b;
          if (w[n*nb + m] > 0.0) {
            sw_ += w[n*nb + m];
            ns_++;
          }
        }
        if (ns_ < 2) continue ;
        for (s = 0; s < S; s++) {
          m = s*nbs + b;
          if (w[n*nb + m] <= 0.0) continue ;
          for (k = 0; k < K1; k++) {
            for (l = 0; l < 4; l++) mj[k*5 + l] += w[n*nb + m]*U[s][k]*cf[4*n + l]*(Hcur[n] + tides[s]);
            mj[k*5 + 4] -= w[n*nb + m]*U[s][k];
          }
        }
        for (s = 0; s < S; s++) {
          m = s*nbs + b;
          if (w[n*nb + m] <= 0.0) continue ;
          for (k = 0; k < K1; k++) {
            for (l = 0; l < 4; l++) {
              Jl[m*nloc + c0 + k*5 + l] = sqrt(w[n*nb + m])*(U[s][k]*cf[4*n + l]*(Hcur[n] + tides[s]) - mj[k*5 + l]/sw_);
            }
            Jl[m*nloc + c0 + k*5 + 4] = sqrt(w[n*nb + m])*(-U[s][k] - mj[k*5 + 4]/sw_);
          }
        }
        for (k = 0; k < K1; k++) {
          for (l = 0; l < 4; l++) loc[c0 + k*5 + l] = b*nb1 + k*nq + cq[4*n + l];
          loc[c0 + k*5 + 4] = b*nb1 + K1*nq + k;
        }
      }
      for (b = 0; b < nbs; b++) {  // bands not seen here still need valid columns (zero)
        int ns_ = 0;
        for (s = 0; s < S; s++) ns_ += w[n*nb + s*nbs + b] > 0.0;
        if (ns_ >= 2) continue ;
        for (k = 0; k < K1; k++) {
          for (l = 0; l < 4; l++) loc[b*K1*5 + k*5 + l] = b*nb1 + k*nq + cq[4*n + l];
          loc[b*K1*5 + k*5 + 4] = b*nb1 + K1*nq + k;
        }
      }
      hh = rH = 0.0;
      for (m = 0; m < nb; m++) {
        res[m] = ra[m] + JH[m]*Hcur[n];
        hh += JH[m]*JH[m];
        rH += JH[m]*res[m];
      }
      for (a1 = 0; a1 < nloc; a1++) {
        b = a1/(5*K1);
        u[a1] = 0.0;
        for (s = 0; s < S; s++) u[a1] += Jl[(s*nbs + b)*nloc + a1]*JH[s*nbs + b];
      }
      // Grid m's row is zero outside its band's columns, so J'J is block diagonal by band; the
      // Schur term couples the bands.
      bool free_h = ! isfinite(Hgiven[n]) && hh > 1.0e-12;
      for (a1 = 0; a1 < nloc; a1++) {
        double gsum = 0.0;
        b = a1/(5*K1);
        for (s = 0; s < S; s++) gsum += Jl[(s*nbs + b)*nloc + a1]*res[s*nbs + b];
        if (free_h) gsum -= u[a1]*rH/hh;
        gr[loc[a1]] += gsum;
        for (a2 = b*5*K1; a2 < (b + 1)*5*K1; a2++) {
          double nsum = 0.0;
          for (s = 0; s < S; s++) nsum += Jl[(s*nbs + b)*nloc + a1]*Jl[(s*nbs + b)*nloc + a2];
          N[loc[a1]*P + loc[a2]] += nsum;
        }
        if (! free_h) continue ;
        for (a2 = 0; a2 < nloc; a2++) N[loc[a1]*P + loc[a2]] -= u[a1]*u[a2]/hh;
      }
    }
    if (it == 0 && prior > 0.0) {
      double tr_eta = 0.0, tr_d = 0.0;
      for (l = 0; l < P; l++) {
        if (l % nb1 < K1*nq) tr_eta += N[l*P + l];
        else tr_d += N[l*P + l];
      }
      for (l = 0; l < P; l++) lam[l] = prior*(l % nb1 < K1*nq ? tr_eta/(nbs*K1*nq) : tr_d/(nbs*K1));
    }
    for (l = 0; l < P; l++) {
      N[l*P + l] += lam[l];
      gr[l] += lam[l]*(th[l] - th1[l]);
    }
    double dmax_n = 0.0;
    for (l = 0; l < P; l++) dmax_n = MAX(dmax_n, N[l*P + l]);
    for (l = 0; l < P; l++) {
      N[l*P + l] += CON_RIDGE*dmax_n;
      dth[l] = -gr[l];
    }
    if (dmax_n <= 0.0 || cholesky_solve(N, dth, P) != 0) {
      printf("  CONSISTENT: the system is singular at iteration %d; stopped\n", it + 1);
      break;
    }

    // Backtrack until the reduced objective falls.
    for (alpha = 1.0; alpha > 1.0e-6; alpha *= 0.5) {
      for (l = 0; l < P; l++) trial[l] = th[l] + alpha*dth[l];
      J1 = 0.0;
      for (n = 0; n < ncell; n++) {
        double H;
        CON_CELL(trial, n);
        hh = rH = 0.0;
        for (m = 0; m < nb; m++) {
          hh += JH[m]*JH[m];
          rH += JH[m]*ra[m];
        }
        H = isfinite(Hgiven[n]) ? Hgiven[n] : (hh > 1.0e-12 ? -rH/hh : NAN);
        if (! isfinite(H)) continue ;
        for (m = 0; m < nb; m++) J1 += (ra[m] + JH[m]*H)*(ra[m] + JH[m]*H);
      }
      for (l = 0; l < P; l++) J1 += lam[l]*(trial[l] - th1[l])*(trial[l] - th1[l]);
      if (J1 < J0) break;
    }
    if (J1 >= J0) {
      printf("  CONSISTENT: converged at iteration %d, objective %.6g\n", it + 1, J0);
      break;
    }
    step = 0.0;
    for (l = 0; l < P; l++) {
      step = MAX(step, fabs(alpha*dth[l]));
      th[l] = trial[l];
    }
    printf("  CONSISTENT: iteration %2d, objective %.6g -> %.6g (step %.3g, largest change %.3g)\n", it + 1, J0, J1,
      alpha, step);
    if (step < 1.0e-9 || (J0 - J1) < 1.0e-12*J0) break;
  }

// Back to each grid's field: g_sb = gbar_b + sum_k U_sk eta_bk, c_sb = cbar_b + sum_k U_sk d_bk.

  for (b = 0; b < nbs; b++) {
    for (s = 0; s < S; s++) {
      m = s*nbs + b;
      af->c[m] = cbar[b];
      for (k = 0; k < K1; k++) af->c[m] += U[s][k]*th[b*nb1 + K1*nq + k];
      for (l = 0; l < nq; l++) {
        af->theta[m*nq + l] = gbar[b*nq + l];
        for (k = 0; k < K1; k++) af->theta[m*nq + l] += U[s][k]*th[b*nb1 + k*nq + l];
      }
    }
  }
  af->consistent = true;
#undef CON_CELL

  free(mask);
  free(cq);
  free(cf);
  free(X);
  free(w0);
  free(w);
  free(Hcur);
  free(Hgiven);
  free(gbar);
  free(cbar);
  free(th);
  free(th1);
  free(lam);
  free(trial);
  free(dth);
  free(gr);
  free(N);
  free(Jl);
  free(u);
  free(loc);
  free(JH);
  free(ra);
}

// ATTENUATION's outputs from stage 1's field: K of every grid, then a depth at every cell from it
// (stage 2), and attenuation.csv; see the header.
static void att_outputs(const att_field *af, const char *prefix, int nb, const int *bands, float ***Lbar,
    const geogrid *g0, const double *lsm, const double *sigma, double deep_k, int n_scenes, const double *tides,
    const double *path, double dmin, double dmax, double sign, const int *cell_i, const int *cell_j,
    const int *status, bool indices) {

  int nr = g0->nrows, nc = g0->ncols, nbs = nb/n_scenes, nq = af->nq, m, s, n, i, j, l, q[4], n_depth = 0,
    n_few = 0;
  double f[4], v, g;
  const double *theta = af->theta, *c = af->c, *support = af->support;
  const bool *fitted = af->fitted;
  const att_nodes nodes = af->nodes;
  float **out = NULL, **depth = NULL, **dsig = NULL, **bottom = NULL, **one[LYZ_MAX_BANDS] = {NULL},
    nodata = g0->nodata_value;
  char name[MAX_ARG_SIZE];
  FILE *fp;

// The attenuation grids, K = g/M.

  for (m = 0; m < nb; m++) {
    if (! fitted[m]) continue ;
    s = m/nbs;
    allocate_float_array_2d(&out, nr, nc);
    for (i = 0; i < nr; i++) {
      for (j = 0; j < nc; j++) {
        out[i][j] = nodata;
        if (! isfinite(Lbar[m][i][j])) continue ;
        node_weights(&nodes, i, j, q, f);
        g = 0.0;
        for (l = 0; l < 4; l++) g += f[l]*theta[m*nq + q[l]];
        out[i][j] = g/path[s];
      }
    }
    if (! att_name(name, prefix, grid_names[bands[m]]) || ! store_grid(name, out, g0)) goto done;
    out = NULL;
  }

// Stage 2: at every cell, X - c + g h = B - g H, a weighted straight line in g.

  allocate_float_array_2d(&depth, nr, nc);
  allocate_float_array_2d(&dsig, nr, nc);
  allocate_float_array_2d(&bottom, nr, nc);
  // CONSISTENT with INDICES: the indices of the one bottom, beta_b rotated with gbar_b (as below).
  int n_one = af->consistent && indices && nbs >= 2 ? nbs - 1 : 0;
  for (l = 0; l < n_one; l++) {
    allocate_float_array_2d(&one[l], nr, nc);
    for (i = 0; i < nr; i++) {
      for (j = 0; j < nc; j++) one[l][i][j] = nodata;
    }
  }
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      depth[i][j] = dsig[i][j] = bottom[i][j] = nodata;
      if (! isfinite(Lbar[0][i][j])) continue ;
      // As stage 1: the grids are chosen and weighted by the observed signal first, then by the
      // fitted one, exp(B + c - g (H + h)).
      double gm[LYZ_MAX_BANDS], ym[LYZ_MAX_BANDS], vm[LYZ_MAX_BANDS], sm[LYZ_MAX_BANDS], wl[LYZ_MAX_BANDS], H = 0.0,
        B = 0.0,
        sw, sgg, sgy, rss, vfit;
      bool use[LYZ_MAX_BANDS], ok2 = false;
      int k = 0, used, pass, changed;
      node_weights(&nodes, i, j, q, f);
      for (m = 0; m < nb; m++) {
        if (! fitted[m]) continue ;
        v = Lbar[m][i][j] - lsm[m];
        if (v <= 0.0) continue ;
        g = 0.0;
        for (l = 0; l < 4; l++) g += f[l]*theta[m*nq + q[l]];
        gm[k] = g;
        ym[k] = log(v) - c[m] + g*tides[m/nbs];
        vm[k] = v;
        sm[k] = sigma[m];
        use[k] = v >= deep_k*sigma[m];
        k++;
      }
      // One bottom for every grid (B, a flat brightness), or with CONSISTENT one per band (beta_b,
      // any colour), shared by the scenes: then only the scenes' differences of g within a band
      // carry the depth. Grid l belongs to group grp[l].
      int grp[LYZ_MAX_BANDS], ngrp = af->consistent ? nbs : 1, ig;
      double gsw[LYZ_MAX_BANDS], gsg[LYZ_MAX_BANDS], gsy[LYZ_MAX_BANDS], beta[LYZ_MAX_BANDS];
      int gn[LYZ_MAX_BANDS];
      for (l = 0, m = 0; m < nb; m++) {
        if (! fitted[m] || Lbar[m][i][j] - lsm[m] <= 0.0) continue ;
        grp[l++] = af->consistent ? m % nbs : 0;
      }
      for (pass = 0; pass < ATT_MAX_PASSES; pass++) {
        sw = sgg = sgy = rss = 0.0;
        used = 0;
        for (ig = 0; ig < ngrp; ig++) {
          gsw[ig] = gsg[ig] = gsy[ig] = 0.0;
          gn[ig] = 0;
        }
        for (l = 0; l < k; l++) {
          if (! use[l]) continue ;
          // The fitted signal is v exp((B - g H) - y).
          vfit = pass == 0 ? vm[l] : vm[l]*exp((beta[grp[l]] - gm[l]*H) - ym[l]);
          wl[l] = 1.0/(sm[l]*sm[l]/(vfit*vfit) + ATT_EPS*ATT_EPS);
          gsw[grp[l]] += wl[l];
          gsg[grp[l]] += wl[l]*gm[l];
          gsy[grp[l]] += wl[l]*ym[l];
          gn[grp[l]]++;
          used++;
        }
        int nseen = 0;
        for (ig = 0; ig < ngrp; ig++) {
          sw += gsw[ig];
          nseen += gn[ig] > 0;
        }
        if (used < 3 || used - nseen < 2) break;
        for (l = 0; l < k; l++) {
          if (! use[l]) continue ;
          sgg += wl[l]*(gm[l] - gsg[grp[l]]/gsw[grp[l]])*(gm[l] - gsg[grp[l]]/gsw[grp[l]]);
          sgy += wl[l]*(gm[l] - gsg[grp[l]]/gsw[grp[l]])*(ym[l] - gsy[grp[l]]/gsw[grp[l]]);
        }
        if (sgg <= 1.0e-12*sw) break;
        double H1 = -sgy/sgg, B1 = 0.0;
        for (ig = 0; ig < ngrp; ig++) {
          if (gn[ig] == 0) continue ;
          beta[ig] = gsy[ig]/gsw[ig] + H1*gsg[ig]/gsw[ig];
          B1 += beta[ig]/nseen;
        }
        for (l = 0; l < k; l++) {
          if (use[l]) rss += wl[l]*(ym[l] - (beta[grp[l]] - gm[l]*H1))*(ym[l] - (beta[grp[l]] - gm[l]*H1));
        }
        H = H1;
        B = B1;
        ok2 = true;
        depth[i][j] = -H;
        dsig[i][j] = used - nseen > 1 ? sqrt(rss/(used - nseen - 1)/sgg) : nodata;
        bottom[i][j] = B;
        changed = 0;
        for (l = 0; l < k; l++) {
          bool now = vm[l]*exp((beta[grp[l]] - gm[l]*H) - ym[l]) >= deep_k*sm[l];
          if (now != use[l]) changed++;
          use[l] = now;
        }
        if (changed == 0) break;
      }
      if (! ok2) {
        n_few++;
        continue ;
      }
      n_depth++;
      for (ig = 0; ig < ngrp && gn[ig] > 0; ig++) ;
      if (n_one > 0 && ig == ngrp) {
        double bb[LYZ_MAX_BANDS], S1, S2, yv;
        int ix, b2;
        for (b2 = 0; b2 < nbs; b2++) {
          bb[b2] = 0.0;
          for (s = 0; s < n_scenes; s++) {
            for (l = 0; l < 4; l++) bb[b2] += f[l]*theta[(s*nbs + b2)*nq + q[l]]/n_scenes;
          }
        }
        for (ix = 0; ix < n_one; ix++) {
          S1 = 0.0;
          for (l = 0; l <= ix; l++) S1 += bb[l]*bb[l];
          S2 = S1 + bb[ix + 1]*bb[ix + 1];
          if (S1 <= 0.0) continue ;
          yv = -sqrt(S1/S2)*beta[ix + 1];
          for (l = 0; l <= ix; l++) yv += bb[ix + 1]*bb[l]/sqrt(S1*S2)*beta[l];
          one[ix][i][j] = yv;
        }
      }
    }
  }

  // In-sample check at the soundings.
  double e, rms = 0.0, bias = 0.0;
  int n_check = 0;
  for (n = 0; n < nsoundings; n++) {
    if (status[n] == 1 || status[n] == 2) continue ;
    v = sign*soundings[n][2];
    if (v < dmin || v > dmax || depth[cell_i[n]][cell_j[n]] == nodata) continue ;
    e = -depth[cell_i[n]][cell_j[n]] - v;
    rms += e*e;
    bias += e;
    n_check++;
  }
  printf("  depth from the attenuation: %d cells, %d with fewer than 3 usable grids; at %d soundings rms %.3f m, "
    "bias %+.3f m\n", n_depth, n_few, n_check, n_check ? sqrt(rms/n_check) : NAN, n_check ? bias/n_check : NAN);

  if (! att_name(name, prefix, "depth") || ! store_grid(name, depth, g0)) goto done;
  depth = NULL;
  if (! att_name(name, prefix, "depth_sigma") || ! store_grid(name, dsig, g0)) goto done;
  dsig = NULL;
  if (! att_name(name, prefix, "bottom") || ! store_grid(name, bottom, g0)) goto done;
  bottom = NULL;
  for (l = 0; l < n_one; l++) {
    char suffix[32];
    snprintf(suffix, sizeof(suffix), "index%d", l + 1);
    if (! att_name(name, prefix, suffix) || ! store_grid(name, one[l], g0)) goto done;
    one[l] = NULL;
  }

// INDICES: Lyzenga 1978's depth-invariant bottom indices for each scene, with the local attenuation.
// Within a scene, X_j = a_j - b_j Z with b_j = g_j(x); the rotation of eqs. 8 and B2-B7 takes the
// bands to Y_N along b (depth) and N - 1 axes orthogonal to it, which depend on the bottom only:
//   Y_i = sum_j A_ij X_j, A_ij = b_i+1 b_j/sqrt(S_i S_i+1) for j <= i, A_i,i+1 = -sqrt(S_i/S_i+1),
// with S_i = sum_k<=i b_k^2. Here b follows the water across the scene, from stage 1.

  for (s = 0; indices && s < n_scenes && nbs >= 2; s++) {
    int ix;
    for (m = s*nbs; m < (s + 1)*nbs; m++) {
      if (! fitted[m]) break;
    }
    if (m < (s + 1)*nbs) {
      printf("  INDICES: scene %d lacks the attenuation of a band; no indices\n", s + 1);
      continue ;
    }
    for (ix = 0; ix < nbs - 1; ix++) {
      char suffix[32];
      allocate_float_array_2d(&out, nr, nc);
      for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
          double b[LYZ_MAX_BANDS], x[LYZ_MAX_BANDS], S1 = 0.0, S2, y = 0.0;
          out[i][j] = nodata;
          if (! isfinite(Lbar[0][i][j])) continue ;
          for (l = 0; l < nbs; l++) {
            m = s*nbs + l;
            v = Lbar[m][i][j] - lsm[m];
            if (v <= 0.0) break;
            x[l] = log(v);
            node_weights(&nodes, i, j, q, f);
            b[l] = 0.0;
            for (n = 0; n < 4; n++) b[l] += f[n]*theta[m*nq + q[n]];
          }
          if (l < nbs) continue ;
          for (l = 0; l <= ix; l++) S1 += b[l]*b[l];
          S2 = S1 + b[ix + 1]*b[ix + 1];
          if (S1 <= 0.0) continue ;
          for (l = 0; l <= ix; l++) y += b[ix + 1]*b[l]/sqrt(S1*S2)*x[l];
          y -= sqrt(S1/S2)*x[ix + 1];
          out[i][j] = y;
        }
      }
      snprintf(suffix, sizeof(suffix), "index%d_s%d", ix + 1, s + 1);
      if (! att_name(name, prefix, suffix) || ! store_grid(name, out, g0)) goto done;
      out = NULL;
    }
  }

// K at every node.

  fp = fopen("attenuation.csv", "w");
  if (fp != NULL) {
    fprintf(fp, "LON,LAT,SOUNDINGS");
    for (m = 0; m < nb; m++) fprintf(fp, ",%s", grid_names[bands[m]]);
    fprintf(fp, "\n");
    for (i = 0; i < nodes.ny; i++) {
      for (j = 0; j < nodes.nx; j++) {
        l = i*nodes.nx + j;
        fprintf(fp, "%.6f,%.6f,%.2f", g0->wlon + j*nodes.dx*g0->cellsize, g0->slat + i*nodes.dy*g0->cellsize, support[l]);
        for (m = 0; m < nb; m++) {
          if (fitted[m]) fprintf(fp, ",%.6g", theta[m*nq + l]/path[m/nbs]);
          else fprintf(fp, ",");
        }
        fprintf(fp, "\n");
      }
    }
    fclose(fp);
  }
done:
  if (out != NULL) free_float_array_2d(out, nr);
  if (depth != NULL) free_float_array_2d(depth, nr);
  if (dsig != NULL) free_float_array_2d(dsig, nr);
  if (bottom != NULL) free_float_array_2d(bottom, nr);
  for (l = 0; l < LYZ_MAX_BANDS; l++) {
    if (one[l] != NULL) free_float_array_2d(one[l], nr);
  }
}

static void usage(void) {
  printf("\nERROR: MODEL lyzenga_log_linear GRIDS /band 1/ ... /band n/ BATHYMETRY /output/ [optional, any order]\n"
    "  LAND /grid/ SHALLOW /grid/ LSM /Lsm 1/ ... /Lsm n/ PERCENTILE /p/ DEEP /k/ RADIUS /r/\n"
    "  DEPTHS /min/ /max/ TIDE /h/ HUBER /k/ BALANCE /bin (m)/ POWER ERROR /output/\n"
    "  ATTENUATION /prefix/ SCENES /n/ TIDES /h_1/ ... SUNZENITH /theta_1/ ... VIEWZENITH /theta_1/ ...\n"
    "  NODES /deg/ SMOOTHING /w/ INDICES CONSISTENT  CONSTRAIN GAIN [BRIGHTNESS]  NORMALISE (with ATTENUATION)\n");
}

static bool parse_number(int *k, double *v) {
  char *end;
  if (*k + 1 >= MAX_ARGS || parsed[*k + 1][0] == '\0') return false;
  *v = strtod(parsed[*k + 1], &end);
  if (*end != '\0') return false;
  (*k)++;
  return true;
}

void run_model_lyzenga(void) {

  int k, m, i, j, ii, jj, n, nb = 0, nr, nc, p, ps_out = -1, ps_error = -1, ps_land = -1, ps_shallow = -1,
    bands[LYZ_MAX_BANDS], radius = depth_radius, n_cells = 0, n_fit = 0, it, *cell_i = NULL, *cell_j = NULL,
    dropped_outside = 0, dropped_masked = 0, dropped_deep = 0, dropped_range = 0, count, *status = NULL, bin;
  int ks;
  double v, lsm[LYZ_MAX_BANDS], sigma[LYZ_MAX_BANDS], pct = 1.0, deep_k = 2.0, dmin = depthmin, dmax = depthmax,
    tide = 0.0, huber = 1.345, balance = 0.0, coef[LYZ_MAX_BANDS + 1] = {0.0},
    R[(LYZ_MAX_BANDS + 1)*(LYZ_MAX_BANDS + 1)], C[(LYZ_MAX_BANDS + 1)*(LYZ_MAX_BANDS + 1)], *A = NULL, *y = NULL,
    *w0 = NULL, *w = NULL, *res = NULL, *tmp = NULL, s, s2, rss, sw, sign = 1.0, mae, bias, ss_tot, mean_y,
    xrow[LYZ_MAX_BANDS + 1], pw_a = 0.0, pw_d = 1.0;
  bool lsm_given = false, error_wanted = false, power = false, power_all = false, ok, scenes_given = false,
    tides_given = false, con_gain = false, con_bright = false, normalise = false, error_noise = false,
    error_gain = false, indices = false, consistent = false;
  double con_prior = 0.0;
  att_field af = {0};
  int n_con = 0, piv[2] = {0, 0}, p_r = 0, freec[LYZ_MAX_BANDS + 1];
  double E[2][LYZ_MAX_BANDS + 1], dcon[2] = {0.0, 0.0}, gbar[LYZ_MAX_BANDS] = {0.0};
  int n_scenes = 1, k_att = -1;
  double pct_lsm[LYZ_MAX_BANDS], att_lsm[LYZ_MAX_BANDS], tides[LYZ_MAX_BANDS], sunz[LYZ_MAX_BANDS] = {0.0},
    viewz[LYZ_MAX_BANDS] = {0.0}, path[LYZ_MAX_BANDS], spacing = 0.05, smooth = 1.0;
  float ***Lbar = NULL, **bathy = NULL, **err = NULL, *vals, nodata, *depths;
  geogrid *g0;
  FILE *fp;

// Read the band grids and the output name.

  if (strcmp(parsed[2], "GRIDS") != 0) {
    usage();
    return ;
  }
  for (k = 3; k < MAX_ARGS && parsed[k][0] != '\0' && strcmp(parsed[k], "BATHYMETRY") != 0; k++) {
    if (nb == LYZ_MAX_BANDS) {
      printf("\nERROR: at most %d bands.\n", LYZ_MAX_BANDS);
      return ;
    }
    if (! parse_old_grid(&bands[nb], k)) return ;
    nb++;
  }
  if (nb == 0 || k >= MAX_ARGS - 1 || strcmp(parsed[k], "BATHYMETRY") != 0 || parsed[k + 1][0] == '\0') {
    usage();
    return ;
  }
  k++;
  int k_out = k++;
  int k_error = -1;

// Options, in any order.

  for (; k < MAX_ARGS && parsed[k][0] != '\0'; k++) {
    if (strcmp(parsed[k], "LAND") == 0 || strcmp(parsed[k], "SHALLOW") == 0) {
      if (! parse_old_grid(strcmp(parsed[k], "LAND") == 0 ? &ps_land : &ps_shallow, k + 1)) return ;
      k++;
    } else if (strcmp(parsed[k], "LSM") == 0) {
      for (m = 0; m < nb; m++) {
        if (! parse_number(&k, &lsm[m])) {
          printf("\nERROR: LSM needs one value per band (%d).\n", nb);
          return ;
        }
      }
      lsm_given = true;
    } else if (strcmp(parsed[k], "PERCENTILE") == 0) {
      if (! parse_number(&k, &pct) || pct <= 0.0 || pct >= 100.0) {
        printf("\nERROR: PERCENTILE /p/ with 0 < p < 100.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "DEEP") == 0) {
      if (! parse_number(&k, &deep_k) || deep_k < 0.0) {
        printf("\nERROR: DEEP /k/ with k >= 0.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "RADIUS") == 0) {
      if (! parse_number(&k, &v) || v < 1.0) {
        printf("\nERROR: RADIUS /r/ with r >= 1.\n");
        return ;
      }
      radius = (int) v;
    } else if (strcmp(parsed[k], "DEPTHS") == 0) {
      if (! parse_number(&k, &dmin) || ! parse_number(&k, &dmax) || dmax <= dmin) {
        printf("\nERROR: DEPTHS /min/ /max/ with min < max.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "TIDE") == 0) {
      if (! parse_number(&k, &tide)) {
        printf("\nERROR: TIDE /h/\n");
        return ;
      }
    } else if (strcmp(parsed[k], "HUBER") == 0) {
      if (! parse_number(&k, &huber) || huber < 0.0) {
        printf("\nERROR: HUBER /k/ with k >= 0.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "BALANCE") == 0) {
      if (! parse_number(&k, &balance) || balance <= 0.0) {
        printf("\nERROR: BALANCE /bin width (m)/ > 0.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "POWER") == 0) {
      power = true;
      if (k + 1 < MAX_ARGS && strcmp(parsed[k + 1], "ALL") == 0) {
        power_all = true;
        k++;
      }
    } else if (strcmp(parsed[k], "ATTENUATION") == 0) {
      if (k + 1 >= MAX_ARGS || parsed[k + 1][0] == '\0') {
        usage();
        return ;
      }
      k_att = ++k;
    } else if (strcmp(parsed[k], "SCENES") == 0) {
      if (! parse_number(&k, &v) || v < 1.0 || v != floor(v) || v > LYZ_MAX_BANDS) {
        printf("\nERROR: SCENES /n/ with n a whole number from 1 to %d.\n", LYZ_MAX_BANDS);
        return ;
      }
      n_scenes = (int) v;
      scenes_given = true;
    } else if (strcmp(parsed[k], "TIDES") == 0 || strcmp(parsed[k], "SUNZENITH") == 0 ||
        strcmp(parsed[k], "VIEWZENITH") == 0) {
      double *dst = parsed[k][0] == 'T' ? tides : parsed[k][0] == 'S' ? sunz : viewz;
      const char *key = parsed[k];
      if (! scenes_given) {
        printf("\nERROR: SCENES /n/ must come before %s.\n", key);
        return ;
      }
      for (m = 0; m < n_scenes; m++) {
        if (! parse_number(&k, &dst[m]) || (dst != tides && (dst[m] < 0.0 || dst[m] >= 90.0))) {
          printf("\nERROR: %s needs one value per scene (%d)%s.\n", key, n_scenes,
            dst == tides ? "" : ", from 0 to 90 degrees");
          return ;
        }
      }
      if (dst == tides) tides_given = true;
    } else if (strcmp(parsed[k], "NODES") == 0) {
      if (! parse_number(&k, &spacing) || spacing <= 0.0) {
        printf("\nERROR: NODES /spacing (degrees)/ > 0.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "SMOOTHING") == 0) {
      if (! parse_number(&k, &smooth) || smooth < 0.0) {
        printf("\nERROR: SMOOTHING /w/ with w >= 0.\n");
        return ;
      }
    } else if (strcmp(parsed[k], "CONSTRAIN") == 0) {
      if (k + 1 >= MAX_ARGS || strcmp(parsed[k + 1], "GAIN") != 0) {
        printf("\nERROR: CONSTRAIN GAIN [BRIGHTNESS]\n");
        return ;
      }
      con_gain = true;
      k++;
      if (k + 1 < MAX_ARGS && strcmp(parsed[k + 1], "BRIGHTNESS") == 0) {
        con_bright = true;
        k++;
      }
    } else if (strcmp(parsed[k], "NORMALISE") == 0) {
      normalise = true;
    } else if (strcmp(parsed[k], "ERROR") == 0) {
      if (k + 1 >= MAX_ARGS || parsed[k + 1][0] == '\0') {
        usage();
        return ;
      }
      error_wanted = true;
      k_error = ++k;
      while (k + 1 < MAX_ARGS && (strcmp(parsed[k + 1], "NOISE") == 0 || strcmp(parsed[k + 1], "GAIN") == 0)) {
        if (parsed[k + 1][0] == 'N') error_noise = true;
        else error_gain = true;
        k++;
      }
    } else if (strcmp(parsed[k], "INDICES") == 0) {
      indices = true;
    } else if (strcmp(parsed[k], "CONSISTENT") == 0) {
      consistent = true;
      if (k + 1 < MAX_ARGS && strcmp(parsed[k + 1], "PRIOR") == 0) {
        k++;
        if (! parse_number(&k, &con_prior) || con_prior < 0.0) {
          printf("\nERROR: CONSISTENT [PRIOR /w/] with w >= 0.\n");
          return ;
        }
      }
    } else {
      printf("\nERROR: MODEL lyzenga_log_linear: unexpected '%s'.\n", parsed[k]);
      usage();
      return ;
    }
  }

// Check the inputs.

  g0 = &gridded_data[bands[0]];
  for (m = 1; m < nb; m++) {
    if (! commensurate_grids(*g0, gridded_data[bands[m]])) {
      printf("\nERROR: input grids are not commensurate.\n");
      return ;
    }
  }
  if ((ps_land >= 0 && ! commensurate_grids(*g0, gridded_data[ps_land])) ||
      (ps_shallow >= 0 && ! commensurate_grids(*g0, gridded_data[ps_shallow]))) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }
  if (nsoundings == 0) {
    printf("\nERROR: no soundings (READ SOUNDINGS).\n");
    return ;
  }
  for (n = 0; n < MAX_GRIDS; n++) {
    if (grid_name_taken(n, parsed[k_out]) || (error_wanted && grid_name_taken(n, parsed[k_error]))) {
      printf("\nERROR: grid '%s' already in use.\n", grid_name_taken(n, parsed[k_out]) ? parsed[k_out] : parsed[k_error]);
      return ;
    }
  }
  if (error_wanted && strcmp(parsed[k_out], parsed[k_error]) == 0) {
    printf("\nERROR: the bathymetry and error grids need different names.\n");
    return ;
  }
  if ((con_gain || normalise || error_gain || indices || consistent) && k_att < 0) {
    printf("\nERROR: CONSTRAIN, NORMALISE, ERROR ... GAIN, INDICES and CONSISTENT need ATTENUATION.\n");
    return ;
  }
  if (k_att >= 0) {
    char name[MAX_ARG_SIZE];
    if (nb % n_scenes != 0) {
      printf("\nERROR: ATTENUATION: %d grids do not divide into %d scenes.\n", nb, n_scenes);
      return ;
    }
    for (m = 0; m < nb; m++) {
      if (! att_name(name, parsed[k_att], grid_names[bands[m]])) return ;
    }
    if (! att_name(name, parsed[k_att], "depth") || ! att_name(name, parsed[k_att], "depth_sigma") ||
        ! att_name(name, parsed[k_att], "bottom")) return ;
    for (m = 0; consistent && indices && m < nb/n_scenes - 1; m++) {
      char suffix[32];
      snprintf(suffix, sizeof(suffix), "index%d", m + 1);
      if (! att_name(name, parsed[k_att], suffix)) return ;
    }
    for (ks = 0; indices && ks < n_scenes; ks++) {
      for (m = 0; m < nb/n_scenes - 1; m++) {
        char suffix[32];
        snprintf(suffix, sizeof(suffix), "index%d_s%d", m + 1, ks + 1);
        if (! att_name(name, parsed[k_att], suffix)) return ;
      }
    }
  }

  nr = g0->nrows;
  nc = g0->ncols;
  nodata = g0->nodata_value;
  p = nb + 1;

// Window means of each band over valid cells (no data outside the grid or at nodata).

  Lbar = malloc(nb*sizeof(float **));
  for (m = 0; m < nb; m++) {
    allocate_float_array_2d(&Lbar[m], nr, nc);
  }
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      ok = (ps_land < 0 || gridded_data[ps_land].array[i][j] == gridded_data[ps_land].nodata_value) &&
        (ps_shallow < 0 || gridded_data[ps_shallow].array[i][j] != gridded_data[ps_shallow].nodata_value);
      for (m = 0; m < nb && ok; m++) {
        if (gridded_data[bands[m]].array[i][j] == gridded_data[bands[m]].nodata_value) ok = false;
      }
      for (m = 0; m < nb; m++) {
        Lbar[m][i][j] = NAN;
        if (! ok) continue ;
        s = 0.0;
        count = 0;
        for (ii = MAX(i - radius + 1, 0); ii <= MIN(i + radius - 1, nr - 1); ii++) {
          for (jj = MAX(j - radius + 1, 0); jj <= MIN(j + radius - 1, nc - 1); jj++) {
            v = gridded_data[bands[m]].array[ii][jj];
            if (v != gridded_data[bands[m]].nodata_value) {
              s += v;
              count++;
            }
          }
        }
        Lbar[m][i][j] = s/count;
      }
      if (ok) n_cells++;
    }
  }
  if (n_cells == 0) {
    printf("\nERROR: no cells to model.\n");
    goto done;
  }

// The deep-water signal and the spread of the darkest cells, per band.

  vals = malloc((size_t) n_cells*sizeof(float));
  for (m = 0; m < nb; m++) {
    n = 0;
    for (i = 0; i < nr; i++) {
      for (j = 0; j < nc; j++) {
        if (isfinite(Lbar[m][i][j])) vals[n++] = Lbar[m][i][j];
      }
    }
    v = percentile(vals, n, pct);
    pct_lsm[m] = v;
    if (! lsm_given) lsm[m] = v;
    // The robust spread (1.4826 MAD) of the cells at or below the percentile.
    count = MAX(1, (int) ceil(pct/100.0*n));
    tmp = malloc((size_t) count*sizeof(double));
    for (i = 0; i < count; i++) tmp[i] = vals[i];
    s = median_double(tmp, count);
    for (i = 0; i < count; i++) tmp[i] = fabs(vals[i] - s);
    sigma[m] = 1.4826*median_double(tmp, count);
    free(tmp);
    tmp = NULL;
  }
  free(vals);

// Soundings: sign, cell, and the fit's rows.

  depths = malloc((size_t) nsoundings*sizeof(float));
  for (n = 0; n < nsoundings; n++) depths[n] = soundings[n][2];
  if (percentile(depths, nsoundings, 50.0) < 0.0) sign = -1.0;
  free(depths);
  if (sign < 0.0) printf("\nThe soundings' median is negative: taken as negative down.\n");

  A = malloc((size_t) nsoundings*p*sizeof(double));
  y = malloc((size_t) nsoundings*sizeof(double));
  w0 = malloc((size_t) nsoundings*sizeof(double));
  w = malloc((size_t) nsoundings*sizeof(double));
  res = malloc((size_t) nsoundings*sizeof(double));
  tmp = malloc((size_t) nsoundings*sizeof(double));
  cell_i = malloc((size_t) nsoundings*sizeof(int));
  cell_j = malloc((size_t) nsoundings*sizeof(int));
  status = malloc((size_t) nsoundings*sizeof(int));

  for (n = 0; n < nsoundings; n++) {
    cell_i[n] = (int) lround((soundings[n][1] - g0->slat)/g0->cellsize);
    cell_j[n] = (int) lround((soundings[n][0] - g0->wlon)/g0->cellsize);
    status[n] = 0;  // 0 used, 1 outside, 2 masked, 3 deep, 4 out of range
    if (cell_i[n] < 0 || cell_i[n] >= nr || cell_j[n] < 0 || cell_j[n] >= nc) {
      status[n] = 1;
      dropped_outside++;
      continue ;
    }
    if (! isfinite(Lbar[0][cell_i[n]][cell_j[n]])) {
      status[n] = 2;
      dropped_masked++;
      continue ;
    }
    A[n_fit*p] = 1.0;
    for (m = 0; m < nb; m++) {
      v = Lbar[m][cell_i[n]][cell_j[n]] - lsm[m];
      if (v <= 0.0 || v < deep_k*sigma[m]) break;
      A[n_fit*p + m + 1] = log(v);
    }
    if (m < nb) {
      status[n] = 3;
      dropped_deep++;
      continue ;
    }
    v = sign*soundings[n][2] + tide;
    if (v < dmin || v > dmax) {
      status[n] = 4;
      dropped_range++;
      continue ;
    }
    y[n_fit] = v;
    tmp[n_fit] = n;  // the sounding of each row
    n_fit++;
  }

  printf("\nLyzenga log-linear model: %d bands, %d modelled cells, radius %d.\n", nb, n_cells, radius);
  for (m = 0; m < nb; m++) {
    printf("  %-20s Lsm %.6g%s, sigma %.3g\n", grid_names[bands[m]], lsm[m], lsm_given ? "" : " (percentile)", sigma[m]);
  }
  printf("  soundings: %d read, %d used; dropped %d outside the grid, %d masked or nodata, %d optically deep, "
    "%d outside %.2f-%.2f m\n", nsoundings, n_fit, dropped_outside, dropped_masked, dropped_deep, dropped_range, dmin, dmax);
  if (n_fit < 2*p) {
    printf("\nERROR: too few soundings to fit %d coefficients.\n", p);
    goto done;
  }

// The attenuation at the soundings (ATTENUATION stage 1), before the fit: CONSTRAIN and NORMALISE
// use it. Lyzenga's model is the linearised X_m = B + c_m - g_m Z (see the header), so its depth
// is a0 + sum a_m (B + c_m) - (sum a_m g_m) Z: unbiased (Lyzenga et al. 2006, eq. 12) only where
// the gain -sum a_m g_m is 1. Least squares meets that only on average over the soundings, and
// less than 1 for noisy bands (regression dilution). CONSTRAIN GAIN imposes it at the soundings'
// mean attenuation, -sum a_m gbar_m = 1, and BRIGHTNESS adds sum a_m = 0, so that a spectrally
// flat change of B leaves the depth unchanged. NORMALISE divides each predictor by the local
// g_m(x), so the gain is -sum a_m at every cell, whatever the attenuation there; CONSTRAIN GAIN
// then imposes sum a_m = -1 (and BRIGHTNESS sum a_m/gbar_m = 0). The constraints are linear
// equalities, solved exactly by eliminating one coefficient each.

  if (k_att >= 0) {
    for (ks = 0; ks < n_scenes; ks++) {
      if (! tides_given) tides[ks] = tide;
      path[ks] = 1.0/cos(asin(sin(sunz[ks]*M_PI/180.0)/N_WATER)) + 1.0/cos(asin(sin(viewz[ks]*M_PI/180.0)/N_WATER));
    }
    for (m = 0; m < nb; m++) att_lsm[m] = lsm_given && lsm[m] > 0.0 ? lsm[m] : pct_lsm[m];
    att_stage1(&af, parsed[k_att], nb, bands, Lbar, g0, att_lsm, sigma, deep_k, n_scenes, tides, path, spacing,
      smooth, huber, dmin, dmax, sign, cell_i, cell_j, status);
    if (consistent) {
      att_consistent(&af, nb, Lbar, g0, att_lsm, sigma, deep_k, n_scenes, tides, huber, con_prior, dmin, dmax, sign,
        cell_i, cell_j, status);
    }
  }
  if (con_gain || normalise) {
    for (m = 0; m < nb; m++) {
      if (! af.fitted[m]) {
        printf("\nERROR: CONSTRAIN/NORMALISE: no attenuation for grid '%s'.\n", grid_names[bands[m]]);
        goto done;
      }
      for (n = 0; n < n_fit; n++) {
        v = att_g(&af, m, cell_i[(int) tmp[n]], cell_j[(int) tmp[n]]);
        if (v <= 0.0) {
          printf("\nERROR: NORMALISE: the attenuation of grid '%s' is not positive at a sounding.\n",
            grid_names[bands[m]]);
          goto done;
        }
        gbar[m] += v/n_fit;
        if (normalise) A[n*p + m + 1] /= v;
      }
    }
  }
  if (con_gain) {
    double Cc[2][LYZ_MAX_BANDS + 1] = {{0.0}};
    // Row 0, the gain; row 1, the brightness (no constraint touches the intercept).
    for (m = 0; m < nb; m++) {
      Cc[0][m + 1] = normalise ? 1.0 : -gbar[m];
      Cc[1][m + 1] = normalise ? 1.0/gbar[m] : 1.0;
    }
    dcon[0] = normalise ? -1.0 : 1.0;
    dcon[1] = 0.0;
    n_con = con_bright ? 2 : 1;
    // Reduced row echelon form, pivoting on the largest entry.
    for (ks = 0; ks < n_con; ks++) {
      int best = 1, r2;
      for (m = 1; m < p; m++) if (fabs(Cc[ks][m]) > fabs(Cc[ks][best])) best = m;
      if (fabs(Cc[ks][best]) < 1.0e-12) {
        printf("\nERROR: CONSTRAIN: the constraints are degenerate.\n");
        goto done;
      }
      piv[ks] = best;
      s = Cc[ks][best];
      for (m = 0; m < p; m++) Cc[ks][m] /= s;
      dcon[ks] /= s;
      for (r2 = 0; r2 < n_con; r2++) {
        if (r2 == ks) continue ;
        s = Cc[r2][best];
        for (m = 0; m < p; m++) Cc[r2][m] -= s*Cc[ks][m];
        dcon[r2] -= s*dcon[ks];
      }
    }
    for (ks = 0; ks < n_con; ks++) {
      for (m = 0; m < p; m++) E[ks][m] = Cc[ks][m];
    }
    // a_piv = d - sum over the free columns of E a: the reduced design and target.
    p_r = 0;
    for (m = 0; m < p; m++) {
      if (m != piv[0] && (n_con < 2 || m != piv[1])) freec[p_r++] = m;
    }
    printf("  constrained: %s%s%s, at the soundings' mean attenuation\n", normalise ? "sum a = -1" : "-sum a gbar = 1",
      con_bright ? (normalise ? ", sum a/gbar = 0" : ", sum a = 0") : "", normalise ? " (predictors / g(x))" : "");
  } else if (normalise) {
    printf("  normalised: predictors divided by the local two-way attenuation g(x)\n");
  }

// Weights: equal, or 1/count in depth bins (BALANCE).

  for (n = 0; n < n_fit; n++) w0[n] = 1.0;
  if (balance > 0.0) {
    int n_bins = (int) floor((dmax - dmin)/balance) + 1, *counts = calloc(n_bins, sizeof(int));
    for (n = 0; n < n_fit; n++) counts[MIN((int) floor((y[n] - dmin)/balance), n_bins - 1)]++;
    for (n = 0; n < n_fit; n++) {
      bin = MIN((int) floor((y[n] - dmin)/balance), n_bins - 1);
      w0[n] = ((double) n_fit)/(counts[bin]*n_bins);
    }
    free(counts);
  }

// Iteratively reweighted least squares (Huber), from the least-squares solution.

  if (n_con == 0) {
    it = huber_irls(A, y, w0, n_fit, p, huber, coef, R, w, res);
  } else {
    double *Ar = malloc((size_t) n_fit*p_r*sizeof(double)), *yr = malloc((size_t) n_fit*sizeof(double)), cr[LYZ_MAX_BANDS + 1];
    for (n = 0; n < n_fit; n++) {
      yr[n] = y[n];
      for (ks = 0; ks < n_con; ks++) yr[n] -= A[n*p + piv[ks]]*dcon[ks];
      for (j = 0; j < p_r; j++) {
        s = A[n*p + freec[j]];
        for (ks = 0; ks < n_con; ks++) s -= A[n*p + piv[ks]]*E[ks][freec[j]];
        Ar[n*p_r + j] = s;
      }
    }
    it = huber_irls(Ar, yr, w0, n_fit, p_r, huber, cr, R, w, res);
    free(Ar);
    free(yr);
    if (it >= 0) {
      for (j = 0; j < p_r; j++) coef[freec[j]] = cr[j];
      for (ks = 0; ks < n_con; ks++) {
        s = dcon[ks];
        for (j = 0; j < p_r; j++) s -= E[ks][freec[j]]*cr[j];
        coef[piv[ks]] = s;
      }
    }
  }
  if (it < 0) {
    printf("\nERROR: the bands are collinear at the soundings; use fewer bands.\n");
    goto done;
  }

// Fit statistics: the weighted residual variance, the coefficients' covariance, and the
// unweighted in-sample errors.

  rss = sw = mae = bias = mean_y = ss_tot = 0.0;
  for (n = 0; n < n_fit; n++) mean_y += y[n]/n_fit;
  for (n = 0; n < n_fit; n++) {
    s = coef[0];
    for (m = 1; m < p; m++) s += coef[m]*A[n*p + m];
    res[n] = y[n] - s;
    rss += w[n]*res[n]*res[n];
    sw += w[n];
    mae += fabs(res[n])/n_fit;
    bias -= res[n]/n_fit;
    ss_tot += (y[n] - mean_y)*(y[n] - mean_y);
  }
  // Scale the weights to sum to n, so s2 is a variance of one sounding of average weight.
  int p_fit = n_con ? p_r : p;
  s2 = rss*n_fit/sw/(n_fit - p_fit);
  covariance_factor(R, p_fit, C);
  for (m = 0; m < p_fit*p_fit; m++) C[m] *= s2*sw/n_fit;
  if (n_con) {
    // a = T a_r + t, with T[free_j][j] = 1 and T[piv_k][j] = -E[k][free_j]: C = T C_r T'.
    double T[LYZ_MAX_BANDS + 1][LYZ_MAX_BANDS + 1] = {{0.0}}, Cr[(LYZ_MAX_BANDS + 1)*(LYZ_MAX_BANDS + 1)];
    for (m = 0; m < p_r*p_r; m++) Cr[m] = C[m];
    for (j = 0; j < p_r; j++) {
      T[freec[j]][j] = 1.0;
      for (ks = 0; ks < n_con; ks++) T[piv[ks]][j] = -E[ks][freec[j]];
    }
    for (ii = 0; ii < p; ii++) {
      for (jj = 0; jj < p; jj++) {
        s = 0.0;
        for (i = 0; i < p_r; i++) {
          for (j = 0; j < p_r; j++) s += T[ii][i]*Cr[i*p_r + j]*T[jj][j];
        }
        C[ii*p + jj] = s;
      }
    }
  }
  double rms = 0.0, rmin = INFINITY, rmax = 0.0;
  for (n = 0; n < n_fit; n++) rms += res[n]*res[n]/n_fit;
  for (m = 0; m < p_fit; m++) {
    rmin = MIN(rmin, fabs(R[m*p_fit + m]));
    rmax = MAX(rmax, fabs(R[m*p_fit + m]));
  }

  printf("  fit: H = %.4f", coef[0]);
  for (m = 1; m < p; m++) printf(" %s %.4f ln(%s - Lsm)", coef[m] < 0.0 ? "-" : "+", fabs(coef[m]), grid_names[bands[m - 1]]);
  printf("\n  standard errors:");
  for (m = 0; m < p; m++) printf(" %.4f", sqrt(C[m*p + m]));
  printf("\n  %s, %d iterations; rms %.3f m, mae %.3f m, bias %+.3f m, R2 %.3f, robust sigma %.3f m; "
    "|R| diagonal ratio %.3g\n", huber > 0.0 ? "Huber IRLS" : "least squares", it, sqrt(rms), mae, bias,
    1.0 - rms*n_fit/MAX(ss_tot, 1.0e-30), sqrt(s2), rmax/rmin);
  if (af.fitted != NULL) {
    // The gain of eq. 12 at the soundings' mean attenuation, and the response to a flat brightness.
    double gain = 0.0, bright = 0.0;
    ok = true;
    for (m = 0; m < nb; m++) {
      if (! af.fitted[m]) ok = false;
      gain -= coef[m + 1]*(normalise ? 1.0 : gbar[m]);
      bright += coef[m + 1]/(normalise ? gbar[m] : 1.0);
    }
    if (ok && (con_gain || normalise)) printf("  gain %.4f (1 is unbiased), brightness response %.4f m per unit ln B\n", gain, bright);
  }

// The power-law rescaling (POWER): ln h = a + d ln z, robustly, over the soundings with h > 0.

  if (power) {
    double *A2 = malloc((size_t) 2*n_fit*sizeof(double)), *y2 = malloc((size_t) n_fit*sizeof(double)),
      *wb = malloc((size_t) n_fit*sizeof(double)), *w2 = malloc((size_t) n_fit*sizeof(double)),
      *r2 = malloc((size_t) n_fit*sizeof(double)), c2[2], R2[4], h, e;
    int n2 = 0, it2;
    for (n = 0; n < n_fit; n++) {
      h = y[n] - res[n];
      if (h <= 0.0 || y[n] <= 0.0) continue ;
      A2[2*n2] = 1.0;
      A2[2*n2 + 1] = log(y[n]);
      y2[n2] = log(h);
      wb[n2] = w0[n];
      n2++;
    }
    it2 = n2 >= 4 ? huber_irls(A2, y2, wb, n2, 2, huber, c2, R2, w2, r2) : -1;
    free(A2);
    free(y2);
    free(wb);
    free(w2);
    free(r2);
    if (it2 < 0 || c2[1] <= 0.0) {
      printf("\nERROR: POWER: the fitted depths do not increase with the soundings' depths.\n");
      goto done;
    }
    pw_a = c2[0];
    pw_d = c2[1];
    double prms = 0.0, pmae = 0.0, pbias = 0.0;
    for (n = 0; n < n_fit; n++) {
      e = power_map(y[n] - res[n], pw_a, pw_d, power_all) - y[n];
      prms += e*e/n_fit;
      pmae += fabs(e)/n_fit;
      pbias += e/n_fit;
    }
    printf("  power law: fit = %.4f depth^%.4f (%d soundings), so depth = (fit/%.4f)^%.4f %s; rms %.3f m, "
      "mae %.3f m, bias %+.3f m\n", exp(pw_a), pw_d, n2, exp(pw_a), 1.0/pw_d, power_all ? "everywhere" : "on the deep side",
      sqrt(prms), pmae, pbias);
    if (! power_all && pw_d != 1.0) printf("  (beyond the crossover at %.2f m)\n", exp(pw_a/(1.0 - pw_d)));
  }

// The bathymetry (and its error), below the datum, negative down.

  // NOISE: the variance of a cell's depth from its bands' noise (Lyzenga 1978, eq. 16).
  double s2_floor = s2;
  #define LYZ_NOISE2(ci, cj, out) do { \
    out = 0.0; \
    for (int mm = 0; mm < nb; mm++) { \
      double vv = Lbar[mm][ci][cj] - lsm[mm], gg = normalise ? att_g(&af, mm, ci, cj) : 1.0; \
      out += coef[mm + 1]*coef[mm + 1]*sigma[mm]*sigma[mm]/(vv*vv*gg*gg); \
    } \
  } while (0)
  if (error_wanted && error_noise) {
    double nbar = 0.0, n2;
    for (n = 0; n < n_fit; n++) {
      LYZ_NOISE2(cell_i[(int) tmp[n]], cell_j[(int) tmp[n]], n2);
      nbar += n2/n_fit;
    }
    s2_floor = MAX(s2 - nbar, 0.0);
    printf("  ERROR NOISE: the soundings' mean noise %.3f m of the fit's %.3f m\n", sqrt(nbar), sqrt(s2));
  }

  allocate_float_array_2d(&bathy, nr, nc);
  if (error_wanted) allocate_float_array_2d(&err, nr, nc);
  int n_deep = 0, n_depth = 0;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      bathy[i][j] = nodata;
      if (error_wanted) err[i][j] = nodata;
      if (! isfinite(Lbar[0][i][j])) continue ;
      xrow[0] = 1.0;
      for (m = 0; m < nb; m++) {
        v = Lbar[m][i][j] - lsm[m];
        if (v <= 0.0 || v < deep_k*sigma[m]) break;
        xrow[m + 1] = log(v);
        if (normalise) {
          double gm = att_g(&af, m, i, j);
          if (gm <= 0.0) break;
          xrow[m + 1] /= gm;
        }
      }
      if (m < nb) {
        n_deep++;
        continue ;
      }
      double h = 0.0, z;
      for (m = 0; m < p; m++) h += coef[m]*xrow[m];
      z = power ? power_map(h, pw_a, pw_d, power_all) : h;
      bathy[i][j] = -(z - tide);
      n_depth++;
      if (error_wanted) {
        double n2 = 0.0, scale;
        if (error_noise) LYZ_NOISE2(i, j, n2);
        s = s2_floor + n2;
        for (ii = 0; ii < p; ii++) {
          for (jj = 0; jj < p; jj++) s += xrow[ii]*C[ii*p + jj]*xrow[jj];
        }
        // Through POWER, an error in h scales by dz/dh = z/(d h).
        scale = power && h > 0.0 && z != h ? z/(pw_d*h) : 1.0;
        if (error_gain) {
          double gain = 0.0, fm;
          for (m = 0; m < nb; m++) {
            fm = Lbar[m][i][j] - att_lsm[m] > 0.0 ? (Lbar[m][i][j] - att_lsm[m])/(Lbar[m][i][j] - lsm[m]) : 0.0;
            fm = MIN(fm, 1.0);
            gain -= coef[m + 1]*fm*(normalise ? 1.0 : att_g(&af, m, i, j));
          }
          scale = MAX(scale, 1.0/MAX(MIN(gain, 1.0), 0.05));
        }
        err[i][j] = sqrt(s)*scale;
      }
    }
  }
  printf("  %d cells with a depth, %d optically deep\n", n_depth, n_deep);

  if (! parse_new_grid(&ps_out, k_out)) goto done;
  gridded_data[ps_out] = *g0;
  gridded_data[ps_out].array = bathy;
  bathy = NULL;
  meminuse += 1.0e-6*((float) nr*nc)*((float) sizeof(float));
  if (error_wanted) {
    if (! parse_new_grid(&ps_error, k_error)) goto done;
    gridded_data[ps_error] = *g0;
    gridded_data[ps_error].array = err;
    err = NULL;
    meminuse += 1.0e-6*((float) nr*nc)*((float) sizeof(float));
  }

// The model as applied, and every sounding.

  fp = fopen("model.bam", "w");
  if (fp != NULL) {
    fprintf(fp, "MODEL=LYZENGA_LOG_LINEAR\nNBANDS=%d\nBANDS=", nb);
    for (m = 0; m < nb; m++) fprintf(fp, "%s%s", grid_names[bands[m]], m < nb - 1 ? "," : "\n");
    fprintf(fp, "LSM=");
    for (m = 0; m < nb; m++) fprintf(fp, "%.8g%s", lsm[m], m < nb - 1 ? "," : "\n");
    fprintf(fp, "SIGMA=");
    for (m = 0; m < nb; m++) fprintf(fp, "%.8g%s", sigma[m], m < nb - 1 ? "," : "\n");
    fprintf(fp, "COEFFICIENTS=");  // H = a0 + sum a_m ln(L_m - Lsm_m)
    for (m = 0; m < p; m++) fprintf(fp, "%.8g%s", coef[m], m < p - 1 ? "," : "\n");
    if (power) {
      fprintf(fp, "POWER=%.8g,%.8g,%s\n", pw_a, pw_d, power_all ? "ALL" : "BEYOND");  // depth = exp((ln H - a)/d)
    } else {
      fprintf(fp, "POWER=none\n");
    }
    fprintf(fp, "CONSTRAIN=%s\nNORMALISE=%s\n", con_gain ? (con_bright ? "GAIN,BRIGHTNESS" : "GAIN") : "none",
      normalise ? "yes" : "no");
    fprintf(fp, "RADIUS=%d\nDEEP=%g\nTIDE=%g\nDEPTHS=%g,%g\nHUBER=%g\nBALANCE=%g\nNSOUNDINGS=%d\nRMS=%.4f\n",
      radius, deep_k, tide, dmin, dmax, huber, balance, n_fit, sqrt(rms));
    fclose(fp);
  }
  fp = fopen("output.csv", "w");
  if (fp != NULL) {
    static const char *why[] = {"used", "outside", "masked", "deep", "range"};
    double *model_at = calloc(nsoundings, sizeof(double)), *weight_at = calloc(nsoundings, sizeof(double));
    for (n = 0; n < n_fit; n++) {
      int q = (int) tmp[n];
      model_at[q] = (power ? power_map(y[n] - res[n], pw_a, pw_d, power_all) : y[n] - res[n]) - tide;
      weight_at[q] = w[n];
    }
    fprintf(fp, "LON,LAT,ROW,COL,DEPTH,MODEL DEPTH,RESIDUAL,WEIGHT,STATUS\n");
    for (n = 0; n < nsoundings; n++) {
      if (status[n] == 0) {
        fprintf(fp, "%.7f,%.7f,%d,%d,%.3f,%.3f,%.3f,%.4f,%s\n", soundings[n][0], soundings[n][1], cell_i[n], cell_j[n],
          sign*soundings[n][2], model_at[n], model_at[n] - sign*soundings[n][2], weight_at[n], why[0]);
      } else {
        fprintf(fp, "%.7f,%.7f,%d,%d,%.3f,,,,%s\n", soundings[n][0], soundings[n][1], cell_i[n], cell_j[n],
          sign*soundings[n][2], why[status[n]]);
      }
    }
    free(model_at);
    free(weight_at);
    fclose(fp);
  }

// The attenuation over the whole grid, and the depth from it (ATTENUATION).

  if (k_att >= 0 && af.fitted != NULL) {
    att_outputs(&af, parsed[k_att], nb, bands, Lbar, g0, att_lsm, sigma, deep_k, n_scenes, tides, path, dmin, dmax,
      sign, cell_i, cell_j, status, indices);
  }

done:
  att_free(&af);
  if (Lbar != NULL) {
    for (m = 0; m < nb; m++) free_float_array_2d(Lbar[m], nr);
    free(Lbar);
  }
  if (bathy != NULL) free_float_array_2d(bathy, nr);
  if (err != NULL) free_float_array_2d(err, nr);
  free(A);
  free(y);
  free(w0);
  free(w);
  free(res);
  free(tmp);
  free(cell_i);
  free(cell_j);
  free(status);
}
