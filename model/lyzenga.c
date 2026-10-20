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
//   LSM /Lsm 1/ ... /Lsm n/      deep-water signal per band (0 gives Lyzenga et al. 2006's ln L)
//   PERCENTILE /p/               else Lsm is the p-th percentile of the window means (default 1)
//   DEEP /k/                     optically deep below k sigma above Lsm (default 2; 0: off)
//   RADIUS /r/                   window radius (default SET DEPTHRADIUS)
//   DEPTHS /min/ /max/           soundings used for the fit (default SET DEPTHMIN, DEPTHMAX)
//   TIDE /h/                     water level above the datum at the imagery's time (m)
//   HUBER /k/                    Huber threshold in robust sigmas (default 1.345; 0: least squares)
//   BALANCE /bin/                weight soundings by 1/count in depth bins of this width (m)
//   POWER [ALL]                  rescale the depths by the inverted power law above
//   ERROR /grid/                 output: the 1-sigma prediction error (m)
//
// Outputs model.bam (the fitted model as applied) and output.csv (every sounding, its cell, the
// model depth and whether it was used).

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lyzenga.h"

#define LYZ_MAX_BANDS 32  // e.g. 3 bands of each of 6 scenes, fitted jointly
#define LYZ_MAX_ITER 100

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

static void usage(void) {
  printf("\nERROR: MODEL lyzenga_log_linear GRIDS /band 1/ ... /band n/ BATHYMETRY /output/ [optional, any order]\n"
    "  LAND /grid/ SHALLOW /grid/ LSM /Lsm 1/ ... /Lsm n/ PERCENTILE /p/ DEEP /k/ RADIUS /r/\n"
    "  DEPTHS /min/ /max/ TIDE /h/ HUBER /k/ BALANCE /bin (m)/ POWER ERROR /output/\n");
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
  double v, lsm[LYZ_MAX_BANDS], sigma[LYZ_MAX_BANDS], pct = 1.0, deep_k = 2.0, dmin = depthmin, dmax = depthmax,
    tide = 0.0, huber = 1.345, balance = 0.0, coef[LYZ_MAX_BANDS + 1] = {0.0},
    R[(LYZ_MAX_BANDS + 1)*(LYZ_MAX_BANDS + 1)], C[(LYZ_MAX_BANDS + 1)*(LYZ_MAX_BANDS + 1)], *A = NULL, *y = NULL,
    *w0 = NULL, *w = NULL, *res = NULL, *tmp = NULL, s, s2, rss, sw, sign = 1.0, mae, bias, ss_tot, mean_y,
    xrow[LYZ_MAX_BANDS + 1], pw_a = 0.0, pw_d = 1.0;
  bool lsm_given = false, error_wanted = false, power = false, power_all = false, ok;
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
    } else if (strcmp(parsed[k], "ERROR") == 0) {
      if (k + 1 >= MAX_ARGS || parsed[k + 1][0] == '\0') {
        usage();
        return ;
      }
      error_wanted = true;
      k_error = ++k;
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

  if ((it = huber_irls(A, y, w0, n_fit, p, huber, coef, R, w, res)) < 0) {
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
  s2 = rss*n_fit/sw/(n_fit - p);
  covariance_factor(R, p, C);
  for (m = 0; m < p*p; m++) C[m] *= s2*sw/n_fit;
  double rms = 0.0, rmin = INFINITY, rmax = 0.0;
  for (n = 0; n < n_fit; n++) rms += res[n]*res[n]/n_fit;
  for (m = 0; m < p; m++) {
    rmin = MIN(rmin, fabs(R[m*p + m]));
    rmax = MAX(rmax, fabs(R[m*p + m]));
  }

  printf("  fit: H = %.4f", coef[0]);
  for (m = 1; m < p; m++) printf(" %s %.4f ln(%s - Lsm)", coef[m] < 0.0 ? "-" : "+", fabs(coef[m]), grid_names[bands[m - 1]]);
  printf("\n  standard errors:");
  for (m = 0; m < p; m++) printf(" %.4f", sqrt(C[m*p + m]));
  printf("\n  %s, %d iterations; rms %.3f m, mae %.3f m, bias %+.3f m, R2 %.3f, robust sigma %.3f m; "
    "|R| diagonal ratio %.3g\n", huber > 0.0 ? "Huber IRLS" : "least squares", it, sqrt(rms), mae, bias,
    1.0 - rms*n_fit/MAX(ss_tot, 1.0e-30), sqrt(s2), rmax/rmin);

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
        s = s2;
        for (ii = 0; ii < p; ii++) {
          for (jj = 0; jj < p; jj++) s += xrow[ii]*C[ii*p + jj]*xrow[jj];
        }
        // Through POWER, an error in h scales by dz/dh = z/(d h).
        err[i][j] = sqrt(s)*(power && h > 0.0 && z != h ? z/(pw_d*h) : 1.0);
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

done:
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
