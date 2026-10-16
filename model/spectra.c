// Spectral data for the semi-analytical model, read from CSV files.
//
// File format (model/data/README.md):
//   - lines starting with # are comments, as are blank lines;
//   - the first other line is the header: wavelength, then one name per column;
//   - each data row starts with a wavelength in nm, strictly increasing down the file;
//   - a bottom library may also hold rows labelled brightness_min and brightness_max: the
//     bounds of each type's reflectance at 550 nm.

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netcdf.h>

#include "common.h"
#include "spectra.h"

#ifndef PHOTIC_DATA_DIR
#define PHOTIC_DATA_DIR "data"
#endif

#define MAX_COLUMNS 16
#define MAX_SRF 16
#define MAX_ROWS 4096
#define LINE_LEN 8192
#define NAME_LEN 64
#define FILE_LEN 1024

// The wavelength at which bottom spectra are normalised, as in the brightness bounds.

#define BOTTOM_NORMALISATION 550.0

static const char *water_names[N_WATER_TABLES] = {"aw", "bbw", "a0", "a1"};
static const char *default_water[] = {"aw_pope_fry_1997.csv", "bbw_morel_1974.csv", "aphi_lee_1994.csv"};

// The bottom when no library is read (READ BOTTOM): bright sand, HydroLight's avg_ooid_sand
// normalised at 550 nm, 400-800 nm every 10 nm (as in model/data/bottom_hydrolight.csv), with
// the brightness bounds of the in situ sediment class.

static const double builtin_sand[41] = {0.634432,0.654965,0.682443,0.712079,0.74123,0.769656,
  0.802245,0.827992,0.854624,0.878036,0.89898,0.924542,0.94899,0.965671,0.982481,1.,1.01237,
  1.02994,1.04658,1.07139,1.08139,1.08884,1.09447,1.09916,1.10557,1.10758,1.101,1.08464,1.08658,
  1.12477,1.17063,1.20125,1.21407,1.22691,1.23186,1.24098,1.24708,1.25257,1.25607,1.26361,1.26789};
static const double builtin_sand_min = 0.15, builtin_sand_max = 0.75;

typedef struct {
  int n;
  double *wl, *v;
  char file[FILE_LEN];
} spectrum;

typedef struct {
  char name[NAME_LEN];
  spectrum shape;
  double min, max;
} bottom_type;

typedef struct {
  char file[FILE_LEN];
  int n_bands, n;
  char names[MAX_COLUMNS][NAME_LEN];
  double *wl, *r[MAX_COLUMNS], centroid[MAX_COLUMNS];
} srf_table;

static spectrum water[N_WATER_TABLES];
static srf_table srfs[MAX_SRF];
static int n_srfs = 0;
static bottom_type bottoms[MAX_BOTTOM_TYPES];
static int n_bottoms = 0;

// A parsed file: column names, data rows and the labelled rows.

typedef struct {
  int n_cols, n_rows, n_labels;
  char names[MAX_COLUMNS][NAME_LEN];
  double *wl, *v[MAX_COLUMNS];
  char labels[4][NAME_LEN];
  double label_v[4][MAX_COLUMNS];
} csv_table;


static char *trim_field(char *s) {
  char *e;
  while (isspace((unsigned char) *s)) s++;
  e = s + strlen(s);
  while (e > s && isspace((unsigned char) e[-1])) *--e = '\0';
  return s;
}

static int split(char *line, char **fields, int max_fields) {
  int n = 0;
  char *p = line, *comma;
  while (n < max_fields) {
    comma = strchr(p, ',');
    if (comma) *comma = '\0';
    fields[n++] = trim_field(p);
    if (!comma) return n;
    p = comma + 1;
  }
  return -1;
}

static int parse_number(const char *s, double *x) {
  char *end;
  if (*s == '\0') return 0;
  *x = strtod(s, &end);
  return *end == '\0';
}

static void free_table(csv_table *t) {
  int k;
  free(t->wl);
  for (k = 0; k < MAX_COLUMNS; k++) free(t->v[k]);
}

static int read_table(const char *file, csv_table *t) {

  FILE *fp;
  char line[LINE_LEN], *fields[MAX_COLUMNS + 1], *s;
  int n, k, line_no = 0, have_header = 0;
  double x;

  memset(t, 0, sizeof(*t));
  fp = fopen(file, "r");
  if (fp == NULL) {
    printf("\nERROR: unable to open '%s'\n", file);
    return 1;
  }
  t->wl = malloc(MAX_ROWS*sizeof(double));

  while (fgets(line, sizeof(line), fp) != NULL) {
    line_no++;
    if (strchr(line, '\n') == NULL && !feof(fp)) {
      printf("\nERROR: %s:%d: line too long\n", file, line_no);
      goto fail;
    }
    s = trim_field(line);
    if (*s == '\0' || *s == '#') continue;
    n = split(s, fields, MAX_COLUMNS + 1);
    if (n < 0) {
      printf("\nERROR: %s:%d: more than %d columns\n", file, line_no, MAX_COLUMNS);
      goto fail;
    }

    // Header.

    if (!have_header) {
      if (strcmp(fields[0], "wavelength") != 0 || n < 2) {
        printf("\nERROR: %s:%d: the header must be 'wavelength,<name>,...'\n", file, line_no);
        goto fail;
      }
      t->n_cols = n - 1;
      for (k = 0; k < t->n_cols; k++) {
        if (fields[k + 1][0] == '\0' || strlen(fields[k + 1]) >= NAME_LEN) {
          printf("\nERROR: %s:%d: bad column name '%s'\n", file, line_no, fields[k + 1]);
          goto fail;
        }
        strcpy(t->names[k], fields[k + 1]);
        t->v[k] = malloc(MAX_ROWS*sizeof(double));
      }
      have_header = 1;
      continue;
    }

    if (n != t->n_cols + 1) {
      printf("\nERROR: %s:%d: %d values, expected %d\n", file, line_no, n, t->n_cols + 1);
      goto fail;
    }

    // Labelled row or data row.

    if (!parse_number(fields[0], &x)) {
      if (t->n_labels == 4 || strlen(fields[0]) >= NAME_LEN) {
        printf("\nERROR: %s:%d: unexpected row '%s'\n", file, line_no, fields[0]);
        goto fail;
      }
      strcpy(t->labels[t->n_labels], fields[0]);
      for (k = 0; k < t->n_cols; k++) {
        if (!parse_number(fields[k + 1], &t->label_v[t->n_labels][k])) {
          printf("\nERROR: %s:%d: bad value '%s'\n", file, line_no, fields[k + 1]);
          goto fail;
        }
      }
      t->n_labels++;
      continue;
    }

    if (t->n_rows == MAX_ROWS) {
      printf("\nERROR: %s: more than %d rows\n", file, MAX_ROWS);
      goto fail;
    }
    if (t->n_rows > 0 && x <= t->wl[t->n_rows - 1]) {
      printf("\nERROR: %s:%d: wavelengths must increase\n", file, line_no);
      goto fail;
    }
    t->wl[t->n_rows] = x;
    for (k = 0; k < t->n_cols; k++) {
      if (!parse_number(fields[k + 1], &t->v[k][t->n_rows])) {
        printf("\nERROR: %s:%d: bad value '%s'\n", file, line_no, fields[k + 1]);
        goto fail;
      }
    }
    t->n_rows++;
  }

  fclose(fp);
  if (!have_header || t->n_rows < 2) {
    printf("\nERROR: %s: needs a header and at least two rows\n", file);
    free_table(t);
    return 1;
  }
  return 0;

fail:
  fclose(fp);
  free_table(t);
  return 1;
}

static void set_spectrum(spectrum *s, const csv_table *t, int col, const char *file) {
  free(s->wl);
  free(s->v);
  s->n = t->n_rows;
  s->wl = malloc(s->n*sizeof(double));
  s->v = malloc(s->n*sizeof(double));
  memcpy(s->wl, t->wl, s->n*sizeof(double));
  memcpy(s->v, t->v[col], s->n*sizeof(double));
  snprintf(s->file, FILE_LEN, "%s", file);
}

static double spectrum_at(const spectrum *s, double lambda) {
  return interp_1d(s->wl, s->v, s->n, lambda);
}


const char *spectra_data_dir(void) {
  const char *dir = getenv("PHOTIC_DATA");
  return (dir && *dir) ? dir : PHOTIC_DATA_DIR;
}


// Read a water file; with only_missing, keep tables already loaded.

static int read_water(const char *file, int only_missing) {

  csv_table t;
  int k, w, found[N_WATER_TABLES] = {0};

  if (read_table(file, &t)) return 1;

  for (k = 0; k < t.n_cols; k++) {
    for (w = 0; w < N_WATER_TABLES; w++) {
      if (strcmp(t.names[k], water_names[w]) == 0) break;
    }
    if (w == N_WATER_TABLES || found[w]) {
      printf("\nERROR: %s: column '%s' is not one of aw, bbw, a0, a1 (or is repeated)\n", file, t.names[k]);
      free_table(&t);
      return 1;
    }
    found[w] = 1;
  }
  if (t.n_labels > 0) {
    printf("\nERROR: %s: labelled row '%s' in a water file\n", file, t.labels[0]);
    free_table(&t);
    return 1;
  }

  for (k = 0; k < t.n_cols; k++) {
    for (w = 0; strcmp(t.names[k], water_names[w]) != 0; w++);
    if (!(only_missing && water[w].n > 0)) set_spectrum(&water[w], &t, k, file);
  }
  free_table(&t);
  return 0;
}

int spectra_read_water(const char *file) {
  return read_water(file, 0);
}


int spectra_read_bottom(const char *file) {

  csv_table t;
  int k, i, j, kmin = -1, kmax = -1;
  double norm;

  if (read_table(file, &t)) return 1;

  if (t.n_cols > MAX_BOTTOM_TYPES) {
    printf("\nERROR: %s: more than %d bottom types\n", file, MAX_BOTTOM_TYPES);
    free_table(&t);
    return 1;
  }
  for (i = 0; i < t.n_labels; i++) {
    if (strcmp(t.labels[i], "brightness_min") == 0 && kmin < 0) kmin = i;
    else if (strcmp(t.labels[i], "brightness_max") == 0 && kmax < 0) kmax = i;
    else {
      printf("\nERROR: %s: unexpected row '%s' (brightness_min, brightness_max)\n", file, t.labels[i]);
      free_table(&t);
      return 1;
    }
  }
  if (t.wl[0] > BOTTOM_NORMALISATION || t.wl[t.n_rows - 1] < BOTTOM_NORMALISATION) {
    printf("\nERROR: %s: the spectra must span %.0f nm\n", file, BOTTOM_NORMALISATION);
    free_table(&t);
    return 1;
  }

  // Check every type before replacing the library.

  for (k = 0; k < t.n_cols; k++) {
    norm = interp_1d(t.wl, t.v[k], t.n_rows, BOTTOM_NORMALISATION);
    for (j = 0; j < t.n_rows; j++) {
      if (!(t.v[k][j] >= 0.0)) break;
    }
    if (!(norm > 0.0) || j < t.n_rows) {
      printf("\nERROR: %s: bottom '%s' must be non-negative and positive at %.0f nm\n", file,
        t.names[k], BOTTOM_NORMALISATION);
      free_table(&t);
      return 1;
    }
    if ((kmin >= 0 && t.label_v[kmin][k] < 0.0) || (kmax >= 0 && t.label_v[kmax][k] <= 0.0) ||
      (kmin >= 0 && kmax >= 0 && t.label_v[kmin][k] > t.label_v[kmax][k])) {
      printf("\nERROR: %s: bad brightness bounds for '%s'\n", file, t.names[k]);
      free_table(&t);
      return 1;
    }
  }

  for (k = 0; k < t.n_cols; k++) {
    snprintf(bottoms[k].name, NAME_LEN, "%s", t.names[k]);
    set_spectrum(&bottoms[k].shape, &t, k, file);
    norm = interp_1d(t.wl, t.v[k], t.n_rows, BOTTOM_NORMALISATION);
    for (j = 0; j < t.n_rows; j++) bottoms[k].shape.v[j] /= norm;
    bottoms[k].min = kmin >= 0 ? t.label_v[kmin][k] : 0.0;
    bottoms[k].max = kmax >= 0 ? t.label_v[kmax][k] : 1.0;
  }
  n_bottoms = t.n_cols;
  free_table(&t);
  return 0;
}


int spectra_load_defaults(void) {

  char path[FILE_LEN];
  int i, w;

  for (i = 0; i < (int) (sizeof(default_water)/sizeof(default_water[0])); i++) {
    for (w = 0; w < N_WATER_TABLES && water[w].n > 0; w++);
    if (w == N_WATER_TABLES) break;
    snprintf(path, sizeof(path), "%s/%s", spectra_data_dir(), default_water[i]);
    if (read_water(path, 1)) return 1;
  }
  for (w = 0; w < N_WATER_TABLES; w++) {
    if (water[w].n == 0) {
      printf("\nERROR: no table of %s (READ WATER)\n", water_names[w]);
      return 1;
    }
  }
  if (n_bottoms == 0) {
    spectrum *s = &bottoms[0].shape;
    s->n = 41;
    s->wl = malloc(s->n*sizeof(double));
    s->v = malloc(s->n*sizeof(double));
    for (i = 0; i < s->n; i++) {
      s->wl[i] = 400.0 + 10.0*i;
      s->v[i] = builtin_sand[i];
    }
    snprintf(s->file, FILE_LEN, "built in");
    snprintf(bottoms[0].name, NAME_LEN, "sand");
    bottoms[0].min = builtin_sand_min;
    bottoms[0].max = builtin_sand_max;
    n_bottoms = 1;
  }
  return 0;
}


void spectra_print(void) {
  int w, k;
  printf("\nWater tables:\n");
  for (w = 0; w < N_WATER_TABLES; w++) {
    if (water[w].n > 0) {
      printf("    %-4s %4d rows, %g-%g nm, %s\n", water_names[w], water[w].n, water[w].wl[0],
        water[w].wl[water[w].n - 1], water[w].file);
    } else {
      printf("    %-4s (not loaded)\n", water_names[w]);
    }
  }
  printf("Bottom library (%d type%s%s%s):\n", n_bottoms, n_bottoms == 1 ? "" : "s", n_bottoms ? ", " : "",
    n_bottoms ? bottoms[0].shape.file : "");
  for (k = 0; k < n_bottoms; k++) {
    printf("    %2d %-20s brightness %.3f-%.3f, 443/550 %.3f, 665/550 %.3f\n", k + 1, bottoms[k].name,
      bottoms[k].min, bottoms[k].max, spectra_bottom(k, 443.0), spectra_bottom(k, 665.0));
  }
}


double spectra_water(int table, double lambda) {
  return spectrum_at(&water[table], lambda);
}

int spectra_n_bottoms(void) {
  return n_bottoms;
}

const char *spectra_bottom_name(int k) {
  return bottoms[k].name;
}

double spectra_bottom(int k, double lambda) {
  return spectrum_at(&bottoms[k].shape, lambda);
}

double spectra_bottom_min(int k) {
  return bottoms[k].min;
}

double spectra_bottom_max(int k) {
  return bottoms[k].max;
}


int spectra_read_srf(const char *file) {

  csv_table t;
  srf_table *f;
  int k, j;
  double sw, swl;

  for (k = 0; k < n_srfs; k++) {
    if (strcmp(srfs[k].file, file) == 0) return k;
  }
  if (n_srfs == MAX_SRF) {
    printf("\nERROR: more than %d spectral response files\n", MAX_SRF);
    return -1;
  }
  if (read_table(file, &t)) return -1;
  if (t.n_labels > 0) {
    printf("\nERROR: %s: labelled row '%s' in a spectral response file\n", file, t.labels[0]);
    free_table(&t);
    return -1;
  }

  f = &srfs[n_srfs];
  for (k = 0; k < t.n_cols; k++) {
    sw = swl = 0.0;
    for (j = 0; j < t.n_rows; j++) {
      if (t.v[k][j] < 0.0) break;
      sw += t.v[k][j];
      swl += t.v[k][j]*t.wl[j];
    }
    if (j < t.n_rows || !(sw > 0.0)) {
      printf("\nERROR: %s: the response of band '%s' must be non-negative and not all zero\n", file, t.names[k]);
      free_table(&t);
      return -1;
    }
    f->centroid[k] = swl/sw;
    strcpy(f->names[k], t.names[k]);
    f->r[k] = malloc(t.n_rows*sizeof(double));
    memcpy(f->r[k], t.v[k], t.n_rows*sizeof(double));
  }
  snprintf(f->file, FILE_LEN, "%s", file);
  f->n_bands = t.n_cols;
  f->n = t.n_rows;
  f->wl = malloc(t.n_rows*sizeof(double));
  memcpy(f->wl, t.wl, t.n_rows*sizeof(double));
  free_table(&t);
  return n_srfs++;
}

int spectra_srf_band(int srf, double lambda) {
  int k, best = -1;
  for (k = 0; k < srfs[srf].n_bands; k++) {
    if (fabs(srfs[srf].centroid[k] - lambda) <= 20.0 &&
      (best < 0 || fabs(srfs[srf].centroid[k] - lambda) < fabs(srfs[srf].centroid[best] - lambda))) {
      best = k;
    }
  }
  return best;
}

const char *spectra_srf_band_name(int srf, int band) {
  return srfs[srf].names[band];
}

double spectra_srf_centroid(int srf, int band) {
  return srfs[srf].centroid[band];
}

// The response-weighted mean of a spectrum over a band (trapezoidal weights on the response's
// wavelength grid).

static double band_average(const spectrum *s, int srf, int band) {
  const srf_table *f = &srfs[srf];
  double w, sw = 0.0, swv = 0.0;
  int j;
  for (j = 0; j < f->n; j++) {
    if (f->r[band][j] == 0.0) continue;
    w = f->r[band][j]*0.5*((j + 1 < f->n ? f->wl[j + 1] : f->wl[j]) - (j > 0 ? f->wl[j - 1] : f->wl[j]));
    sw += w;
    swv += w*spectrum_at(s, f->wl[j]);
  }
  return swv/sw;
}

double spectra_water_srf(int table, int srf, int band) {
  return band_average(&water[table], srf, band);
}

double spectra_bottom_srf(int k, int srf, int band) {
  return band_average(&bottoms[k].shape, srf, band);
}


// Bottom libraries from NetCDF files of individual spectra (bottom_spectra_insitu.nc, written by
// validation/bottom_library.py pack). Each spectrum is smoothed with an 11 nm boxcar onto 410-700
// nm every 5 nm, as bottom_library.py build does, so the result matches bottom_measured.csv.

#define NC_LAM0 410.0
#define NC_DLAM 5.0
#define NC_NLAM 59
#define NC_SMOOTH 5.0

static int compare_double(const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

// Linear-interpolated percentile (numpy's default), sorting v in place.

static double percentile(double *v, int n, double p) {
  double pos = p/100.0*(n - 1);
  int i = (int) pos;
  qsort(v, n, sizeof(double), compare_double);
  return i + 1 < n ? v[i] + (pos - i)*(v[i + 1] - v[i]) : v[n - 1];
}

static int nc_fail(const char *file, int status) {
  printf("\nERROR: %s: %s\n", file, nc_strerror(status));
  return 1;
}

int spectra_read_bottom_nc(const char *file, char **classes, int n_classes, int all_halves) {

  int ncid, status, dim, var, nw, ns, i, j, k, m, n_types = 0, n_lib, n_all, *member, *half_k;
  size_t len, nw_t, ns_t;
  float *wl, *R;
  signed char *half;
  char **cls, *order = NULL, *tok, names[MAX_BOTTOM_TYPES][NAME_LEN];
  double *smooth, *shape, *col, sum, x;
  int count;

  if ((status = nc_open(file, NC_NOWRITE, &ncid)) != NC_NOERR) return nc_fail(file, status);
  if ((status = nc_inq_dimid(ncid, "wavelength", &dim)) != NC_NOERR ||
    (status = nc_inq_dimlen(ncid, dim, &nw_t)) != NC_NOERR ||
    (status = nc_inq_dimid(ncid, "spectrum", &dim)) != NC_NOERR ||
    (status = nc_inq_dimlen(ncid, dim, &ns_t)) != NC_NOERR) {
    nc_close(ncid);
    return nc_fail(file, status);
  }
  nw = (int) nw_t;
  ns = (int) ns_t;
  wl = malloc(nw*sizeof(float));
  R = malloc((size_t) ns*nw*sizeof(float));
  half = malloc(ns*sizeof(signed char));
  cls = malloc(ns*sizeof(char *));
  if ((status = nc_inq_varid(ncid, "wavelength", &var)) != NC_NOERR || (status = nc_get_var_float(ncid, var, wl)) != NC_NOERR ||
    (status = nc_inq_varid(ncid, "reflectance", &var)) != NC_NOERR || (status = nc_get_var_float(ncid, var, R)) != NC_NOERR ||
    (status = nc_inq_varid(ncid, "half", &var)) != NC_NOERR || (status = nc_get_var_schar(ncid, var, half)) != NC_NOERR ||
    (status = nc_inq_varid(ncid, "class", &var)) != NC_NOERR || (status = nc_get_var_string(ncid, var, cls)) != NC_NOERR) {
    free(wl); free(R); free(half); free(cls);
    nc_close(ncid);
    return nc_fail(file, status);
  }

  // The types: the classes asked for, or the file's order (global attribute classes).

  if (n_classes > 0) {
    for (k = 0; k < n_classes && n_types < MAX_BOTTOM_TYPES; k++) {
      snprintf(names[n_types++], NAME_LEN, "%s", classes[k]);
    }
  } else if (nc_inq_attlen(ncid, NC_GLOBAL, "classes", &len) == NC_NOERR) {
    order = calloc(len + 1, 1);
    nc_get_att_text(ncid, NC_GLOBAL, "classes", order);
    for (tok = strtok(order, " "); tok && n_types < MAX_BOTTOM_TYPES; tok = strtok(NULL, " ")) {
      snprintf(names[n_types++], NAME_LEN, "%s", tok);
    }
    free(order);
  }
  nc_close(ncid);
  if (n_types == 0 || n_classes > MAX_BOTTOM_TYPES) {
    printf("\nERROR: %s: no classes given (CLASSES) or in the file, or more than %d\n", file, MAX_BOTTOM_TYPES);
    goto fail;
  }

  // Smooth every spectrum onto the library grid.

  smooth = malloc((size_t) ns*NC_NLAM*sizeof(double));
  for (i = 0; i < ns; i++) {
    for (j = 0; j < NC_NLAM; j++) {
      x = NC_LAM0 + NC_DLAM*j;
      sum = 0.0;
      count = 0;
      for (m = 0; m < nw; m++) {
        if (fabs(wl[m] - x) <= NC_SMOOTH && isfinite(R[(size_t) i*nw + m])) {
          sum += R[(size_t) i*nw + m];
          count++;
        }
      }
      smooth[(size_t) i*NC_NLAM + j] = count >= 3 ? sum/count : NAN;
    }
  }

  // Each type from its class's spectra; check them all before replacing the library.

  member = malloc(ns*sizeof(int));
  half_k = malloc(ns*sizeof(int));
  shape = malloc((size_t) MAX_BOTTOM_TYPES*NC_NLAM*sizeof(double));
  col = malloc(ns*sizeof(double));
  {
    double lo[MAX_BOTTOM_TYPES], hi[MAX_BOTTOM_TYPES];
    int bad = 0;
    for (k = 0; k < n_types && !bad; k++) {
      n_all = 0;
      for (i = 0; i < ns; i++) {
        if (cls[i] && strcmp(cls[i], names[k]) == 0) {
          for (j = 0; j < NC_NLAM && isfinite(smooth[(size_t) i*NC_NLAM + j]); j++);
          if (j < NC_NLAM) continue;  // incomplete over 410-700 nm
          member[n_all] = i;
          half_k[n_all++] = half[i];
        }
      }
      n_lib = 0;
      for (m = 0; m < n_all; m++) {
        if (all_halves || half_k[m] == 0) n_lib++;
      }
      if (n_lib < 1) {
        printf("\nERROR: %s: no spectra of class '%s'\n", file, names[k]);
        bad = 1;
        break;
      }

      // Brightness bounds over the whole class.

      for (m = 0; m < n_all; m++) col[m] = smooth[(size_t) member[m]*NC_NLAM + (int) ((550.0 - NC_LAM0)/NC_DLAM)];
      lo[k] = percentile(col, n_all, 5.0);
      hi[k] = percentile(col, n_all, 95.0);

      // The median shape (photic normalises at 550 nm, so the median reflectance there, which
      // bottom_measured.csv multiplies in, does not matter).

      for (j = 0; j < NC_NLAM; j++) {
        n_lib = 0;
        for (m = 0; m < n_all; m++) {
          if (all_halves || half_k[m] == 0) {
            i = member[m];
            col[n_lib] = smooth[(size_t) i*NC_NLAM + j]/smooth[(size_t) i*NC_NLAM + (int) ((550.0 - NC_LAM0)/NC_DLAM)];
            n_lib++;
          }
        }
        shape[k*NC_NLAM + j] = percentile(col, n_lib, 50.0);
      }
    }
    if (!bad) {
      for (k = 0; k < n_types; k++) {
        spectrum *s = &bottoms[k].shape;
        free(s->wl);
        free(s->v);
        s->n = NC_NLAM;
        s->wl = malloc(NC_NLAM*sizeof(double));
        s->v = malloc(NC_NLAM*sizeof(double));
        for (j = 0; j < NC_NLAM; j++) {
          s->wl[j] = NC_LAM0 + NC_DLAM*j;
          s->v[j] = shape[k*NC_NLAM + j];  // the median shape is 1 at 550 nm
        }
        snprintf(s->file, FILE_LEN, "%s (%s)", file, all_halves ? "all spectra" : "library half");
        snprintf(bottoms[k].name, NAME_LEN, "%s", names[k]);
        bottoms[k].min = lo[k];
        bottoms[k].max = hi[k];
      }
      n_bottoms = n_types;
    }
    free(member); free(half_k); free(shape); free(col); free(smooth);
    nc_free_string(ns, cls);
    free(wl); free(R); free(half); free(cls);
    return bad;
  }

fail:
  nc_free_string(ns, cls);
  free(wl); free(R); free(half); free(cls);
  return 1;
}
