//
//    synth_gen.c
//

/* Synthetic test data for the semi-analytical model. */

// Writes Landsat-8-like band grids (443, 482, 561, 655 nm) for two scenes over a known
// depth field, using photic's own forward model (samodel_Rrs) and IOP/bottom tables,
// so MODEL Lee_Semi_Analytical should recover the depths. Also writes truth_H.asc.
//
// Usage: synth_gen /output directory/

#include "common.h"
#include "samodel.h"

#define NROWS 16
#define NCOLS 16
#define NSCENES 2
#define NBANDS 4

extern int n_ref_wlens;
extern double ref_wlens[], aw_Pope_Fry1997[], bbw_Morel_1974[], a0_Lee[], a1_Lee[], bottom_type_sand[];

static const double wlon = 115.0, slat = -21.0, cellsize = 0.0003;
static const double wavelengths[NBANDS] = {443.0, 482.0, 561.0, 655.0};

// Per scene: P, G, X, tide (m), subsurface sun and view zenith angles (degrees). The water
// must be turbid: samodel_error() heavily penalises K_min below ~0.2-0.28 /m at depths
// under 5 m, so clear water is pushed deep. These IOPs give K_min ~0.30 and ~0.35 /m.
static const double P_s[NSCENES] = {0.30, 0.35};
static const double G_s[NSCENES] = {0.50, 0.60};
static const double X_s[NSCENES] = {0.08, 0.09};
static const double tide_s[NSCENES] = {0.0, 0.4};
static const double theta_sun_s[NSCENES] = {30.0, 40.0};
static const double theta_view_s[NSCENES] = {0.0, 0.0};

// Sand bottom with albedo B.
static const double B_true = 0.6;

static double truth_depth(int i, int j) {
  return 1.5 + 0.35*i + 0.1*j;
}

static void write_grid(const char *dir, const char *name, float **grid) {

  char file[MAX_STRING_LEN];
  FILE *fp;
  int i, j;

  snprintf(file, sizeof(file), "%s/%s.asc", dir, name);
  fp = fopen(file, "w");
  if (fp == NULL) {
    printf("\nERROR: unable to write '%s'\n", file);
    exit(1);
  }

  fprintf(fp, "ncols %d\nnrows %d\nxllcorner %.6f\nyllcorner %.6f\ncellsize %.6f\nNODATA_value -9999\n",
    NCOLS, NROWS, wlon, slat, cellsize);

  // ESRI ascii grids are north row first; row 0 is south in photic.

  for (i = NROWS - 1; i >= 0; i--) {
    for (j = 0; j < NCOLS; j++)
      fprintf(fp, "%.8f ", grid[i][j]);
    fprintf(fp, "\n");
  }

  fclose(fp);
}

int main(int argc, char **argv) {

  model_data md;
  double ***Rrs_modelled, ***rrs_modelled, ***rrs_bottom, ***rrs_dp, ***rho, ***bottom_reflectance,
    **K, **a_0, **a_1, **a_w, **b_bw, **wl, **Rrs440, **Rrs490, **Rrs750,
    h_tide[NSCENES], sec_sun[NSCENES], sec_view[NSCENES], B[1] = {B_true}, q[1] = {1.0},
    Rrs[NBANDS], R440, R490, change;
  float ***bands, **truth;
  char name[64];
  int i, j, s, b, it;

  if (argc != 2) {
    printf("\nUsage: synth_gen /output directory/\n\n");
    return 1;
  }

  // Model data for one region.

  memset(&md, 0, sizeof(md));
  md.n_scenes = NSCENES;
  md.n_bottoms = 1;
  md.bottom_type_indexes[0] = 0;

  allocate_double_array_2d(&a_0, NSCENES, NBANDS);
  allocate_double_array_2d(&a_1, NSCENES, NBANDS);
  allocate_double_array_2d(&a_w, NSCENES, NBANDS);
  allocate_double_array_2d(&b_bw, NSCENES, NBANDS);
  allocate_double_array_2d(&wl, NSCENES, NBANDS);
  allocate_double_array_2d(&K, NSCENES, NBANDS);
  allocate_double_array_2d(&Rrs440, 1, NSCENES);
  allocate_double_array_2d(&Rrs490, 1, NSCENES);
  allocate_double_array_2d(&Rrs750, 1, NSCENES);
  allocate_double_array_3d(&bottom_reflectance, 1, NSCENES, NBANDS);
  allocate_double_array_3d(&rho, 1, NSCENES, NBANDS);
  allocate_double_array_3d(&rrs_dp, 1, NSCENES, NBANDS);
  allocate_double_array_3d(&rrs_bottom, 1, NSCENES, NBANDS);
  allocate_double_array_3d(&rrs_modelled, 1, NSCENES, NBANDS);
  allocate_double_array_3d(&Rrs_modelled, 1, NSCENES, NBANDS);

  // As in samodel(): tables interpolated at the band wavelengths.

  for (s = 0; s < NSCENES; s++) {
    h_tide[s] = tide_s[s];
    sec_sun[s] = 1.0/cos(theta_sun_s[s]*PI/180.0);
    sec_view[s] = 1.0/cos(theta_view_s[s]*PI/180.0);
    Rrs750[0][s] = 0.0;
    for (b = 0; b < NBANDS; b++) {
      wl[s][b] = wavelengths[b];
      a_0[s][b] = interp_1d(ref_wlens, a0_Lee, n_ref_wlens, wavelengths[b]);
      a_1[s][b] = interp_1d(ref_wlens, a1_Lee, n_ref_wlens, wavelengths[b]);
      a_w[s][b] = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, wavelengths[b]);
      b_bw[s][b] = interp_1d(ref_wlens, bbw_Morel_1974, n_ref_wlens, wavelengths[b]);
      bottom_reflectance[0][s][b] = interp_1d(ref_wlens, bottom_type_sand, n_ref_wlens, wavelengths[b]);
    }
  }

  md.a_0 = a_0; md.a_1 = a_1; md.a_w = a_w; md.b_bw = b_bw; md.wavelengths = wl; md.K = K;
  md.Rrs440 = Rrs440; md.Rrs490 = Rrs490; md.Rrs750 = Rrs750;
  md.H_tide = h_tide; md.sec_theta_sun = sec_sun; md.sec_theta_view = sec_view;
  md.bottom_reflectance = bottom_reflectance; md.rho = rho; md.rrs_dp = rrs_dp;
  md.rrs_bottom = rrs_bottom; md.rrs_modelled = rrs_modelled; md.Rrs_modelled = Rrs_modelled;

  allocate_float_array_3d(&bands, NSCENES*NBANDS, NROWS, NCOLS);
  allocate_float_array_2d(&truth, NROWS, NCOLS);

  for (i = 0; i < NROWS; i++) {
    for (j = 0; j < NCOLS; j++) {

      truth[i][j] = truth_depth(i, j);

      for (s = 0; s < NSCENES; s++) {

        // The backscattering slope depends on Rrs(440)/Rrs(490), which the model
        // interpolates from the bands themselves, so iterate to a fixed point.

        Rrs440[0][s] = 0.01;
        Rrs490[0][s] = 0.01;

        for (it = 0; it < 100; it++) {
          for (b = 0; b < NBANDS; b++) {
            samodel_Rrs(P_s[s], G_s[s], X_s[s], B, q, truth[i][j], 0.0, &md, s, b, 0);
            Rrs[b] = Rrs_modelled[0][s][b];
          }
          R440 = interp_1d(wl[s], Rrs, NBANDS, 440.0);
          R490 = interp_1d(wl[s], Rrs, NBANDS, 490.0);
          change = fabs(R440 - Rrs440[0][s]) + fabs(R490 - Rrs490[0][s]);
          Rrs440[0][s] = R440;
          Rrs490[0][s] = R490;
          if (change < 1.0e-12)
            break ;
        }

        for (b = 0; b < NBANDS; b++)
          bands[s*NBANDS + b][i][j] = (float) Rrs[b];

        if (i == 0 && j == 0) {
          double K_min = K[s][0];
          for (b = 1; b < NBANDS; b++)
            K_min = MIN(K_min, K[s][b]);
          printf("scene %d: K_min = %.3f /m, Rrs = %.4f %.4f %.4f %.4f\n", s + 1, K_min, Rrs[0], Rrs[1], Rrs[2], Rrs[3]);
        }
      }
    }
  }

  for (s = 0; s < NSCENES; s++) {
    for (b = 0; b < NBANDS; b++) {
      snprintf(name, sizeof(name), "s%d_b%d", s + 1, (int) wavelengths[b]);
      write_grid(argv[1], name, bands[s*NBANDS + b]);
    }
  }
  write_grid(argv[1], "truth_H", truth);

  return 0;
}
