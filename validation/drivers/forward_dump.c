//
//    forward_dump.c
//

/* Evaluate photic's semi-analytical forward model for a list of cases. */

// Reads whitespace-separated cases from stdin, one per line:
//
//   id wavelength P G X B bottom_index H tide theta_sun theta_view Rrs440 Rrs490 Rrs750
//
// (angles are subsurface, in degrees; bottom_index selects a type of photic's default bottom
// library, model/data/bottom_hydrolight.csv: 0 sand, 1 seagrass, 2 coral, 3 macrophytes,
// 4 dark sediment, 5 coral sand, 6 green algae, 7 red algae). For each case, sets up model_data the way samodel()
// does (tables interpolated at the wavelength), calls samodel_Rrs() and writes a CSV
// row with the interpolated table values and the model outputs.

#include "common.h"
#include "samodel.h"
#include "spectra.h"

// The spectral data the benchmark's truth was built with (spectra_inputs.py): the legacy water
// tables and the bottom library $PHOTIC_BENCH_BOTTOM or bottom_hydrolight.csv.

static int read_bench_spectra(void) {
  static const char *water[] = {"aw_pope_fry_1997.csv", "bbw_morel_1974.csv", "aphi_lee_1994.csv"};
  char path[1024];
  const char *env = getenv("PHOTIC_BENCH_BOTTOM");
  int i;
  for (i = 0; i < 3; i++) {
    snprintf(path, sizeof(path), "%s/%s", spectra_data_dir(), water[i]);
    if (spectra_read_water(path) != 0) return 1;
  }
  if (env && *env) return spectra_read_bottom(env);
  snprintf(path, sizeof(path), "%s/bottom_hydrolight.csv", spectra_data_dir());
  return spectra_read_bottom(path);
}


int main(void) {

  model_data md;
  double **a_0, **a_1, **a_w, **b_bw, **wl, **K, **Rrs440, **Rrs490, **Rrs750,
    ***bottom_reflectance, ***rho, ***rrs_dp, ***rrs_bottom, ***rrs_modelled, ***Rrs_modelled,
    h_tide[1], sec_sun[1], sec_view[1], B[1], q[1] = {1.0},
    lambda, P, G, X, H, tide, theta_sun, theta_view, r440, r490, r750;
  long id;
  int bottom;
  char line[1024];

  if (read_bench_spectra() != 0 || spectra_load_defaults() != 0) return 1;

  memset(&md, 0, sizeof(md));
  md.n_scenes = 1;
  md.n_bottoms = 1;
  md.bottom_type_indexes[0] = 0;

  allocate_double_array_2d(&a_0, 1, 1);
  allocate_double_array_2d(&a_1, 1, 1);
  allocate_double_array_2d(&a_w, 1, 1);
  allocate_double_array_2d(&b_bw, 1, 1);
  allocate_double_array_2d(&wl, 1, 1);
  allocate_double_array_2d(&K, 1, 1);
  allocate_double_array_2d(&Rrs440, 1, 1);
  allocate_double_array_2d(&Rrs490, 1, 1);
  allocate_double_array_2d(&Rrs750, 1, 1);
  allocate_double_array_3d(&bottom_reflectance, 1, 1, 1);
  allocate_double_array_3d(&rho, 1, 1, 1);
  allocate_double_array_3d(&rrs_dp, 1, 1, 1);
  allocate_double_array_3d(&rrs_bottom, 1, 1, 1);
  allocate_double_array_3d(&rrs_modelled, 1, 1, 1);
  allocate_double_array_3d(&Rrs_modelled, 1, 1, 1);

  md.a_0 = a_0; md.a_1 = a_1; md.a_w = a_w; md.b_bw = b_bw; md.wavelengths = wl; md.K = K;
  md.Rrs440 = Rrs440; md.Rrs490 = Rrs490; md.Rrs750 = Rrs750;
  md.H_tide = h_tide; md.sec_theta_sun = sec_sun; md.sec_theta_view = sec_view;
  md.bottom_reflectance = bottom_reflectance; md.rho = rho; md.rrs_dp = rrs_dp;
  md.rrs_bottom = rrs_bottom; md.rrs_modelled = rrs_modelled; md.Rrs_modelled = Rrs_modelled;

  printf("id,wavelength,a_0,a_1,a_w,b_bw,rho_table,K,rrs_dp,rrs_B,rrs,Rrs,Y,y_model\n");

  while (fgets(line, sizeof(line), stdin) != NULL) {

    if (sscanf(line, "%ld %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf",
        &id, &lambda, &P, &G, &X, &B[0], &bottom, &H, &tide, &theta_sun, &theta_view,
        &r440, &r490, &r750) != 14)
      continue ;

    if (bottom < 0 || bottom >= spectra_n_bottoms()) {
      printf("\nERROR: bottom index %d out of range.\n", bottom);
      return 1;
    }

    // As in samodel().

    wl[0][0] = lambda;
    a_0[0][0] = spectra_water(WATER_A0, lambda);
    a_1[0][0] = spectra_water(WATER_A1, lambda);
    a_w[0][0] = spectra_water(WATER_AW, lambda);
    b_bw[0][0] = spectra_water(WATER_BBW, lambda);
    bottom_reflectance[0][0][0] = spectra_bottom(bottom, lambda);
    h_tide[0] = tide;
    sec_sun[0] = 1.0/cos(theta_sun*PI/180.0);
    sec_view[0] = 1.0/cos(theta_view*PI/180.0);
    Rrs440[0][0] = r440;
    Rrs490[0][0] = r490;
    Rrs750[0][0] = r750;

    samodel_Rrs(P, G, X, B, q, H, 0.0, &md, 0, 0, 0);

    printf("%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",
      id, lambda, a_0[0][0], a_1[0][0], a_w[0][0], b_bw[0][0], bottom_reflectance[0][0][0],
      K[0][0], rrs_dp[0][0][0], rrs_bottom[0][0][0], rrs_modelled[0][0][0], Rrs_modelled[0][0][0],
      samodel_backscatter_slope(P, G, X, &md, 0, 0), samodel_y_model());
  }

  return 0;
}
