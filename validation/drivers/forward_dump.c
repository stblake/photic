//
//    forward_dump.c
//

/* Evaluate photic's semi-analytical forward model for a list of cases. */

// Reads whitespace-separated cases from stdin, one per line:
//
//   id wavelength P G X B bottom_index H tide theta_sun theta_view Rrs440 Rrs490 Rrs750
//
// (angles are subsurface, in degrees; bottom_index selects photic's bottom table:
// 0 sand, 1 seagrass, 2 coral, 3 macrophytes, 4 dark sediment, 5 coral sand,
// 6 green algae, 7 red algae). For each case, sets up model_data the way samodel()
// does (tables interpolated at the wavelength), calls samodel_Rrs() and writes a CSV
// row with the interpolated table values and the model outputs.

#include "common.h"
#include "samodel.h"

extern int n_ref_wlens;
extern double ref_wlens[], aw_Pope_Fry1997[], bbw_Morel_1974[], a0_Lee[], a1_Lee[],
  bottom_type_sand[], bottom_type_seagrass[], bottom_type_coral[], bottom_type_macrophytes[],
  bottom_type_dark_sediment[], bottom_type_coral_sand[], bottom_type_green_algae[],
  bottom_type_red_algae[];

int main(void) {

  double *bottom_tables[8] = {bottom_type_sand, bottom_type_seagrass, bottom_type_coral,
    bottom_type_macrophytes, bottom_type_dark_sediment, bottom_type_coral_sand,
    bottom_type_green_algae, bottom_type_red_algae};
  model_data md;
  double **a_0, **a_1, **a_w, **b_bw, **wl, **K, **Rrs440, **Rrs490, **Rrs750,
    ***bottom_reflectance, ***rho, ***rrs_dp, ***rrs_bottom, ***rrs_modelled, ***Rrs_modelled,
    h_tide[1], sec_sun[1], sec_view[1], B[1], q[1] = {1.0},
    lambda, P, G, X, H, tide, theta_sun, theta_view, r440, r490, r750;
  long id;
  int bottom;
  char line[1024];

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

  printf("id,wavelength,a_0,a_1,a_w,b_bw,rho_table,K,rrs_dp,rrs_B,rrs,Rrs\n");

  while (fgets(line, sizeof(line), stdin) != NULL) {

    if (sscanf(line, "%ld %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf",
        &id, &lambda, &P, &G, &X, &B[0], &bottom, &H, &tide, &theta_sun, &theta_view,
        &r440, &r490, &r750) != 14)
      continue ;

    if (bottom < 0 || bottom > 7) {
      printf("\nERROR: bottom index %d out of range.\n", bottom);
      return 1;
    }

    // As in samodel().

    wl[0][0] = lambda;
    a_0[0][0] = interp_1d(ref_wlens, a0_Lee, n_ref_wlens, lambda);
    a_1[0][0] = interp_1d(ref_wlens, a1_Lee, n_ref_wlens, lambda);
    a_w[0][0] = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, lambda);
    b_bw[0][0] = interp_1d(ref_wlens, bbw_Morel_1974, n_ref_wlens, lambda);
    bottom_reflectance[0][0][0] = interp_1d(ref_wlens, bottom_tables[bottom], n_ref_wlens, lambda);
    h_tide[0] = tide;
    sec_sun[0] = 1.0/cos(theta_sun*PI/180.0);
    sec_view[0] = 1.0/cos(theta_view*PI/180.0);
    Rrs440[0][0] = r440;
    Rrs490[0][0] = r490;
    Rrs750[0][0] = r750;

    samodel_Rrs(P, G, X, B, q, H, 0.0, &md, 0, 0, 0);

    printf("%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
      id, lambda, a_0[0][0], a_1[0][0], a_w[0][0], b_bw[0][0], bottom_reflectance[0][0][0],
      K[0][0], rrs_dp[0][0][0], rrs_bottom[0][0][0], rrs_modelled[0][0][0], Rrs_modelled[0][0][0]);
  }

  return 0;
}
