//
//    invert_lab.c
//

/* Run photic's inversion on a single pixel, for diagnosing the SA model (EXPERIMENTS.md). */

// Builds model_data for one pixel and its spatial window the way samodel() does, then
// runs commands. The setup block comes first, one item per line:
//
//   SETUP n_scenes n_bands n_regions origin n_bottoms
//   SCENE k theta_sun theta_view tide w_1 .. w_n      (subsurface angles in degrees)
//   RRS region scene Rrs_1 .. Rrs_n                   (one line per region and scene)
//   END
//
// Commands follow, one per line. Parameter vectors use photic's internal (transformed)
// layout; see samodel_n_params() and samodel_pack() in samodel.c.
//
//   EVAL p_1 .. p_m
//       prints the total error, its components, K_min and the modelled Rrs.
//   NM reqmin konvge kcount  p_1 .. p_m  s_1 .. s_m
//       Nelder-Mead from p with steps s; a zero step holds that parameter fixed.
//       Prints the result vector, error, evaluations and ifault.
//   PHOTIC h_empirical
//       samodel_optimise() exactly as samodel() calls it for one pixel, with no lookup
//       table and no carry-over from a previous pixel. h_empirical < 0 means no prior
//       (h_empirical is then 0, as with a fresh allocation). Prints depth, albedo, K_min,
//       errors, P, G, X per scene and the bottom type percentages.
//   START h [albedo]
//       photic's starting vector and steps for depth h (samodel_initial).
//   PACK H.. albedo.. q.. P G X..
//       the parameter vector for physical values (samodel_pack).
//
// Each command's output ends with a line DONE.

#include "common.h"
#include "samodel.h"

extern int n_ref_wlens;
extern double ref_wlens[], aw_Pope_Fry1997[], bbw_Morel_1974[], a0_Lee[], a1_Lee[],
  bottom_type_sand[], bottom_type_seagrass[], bottom_type_coral[], bottom_type_macrophytes[],
  bottom_type_dark_sediment[], bottom_type_coral_sand[], bottom_type_green_algae[],
  bottom_type_red_algae[];

#define MAX_P 512
#define N_TYPES 8

static model_data *lab_md;
static double lab_full[MAX_P];
static int lab_free[MAX_P], lab_n_free;

// Error as a function of the free parameters only.

static double free_error(double *x, model_data *md) {
  int i;
  for (i = 0; i < lab_n_free; i++) {
    lab_full[lab_free[i]] = x[i];
  }
  return samodel_error(lab_full, lab_md);
}

static int read_doubles(char **s, double *v, int n) {
  int i, used;
  for (i = 0; i < n; i++) {
    if (sscanf(*s, "%lf%n", &v[i], &used) != 1) return 0;
    *s += used;
  }
  return 1;
}

static void print_errors(model_data *md) {
  int k_scene, k_band, k_region;
  double K_min = 1.0e4;
  for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
    for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
      if (md->K[k_scene][k_band] > 0.0 && md->K[k_scene][k_band] < K_min) K_min = md->K[k_scene][k_band];
    }
  }
  printf("model_error %.10g Rrs_error %.10g depth_error %.10g bottom_error %.10g K_error %.10g K_min %.10g\n",
    md->model_error, md->Rrs_error, md->depth_error, md->bottom_error, md->K_error, K_min);
  for (k_region = 0; k_region < md->n_regions; k_region++) {
    for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
      printf("Rrs_modelled %d %d", k_region, k_scene);
      for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
        printf(" %.10g", md->Rrs_modelled[k_region][k_scene][k_band]);
      }
      printf("\n");
    }
  }
}

int main(void) {

  double *bottom_tables[N_TYPES] = {bottom_type_sand, bottom_type_seagrass, bottom_type_coral,
    bottom_type_macrophytes, bottom_type_dark_sediment, bottom_type_coral_sand,
    bottom_type_green_algae, bottom_type_red_algae};
  model_data *md;
  int n_scenes, n_bands, n_regions, origin, n_bottoms, k, k_scene, k_band, k_region, i,
    n_params, used, konvge, kcount, icount, numres, ifault;
  double theta_sun, theta_view, tide, v[MAX_P], step[MAX_P], x[MAX_P], xmin[MAX_P], xstep[MAX_P],
    reqmin, error, h;
  char line[65536], word[32], *s;

  // Setup.

  if (fgets(line, sizeof(line), stdin) == NULL ||
      sscanf(line, "SETUP %d %d %d %d %d", &n_scenes, &n_bands, &n_regions, &origin, &n_bottoms) != 5) {
    printf("ERROR: expected SETUP\n");
    return 1;
  }

  md = calloc(1, sizeof(model_data));
  lab_md = md;
  md->n_scenes = n_scenes;
  md->n_regions = n_regions;
  md->max_n_regions = n_regions;
  md->max_n_bands = n_bands;
  md->origin = origin;
  md->n_bottoms = n_bottoms;
  md->n_bands = malloc(n_scenes*sizeof(int));
  md->n_raw_bands = malloc(n_scenes*sizeof(int));
  for (k_scene = 0; k_scene < n_scenes; k_scene++) {
    md->n_bands[k_scene] = n_bands;
    md->n_raw_bands[k_scene] = n_bands;
  }

  allocate_double_array_3d(&md->Rrs_measured, n_regions, n_scenes, n_bands);
  allocate_double_array_3d(&md->Rrs_modelled, n_regions, n_scenes, n_bands);
  allocate_double_array_3d(&md->rrs_modelled, n_regions, n_scenes, n_bands);
  allocate_double_array_3d(&md->rrs_bottom, n_regions, n_scenes, n_bands);
  allocate_double_array_3d(&md->rrs_dp, n_regions, n_scenes, n_bands);
  allocate_double_array_3d(&md->rho, n_regions, n_scenes, n_bands);
  allocate_double_array_3d(&md->bottom_reflectance, N_TYPES, n_scenes, n_bands);
  allocate_double_array_2d(&md->wavelengths, n_scenes, n_bands);
  allocate_double_array_2d(&md->raw_wavelengths, n_scenes, n_bands);
  allocate_double_array_2d(&md->a_0, n_scenes, n_bands);
  allocate_double_array_2d(&md->a_1, n_scenes, n_bands);
  allocate_double_array_2d(&md->a_w, n_scenes, n_bands);
  allocate_double_array_2d(&md->b_bw, n_scenes, n_bands);
  allocate_double_array_2d(&md->K, n_scenes, n_bands);
  allocate_double_array_2d(&md->rrs_noise, n_scenes, n_bands);
  allocate_double_array_2d(&md->Rrs440, n_regions, n_scenes);
  allocate_double_array_2d(&md->Rrs490, n_regions, n_scenes);
  allocate_double_array_2d(&md->Rrs550, n_regions, n_scenes);
  allocate_double_array_2d(&md->Rrs640, n_regions, n_scenes);
  allocate_double_array_2d(&md->Rrs750, n_regions, n_scenes);
  md->H_tide = calloc(n_scenes, sizeof(double));
  md->theta_sun = calloc(n_scenes, sizeof(double));
  md->theta_view = calloc(n_scenes, sizeof(double));
  md->sec_theta_sun = calloc(n_scenes, sizeof(double));
  md->sec_theta_view = calloc(n_scenes, sizeof(double));
  md->P = calloc(n_scenes, sizeof(double));
  md->G = calloc(n_scenes, sizeof(double));
  md->X = calloc(n_scenes, sizeof(double));
  md->D = calloc(n_scenes, sizeof(double));
  md->prev = calloc(4*n_scenes, sizeof(double));
  md->a_w640 = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, 640.0);

  for (k = 0; k < N_TYPES; k++) {
    strcpy(md->B_type_name[k], "type");
  }

  while (fgets(line, sizeof(line), stdin) != NULL) {
    if (strncmp(line, "END", 3) == 0) break;
    s = line;
    if (sscanf(s, "%31s%n", word, &used) != 1) continue;
    s += used;
    if (strcmp(word, "SCENE") == 0) {
      if (sscanf(s, "%d %lf %lf %lf%n", &k_scene, &theta_sun, &theta_view, &tide, &used) != 4 ||
          k_scene < 0 || k_scene >= n_scenes) {
        printf("ERROR: bad SCENE\n");
        return 1;
      }
      s += used;
      if (! read_doubles(&s, v, n_bands)) { printf("ERROR: bad SCENE wavelengths\n"); return 1; }
      md->H_tide[k_scene] = tide;
      md->theta_sun[k_scene] = theta_sun*PI/180.0;
      md->theta_view[k_scene] = theta_view*PI/180.0;
      md->sec_theta_sun[k_scene] = 1.0/cos(md->theta_sun[k_scene]);
      md->sec_theta_view[k_scene] = 1.0/cos(md->theta_view[k_scene]);
      for (k_band = 0; k_band < n_bands; k_band++) {
        md->wavelengths[k_scene][k_band] = v[k_band];
        md->raw_wavelengths[k_scene][k_band] = v[k_band];
        md->a_0[k_scene][k_band] = interp_1d(ref_wlens, a0_Lee, n_ref_wlens, v[k_band]);
        md->a_1[k_scene][k_band] = interp_1d(ref_wlens, a1_Lee, n_ref_wlens, v[k_band]);
        md->a_w[k_scene][k_band] = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, v[k_band]);
        md->b_bw[k_scene][k_band] = interp_1d(ref_wlens, bbw_Morel_1974, n_ref_wlens, v[k_band]);
        for (k = 0; k < N_TYPES; k++) {
          md->bottom_reflectance[k][k_scene][k_band] = interp_1d(ref_wlens, bottom_tables[k], n_ref_wlens, v[k_band]);
        }
      }
    } else if (strcmp(word, "RRS") == 0) {
      if (sscanf(s, "%d %d%n", &k_region, &k_scene, &used) != 2 ||
          k_region < 0 || k_region >= n_regions || k_scene < 0 || k_scene >= n_scenes) {
        printf("ERROR: bad RRS\n");
        return 1;
      }
      s += used;
      if (! read_doubles(&s, v, n_bands)) { printf("ERROR: bad RRS values\n"); return 1; }
      for (k_band = 0; k_band < n_bands; k_band++) {
        md->Rrs_measured[k_region][k_scene][k_band] = v[k_band];
      }
    }
  }

  // As at the start of samodel_optimise().

  for (k_region = 0; k_region < n_regions; k_region++) {
    for (k_scene = 0; k_scene < n_scenes; k_scene++) {
      md->Rrs440[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], n_bands, 440.0);
      md->Rrs490[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], n_bands, 490.0);
      md->Rrs550[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], n_bands, 550.0);
      md->Rrs640[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], n_bands, 640.0);
      md->Rrs750[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], n_bands, 750.0);
      if (md->Rrs440[k_region][k_scene] < 0.0) md->Rrs440[k_region][k_scene] = 0.0001;
    }
  }

  for (k = 0; k < n_bottoms; k++) {
    md->bottom_type_indexes[k] = k;
  }
  n_params = samodel_n_params(md);
  md->n_params = n_params;
  if (n_params > MAX_P) { printf("ERROR: too many parameters\n"); return 1; }

  // Commands.

  while (fgets(line, sizeof(line), stdin) != NULL) {
    s = line;
    if (sscanf(s, "%31s%n", word, &used) != 1) continue;
    s += used;

    if (strcmp(word, "EVAL") == 0) {
      if (! read_doubles(&s, v, n_params)) { printf("ERROR: bad EVAL\n"); continue; }
      md->n_bottoms = n_bottoms;
      md->n_params = n_params;
      error = samodel_error(v, md);
      print_errors(md);

    } else if (strcmp(word, "NM") == 0) {
      if (sscanf(s, "%lf %d %d%n", &reqmin, &konvge, &kcount, &used) != 3) { printf("ERROR: bad NM\n"); continue; }
      s += used;
      if (! read_doubles(&s, v, n_params) || ! read_doubles(&s, step, n_params)) { printf("ERROR: bad NM\n"); continue; }
      md->n_bottoms = n_bottoms;
      md->n_params = n_params;
      lab_n_free = 0;
      for (i = 0; i < n_params; i++) {
        lab_full[i] = v[i];
        if (step[i] != 0.0) {
          x[lab_n_free] = v[i];
          xstep[lab_n_free] = step[i];
          lab_free[lab_n_free++] = i;
        }
      }
      icount = numres = ifault = 0;
      if (lab_n_free > 0) {
        nelmin(free_error, md, lab_n_free, x, xmin, &error, reqmin, xstep, konvge, kcount,
          &icount, &numres, &ifault);
        for (i = 0; i < lab_n_free; i++) lab_full[lab_free[i]] = xmin[i];
      }
      error = samodel_error(lab_full, md);
      printf("NM error %.10g icount %d numres %d ifault %d\nparams", error, icount, numres, ifault);
      for (i = 0; i < n_params; i++) printf(" %.10g", lab_full[i]);
      printf("\n");
      printf("depths");
      for (i = 0; i < n_regions; i++) printf(" %.10g", samodel_depth(lab_full, i));
      printf("\n");
      print_errors(md);

    } else if (strcmp(word, "START") == 0) {

      // photic's starting point and steps for depth h (and optionally albedo).

      double albedo = 0.3;
      if (sscanf(s, "%lf %lf", &h, &albedo) < 1) { printf("ERROR: bad START\n"); continue; }
      md->n_bottoms = n_bottoms;
      samodel_initial(md, h, albedo, v, step);
      printf("START");
      for (i = 0; i < n_params; i++) printf(" %.10g", v[i]);
      printf("\nSTEP");
      for (i = 0; i < n_params; i++) printf(" %.10g", step[i]);
      printf("\n");

    } else if (strcmp(word, "PACK") == 0) {

      // Physical values to a parameter vector:
      // PACK H_1..H_nr  albedo_1..albedo_nr  q_11..q_nr,nb  P_1 G_1 X_1 .. per scene

      double pH[MAX_P], pA[MAX_P], pq[MAX_P], pP[MAX_P], pG[MAX_P], pX[MAX_P], pgx[3*MAX_P];
      md->n_bottoms = n_bottoms;
      if (! read_doubles(&s, pH, n_regions) || ! read_doubles(&s, pA, n_regions) ||
          ! read_doubles(&s, pq, n_regions*n_bottoms) || ! read_doubles(&s, pgx, 3*n_scenes)) {
        printf("ERROR: bad PACK\n");
        continue;
      }
      for (k_scene = 0; k_scene < n_scenes; k_scene++) {
        pP[k_scene] = pgx[3*k_scene];
        pG[k_scene] = pgx[3*k_scene + 1];
        pX[k_scene] = pgx[3*k_scene + 2];
      }
      samodel_pack(md, pH, pA, pq, pP, pG, pX, v);
      printf("PACKED");
      for (i = 0; i < n_params; i++) printf(" %.10g", v[i]);
      printf("\n");

    } else if (strcmp(word, "PHOTIC") == 0) {
      if (sscanf(s, "%lf", &h) != 1) { printf("ERROR: bad PHOTIC\n"); continue; }
      md->n_bottoms = n_bottoms;
      md->start_at_previous = false;
      md->depth_prev = 0.0;
      if (h < 0.0) {
        md->empirical_depth_present = false;
        md->h_empirical = 0.0;
      } else {
        md->empirical_depth_present = true;
        md->h_empirical = h;
      }
      samodel_optimise(md);
      printf("PHOTIC depth %.10g albedo %.10g K_min %.10g iterations %d converged %d n_params %d\n",
        md->depth, md->bottom_albedo, md->K_min, md->n_iterations, md->converged, md->n_params);
      printf("PGX");
      for (k_scene = 0; k_scene < n_scenes; k_scene++) printf(" %.10g %.10g %.10g", md->P[k_scene], md->G[k_scene], md->X[k_scene]);
      printf("\nbottom_percent");
      for (k = 0; k < N_TYPES; k++) printf(" %.4g", md->B_type_percent[k]);
      printf("\n");
      print_errors(md);
    }
    printf("DONE\n");
    fflush(stdout);
  }

  return 0;
}
