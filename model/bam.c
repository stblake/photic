/*
 *
 *
 *                              P H O T I C
 *                              -----------
 *
 *
 *    A Multi-Spectral, Semi-Analytical, Shallow Water, Bathymetry Model.
 * 
 *
 *
 *
 * 
 */

/* Written by Sam Blake. */

/* Started on 27 March 2016. */

/*
  Updates
  -------
    1.01, 20/4/2017 - Fixed a bug in the estimation of attenuation coefficients. Included 
            a correlation plot with a line of best fit in COMPUTE K. 

    1.02, 22/4/2017 - Removed the rescaling of pan-sharpened imagery. 
    
    1.03, 24/4/2017 - Complete rewrite of COMPUTE LSB, now called via COMPUTE BOTTOM. We 
            now use the soundings and rewrite L_i = Ls_i + Lb_i*e^{-\alpha*d} as 
            Lb_i = e^{\alpha*d}*(L_i - Ls_i) and compute Lb_i at each of the soundings in 
            optically shallow water. Finally Lb_i is computed via the mean of 
            Lb_i at each sounding and the standard deviation of Lb_i is also computed. 

    1.04, 28/4/2017 - Code for computing the water-type and corresponding consistently 
            estimating attenuation coefficients. Water type is stored as water_type, 
            where 0 - 1 corresponds to OI, 1 - 2 corresponds to O2, and so on 
            to 9 - 10 corresponding to C9. 

    1.05, 3/5/2017 - Moved away from bicubic interpolation in favour of bilinear 
            interpolation to preserve as many optically shallow points as possible. 

    1.06, 11/5/2017 - Fixed division by zero bug in computing mean relative error. 

    1.07, 16/5/2017 - Fixed a major bug in the computation of depths whereby the derived 
            coefficient of the logarithmically transformed radiance was mismatches under 
            certain conditions. This occured in the routine zdepth_lyzenga_log_linear. 

    1.08, 17/5/2017 - Added LSB_SIGMA flag. 

    1.09, 3/6/2017 - Updated Land/Sea masking algorithm. It now takes in the green and NIR 
            bands and computes the so-called "Modified Normalized Difference Water Index". 

    1.10, 7/6/2017 - Included the regression line slope as a weighted input to the 
            model error function. 

    1.11, 13/6/2017 - Increased SA_MAX_ITERATIONS, modified the restarting mechanism 
            within the simulating annealing code to restart to the best solution instead 
            of a randomly generated solution. Rewrote pertubate_vec to discretize the 
            coefficients. Added a non-linear terms to the depth calculation. Modified the 
            constraints on the logarithmic coefficients to demand the first is 
            negative and the rest are positive. 

    1.12, 16/6/2017 - Updated PROCESS POLYGONS to exclude either interior or exterior points. 

    1.13, 19/6/2017 - Implemented estimation of diffuse attenuation coefficient per the 
            Mueller algorithm. This is given by Kd_490_Mueller. Also implemented estimation 
            of chlorophyll-a content per Morel, 1988. 

    1.14, 24/6/2017 - Implemented Lee's semi-analytical model for diffuse attenuation 
            coefficients and Secchi disk depth. Implemented in MODEL Lee_Kd_LS8 and 
            MODEL Lee_Secchi_LS8. Also added a pseudo-Secchi disk depth calculation 
            MODEL Kd_490_Mueller_Secchi - which is 1.7 over the Kd_490 attenuation coefficient. 

    1.50, 16/7/2017 - Implemented Lee's semi-analytical model. A large number of additional 
            improvements. All legacy functionality is still present. 

    1.60, 26/7/2017 - Implemented SCENE encapsulation and recoded the semi-analytical model 
            to use multiple scenes. 

    1.62, 20/8/2017 - Lots of improvements to the SA model. We have support for processing 
            multiple scenes. Started implementing an experimental depth calibration based 
            on the model derived by Maritorena, 1994.

    1.70, 2/9/2017 - Successfully tested hybrid model. Polcyn's log ratio model is used to 
            normalise the depths in semi-analytical model. This allows the semi-analytical 
            model to retrieve stable depths without a large number of spectral 
            bands (such as in Sentinel-2). 

    1.71, 4/9/2017 - Incorporated tidal heights into the semi-analytical model. 

    1.72, 6/9/2017 - Moved the Polcyn log ratio to use rrs instead of Rrs. 

    1.73, 11/9/2017 - Fixed a bug in PROCESS UNITE for METHOD either MIN or MAX. 

    1.74, 17/12/2017 - Added additional constraints on the Lyzenga log coefficients. 

    1.75, 22/1/2018 - Moved all the PROCESS commands into process.c. 

    1.80, 28/1/2018 - Some updates to the semi-analytical model. Added code to interpolate
            input spectral bands to a set of user-specified spectral bands. The idea is to 
            increase the convergence reliability of the optimisation scheme. This is still
            highly experimental. Still trying to reconcile the SA depths with log-ratio 
            depths. 

    1.81, 30/1/2018 - Fixed a bug in the linear regression slope error estimate in the 
            empirical model. 

    1.82, 28/2/2018 - Improved the model depth error estimate in the empirical model.

    1.83, 23/10/2018 - Cleaned up some code in samodel prior to refactoring. 

    1.84, 15/11/2018 - 

    1.85, 17/11/2018 - INTERPOLATE now takes a grid as the REGION spec. In such cases it will
            interpolate onto the region of the specified grid. Added WEIGHTS option to PROCESS 
            UNITE. 

    1.86, 4/12/2018 - Fixed a bug in INTERPOLATE when REGION is not a grid. 

    1.87, 6/12/2018 - Implemented 2D histogram for attenutation coefficient estimation. 

    1.88, 20/02/2019 - Implemented COMPARE. 

    1.89, 3/3/2019 - Improved graphics in COMPUTE K, and fixed a bug in linear_fit and linear_fit2. 

    1.90, 25/3/2019 - Added a variant of the semi-analytical model where R_deep and K are 
            stastically estimated from the imagery prior to inverting the RTE. This model produces 
            the depth, z, the one sigma depth error, z_sigma, and the bottom type. 

    1.95, 26/6/2019 - Greatly rewritten and designed the semi-analytical model. 
*/

/*

  TODOs: 
          * Store arrays (internally within memory of BAM) as short 
              ints with scale_factor and add_offset. This will decrease RAM
              consumption by 50%. Then arrays are accessed with something like: 

              #define array_elem(arrno, i, j) \ 
                (gridded_data[arrno].scale_factor*((float) gridded_data[arrno].array[i][j]) + gridded_data[arrno].add_offset)

              Furthermore, they could be stored at 8-bit ints (unsigned char or uint8_t with -std=c99) with further 
              significant reductions in memory. 

*/

#include "bam.h"

#define BAM_VERSION "1.95"


void interpreter();
int parse_input(char line[MAX_LINE]);
void tokenize(char line[MAX_LINE]);
void run_help();
void run_command();
void run_read();
void run_write();
void run_write_grid();
void run_write_line();
void run_write_points();
void run_set();
void run_print();
void run_plot();
void run_model();
void run_history();
void run_copy();
void run_scene();
void run_copy_lsb();
void run_copy_lsm();
void run_read_commands();
void run_read_grid();
void run_read_soundings();
void run_read_landsat();
void run_read_bottom();
void run_delete();
void run_plot_grid();
void run_plot_line();
void run_plot_logscatter();
void run_plot_qa();
void run_plot_correlation();
void run_model_lyzenga1985();
void run_model_stumpf2003();
void run_model_spline();
void run_model_mueller();
void run_model_mueller_secchi();
void run_model_lee_semi_analytical();
void run_model_morel();
void run_model_lee_kd();
void run_model_lee_zsd();
void random_search();
void copy_vec(float *source, float *dest, int len);
void pertubate_vec(float *vec, int len, float scale, float scale0);
float model_error();
void pertubate_vec2(float *vec, int len, float scale);
int run_input_command_file();
void run_empirical_model();
void simulated_annealing(int nparams);
float zdepth_blake(int i, int j);
float zdepth_polcyn_log_ratio(int i, int j);
void lyzenga_log_linear();
float zdepth_lyzenga_log_linear(int i, int j);
void polcyn_log_ratio();
float bottom_type(int i, int j);
float geodistance(float lon1, float lat1, float lon2, float lat2);
void blake_local();
float radiance_min_cutoff(float **array, int nrows, int ncols, float spval,
    int i_min, int j_min, int i_max, int j_max);
void run_compute();
void run_compute_area();
void run_compute_lsm();
void run_compute_lsb();
void run_compute_bottom();
void run_compute_k();
void run_compute_percentile();
void compute_weights(float **weights, float **soundings, int nsoundings);
bool point_in_polygon(float *x, float *y, int n, float pt_x, float pt_y);
void run_compare();
void run_combine();
// void run_create_empirical();


#define OK 1000 
#define COMMENT 1001
#define EXIT 1002
#define READ 1003
#define WRITE 1004
#define STOP 1005

int main(int argc, char **argv) {

  int n, ret_code;

  printf("\n\n                                   --    P H O T I C    -- \n\n\n");
  printf("      A Multi-Spectral, Multi-Temporal, Semi-Analytical, Shallow Water, Bathymetry Model.\n\n");
  printf("Version %s (%s)\n\n", BAM_VERSION, __TIMESTAMP__);

  printf("\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR \n\
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, \n\
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT \n\
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE \n\
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, \n\
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER \n\
DEALINGS IN THE SOFTWARE.\n\n");

// Check program expiry dates. 

  // TODO

// Set model defaults. 

  for (n = 0; n < MAX_SCENES; n++) {
    allocated_scenes[n] = false;
  }

  for (n = 0; n < MAX_GRIDS; n++) {
    allocated_grids[n] = false;
    thetas[n] = 0.0;
    Lsm[n] = 0.0;
    Lsm_sigma[n] = 0.0;
    Lsb[n] = 0.0;
    Lsb_sigma[n] = 0.0;
    Alpha[n] = 0.0;
    Alpha_sigma[n] = 0.0; 
    Rrs_sigma[n] = 0.0;
    Hj[n]  = 0.0;
    NAlpha[n] = 0;
    model_wavelengths[n] = 0.0;
    radiance_mult_band[n]    = BAM_NODATA;
    radiance_add_band[n]     = BAM_NODATA;
    quantize_cal_min_band[n] = BAM_NODATA;
    quantize_cal_max_band[n] = BAM_NODATA;
    radiance_minimum_band[n] = BAM_NODATA;
    radiance_maximum_band[n] = BAM_NODATA;
    K_ratios[n][0] = 0.0;
    K_ratios[n][1] = 0.0;
    K_ratios[n][2] = 0.0;
    K_ratios[n][3] = 0.0;
  }

  for (n = 0; n < MAX_SCENES; n++) {
    allocated_scenes[n] = false;
  }

  zdepth = NULL;

// Set previous inputs. 

  for (n = 0; n < MAX_STRING_LEN; n++) {
    prev1[n] = '\0';
    prev2[n] = '\0';
    prev3[n] = '\0';
    prev4[n] = '\0';
    prev5[n] = '\0';
    prev6[n] = '\0';
    prev7[n] = '\0';
    prev8[n] = '\0';
    prev9[n] = '\0';
  }

// Read command line inputs 

  if (argc == 2) {
  	strcpy(command_file, argv[1]);
  	ret_code = run_input_command_file();
    if (ret_code == EXIT) {
      return 0;
    }
  } else {
  	for (n = 0; n < argc; n++) {
  		// ...
  	}
  }

// Call the interpreter. 

  interpreter();

	return 0;
}



int run_input_command_file() {

	char line[MAX_LINE] = "READ COMMANDS ";
	int n, code;

// Check input file is present.

	if (! file_exists(command_file)) {
		printf("\nERROR: missing file '%s'\n", command_file);
  	return 0;
	}

	strcat(line, trim(command_file));
	code = parse_input(line);
	run_command();

  return code;
}

void interpreter() {
  
  char line[MAX_LINE], *raw_line;
  int code, n;
  FILE *logfile;

// Load the history at startup.

  linenoiseHistoryLoad("history.bam"); 

// Main loop. 

  while (true) {

    // Clear input line.
    for (n = 0; n < MAX_LINE; n++)
      line[n] = '\0';

  	// Read input from user. 

    raw_line = linenoise(">> ");

    if (strlen(raw_line) > MAX_LINE) {
      printf("\nERROR: input commands too long.\n");
      continue ;
    }

    linenoiseHistoryAdd(raw_line);
    linenoiseHistorySave("history.bam");

    strcpy(line, raw_line);
    trim(line);

    free(raw_line);

    // Check for %, %n.

    if (line[0] == '%') {
      if (line[1] == '\0' || line[1] == '1')
        strcpy(line, prev1);
      else if (line[1] == '2')
        strcpy(line, prev2);
      else if (line[1] == '3')
        strcpy(line, prev3);
      else if (line[1] == '4')
        strcpy(line, prev4);
      else if (line[1] == '5')
        strcpy(line, prev5);
      else if (line[1] == '6')
        strcpy(line, prev6);
      else if (line[1] == '7')
        strcpy(line, prev7);
      else if (line[1] == '8')
        strcpy(line, prev8);
      else if (line[1] == '9')
        strcpy(line, prev9);
    }

    strcpy(prev9, prev8);
    strcpy(prev8, prev7);
    strcpy(prev7, prev6);
    strcpy(prev6, prev5);
    strcpy(prev5, prev4);
    strcpy(prev4, prev3);
    strcpy(prev3, prev2);
    strcpy(prev2, prev1);
    strcpy(prev1, line);

    // Parse input.

  	code = parse_input(line);

  	switch (code) {
  		case COMMENT:
  			break ;
  		case EXIT:
  			printf("\nAre you sure you want to exit? (yes/no) ");
  			fgets(line, sizeof(line), stdin);
  			if (strncmp(line, "yes", 3) == 0) {
  				goto finish; 
        }
  			break ;
  		default:
  			run_command();
  	}

    fflush(stdout);
  }

finish: ;

// Clean up memory.

	for (n = 0; n < MAX_GRIDS; n++) {
		if (allocated_grids[n]) {
			free_float_array_2d(gridded_data[n].array, gridded_data[n].nrows);
		}
	}

}


void run_command() {

// Top level commands. 

	if (strcmp(parsed[0], "READ") == 0) {
		run_read();
	} else if (strcmp(parsed[0], "WRITE") == 0) {
		run_write();
	} else if (strcmp(parsed[0], "SET") == 0) {
		run_set();
	} else if (strcmp(parsed[0], "PRINT") == 0) {
		run_print();
	} else if (strcmp(parsed[0], "PLOT") == 0) {
		run_plot();
	} else if (strcmp(parsed[0], "PROCESS") == 0) {
		run_process();
	} else if (strcmp(parsed[0], "MODEL") == 0) {
		run_model();
	} else if (strcmp(parsed[0], "REFINE") == 0){
		run_refine();
	} else if (strcmp(parsed[0], "HELP") == 0) {
		run_help();
	} else if (strcmp(parsed[0], "DELETE") == 0) {
		run_delete();
  } else if (strcmp(parsed[0], "HISTORY") == 0) {
    run_history();
  } else if (strcmp(parsed[0], "COMPUTE") == 0) {
    run_compute();
  } else if (strcmp(parsed[0], "COPY") == 0) {
    run_copy();
  } else if (strcmp(parsed[0], "SCENE") == 0) {
    run_scene();
  } else if (strcmp(parsed[0], "COMPARE") == 0) {
    run_compare();
  } else if (strcmp(parsed[0], "COMBINE") == 0) {
    run_combine();
  } else if (strcmp(parsed[0], "ALIGN") == 0) {
    run_align();
  } else if (strcmp(parsed[0], "CREATE") == 0) {
    // run_create_empirical();
	} else if (strcmp(parsed[0], "") != 0) {
		printf("\nERROR: unknown command '%s'\n", parsed[0]);
	}
}



//
//  CREATE 
//

// CREATE /blue/ /green/ /red/ /bathy/

// M1 = max(Rblue/Rgreen, Rgreen/Rblue)
// M2 = max(Rblue/Rred, Rred/Rblue)

// z == c[0] log(M1) + c[1] log(M2) + c[2] log(M3)

#if 0
void run_create_empirical() {



}
#endif


//
//  COMBINE
//

// COMBINE IN /scene 1/ /scene 2/ ... OUT /output scene/ [optional] INDEXES /blue band index/
//    /green band index/ /red band index/

void run_combine() {

  int i, j, k, n, m, nscenes, nbands, nrows, ncols, ps_output_scene, ps_in[MAX_SCENES], 
    ps_out, best_scene, blue_band_index, green_band_index, red_band_index, src_index, dest_index;
  float Rsigma, Ksigma, Kmin, H_tide_mean, Rinf_mean, nobs, ***arrays, r, lowest, R_blue, R_green, R_red;
  bool success, nodata;
  char band_name[32];

  // Defaults. 

  blue_band_index  = 1;
  green_band_index = 2;
  red_band_index   = 3;

  // Read input scenes. 

  if (strcmp(parsed[1], "IN") != 0) {
    printf("\nERROR: COMBINE IN /scene 1/ /scene 2/ ... OUT /output scene/ [optional] INDEXES /blue band index/ \
/green band index/ /red band index/ \n");
    return ;
  } 

  nscenes = 0;
  k = 2;

  while (k < MAX_SCENES && strcmp(parsed[k], "OUT") != 0) {

    success = parse_old_scene(&ps_in[nscenes++], k++);
    printf("\nps_in[%d] = %d\n", nscenes-1, ps_in[nscenes-1]);
    if (! success) return ;
  }

  printf("\nn_scenes = %d", nscenes);

  // Check spectral resolution is commensurate. 

  for (n = 0; n < nscenes - 1; n++) {
    for (m = n + 1; m < nscenes; m++) {
      if (scene_data[ps_in[n]].n_bands != scene_data[ps_in[m]].n_bands) {
        printf("\n\nERROR: input scenes have different spectral resolutions.\n\n");
        return ;
      }
    }
  }  

  // Check input grids are commensurate. 

  for (n = 0; n < nscenes; n++) {
    for (m = 0; m < nscenes; m++) {
      if (! commensurate_grids(
          gridded_data[scene_data[ps_in[n]].band_indexes[0]], 
          gridded_data[scene_data[ps_in[m]].band_indexes[0]])) {
        printf("\nERROR: scenes are not spatially commensurate.\n");
        return ; 
      }
    }
  }

  // Create output scene. 

  if (strcmp(parsed[k], "OUT") != 0) {
    printf("\nERROR: COMBINE IN /scene 1/ /scene 2/ ... OUT /output scene/ [optional] INDEXES /blue band index/ \
/green band index/ /red band index/ \n");
    return ;
  } 

  success = parse_new_scene(&ps_output_scene, ++k);
  if (! success) return ;

  strcpy(scene_data[ps_output_scene].scene_name, parsed[k]);

  // [optional] INDEXES. 

  if (strcmp(parsed[++k], "INDEXES") == 0) {
    blue_band_index  = atoi(parsed[++k]);
    green_band_index = atoi(parsed[++k]);
    red_band_index   = atoi(parsed[++k]);
  }

  scene_data[ps_output_scene].n_bands = scene_data[ps_in[0]].n_bands;
  nbands = scene_data[ps_output_scene].n_bands;
  printf("\nn_bands = %d\n", nbands);

  scene_data[ps_output_scene].nrows = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nrows;
  scene_data[ps_output_scene].ncols = gridded_data[scene_data[ps_in[0]].band_indexes[0]].ncols;
  nrows = scene_data[ps_output_scene].nrows;
  ncols = scene_data[ps_output_scene].ncols;
  printf("\nnrows,ncols = %d,%d\n", nrows, ncols);

  printf("\nwavelengths = ");
  for (n = 0; n < nbands; n++) {
    scene_data[ps_output_scene].wavelengths[n] = scene_data[ps_in[0]].wavelengths[n];
    printf("%d ", scene_data[ps_output_scene].wavelengths[n]);
  }
  printf("\n");

  // theta_v/theta_w

  scene_data[ps_output_scene].theta_v = 0.0;
  scene_data[ps_output_scene].theta_w = 0.0;
  nobs = 0.0;
  for (m = 0; m < nscenes; m++) {
    scene_data[ps_output_scene].theta_v += scene_data[ps_in[0]].theta_v;
    scene_data[ps_output_scene].theta_w += scene_data[ps_in[0]].theta_w;
    nobs += 1.0;
  }
  scene_data[ps_output_scene].theta_v /= nobs;
  scene_data[ps_output_scene].theta_w /= nobs;

  // R_inf

    printf("\nR_inf = ");
    for (n = 0; n < nbands; n++) {
      for (m = 0; m < nscenes; m++) {
        Rinf_mean = 0.0;
        nobs = 0.0;
        if (! approx_equal(scene_data[ps_in[m]].R_inf[n], 0.0, 1.0e-6)) {
          Rinf_mean += scene_data[ps_in[m]].R_inf[n];
          nobs += 1.0;
        }
      }

      if (nobs > 0.5) {
        scene_data[ps_output_scene].R_inf[n] = Rinf_mean/nobs;
      } else {
        scene_data[ps_output_scene].R_inf[n] = 0.0;
      }
      printf("%f ", scene_data[ps_output_scene].R_inf[n]);
    }
    printf("\n");

  // H_tide. 

  H_tide_mean = 0.0;
  nobs = 0.0;
  for (m = 0; m < nscenes; m++) {
    if (! approx_equal(scene_data[ps_in[m]].H_tide, 0.0, 1.0e-6)) {
      H_tide_mean += scene_data[ps_in[m]].H_tide;
      nobs += 1.0;
    }
  }

  if (nobs > 0.5) {
    scene_data[ps_output_scene].H_tide = H_tide_mean/nobs;
  } else {
    scene_data[ps_output_scene].H_tide = 0.0;
  }
  printf("\n\nH_tide = %f", scene_data[ps_output_scene].H_tide);

  // K. 

  printf("\n\nK = ");
  for (n = 0; n < nbands; n++) {
    Kmin = 100.0;

    for (m = 0; m < nscenes; m++) {
      if (! approx_equal(scene_data[ps_in[m]].K[n], 0.0, 1.0e-6)) {
        if (scene_data[ps_in[m]].K[n] < Kmin) {
          Kmin = scene_data[ps_in[m]].K[n];
        }
      }
    }

    scene_data[ps_output_scene].K[n] = Kmin;
    printf("%f ", scene_data[ps_output_scene].K[n]);
  }

  // R_sigma.

  printf("\n\nR_sigma = ");
  for (n = 0; n < nbands; n++) {
    Rsigma = 0.0;
    nobs = 0.0;

    for (m = 0; m < nscenes; m++) {
      if (! approx_equal(scene_data[ps_in[m]].R_sigma[n], 0.0, 1.0e-6)) {
        Rsigma += scene_data[ps_in[m]].R_sigma[n];
        nobs += 1.0;
      }
    }

    if (nobs > 0.5) {
      scene_data[ps_output_scene].R_sigma[n] = Rsigma/nobs;
    } else {
      scene_data[ps_output_scene].R_sigma[n] = 0.0;
    }
    printf("%f ", scene_data[ps_output_scene].R_sigma[n]);
  }

  // K_sigma.

  printf("\n\nK_sigma = ");
  for (n = 0; n < nbands; n++) {
    Ksigma = 0.0;
    nobs = 0.0;

    for (m = 0; m < nscenes; m++) {
      if (! approx_equal(scene_data[ps_in[m]].K_sigma[n], 0.0, 1.0e-6)) {
        Ksigma += scene_data[ps_in[m]].K_sigma[n];
        nobs += 1.0;
      }
    }

    if (nobs > 0.5) {
      scene_data[ps_output_scene].K_sigma[n] = Ksigma/nobs;
    } else {
      scene_data[ps_output_scene].K_sigma[n] = 0.0;
    }
    printf("%f ", scene_data[ps_output_scene].K_sigma[n]);
  }
  printf("\n");

  // Create new spectral bands. 

  arrays = malloc(nbands*sizeof(float**));

  for (i = 0; i < nbands; i++) {

    sprintf(band_name, "combined_b%d", scene_data[ps_output_scene].wavelengths[i]);
    printf("\n%s\n", band_name);

    // Check output grid name is not in use. 

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], band_name) == 0) {
        printf("\nERROR: grid '%s' already in use.\n", band_name);
        return ;
      }
    }

    // Allocate grid into the first empty slot. 

    for (n = 0; n < MAX_GRIDS; n++) {
      if (! allocated_grids[n]) {
        ps_out = n;
        break;
      }
    }

    // Copy grid name. 

    strcpy(grid_names[ps_out], band_name);

    // This grid is now in use. 

    allocated_grids[ps_out] = true;

    // Allocate memory. 

    allocate_float_array_2d(&(arrays[i]), nrows, ncols);

    // Update memory usage (in units of megabytes). 

    meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

    // Update number of bands in use. 

    n_grids_in_use++;

    gridded_data[ps_out].array = arrays[i];
    gridded_data[ps_out].ncols = ncols; 
    gridded_data[ps_out].nrows = nrows; 
    gridded_data[ps_out].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon; 
    gridded_data[ps_out].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
    gridded_data[ps_out].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
    gridded_data[ps_out].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
    gridded_data[ps_out].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
    gridded_data[ps_out].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value;

    scene_data[ps_output_scene].band_indexes[i] = ps_out;
  }  

  // Combine scenes. The scene with largest (R_green - R_red)/(R_blue - R_red) is output. 

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {

      // Check for nodata in each band of each scene. 

      nodata = false;
      for (k = 0; k < nscenes; k++) {
        if (approx_equal(
              gridded_data[scene_data[ps_in[k]].band_indexes[blue_band_index]].array[i][j],
              gridded_data[scene_data[ps_in[k]].band_indexes[blue_band_index]].nodata_value, 1.0e-8) || 
            approx_equal(
              gridded_data[scene_data[ps_in[k]].band_indexes[green_band_index]].array[i][j],
              gridded_data[scene_data[ps_in[k]].band_indexes[green_band_index]].nodata_value, 1.0e-8) || 
            approx_equal(
              gridded_data[scene_data[ps_in[k]].band_indexes[red_band_index]].array[i][j],
              gridded_data[scene_data[ps_in[k]].band_indexes[red_band_index]].nodata_value, 1.0e-8)) {
          nodata = true;
          break ;
        }
      }

      if (nodata) {
        for (k = 0; k < nbands; k++) {
          dest_index = scene_data[ps_output_scene].band_indexes[k];
          gridded_data[dest_index].array[i][j] = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value;
        }
        continue ;
      }

      best_scene = -1;
      lowest = 1000000.0;

      for (k = 0; k < nscenes; k++) {

        // R_blue  = gridded_data[scene_data[ps_in[k]].band_indexes[blue_band_index]].array[i][j];
        R_green = gridded_data[scene_data[ps_in[k]].band_indexes[green_band_index]].array[i][j];
        // R_red   = gridded_data[scene_data[ps_in[k]].band_indexes[red_band_index]].array[i][j];

        r = R_blue/(1.0e-6 + R_green);

        if (r < lowest) {
          lowest = r;
          best_scene = ps_in[k];
        }
      }

      for (k = 0; k < nbands; k++) {
        src_index = scene_data[best_scene].band_indexes[k];
        dest_index = scene_data[ps_output_scene].band_indexes[k];
        gridded_data[dest_index].array[i][j] = gridded_data[src_index].array[i][j];
      }
    }
  }
}

//
//  COMPARE
//

// COMPARE GRID /grid name/ SOUNDINGS /sounding file/ [optional] CHARTDATUM /chart datum for file/ 
//    DEPTHRANGE /depth min/ /depth max/ SPVAL /no data value/ TABULATE NEGATEGRID NEGATESOUNDINGS

void run_compare() {

  bool present, tabulate, success, negate_grid, negate_soundings;
  int i, j, k, n, nrows, ncols, nobs, ps_grid, ps_grid2, n_files, n_soundings, n_valid_soundings;
  float lon, lat, *x, *y, *z_soundings, *z_model, *valid_model_depths, *valid_sounding_depths, 
    *rel_error_vec, *abs_error_vec, median_rel_error, median_abs_error,
    min_depth, max_depth, chart_datum, spval, mean_rel_error, mean_abs_error, 
    Xmin, Xmax, Ymin, Ymax, a, b, d, r, error, sign, depth_error, n_within, n_total; 

  // Read grid name. 

  if (strcmp(parsed[1], "GRID") != 0) {
    printf("\nERROR: COMPARE GRID /grid name/ SOUNDINGS /sounding file/ [optional]\
CHARTDATUM /chart datum for file/ DEPTHRANGE /depth min/ /depth max/ SPVAL /sounding no data value/ TABULATE NEGATEGRID NEGATESOUNDINGS\n");
    return ;
  } 

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[2]) == 0) {
      present = true;
      ps_grid = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
    return ;
  }

  // Read sounding files.

  if (strcmp(parsed[3], "SOUNDINGS") != 0) {
    printf("\nERROR: COMPARE GRID /grid name/ SOUNDINGS /sounding file/ [optional]\
CHARTDATUM /chart datum for file/ DEPTHRANGE /depth min/ /depth max/ SPVAL /sounding no data value/ TABULATE NEGATEGRID NEGATESOUNDINGS\n");
    return ;
  } 

// Read sounding data. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[4]) == 0) {
      present = true;
      ps_grid2 = n;
      break;
    }
  }

  if (present) {
    nrows = gridded_data[ps_grid2].nrows; 
    ncols = gridded_data[ps_grid2].ncols; 
    x = malloc(nrows*ncols*sizeof(float));
    y = malloc(nrows*ncols*sizeof(float));
    z_soundings = malloc(nrows*ncols*sizeof(float));

    k = 0; 
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        if (! approx_equal(gridded_data[ps_grid2].array[i][j], gridded_data[ps_grid2].nodata_value, 1.0e-4)) {
          x[k] = gridded_data[ps_grid2].wlon + ((float) j)*gridded_data[ps_grid2].cellsize;
          y[k] = gridded_data[ps_grid2].slat + ((float) i)*gridded_data[ps_grid2].cellsize;
          z_soundings[k] = fabs(gridded_data[ps_grid2].array[i][j]);
          k++;
        }
      }
    }
    n_soundings = k;
    printf("\nn_soundings = %d\n", k);
  } else if (! file_exists(parsed[4])) {
    printf("\nERROR: File %s does not exist!\n\n", parsed[4]);
    return ;
  } else {
    read_xyz(parsed[4], &x, &y, &z_soundings, &n_soundings);
  }

// Read CHARTDATUM. 

  k = 5; 

  if (strcmp(parsed[k], "CHARTDATUM") == 0) {
    chart_datum = atof(parsed[++k]);
    k++; 
  } else {
    chart_datum = 0.0;
  }

// Read DEPTHRANGE. 

  if (strcmp(parsed[k], "DEPTHRANGE") == 0) {
    min_depth = atof(parsed[++k]);
    max_depth = atof(parsed[++k]);
    k++;
    printf("\nmin/max depth = %.2f, %.2f\n", min_depth, max_depth);
  } else {
    min_depth = depthmin;
    max_depth = depthmax;
  }

// Read no data value. 

  if (strcmp(parsed[k], "SPVAL") == 0) {
    spval = atof(parsed[++k]);
    k++; 
  } else {
    spval = 0.0;
  }

// Read VERBOSE. 

  if (strcmp(parsed[k], "TABULATE") == 0) {
    tabulate = true; 
    k++;
  } else {
    tabulate = false;
  }

// Read NEGATEGRID. 

  if (strcmp(parsed[k], "NEGATEGRID") == 0) {
    negate_grid = true; 
    k++;
  } else {
    negate_grid = false;
  }

// Read NEGATESOUNDINGS. 

  if (strcmp(parsed[k], "NEGATESOUNDINGS") == 0) {
    negate_soundings = true; 
    k++;
  } else {
    negate_soundings = false;
  }

// Correct for chart datum. 

  for (i = 0; i < n_soundings; i++) {
    if (! approx_equal(z_soundings[i], spval, 1.0e-4)) {
      z_soundings[i] += chart_datum;
    }
  }


// Compute model depths at the soundings. 

  mean_rel_error = 0.0;
  mean_abs_error = 0.0;
  nobs = 0; 

  z_model = malloc(n_soundings*sizeof(float));

  rel_error_vec = malloc(n_soundings*sizeof(float));
  abs_error_vec = malloc(n_soundings*sizeof(float));

  for (i = 0; i < n_soundings; i++) {

    z_soundings[i] = fabs(z_soundings[i]);

    // Sounding is no data value. 

    if (approx_equal(z_soundings[i], spval, 1.0e-4)) {
      z_model[i] = gridded_data[ps_grid].nodata_value; 
      continue ;
    }

    // Sounding is too deep. 

    sign = (z_soundings[i] > 0) ? 1.0 : -1.0;

    if (negate_soundings) {
      z_soundings[i] *= sign;
    }

    if (fabs(z_soundings[i]) > max_depth || fabs(z_soundings[i]) < min_depth) {
      z_model[i] = gridded_data[ps_grid].nodata_value;
      continue ;
    }

    lon = x[i];
    lat = y[i];

    if (lon > gridded_data[ps_grid].wlon && 
        lon < gridded_data[ps_grid].elon && 
        lat > gridded_data[ps_grid].slat && 
        lat < gridded_data[ps_grid].nlat) {

      z_model[i] = fabs(interp_bilinear(gridded_data[ps_grid].array, 
          gridded_data[ps_grid].ncols,
          gridded_data[ps_grid].nrows, 
          gridded_data[ps_grid].wlon, 
          gridded_data[ps_grid].slat, 
          gridded_data[ps_grid].cellsize, 
          gridded_data[ps_grid].nodata_value, 
          lon, 
          lat));

    } else {
      z_model[i] = gridded_data[ps_grid].nodata_value;
    }

    // printf("\nz_model[i],z_model[i] = %.3f, %.3f", z_model[i], z_soundings[i]);

    if (! approx_equal(z_model[i], gridded_data[ps_grid].nodata_value, 1.0e-4)) {

      if (negate_grid) {
        z_model[i] *= -1.0; 
      }

      rel_error_vec[nobs] = 100.0*fabs(fabs(z_model[i]) - fabs(z_soundings[i]))/fabs(z_soundings[i]);
      abs_error_vec[nobs] = fabs(fabs(z_model[i]) - fabs(z_soundings[i]));

      mean_rel_error += rel_error_vec[nobs];
      mean_abs_error += abs_error_vec[nobs];
      nobs += 1;
    }
  }

  if (nobs > 0) {
    mean_rel_error /= (float) nobs;
    mean_abs_error /= (float) nobs; 
  }

  qsort(rel_error_vec, nobs, sizeof(float), compare_floats);
  median_rel_error = rel_error_vec[nobs/2];

  qsort(abs_error_vec, nobs, sizeof(float), compare_floats);
  median_abs_error = abs_error_vec[nobs/2];

  free(rel_error_vec);
  free(abs_error_vec);

  printf("\n N = %d", (int) nobs);

  printf("\n\nMEAN ABSOLUTE ERROR   = %.2f  (m)", mean_abs_error);
  printf("\nMEAN RELATIVE ERROR   = %.2f (%%)", mean_rel_error);

  printf("\n\nMEDIAN ABSOLUTE ERROR = %.2f  (m)", median_abs_error);
  printf("\nMEDIAN RELATIVE ERROR = %.2f (%%)\n\n", median_rel_error);

// Tabulate results at each sounding. 

  if (tabulate) {
    printf("\n         LON         |          LAT        |  SOUNDING DEPTH  |    MODEL DEPTH   |    ABS ERROR (m)    |    REL %% ERROR      ");
    printf("\n---------------------------------------------------------------------------------------------------------------------------------------");
    for (n = 0; n < n_soundings; n++) {
      if (approx_equal(z_model[n], gridded_data[ps_grid].nodata_value, 1.0e-3)) {
        continue ;
      }
      printf("\n %14.6f      | %14.6f      |    %8.2f      |    %8.2f      |    %8.2f         |    %8.2f         ", 
          x[n], y[n], z_soundings[n], z_model[n], 
          fabs(fabs(z_model[n]) - fabs(z_soundings[n])), 
          100.0*fabs(z_model[n] - fabs(z_soundings[n]))/fabs(z_soundings[n]));
    }
  } 

  // Display some simple summary statistics. 

  float ref_errors[] = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0};

  printf("\nSUMMARY OF ABSOLUTE ERRORS:\n\n");

  for (i = 0; i < 6; i++) {
    depth_error = ref_errors[i];
    n_total  = 0.0;
    n_within = 0.0; 
    for (n = 0; n < n_soundings; n++) {
      if (! approx_equal(z_model[n], gridded_data[ps_grid].nodata_value, 1.0e-3)) {
        n_total += 1.0;
        if (fabs(fabs(z_model[n]) - fabs(z_soundings[n])) < depth_error + epsilon) {
          n_within += 1.0;
        }
      }
    }

    printf("%5.2f%% within %5.2f (m)\n", 100.0*n_within/n_total, depth_error);
  }

  printf("\n");

  float ref_errors2[] = {2.0, 5.0, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0};

  printf("\nSUMMARY OF RELATIVE %% ERRORS:\n\n");

  for (i = 0; i < 8; i++) {
    depth_error = ref_errors2[i];
    n_total = 0.0;
    n_within = 0.0; 
    for (n = 0; n < n_soundings; n++) {
      if (! approx_equal(z_model[n], gridded_data[ps_grid].nodata_value, 1.0e-3)) {
        n_total += 1.0;
        if (100.0*(fabs(fabs(z_model[n]) - fabs(z_soundings[n]))/fabs(z_soundings[n])) < depth_error + epsilon) {
          n_within += 1.0;
        }
      }
    }

    printf("%5.2f%% within %5.2f (rel %% error)\n", 100.0*n_within/n_total, depth_error);
  }

  printf("\n");

  printf("\nSUMMARY OF RELATIVE %% ERRORS +/- 0.5m:\n\n");

  for (i = 0; i < 8; i++) {
    depth_error = ref_errors2[i];
    n_total = 0.0;
    n_within = 0.0; 
    for (n = 0; n < n_soundings; n++) {
      if (! approx_equal(z_model[n], gridded_data[ps_grid].nodata_value, 1.0e-3)) {
        n_total += 1.0;
        for (d = 0.0; d < 1.01; d += 0.01) {
          if (100.0*(fabs(d - 0.5 + fabs(z_model[n]) - fabs(z_soundings[n]))/fabs(z_soundings[n])) < depth_error + epsilon) {
            n_within += 1.0;
            break ;
          }
        }
      }
    }

    printf("%5.2f%% within %5.2f (rel %% error) +/- 0.5m\n", 100.0*n_within/n_total, depth_error);
  }

  printf("\n");

  printf("\nSUMMARY OF RELATIVE %% ERRORS +/- 1.0m:\n\n");

  for (i = 0; i < 8; i++) {
    depth_error = ref_errors2[i];
    n_total = 0.0;
    n_within = 0.0; 
    for (n = 0; n < n_soundings; n++) {
      if (! approx_equal(z_model[n], gridded_data[ps_grid].nodata_value, 1.0e-3)) {
        n_total += 1.0;
        for (d = 0.0; d < 2.01; d += 0.01) {
          if (100.0*(fabs(d - 1.0 + fabs(z_model[n]) - fabs(z_soundings[n]))/fabs(z_soundings[n])) < depth_error + epsilon) {
            n_within += 1.0;
            break ;
          }
        }
      }
    }

    printf("%5.2f%% within %5.2f (rel %% error) +/- 1.0m\n", 100.0*n_within/n_total, depth_error);
  }

  printf("\n");

  // Correlation plot of model depth vs sounding depth. 

#if PGPLOT
  if (nobs > 0.0) {

    char plot_title[] = "Correlation Plot", 
        xlabel[] = "single-beam SONAR depths (m)", // xlabel[] = "various nautical chart depths (m)", xlabel[] = "LIDAR depths (m)"
        ylabel[] = "\\fi photic\\fn model-derived depths (m)";

    valid_model_depths = malloc(n_soundings*sizeof(float));
    valid_sounding_depths = malloc(n_soundings*sizeof(float));

    n_valid_soundings = 0;

    for (n = 0; n < n_soundings; n++) {
      if (! approx_equal(z_model[n], gridded_data[ps_grid].nodata_value, 1.0e-4) &&  
          ! approx_equal(z_soundings[n], spval, 1.0e-4) && 
          fabs(z_model[n]) < fabs(depthmax) && fabs(z_soundings[n]) < fabs(depthmax)) {
        valid_model_depths[n_valid_soundings]    = fabs(z_model[n]);
        valid_sounding_depths[n_valid_soundings] = fabs(z_soundings[n]);
        n_valid_soundings++;
      }
    }

    error = cpgopen("?");
    if (error < 1) {
      printf("\nERROR: cannot load graphics library.\n");
      return ;
    }


    // Set line width. 

    cpgslw(line_width);

    // Set the page size and aspect ratio to 1.0. 

    cpgpap(pagesize, 1.0);

    // Background and text colours. 

    if (background == BLACK) {
      cpgscr(0, 0.0, 0.0, 0.0);
      cpgscr(1, 1.0, 1.0, 1.0);
    } else {
      cpgscr(0, 1.0, 1.0, 1.0);
      cpgscr(1, 0.0, 0.0, 0.0);
    }

    // Plot the line of best fit.

    // Xmin = vec_min(valid_sounding_depths, n_valid_soundings);
    // Xmax = vec_max(valid_sounding_depths, n_valid_soundings);
    // Ymin = vec_min(valid_model_depths, n_valid_soundings);
    // Ymax = vec_max(valid_model_depths, n_valid_soundings);
    Ymin = Xmin; 
    Ymax = Xmax;

    Xmin = 0.0; 
    Xmax = fabs(depthmax);
    // Xmax = MAX(vec_max(valid_sounding_depths, n_valid_soundings), vec_max(valid_model_depths, n_valid_soundings)); 
    Ymin = 0.0; 
    Ymax = fabs(depthmax);
    // Ymax = MAX(vec_max(valid_sounding_depths, n_valid_soundings), vec_max(valid_model_depths, n_valid_soundings)); 

#if 0
    plot_scatter(plot_title, xlabel, ylabel, 2, 16, 
      n_valid_soundings, valid_sounding_depths, valid_model_depths, line_width,
      Xmin, Xmax, Xmin, Xmax, true);
#endif

#if 1
    plot_bidimensional_histogram(plot_title, xlabel, ylabel,  
      n_valid_soundings, valid_sounding_depths, valid_model_depths, line_width,
      Xmin, Xmax, Xmin, Xmax, true, WOR, 200, 200); 
#endif
    // Compute the lines of best fit and correlation coefficient.

    success = linear_fit(valid_sounding_depths, valid_model_depths, n_valid_soundings, &a, &b, &r);

    printf("\nm, b, r, r^2 = %.3f, %.3f, %.3f, %.3f\n", a, b, r, pow(r, 2));

    cpgsls(2);
    cpgsci(8);
    cpgmove(Xmin, a*Xmin + b);
    cpgdraw(Xmax, a*Xmax + b);
    cpgsls(1);

    // Plot the line y == x.

    cpgsls(1);
    cpgsci(15);
    cpgmove(Xmin, Xmin);
    cpgdraw(Xmax, Xmax);
        
    cpgsci(1);

    // Write r^2, slope, mean and relative error to screen. 

    char n_label[64], r_label[64], slope_label[64], abs_err_label[64], rel_err_label[64];
    sprintf(n_label,       "N               = %d", n_valid_soundings);
    sprintf(r_label,       "R\\u2\\d              = %.2f", pow(r, 2));
    sprintf(slope_label,   "slope           = %.2f", a);
    sprintf(abs_err_label, "mean abs error = %.2f (m)",  mean_abs_error);
    sprintf(rel_err_label, "mean rel error  = %.2f (%%)", mean_rel_error);

    sprintf(abs_err_label, "median abs error = %.2f (m)",  median_abs_error);
    sprintf(rel_err_label, "median rel error  = %.2f (%%)", median_rel_error);

    cpgslw(line_width);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.025*(Ymax - Ymin), n_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.075*(Ymax - Ymin), r_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.125*(Ymax - Ymin), slope_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.175*(Ymax - Ymin), abs_err_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.225*(Ymax - Ymin), rel_err_label);

    free(valid_model_depths);
    free(valid_sounding_depths);

    // Close the graphics device. 

    cpgclos();
  }
#endif

  free(x);
  free(y);
  free(z_soundings);
  free(z_model);
}



//
//  SCENE
//

// SCENE NAME /scene name/ BANDS /band 1/ /band 2/ ... WAVELENGTHS /lambda_1/ /lambda_2/ ... 
//    THETAVIEW /subsurface viewing angle from nadir/ THETASUN /subsurface solar zenith angle/ ...
//    RINF /R_inf band 1/ /R_inf band 2/ ... HTIDE /tidal offset to MSL/ K /K band 1/ 
//    /K band 2/ ... RSIGMA /R sigma band 1/ /R sigma band 2/ ... KSIGMA /K sigma band 1/ 
//    /K sigma band 2/ ... [optional] P /P grid/ G /G grid/ X /X grid/

void run_scene() {

  int k, n, m, ps, ps_scene, nbands, ps_p = 0, ps_g = 0, ps_x = 0;
  bool success, p_present = false, g_present = false, x_present = false;

  // Read scene name. 

  if (strcmp(parsed[1], "NAME") != 0) {
    printf("\nERROR: SCENE NAME /scene name/ BANDS /spectral band 1/ /spectral band 2/ ... WAVELENGTHS /lambda_1/ /lambda_2/ ...\
THETAVIEW /subsurface viewing angle from nadir/ THETASUN /subsurface solar zenith angle/ \
RINF /R_inf band 1/ /R_inf band 2/ ... HTIDE /tidal offset to MSL/ K /K band 1/ \
/K band 2/ ... RSIGMA /R sigma band 1/ /R sigma band 2/ ... KSIGMA /K sigma band 1/ /K sigma band 2/ ... \
[optional] P /P grid/ G /G grid/ X /X grid/ \n");
    return ;
  }

  success = parse_new_scene(&ps_scene, 2);
  if (! success) return ;

  strcpy(scene_data[ps_scene].scene_name, parsed[2]);

  // Read spectral bands.

  if (strcmp(parsed[3], "BANDS") != 0) {
    printf("\nERROR: SCENE NAME /scene name/ BANDS /spectral band 1/ /spectral band 2/ ... WAVELENGTHS /lambda_1/ /lambda_2/ ...\
THETAVIEW /subsurface viewing angle from nadir/ THETASUN /subsurface solar zenith angle/ \
RINF /R_inf band 1/ /R_inf band 2/ ... HTIDE /tidal offset to MSL/ K /K band 1/ \
/K band 2/ ... RSIGMA /R sigma band 1/ /R sigma band 2/ ... KSIGMA /K sigma band 1/ /K sigma band 2/ ... \
[optional] P /P grid/ G /G grid/ X /X grid/ \n");
    return ;
  }

  nbands = 0;
  k = 4;

  while (k < MAX_GRIDS && strcmp(parsed[k], "WAVELENGTHS") != 0) {

    success = parse_old_grid(&(scene_data[ps_scene].band_indexes[nbands++]), k++);
    if (! success) return ;
  }

  // Check input grids are commensurate. 

  for (n = 0; n < nbands; n++) {
    for (m = 0; m < nbands; m++) {
      if (! commensurate_grids(
          gridded_data[scene_data[ps_scene].band_indexes[n]], 
          gridded_data[scene_data[ps_scene].band_indexes[m]])) {
        printf("\nERROR: input grids are not commensurate.\n");
        return ; 
      }
    }
  }

  // Set number of spectral bands. 

  scene_data[ps_scene].n_bands = nbands;

  // Read wavelengths. 

  if (strcmp(parsed[k], "WAVELENGTHS") != 0) {
    printf("\nERROR: SCENE NAME /scene name/ BANDS /spectral band 1/ /spectral band 2/ ... WAVELENGTHS /lambda_1/ /lambda_2/ ...\
THETAVIEW /subsurface viewing angle from nadir/ THETASUN /subsurface solar zenith angle/ \
RINF /R_inf band 1/ /R_inf band 2/ ... HTIDE /tidal offset to MSL/ K /K band 1/ \
/K band 2/ ... RSIGMA /R sigma band 1/ /R sigma band 2/ ... KSIGMA /K sigma band 1/ /K sigma band 2/ ... \
[optional] P /P grid/ G /G grid/ X /X grid/\n");
    return ;
  }

  for (n = 0; n < MAX_GRIDS; n++) {
    scene_data[ps_scene].wavelengths[n] = 0.0;
  }  

  for (n = 0; n < nbands; n++) {
    scene_data[ps_scene].wavelengths[n] = atof(parsed[++k]);
  }

  // Read theta_v / theta_view. 

  k++;
  if (! (strcmp(parsed[k], "THETAV") == 0 || strcmp(parsed[k], "THETAVIEW") == 0)) {
    printf("\nERROR: SCENE NAME /scene name/ BANDS /spectral band 1/ /spectral band 2/ ... WAVELENGTHS /lambda_1/ /lambda_2/ ...\
THETAVIEW /subsurface viewing angle from nadir/ THETASUN /subsurface solar zenith angle/ \
RINF /R_inf band 1/ /R_inf band 2/ ... HTIDE /tidal offset to MSL/ K /K band 1/ \
/K band 2/ ... RSIGMA /R sigma band 1/ /R sigma band 2/ ... KSIGMA /K sigma band 1/ /K sigma band 2/ ... \
[optional] P /P grid/ G /G grid/ X /X grid/ \n");
    return ;
  }

  scene_data[ps_scene].theta_v = atof(parsed[++k]);

  // Read theta_w / theta_sun. 

  k++;
  if (! (strcmp(parsed[k], "THETAW") == 0|| strcmp(parsed[k], "THETASUN") == 0)) {
    printf("\nERROR: SCENE NAME /scene name/ BANDS /spectral band 1/ /spectral band 2/ ... WAVELENGTHS /lambda_1/ /lambda_2/ ...\
THETAVIEW /subsurface viewing angle from nadir/ THETASUN /subsurface solar zenith angle/ \
RINF /R_inf band 1/ /R_inf band 2/ ... HTIDE /tidal offset to MSL/ K /K band 1/ \
/K band 2/ ... RSIGMA /R sigma band 1/ /R sigma band 2/ ... KSIGMA /K sigma band 1/ /K sigma band 2/ ... \
[optional] P /P grid/ G /G grid/ X /X grid/ \n");
    return ;
  }

  scene_data[ps_scene].theta_w = atof(parsed[++k]);

  // [optional] Read R_inf. 

  if (strcmp(parsed[++k], "RINF") == 0) {
    printf("\nRINF = ");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].R_inf[n] = atof(parsed[++k]);
      printf("%f ", scene_data[ps_scene].R_inf[n]);
    }
    printf("\n");
    k++;
  } else {
    for (n = 0; n < MAX_GRIDS; n++) {
      scene_data[ps_scene].R_inf[n] = 0.0;
    }
    printf("\n\nWARNING: RINF not specified! Using modelled values.\n");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].R_inf[n] = Lsm[scene_data[ps_scene].band_indexes[n]];
      printf("%f ", scene_data[ps_scene].R_inf[n]);
    }
    printf("\n");
  }

// [optional] Read HTIDE. 

  if (strcmp(parsed[k], "HTIDE") == 0) {
    scene_data[ps_scene].H_tide = atof(parsed[++k]);
    printf("\nHTIDE = %f\n", scene_data[ps_scene].H_tide);
    k++;
  } else {
    scene_data[ps_scene].H_tide = 0.0;
  }

// [optional] Read K. 

  if (strcmp(parsed[k], "K") == 0) {
    printf("\nK = ");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].K[n] = atof(parsed[++k]);
      printf("%.4f ", scene_data[ps_scene].K[n]);
    }
    printf("\n");
    k++;
  } else {
    for (n = 0; n < MAX_GRIDS; n++) {
      scene_data[ps_scene].K[n] = 0.0;
    }
    printf("\n\nWARNING: K not specified! Using modelled values.\n");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].K[n] = Alpha[scene_data[ps_scene].band_indexes[n]];
      printf("%f ", scene_data[ps_scene].K[n]);
    }
    printf("\n");
  }

// [optional] Read RSIGMA. 

  if (strcmp(parsed[k], "RSIGMA") == 0) {
    printf("\nRSIGMA = ");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].R_sigma[n] = atof(parsed[++k]);
      printf("%.6f ", scene_data[ps_scene].R_sigma[n]);
    }
    printf("\n");
    k++;
  } else {
    for (n = 0; n < MAX_GRIDS; n++) {
      scene_data[ps_scene].R_sigma[n] = 0.0;
    }
    printf("\n\nWARNING: RSIGMA not specified! Using modelled values.\n");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].R_sigma[n] = Lsm_sigma[scene_data[ps_scene].band_indexes[n]];
      printf("%f ", scene_data[ps_scene].R_sigma[n]);
    }
    printf("\n");
  }

// [optional] Read KSIGMA. 

  if (strcmp(parsed[k], "KSIGMA") == 0) {
    printf("\nKSIGMA = ");
    for (n = 0; n < nbands; n++) {
      scene_data[ps_scene].K_sigma[n] = atof(parsed[++k]);
      printf("%.6f ", scene_data[ps_scene].K_sigma[n]);
    }
    printf("\n");
    k++;
  } else {
    printf("\n\nWARNING: KSIGMA not specified!\n");
    for (n = 0; n < MAX_GRIDS; n++) {
      scene_data[ps_scene].K_sigma[n] = scene_data[ps_scene].K[n]/10.0; // default is 10% error
    }
  }

// [optional] Read P /P grid/ G /G grid/ X /X grid/

  if (strcmp(parsed[k], "P") == 0) { 
    p_present = true;
    success = parse_old_grid(&ps_p, ++k);
    if (! success) return ;
    k++;
  }

  if (strcmp(parsed[k], "G") == 0) { 
    g_present = true;
    success = parse_old_grid(&ps_g, ++k);
    if (! success) return ;
    k++;
  }

  if (strcmp(parsed[k], "X") == 0) { 
    x_present = true;
    success = parse_old_grid(&ps_x, ++k);
    if (! success) return ;
    k++;
  }

  scene_data[ps_scene].pgx_present = p_present && g_present && x_present;
  printf("\npgx_present = %s\n\n", scene_data[ps_scene].pgx_present ? "true" : "false");
  scene_data[ps_scene].ps_p = ps_p;
  scene_data[ps_scene].ps_g = ps_g;
  scene_data[ps_scene].ps_x = ps_x;

// Store grid dimensions. 

  scene_data[ps_scene].nrows = gridded_data[scene_data[ps_scene].band_indexes[0]].nrows;
  scene_data[ps_scene].ncols = gridded_data[scene_data[ps_scene].band_indexes[0]].ncols;
  printf("\nnrows,ncols = %d,%d\n", scene_data[ps_scene].nrows, scene_data[ps_scene].ncols);

}


//
//    COPY
//

// COPY LSB /source grid/ /destination grid/

// COPY LSM /source grid/ /destination grid/

void run_copy() {

  if (strcmp(parsed[1], "LSB") == 0) {
    run_copy_lsb();
  } else if (strcmp(parsed[1], "LSM") == 0) {
    run_copy_lsm();
  } else {
    printf("\nERROR: unknown command '%s'\n", parsed[0]);
  }
}


void run_copy_lsm() {

  int n, ps_src, ps_dest;
  bool success;

// Read in source grid. 

  success = parse_old_grid(&ps_src, 2);
  if (! success) return ;

// Read in destination grid.

  success = parse_old_grid(&ps_dest, 3);
  if (! success) return ;

// Copy Lsb and Lsb_sigma. 

  if (Lsm[ps_src] == 0.0) {
    printf("\nWARNING! Lsm is zero for grid '%s'\n", parsed[2]);
  }

  Lsm[ps_dest] = Lsm[ps_src];
  Lsm_sigma[ps_dest] = Lsm_sigma[ps_src];
}


void run_copy_lsb() {

  int n, ps_src, ps_dest;
  bool success;

// Read in source grid. 

  success = parse_old_grid(&ps_src, 2);
  if (! success) return ;

// Read in destination grid.

  success = parse_old_grid(&ps_dest, 3);
  if (! success) return ;

// Copy Lsb and Lsb_sigma. 

  if (Lsb[ps_src] == 0.0) {
    printf("\nWARNING! Lsb is zero for grid '%s'\n", parsed[2]);
  }

  Lsb[ps_dest] = Lsb[ps_src];
  Lsb_sigma[ps_dest] = Lsb_sigma[ps_src];
}

//
//    COMPUTE
//

// COMPUTE AREA WET|DRY|TOTAL /grid name/ [optional] COST /cost per km^2/

// COMPUTE LSM /spectral grid name/ [optional] REGION /wlon/ /slat/ /elon/ /nlat/ DEEP /lon/ /lat/

// COMPUTE LSB /spectral grid name/ LAND /land grid name/ SHALLOW /shallow grid name/ ...
//    [optional] SPREAD /spreading coefficient/

// COMPUTE K /spectral grid name 1/ /spectral grid name 2/ WAVELENGTHS /wavelength of band 1/ ...
//    /wavelength of band 2/ COORDINATES /x1/ /y1/ /x2/ /y2/ [optional] LSM ...
//    /Lsm of spectral grid 1/ /Lsm of spectral grid 2/

// COMPUTE BOTTOM IN /spectral grid name/ [optional] LSM /Lsm for the spectral grid/ ... 
//    K /attenuation coefficient for the spectral grid/ DEPTHMAX /maximum depth/

void run_compute() {
  if (strcmp(parsed[1], "AREA") == 0) {
    run_compute_area();
  } else if (strcmp(parsed[1], "LSM") == 0) {
    run_compute_lsm();
  } else if (strcmp(parsed[1], "LSB") == 0) {
    run_compute_lsb();
  } else if (strcmp(parsed[1], "BOTTOM") == 0) {
    run_compute_bottom();
  } else if (strcmp(parsed[1], "K") == 0) {
    run_compute_k();
  } else if (strcmp(parsed[1], "PERCENTILE") == 0) {
    run_compute_percentile();
  }
}




// COMPUTE BOTTOM IN /spectral grid 1/ /spectral grid 2/ OUT /bottom delineation grid/ 
//    K /K for spectral grid 1/ /K for spectral grid 2/ [optional] SCALE /multiplicative scale/ 
//    OFFSET /additive offset/ 


void run_compute_bottom() {

  int n, i, j, k, ps_in_1, ps_in_2, ps_bottom, nrows, ncols;
  float k1, k2, scale, offset, **array;
  bool success;

// Read input spectral grids. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: COMPUTE BOTTOM IN /spectral grid 1/ /spectral grid 2/ OUT /bottom delineation grid/ \
K /K for spectral grid 1/ /K for spectral grid 2/ [optional] SCALE /multiplicative scale/ OFFSET /additive offset/ \n");
    return ;
  }

  success = parse_old_grid(&ps_in_1, 3);
  if (! success) return ;

  success = parse_old_grid(&ps_in_2, 4);
  if (! success) return ;

// Check grids are commensurate. 

  if (gridded_data[ps_in_1].nrows != gridded_data[ps_in_2].nrows || 
    gridded_data[ps_in_1].ncols != gridded_data[ps_in_2].ncols || 
    fabs(gridded_data[ps_in_1].slat - gridded_data[ps_in_2].slat) > extent_eps || 
    fabs(gridded_data[ps_in_1].nlat - gridded_data[ps_in_2].nlat) > extent_eps || 
    fabs(gridded_data[ps_in_1].wlon - gridded_data[ps_in_2].wlon) > extent_eps || 
    fabs(gridded_data[ps_in_1].elon - gridded_data[ps_in_2].elon) > extent_eps || 
    fabs(gridded_data[ps_in_1].cellsize - gridded_data[ps_in_2].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Read output grid. 

  if (strcmp(parsed[5], "OUT") != 0) {
    printf("\nERROR: COMPUTE BOTTOM IN /spectral grid 1/ /spectral grid 2/ OUT /bottom delineation grid/ \
K /K for spectral grid 1/ /K for spectral grid 2/ [optional] SCALE /multiplicative scale/ OFFSET /additive offset/\n");
    return ;
  }

  success = parse_new_grid(&ps_bottom, 6);
  if (! success) return ;

// [optional] Read K. 

  k = 7;

  if (strcmp(parsed[k], "K") == 0) {

    k1 = atof(parsed[++k]);
    k2 = atof(parsed[++k]);
    k++;
  } else {
    if (approx_equal(Alpha[ps_in_1], 0.0, 1.0e-6) || approx_equal(Alpha[ps_in_2], 0.0, 1.0e-6)) {
      printf("\n\nERROR: undefined attenuation coefficients.\n");
      return ;
    } else {
      k1 = Alpha[ps_in_1];
      k2 = Alpha[ps_in_2]; 
    }
  }

  printf("\nk1,k2 = %.5f, %.5f\n", k1,k2);

  // Read normalise.
  
  bool normalise = false;

  if (strcmp(parsed[k], "NORMALISE") == 0) {
    normalise = true;
    printf("\nnormalise = true\n");
    k++;
  }

  // Read scale and offset. 

  if (strcmp(parsed[k], "SCALE") == 0) {
    scale = atof(parsed[++k]);
    k++;
  } else {
    scale = 1.0; 
  }

  if (strcmp(parsed[k], "OFFSET") == 0) {
    offset = atof(parsed[++k]);
    k++;
  } else {
    offset = 0.0; 
  }

  // Check deep water signal has been computed. 

  if (approx_equal(Lsm[ps_in_1], 0.0, 1.0e-6) || approx_equal(Lsm[ps_in_2], 0.0, 1.0e-6)) {
    printf("\n\nERROR: deep water signal not present.\n");
    return ;
  }

  // Allocate output grid. 

  nrows = gridded_data[ps_in_1].nrows;
  ncols = gridded_data[ps_in_1].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  // Compute bottom delineation. 

  gridded_data[ps_bottom].array = array;
  gridded_data[ps_bottom].ncols = ncols; 
  gridded_data[ps_bottom].nrows = nrows; 
  gridded_data[ps_bottom].wlon = gridded_data[ps_in_1].wlon; 
  gridded_data[ps_bottom].slat = gridded_data[ps_in_1].slat;
  gridded_data[ps_bottom].elon = gridded_data[ps_in_1].elon;
  gridded_data[ps_bottom].nlat = gridded_data[ps_in_1].nlat;
  gridded_data[ps_bottom].cellsize = gridded_data[ps_in_1].cellsize;
  gridded_data[ps_bottom].nodata_value = gridded_data[ps_in_1].nodata_value;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(gridded_data[ps_in_1].array[i][j], gridded_data[ps_in_1].nodata_value, 1.0e-6) || 
          approx_equal(gridded_data[ps_in_2].array[i][j], gridded_data[ps_in_2].nodata_value, 1.0e-6)) {

        array[i][j] = gridded_data[ps_in_1].nodata_value;

      } else {

        array[i][j] = k2*log(gridded_data[ps_in_1].array[i][j] - Lsm[ps_in_1] + log_epsilon) - k1*log(gridded_data[ps_in_2].array[i][j] - Lsm[ps_in_2] + log_epsilon);
        array[i][j] /= sqrt(k1*k1 + k2*k2);

        if (array[i][j] != array[i][j]) {
          array[i][j] = gridded_data[ps_in_1].nodata_value;
        }
      }
    }
  }


  // Post process bottom reflectance. 

  if (normalise) {

    int ii, jj, iii, jjj, wx = 32, wy = 32;
    float nobs, background, **array_copy;

    allocate_float_array_2d(&array_copy, nrows, ncols);
    copy_float_array_2d(array, array_copy, nrows, ncols);

    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        if (! approx_equal(array_copy[i][j], gridded_data[ps_bottom].nodata_value, 1.0e-6)) {

          nobs = 0.0;
          background = 0.0;

          for (ii = -wy; ii < wy; ii++) {

            iii = i + ii;
            if (iii < 0) {
              iii = 0;
            } else if (iii > nrows - 1) {
              iii = nrows - 1;
            }

            for (jj = -wx; jj < wx; jj++) {

              jjj = j + jj;
              if (jjj < 0) {
                jjj = 0;
              } else if (jjj > ncols - 1) {
                jjj = ncols - 1;
              }

              if (! approx_equal(array_copy[iii][jjj], gridded_data[ps_bottom].nodata_value, 1.0e-4)) {
                if (fabs(array_copy[iii][jjj]) < 4.0) {
                  background += fabs(array_copy[iii][jjj]);
                  nobs += 1.0;
                }
              }
            }
          }

          if (nobs > 0.5) {
            background /= nobs;
            array[i][j] -= background;
          }
        }
      }
    }    

    free_float_array_2d(array_copy, nrows);

#if 0
    float B_mean = array_mean(array, nrows, ncols, gridded_data[ps_bottom].nodata_value);
    float B_stddev = array_stddev(array, nrows, ncols, gridded_data[ps_bottom].nodata_value, B_mean);

    printf("\nB_mean   = %f", B_mean);
    printf("\nB_stddev = %f\n\n", B_stddev);

    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        if (! approx_equal(array[i][j], gridded_data[ps_bottom].nodata_value, 1.0e-6)) {
          array[i][j] = (array[i][j] - B_mean)/B_stddev;
        }
      }
    }
#endif
  }

}


// COMPUTE PERCENTILE IN /input grid name/ PERCENTILE /n-th percentile/

void run_compute_percentile() {

  int n, i, j, k, ps_in, len, indx, nrows, ncols;
  float *vec, result, spval, percentile;
  bool success;

// Read input grid name. 

  if (strcmp(parsed[2], "IN") != 0 ) {
    printf("\nERROR: COMPUTE PERCENTILE IN /input grid name/ PERCENTILE /n-th percentile/\n");
    return ;
  }

  success = parse_old_grid(&ps_in, 3);
  if (! success) return ;

// Read n-th percentile.

  if (strcmp(parsed[4], "PERCENTILE") != 0 ) {
    printf("\nERROR: COMPUTE PERCENTILE IN /input grid name/ PERCENTILE /n-th percentile/\n");
    return ;
  }

  percentile = atof(parsed[5]);

// Compute the n-th percentile. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  spval = gridded_data[ps_in].nodata_value;

  len = nrows*ncols;
  vec = malloc(len*sizeof(float));

  k = 0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (gridded_data[ps_in].array[i][j] != spval) {
        vec[k++] = gridded_data[ps_in].array[i][j];
      }
    }
  }

  qsort(vec, k, sizeof(float), compare_floats);

  indx = 0.01*percentile*((float) k);
  result = vec[indx];

//  printf("\nk = %d, indx = %d, result = %f\n", k, indx, result);

  free(vec);

  printf("\n%.2f-th percentile = %.3f\n\n", percentile, result);
}


// COMPUTE K /spectral grid name 1/ /spectral grid name 2/ ... WAVELENGTHS /wavelength of band 1/ ...
//    /wavelength of band 2/ ... COORDINATES||REGION /x1/ /y1/ /x2/ /y2/ [optional] LSM ...
//    /Lsm of spectral grid 1/ /Lsm of spectral grid 2/

void run_compute_k() {

  int n, i, j, k, ii, jj, np, nrows, ncols, ps_i, ps_j, npoints, nspec, 
    spectral_indexes[MAX_GRIDS], ni_bins = 512, nj_bins = 512, n_bpl_pts; 
  float wavelengths[MAX_GRIDS], Alpha_computed[MAX_GRIDS], 
    wlen_i, wlen_j, lsm_i, lsm_j, x1, y1, x2, y2, manual_ratio = 0.0, 
    dist, grad, del, x, y, rad, ki, kj, *Li, *Lj, *Ri, *Rj, m, c, r,
    error, xmin, xmax, ymin, ymax, w, water_type, *bpl_i, *bpl_j, *bpl_xi, *bpl_xj, *bpl_ratio, 
    bp_max, bin_min, bin_max, bp_Ri, bp_Rj, bp_xi, bp_xj, xi, xj, ratio, i_delta, 
    raw_xmin, raw_xmax, raw_ymin, raw_ymax, threshold = 0.0, m_all, c_all, r_all;
  bool present, success, line = false, region = false, region_all = false;
  char xlabel[128], ylabel[128], water_type_str[32];

// Read in first spectral grid. 

  nspec = 0;
  k = 2;

  while (true) {

    if (k == MAX_ARGS || nspec >= MAX_GRIDS || strcmp(parsed[k], "WAVELENGTHS") == 0)
      break ;

    present = false;

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        spectral_indexes[nspec++] = n;
        present = true;
        break ;
      }
    }

    if (! present) {
      printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
      return ;
    }

    k++;
  }

  if (nspec < 2) {
    printf("\nERROR: COMPUTE K /spectral grid name 1/ /spectral grid name 2/ ... WAVELENGTHS /wavelength of band 1/ \
/wavelength of band 2/ ... COORDINATES /x1/ /y1/ /x2/ /y2/ [optional] LSM /Lsm of spectral grid 1/ /Lsm of spectral grid 2/\n");
    return ;
  }

  ps_i = spectral_indexes[0];
  ps_j = spectral_indexes[1];

// Check grids are on the same extents. 

  for (n = 1; n < nspec; n++) {
    if (gridded_data[spectral_indexes[0]].nrows != gridded_data[spectral_indexes[n]].nrows || 
      gridded_data[spectral_indexes[0]].ncols != gridded_data[spectral_indexes[n]].ncols || 
      fabs(gridded_data[spectral_indexes[0]].slat - gridded_data[spectral_indexes[n]].slat) > extent_eps || 
      fabs(gridded_data[spectral_indexes[0]].nlat - gridded_data[spectral_indexes[n]].nlat) > extent_eps || 
      fabs(gridded_data[spectral_indexes[0]].wlon - gridded_data[spectral_indexes[n]].wlon) > extent_eps || 
      fabs(gridded_data[spectral_indexes[0]].elon - gridded_data[spectral_indexes[n]].elon) > extent_eps || 
      fabs(gridded_data[spectral_indexes[0]].cellsize - gridded_data[spectral_indexes[n]].cellsize) > extent_eps) {
      printf("\nERROR: input grids are not commensurate.\n");
      return ;
    }
  }

// Read WAVELENGTHS.

  if (strcmp(parsed[k], "WAVELENGTHS") != 0) {
    printf("\nERROR: COMPUTE K /spectral grid name 1/ /spectral grid name 2/ ... WAVELENGTHS /wavelength of band 1/ \
/wavelength of band 2/ ... COORDINATES /x1/ /y1/ /x2/ /y2/ [optional] LSM /Lsm of spectral grid 1/ /Lsm of spectral grid 2/\n");
    return ;
  }

  for (n = 0; n < nspec; n++) {
    wavelengths[n] = atof(parsed[++k]);
  }

  k++;

  wlen_i = wavelengths[0];
  wlen_j = wavelengths[1];

  // printf("\nwlen_i, wlen_j = %.1f, %.1f\n", wlen_i, wlen_j);

// Read COORDINATES. 

  if (strcmp(parsed[k], "COORDINATES") == 0) {
    line = true;
  } else if (strcmp(parsed[k], "REGION") == 0) {
    region = true;
    if (strcmp(parsed[k + 1], "ALL") == 0) {
      region_all = true;
      k += 2;
    }
  } else {
    printf("\nERROR: COMPUTE K /spectral grid name 1/ /spectral grid name 2/ ... WAVELENGTHS /wavelength of band 1/ \
/wavelength of band 2/ ... COORDINATES /x1/ /y1/ /x2/ /y2/ [optional] LSM /Lsm of spectral grid 1/ /Lsm of spectral grid 2/\n");
    return ;
  }

  if (region_all) {
    x1 = gridded_data[ps_i].wlon;
    y1 = gridded_data[ps_i].slat;
    x2 = gridded_data[ps_i].elon;
    y2 = gridded_data[ps_i].nlat;
  } else {
    x1 = atof(parsed[++k]);
    y1 = atof(parsed[++k]);
    x2 = atof(parsed[++k]);
    y2 = atof(parsed[++k]);
    k++;
  }

// [optional] Read LSM. 

  if (strcmp(parsed[k], "LSM") == 0) {
    lsm_i = atof(parsed[++k]);
    lsm_j = atof(parsed[++k]);
    k++;
  } else {
    if (Lsm[ps_i] != 0.0) {
      lsm_i = Lsm[ps_i]; // Read from Lsm table. 
      printf("\nlsm_i = %f\n", lsm_i);
    } else {
      printf("\nERROR: Lsm not specified.\n");
      return ;
    }
    
    if (Lsm[ps_j] != 0.0) {
      lsm_j = Lsm[ps_j]; // Read from Lsm table.
      printf("\nlsm_j = %f\n", lsm_j); 
    } else {
      printf("\nERROR: Lsm not specified.\n");
      return ;
    }
  }

// [optional] Read PALETTE 

  palette = GRAYSCALE; // Default. 

  if (strcmp(parsed[k], "PALETTE") == 0) {
    if (strcmp(parsed[++k], "jet") == 0)
      palette = JET;
    else if (strcmp(parsed[k], "inverse_jet") == 0)
      palette = INVERSEJET;
    else if (strcmp(parsed[k], "grayscale") == 0)
      palette = GRAYSCALE;
    else if (strcmp(parsed[k], "ocean") == 0)
      palette = OCEAN;
    else if (strcmp(parsed[k], "gmt") == 0)
      palette = GMT;
    else if (strcmp(parsed[k], "rygb") == 0)
      palette = RYGB;
    else if (strcmp(parsed[k], "wor") == 0)
      palette = WOR;
    else if (strcmp(parsed[k], "viridis") == 0)
      palette = VIRIDIS;
    else if (strcmp(parsed[k], "plasma") == 0)
      palette = PLASMA;
    else {
      printf("\nERROR: unknown colour palette '%s'\n", parsed[k]);
    }
    k++;
  }

// [optional] NBINS.

  if (strcmp(parsed[k], "NBINS") == 0) {
    ni_bins = atoi(parsed[++k]);
    nj_bins = ni_bins;
    k++;
  } 


// [optional] THRESHOLD.

  if (strcmp(parsed[k], "THRESHOLD") == 0) {
    threshold = atof(parsed[++k]);
    k++;
  } 

  printf("\nthreshold = %f\n", threshold);

// [optional] RATIO.

  if (strcmp(parsed[k], "RATIO") == 0) {
    manual_ratio = atof(parsed[++k]);
    k++;
  } 

  printf("\nratio = %f\n", manual_ratio);

// Compute points on the line. 

  if (line) {
    
    // LINE (COORDINATES)

    dist = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    npoints = floor(dist/gridded_data[ps_i].cellsize);

    Ri = malloc(npoints*sizeof(float));
    Rj = malloc(npoints*sizeof(float));
    Li = malloc(npoints*sizeof(float));
    Lj = malloc(npoints*sizeof(float));

    grad = (y2 - y1)/(x2 - x1);

    for (k = 0; k < npoints; k++) {

      // Compute the point on the line. 

      del = ((float) k)/((float) npoints);
      x = x1 + del*(x2 - x1);
      y = y1 + grad*del*(x2 - x1);

      // Compute grid coordinates of point.

      i = round( ((float) gridded_data[ps_i].nrows)*(y - gridded_data[ps_i].slat)/
          (gridded_data[ps_i].nlat - gridded_data[ps_i].slat) );
      j = round( ((float) gridded_data[ps_i].ncols)*(x - gridded_data[ps_i].wlon)/
          (gridded_data[ps_i].elon - gridded_data[ps_i].wlon) );

      // Radiance at i,j for Li. 

      rad = 0.0; 
      np = 0;

      for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
        for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
          if (gridded_data[ps_i].array[i + ii][j + jj] != gridded_data[ps_i].nodata_value) {
            rad += gridded_data[ps_i].array[i + ii][j + jj];
            np++;
          }
        }
      }

      if (np == 0)
        Ri[k] = 0.0;
      else
        Ri[k] = rad/((float) np);

      // Radiance at i,j for Lj. 

      rad = 0.0; 
      np = 0;

      for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
        for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
          if (gridded_data[ps_j].array[i + ii][j + jj] != gridded_data[ps_j].nodata_value) {
            rad += gridded_data[ps_j].array[i + ii][j + jj];
            np++;
          }
        }
      }

      if (np == 0)
        Rj[k] = 0.0;
      else
        Rj[k] = rad/((float) np);
    }
  } else {

    // REGION

    ncols = floor((x2 - x1)/gridded_data[ps_i].cellsize);
    nrows = floor((y2 - y1)/gridded_data[ps_i].cellsize);
    printf("\nnrows,ncols = %d, %d\n", nrows, ncols);
    npoints = nrows*ncols;

    Ri = malloc(npoints*sizeof(float));
    Rj = malloc(npoints*sizeof(float));
    Li = malloc(npoints*sizeof(float));
    Lj = malloc(npoints*sizeof(float));

    k = 0; 
    for (i = 0; i < nrows; i++) {
      y = y1 + ((float) i)*gridded_data[ps_i].cellsize;
      for (j = 0; j < ncols; j++) {
        x = x1 + ((float) j)*gridded_data[ps_i].cellsize;

        if (! (x > gridded_data[ps_i].wlon && x < gridded_data[ps_i].elon && 
          y > gridded_data[ps_i].slat && y < gridded_data[ps_i].nlat)) {
          continue ;
        }

        rad = interp_bilinear(gridded_data[ps_i].array, 
              gridded_data[ps_i].ncols,
              gridded_data[ps_i].nrows, 
              gridded_data[ps_i].wlon, 
              gridded_data[ps_i].slat, 
              gridded_data[ps_i].cellsize, 
              gridded_data[ps_i].nodata_value, 
              x, 
              y);
        
        if (! approx_equal(rad, gridded_data[ps_i].nodata_value, 1.0e-4)) {
          Ri[k] = rad;
        } else {
          Ri[k] = 0.0;
        }

        rad = interp_bilinear(gridded_data[ps_j].array, 
              gridded_data[ps_j].ncols,
              gridded_data[ps_j].nrows, 
              gridded_data[ps_j].wlon, 
              gridded_data[ps_j].slat, 
              gridded_data[ps_j].cellsize, 
              gridded_data[ps_j].nodata_value, 
              x, 
              y);

        if (! approx_equal(rad, gridded_data[ps_j].nodata_value, 1.0e-4)) {
          Rj[k] = rad;
        } else {
          Rj[k] = 0.0;
        }

        k++;
      }
    }
  }

// Compute approximation to attenuation coefficients via two way linear 
// interpolation from Jerlov's table of attenuation coefficients. 

  raw_xmin = vec_min2(Ri, npoints, gridded_data[ps_i].nodata_value);
  raw_xmax = vec_max2(Ri, npoints, gridded_data[ps_i].nodata_value);
  raw_ymin = vec_min2(Rj, npoints, gridded_data[ps_j].nodata_value);
  raw_ymax = vec_max2(Rj, npoints, gridded_data[ps_j].nodata_value);

  if (region) {

    // Estimate BPL. 

    bpl_i = malloc(ni_bins*sizeof(float));
    bpl_j = malloc(ni_bins*sizeof(float));

    bpl_xi = malloc(ni_bins*sizeof(float));
    bpl_xj = malloc(ni_bins*sizeof(float));

    bpl_ratio = malloc(ni_bins*sizeof(float));

    i_delta = (raw_xmax - raw_xmin)/((float) ni_bins);

    n_bpl_pts = 0;

    for (n = 0; n < ni_bins; n++) {

      if (raw_xmin + ((float) n)*i_delta - lsm_i < 1.0) {
        continue ;
      }

      bp_max = 0.0;
      bin_min = log(raw_xmin + ((float) n)*i_delta - lsm_i);
      bin_max = log(raw_xmin + ((float) n + 1)*i_delta - lsm_i);

      for (k = 0; k < npoints; k++) {

        if (Ri[k] < lsm_i + 2.0 || Rj[k] < lsm_j + 2.0) {
          continue ;
        }

        xi = log(Ri[k] - lsm_i);

        if (xi > threshold && bin_min < xi && xi < bin_max) {
          xj = log(Rj[k] - lsm_j);
          ratio = xj/xi;
          if (ratio > bp_max) {
            bp_max = ratio;
            bp_Ri = Ri[k];
            bp_Rj = Rj[k];
            bp_xi = xi;
            bp_xj = xj;
          }
        }
      }

      if (bp_max > epsilon) {
        bpl_i[n_bpl_pts] = bp_Ri;
        bpl_j[n_bpl_pts] = bp_Rj;
        bpl_ratio[n_bpl_pts] = bp_xj/bp_xi;
        bpl_xi[n_bpl_pts] = bp_xi;
        bpl_xj[n_bpl_pts++] = bp_xj;
        // printf("\n%.4f %.4f %.4f", bp_Ri, bp_Rj, bp_Rj/bp_Ri);
      }
    }

    // Compute K model. 

#if 0
    float mean, stddev, bpl_i_endpoints[32], bpl_j_endpoints[32], min_bpl_slope, max_bpl_slope, min_bpl_ratio, max_bpl_ratio; 

    for (n = 0; n < 32; n++) {
      bpl_i_endpoints[n] = bpl_i[n];
      bpl_j_endpoints[n] = bpl_j[n];
    }

    success = jerlov(wlen_i, wlen_j, lsm_i, lsm_j, bpl_i_endpoints, bpl_j_endpoints, 32, &ki, &kj, &m, &c, &r, &water_type, 0.0);
    min_bpl_slope = m;

    for (n = 0; n < 32; n++) {
      bpl_i_endpoints[n] = bpl_i[n_bpl_pts - n - 1];
      bpl_j_endpoints[n] = bpl_j[n_bpl_pts - n - 1];
    }

    success = jerlov(wlen_i, wlen_j, lsm_i, lsm_j, bpl_i_endpoints, bpl_j_endpoints, 32, &ki, &kj, &m, &c, &r, &water_type, 0.0);
    max_bpl_slope = m;

    mean = vec_mean(bpl_ratio, n_bpl_pts);
    stddev = vec_stddev(bpl_ratio, n_bpl_pts, 0.0, mean);
    min_bpl_ratio = mean - 2.0*stddev;
    max_bpl_ratio = mean + 2.0*stddev;

    K_ratios[ps_i][0] = min_bpl_ratio;
    K_ratios[ps_j][0] = min_bpl_ratio;
    K_ratios[ps_i][1] = max_bpl_ratio;
    K_ratios[ps_j][1] = max_bpl_ratio;
    K_ratios[ps_i][2] = min_bpl_slope;
    K_ratios[ps_j][2] = min_bpl_slope;
    K_ratios[ps_i][3] = max_bpl_slope;
    K_ratios[ps_j][3] = max_bpl_slope;
#endif
    printf("\nn_bpl_pts = %d\n", n_bpl_pts);
    success = jerlov(wlen_i, wlen_j, lsm_i, lsm_j, bpl_i, bpl_j, n_bpl_pts, &ki, &kj, &m, &c, &r, &water_type, manual_ratio);
  } else {
    success = jerlov(wlen_i, wlen_j, lsm_i, lsm_j, Ri, Rj, npoints, &ki, &kj, &m, &c, &r, &water_type, manual_ratio);
  }

  if (success) {

    // Compute (possibly updated) Jerlov water type. 

    jerlov_water_type = water_type;
    jerlov_water_type_str(jerlov_water_type, water_type_str, 32);    
    printf("\nWater type is %s\n", water_type_str);

    for (n = 0; n < MAX_GRIDS; n++) {
      Alpha_computed[n] = 0.0; 
    }

    compute_k_from_jerlov(jerlov_water_type, Alpha_computed, spectral_indexes, wavelengths, nspec);

    n_water_type_estimates++;

    if (n_water_type_estimates == 0) {
      jerlov_water_type = water_type;
      for (n = 0; n < nspec; n++) {
        Alpha[spectral_indexes[n]] = Alpha_computed[spectral_indexes[n]];
      }
    } else {
      w = (float) n_water_type_estimates;
      jerlov_water_type = (jerlov_water_type*(w - 1.0) + water_type)/w;
      for (n = 0; n < nspec; n++) {
        Alpha[spectral_indexes[n]] = (Alpha[spectral_indexes[n]]*(w - 1.0) + Alpha_computed[spectral_indexes[n]])/w;
      }
    }

    printf("\n");
    printf("    BAND NAME        K    \n");
    printf("--------------------------\n");
    for (n = 0; n < MAX_GRIDS; n++) {
      if (allocated_grids[n] && Alpha_computed[n] != 0.0) {
        printf("%12.12s        %.3f\n", grid_names[n], Alpha_computed[n]);
      }
    }
    printf("\n");

    // Write K to file k.bam. 

    FILE *fp;
    fp = fopen("k.bam", "w+");
    fprintf(fp, "// BAM ESTIMATES OF ATTENUATION COEFFICIENTS \n");

    for (n = 0; n < nspec; n++) {
      fprintf(fp, "SET K %s %.3f \n", grid_names[spectral_indexes[n]], Alpha[spectral_indexes[n]]);
    }
    fprintf(fp, "//\n");
    fclose(fp);
  }

#if PGPLOT
// Load PGPLOT. 

  error = cpgopen("?");
  if (error < 1) {
    printf("\nERROR: cannot load graphics library.\n");
    return ;
  }

// Set the page size and aspect ratio to 1.0. 

  cpgpap(pagesize, 1.0);

// Set the line width. 

  cpgslw(line_width);

// Background and text colours. 

  if (background == BLACK) {
    cpgscr(0, 0.0, 0.0, 0.0);
    cpgscr(1, 1.0, 1.0, 1.0);
  } else {
    cpgscr(0, 1.0, 1.0, 1.0);
    cpgscr(1, 0.0, 0.0, 0.0);
  }

//  Compute the logarithm of the radiances.

  int n_shallow = 0;

  for (n = 0; n < npoints; n++) {
    if (Ri[n] > lsm_i + 1.0 && Rj[n] > lsm_j + 1.0) {
      
      Li[n_shallow] = log(Ri[n] - lsm_i);
      Lj[n_shallow] = log(Rj[n] - lsm_j);
      if (region && (Li[n_shallow] < threshold || Lj[n_shallow] < threshold)) {
        Li[n_shallow] = 0.0; 
        Lj[n_shallow] = 0.0;
      }

      Ri[n_shallow] = Ri[n];
      Rj[n_shallow] = Rj[n];
      n_shallow++;
    }
  }

//  Extents of the plotting window. 

  xmin = vec_min(Li, n_shallow);
  xmax = vec_max(Li, n_shallow);
  ymin = vec_min(Lj, n_shallow);
  ymax = vec_max(Lj, n_shallow);

  raw_xmin = vec_min(Ri, n_shallow);
  raw_xmax = vec_max(Ri, n_shallow);
  raw_ymin = vec_min(Rj, n_shallow);
  raw_ymax = vec_max(Rj, n_shallow);

//  Frame labels.

  strcpy(xlabel, "log(R\\d");
  strcat(xlabel, parsed[2]);
  strcat(xlabel, "\\u - R\\d\\(2270)\\d");
  strcat(xlabel, parsed[2]);
  strcat(xlabel, "\\u\\u)");

  strcpy(ylabel, "log(R\\d");
  strcat(ylabel, parsed[3]);
  strcat(ylabel, "\\u - R\\d\\(2270)\\d");
  strcat(ylabel, parsed[3]);
  strcat(ylabel, "\\u\\u)");

//  Call the graphics routine. 

  if (line) {
    plot_scatter("ATTENUATION COEFFICIENT SCATTER PLOT", 
      xlabel, ylabel, 
      4, // colour
      16, // marker style
      n_shallow, Li, Lj, line_width,
      xmin, xmax, ymin, ymax, true);
  } else {
    plot_log_bidimensional_histogram("ATTENUATION COEFFICIENT 2D HISTOGRAM", 
      xlabel, ylabel, n_shallow, Ri, Rj, lsm_i, lsm_j, line_width,
      raw_xmin, raw_xmax, raw_ymin, raw_ymax, true, palette, nj_bins, ni_bins);

    // plot BPL. 
    
    cpgsci(4);

    cpgbbuf();

    for (i = 0; i < n_bpl_pts; i++) {
      cpgpt1(bpl_xi[i], bpl_xj[i], 16);
    }

    cpgebuf();

    cpgsci(1);
  }

//  Plot regression line.

  cpgslw(3*line_width);
  cpgsls(2);
  cpgsci(8);
  cpgmove(xmin, m*xmin + c);
  cpgdraw(xmax, m*xmax + c);
  cpgsls(1);
  cpgslw(line_width);

// Compute and plot the line of best fit for all data. 

  if (region) {

    success = linear_fit2(Li, Lj, n_shallow, &m_all, &c_all, &r_all, 0.0);

    printf("\nbest fit over entire region m, c, r**2 = %.3f, %.3f, %.3f\n\n\n", m_all, c_all, pow(r_all, 2));

    cpgslw(3*line_width);
    cpgsls(2);
    cpgsci(10);
    cpgmove(xmin, m_all*xmin + c_all);
    cpgdraw(xmax, m_all*xmax + c_all);
    cpgsls(1);
    cpgslw(line_width);
  }

//  Write data to plotting window.

  char r_label[64], slope_label[64], jerlov_k_label[64];
  char r_rhs_label[64], slope_rhs_label[64], jerlov_k_rhs_label[64];

  sprintf(r_label,       "R\\u2\\d ");
  sprintf(slope_label,   "slope ");

  sprintf(r_rhs_label,         " = %.3f", pow(r, 2));
  sprintf(slope_rhs_label,     " = %.3f", m);

  cpgsci(1);
  cpgtext(xmin + 0.025*(xmax - xmin), ymax - 0.025*(ymax - ymin), r_label);
  cpgtext(xmin + 0.025*(xmax - xmin), ymax - 0.07*(ymax - ymin), slope_label);
  cpgtext(xmin + 0.14*(xmax - xmin), ymax - 0.025*(ymax - ymin), r_rhs_label);
  cpgtext(xmin + 0.14*(xmax - xmin), ymax - 0.07*(ymax - ymin), slope_rhs_label);

  if (success) {
    float offset = 0.115; 
    float step = 0.045;
    for (n = 0; n < nspec; n++) {
      sprintf(jerlov_k_label,"K\\d%s\\u ", grid_names[spectral_indexes[n]]);
      sprintf(jerlov_k_rhs_label, " = %.3f", Alpha_computed[spectral_indexes[n]]);

      cpgtext(xmin + 0.025*(xmax - xmin), ymax - (offset + ((float) n)*step)*(ymax - ymin), jerlov_k_label);

      cpgtext(xmin + 0.14*(xmax - xmin), ymax - (offset + ((float) n)*step)*(ymax - ymin), jerlov_k_rhs_label);
    }
  }

// Close the graphics device. 

  cpgclos();
#endif

  free(Ri);
  free(Rj);
  free(Li);
  free(Lj); 
  free(bpl_i); 
  free(bpl_j);
  free(bpl_xi);
  free(bpl_xj);
  free(bpl_ratio);
}


// COMPUTE LSM /spectral grid name/ [optional] REGION /wlon/ /slat/ /elon/ /nlat/ DEEP /lon/ /lat/

void run_compute_lsm() {

  int n, i, j, k, ps_in, nrows, ncols;
  float cellsize, region_wlon, region_slat, region_elon, region_nlat, wlon, slat, lon, lat, londeep, latdeep, spval;
  bool success, region_present = false, deep_present = false;
  float lsm, A = 0.0, B = 0.0, mean_est;
  double nobs; 

// Check input spectral grid is present. 

  success = parse_old_grid(&ps_in, 2);
  if (! success) return ;

  k = 3;

// Read REGION. [optional]

  if (strcmp(parsed[k], "REGION") == 0) {
    region_present = true;
    region_wlon = atof(parsed[++k]);
    region_slat = atof(parsed[++k]);
    region_elon = atof(parsed[++k]);
    region_nlat = atof(parsed[++k]);
    k++;
  }

// Read DEEP. [optional]

  if (strcmp(parsed[k], "DEEP") == 0) {
    deep_present = true;
    londeep = atof(parsed[++k]);
    latdeep = atof(parsed[++k]);
    k++;
  }

// Check we have sufficient info to compute Lsm. 

  if (!region_present && ! deep_present && deep_lon == 0.0 && deep_lat == 0.0) {
    printf("\nERROR: must specify REGION, DEEP, or SET DEEP to the location of optically deep water\n");
    return ;
  }

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  wlon  = gridded_data[ps_in].wlon;
  slat  = gridded_data[ps_in].slat;
  cellsize = gridded_data[ps_in].cellsize; 
  spval = gridded_data[ps_in].nodata_value;

  if (region_present) {
    // Compute Lsm within the user-defined region. We define Lsm to be the mean 
    // of the records within the user-defined region. 

    nobs = 0.0;
    lsm = 0.0;

    // Crude estimate of the mean. This is used to help prevent "catastrophic cancellation". See 
    // "Accuracy and Reliability in Scientific Computing" for more information. 

    mean_est = 0.0;

    lon = (region_wlon + region_elon)/2.0;
    lat = (region_slat + region_nlat)/2.0;
    i = round((lat - slat)/cellsize);
    j = round((lon - wlon)/cellsize);
    mean_est += gridded_data[ps_in].array[i][j];

    lon = region_wlon;
    lat = region_slat;
    i = round((lat - slat)/cellsize);
    j = round((lon - wlon)/cellsize);
    mean_est += gridded_data[ps_in].array[i][j];

    lon = region_wlon;
    lat = region_nlat;
    i = round((lat - slat)/cellsize);
    j = round((lon - wlon)/cellsize);
    mean_est += gridded_data[ps_in].array[i][j];

    lon = region_elon;
    lat = region_slat;
    i = round((lat - slat)/cellsize);
    j = round((lon - wlon)/cellsize);
    mean_est += gridded_data[ps_in].array[i][j];

    lon = region_elon;
    lat = region_nlat;
    i = round((lat - slat)/cellsize);
    j = round((lon - wlon)/cellsize);
    mean_est += gridded_data[ps_in].array[i][j];

    mean_est /= 5.0;

    // printf("\nmean_est = %f", mean_est);

    for (i = 0; i < nrows; i++) {
      lat = slat + i*cellsize;
      for (j = 0; j < ncols; j++) {
        lon = wlon + j*cellsize;
        if (lon > region_wlon && lon < region_elon && lat > region_slat && lat < region_nlat && 
            gridded_data[ps_in].array[i][j] != gridded_data[ps_in].nodata_value) {
          A    += ((double) gridded_data[ps_in].array[i][j] - mean_est );
          B    += ((double) pow(gridded_data[ps_in].array[i][j] - mean_est, 2) );
          lsm  += ((double) gridded_data[ps_in].array[i][j] );
          nobs += 1.0;
        }
      }
    }

    if (nobs < 1.0) {
      lsm = 0.0;
      Lsm_sigma[ps_in] = 0.0;
    } else {
      lsm /= nobs;
      Lsm_sigma[ps_in] = ( (float) sqrt((B - pow(A, 2)/nobs)/nobs) );
    }
  } else {
    if (! deep_present) {
      londeep = deep_lon; 
      latdeep = deep_lat;
    }
    i = floor((lon - wlon)/cellsize);
    j = floor((lat - slat)/cellsize);
    if (i < 0 || i >= ncols || j < 0 || j >= nrows) {
      printf("ERROR: location of optically deep water outside grid extents.");
      return ;
    }
    lsm = interp_bicubic(
            gridded_data[ps_in].array,
            ncols, nrows, wlon, slat, cellsize, spval,
            lon, lat);
  }

  // Store result. 

  Lsm[ps_in] = ((float) lsm);
  printf("\n    R_deep       = %13.6f", Lsm[ps_in]);
  printf("\n    R_deep sigma = %13.6f\n", Lsm_sigma[ps_in]);
  printf("\n    R_deep  + 1 sigma = %13.6f", Lsm[ps_in] + Lsm_sigma[ps_in]);
  printf("\n    R_deep  + 2 sigma = %13.6f", Lsm[ps_in] + 2.0*Lsm_sigma[ps_in]);
  printf("\n    R_deep  + 3 sigma = %13.6f\n\n", Lsm[ps_in] + 3.0*Lsm_sigma[ps_in]);
}


// COMPUTE LSB /spectral grid name/ LAND /land grid name/ SHALLOW /shallow grid name/ ...
//    [optional] SPREAD /spreading coefficient/ REGION /wlon/ /slat/ /elon/ /nlat/

void run_compute_lsb() {

  int n, k, ps_in, land_index, shallow_index; 
  float lsb, spread = 0.0, sigma;
  bool success, present = false, land_present = false, shallow_present = false; 

// Check grid name is present. 

  success = parse_old_grid(&ps_in, 2);
  if (! success) return ;

  k = 3;

// Read LAND grid name. 

  if (strcmp(parsed[k], "LAND") != 0) {
    printf("\nERROR: COMPUTE LSB /spectral grid name/ LAND /land grid name/ SHALLOW \
/shallow grid name/ [optional] SPREAD /spreading coefficient/\n");
  } else {

// Check grid name is present. 

    k++;

    success = parse_old_grid(&land_index, k);
    if (! success) {
      return ;
    } else {
      land_present = true;
    }

    k++;
  }

// Read SHALLOW grid name. 

  if (strcmp(parsed[k], "SHALLOW") != 0) {
    printf("\nERROR: COMPUTE LSB /spectral grid name/ LAND /land grid name/ SHALLOW \
/shallow grid name/ [optional] SPREAD /spreading coefficient/\n");
  } else {

    // Check grid name is present. 
    k++;

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        shallow_index = n;
        shallow_present = true;
        break ;
      }
    }

    if (! shallow_present) {
      printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
      return ;
    }
    k++;
  }  

// Read SPREADing coefficient. 

  if (strcmp(parsed[k], "SPREAD") == 0) {
    spread = atof(parsed[++k]);
    k++;
  }

// Compute Lsb. 

  lsb = Lsb_shallow(
    gridded_data[ps_in].array, 
    gridded_data[ps_in].nrows, 
    gridded_data[ps_in].ncols, 
    gridded_data[ps_in].nodata_value, 
    gridded_data[shallow_index].array, 
    gridded_data[shallow_index].nodata_value, 
    gridded_data[land_index].array, 
    gridded_data[land_index].nodata_value, 
    spread, &sigma);

// Store result. 

  Lsb[ps_in] = lsb;
  Lsb_sigma[ps_in] = sigma;

  printf("\n    Lsb       = %10.2f\n\n", lsb);
}


// COMPUTE AREA WET|DRY|TOTAL /grid name/ [optional] COST /cost per km^2/

void run_compute_area() {
  int n, i, j, ps_in, nrows, ncols;
  float cost, wlon, elon, slat, nlat, midlat, midlon, scalex, scaley;
  bool present = false, wet = false, dry = false, total = false, 
    cost_present = false;
  double mapscale, cell_area, area = 0.0, ncells = 0.0;;

  // WET, DRY, or TOTAL

  if (strcmp(parsed[2], "WET") == 0) {
    wet = true;
  } else if (strcmp(parsed[2], "DRY") == 0) {
    dry = true;
  } else if (strcmp(parsed[2], "TOTAL") == 0) {
    total = true;
  } else {
    printf("\nERROR: COMPUTE AREA WET|DRY|TOTAL /grid name/ [optional] COST /cost per km^2/\n");
    return ;
  }

  // Check grid name is present. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      ps_in = n;
      present = true;
      break ;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

  if (strcmp(parsed[4], "COST") == 0) {
    cost_present = true;
    cost = atof(parsed[5]);
  }

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  wlon = gridded_data[ps_in].wlon;
  elon = gridded_data[ps_in].elon;
  slat = gridded_data[ps_in].slat;
  nlat = gridded_data[ps_in].nlat;

  midlat = 0.5*(slat + nlat);
  scalex = geodistance(wlon, midlat, elon, midlat);
  scalex /= 1.e3;

  // printf("\nscalex = %f\n", scalex);

  midlon = 0.5*(wlon + elon);
  scaley = geodistance(midlon, slat, midlon, nlat);
  scaley /= 1.e3;

  // printf("\nscaley = %f\n", scaley);

  mapscale = 0.5*(scalex/((float) ncols) + scaley/((float) nrows));

  // printf("\nmapscale = %f\n", mapscale);

  cell_area = pow(mapscale, 2);
  // printf("\ncell_area = %f\n", cell_area);

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {

      if (gridded_data[ps_in].array[i][j] != gridded_data[ps_in].nodata_value) {
        if (wet && gridded_data[ps_in].array[i][j] < -0.01) {
          ncells++;
        } else if (dry && gridded_data[ps_in].array[i][j] > 0.01) {
          ncells++;
        } else if (total) {
          ncells++;
        }
      }
    }
  }

  printf("\n    AREA = %12.2f (km^2)\n", ncells*cell_area);

  if (cost_present) {
    printf("\n    COST = %12.2f ($)\n", cost*ncells*cell_area);
  }
}


//
//    HISTORY
//

void run_history() {

  printf("\n%%1: %s", trim(prev1));
  printf("\n%%2: %s", trim(prev2));
  printf("\n%%3: %s", trim(prev3));
  printf("\n%%4: %s", trim(prev4));
  printf("\n%%5: %s", trim(prev5));
  printf("\n%%6: %s", trim(prev6));
  printf("\n%%7: %s", trim(prev7));
  printf("\n%%8: %s", trim(prev8));
  printf("\n%%9: %s", trim(prev9));

}

//
//    MODEL
//

// MODEL /model name/ GRIDS /band A/ /band B/ /.../ LAND /land grid name/ ... 
//    SHALLOW /shallow water grid name/ BATHYMETRY /output bathymetry grid name/ 
//    [optional] LINE  /xa/ /ya/ /xb/ /yb/ WEIGHTS /sounding weight/ /line weight/ 
//    BOTTOM /output bottom type grid name/ ERROR /depth error estimate/
//    LSM /deep water radiances/ LSB /radiances at zero depth/ ALPHA /attenuation coefficients/
//    BOUNDARY /boundary grid/

void run_model() {
  if (strcmp(parsed[1], "lyzenga_log_linear") == 0 || 
      strcmp(parsed[1], "polcyn_log_ratio") == 0 || 
      strcmp(parsed[1], "stumpf_log_ratio") == 0 || 
      strcmp(parsed[1], "blake_local") == 0) {
    run_empirical_model();
  } else if (strcmp(parsed[1], "Kd_490_Mueller") == 0) {
    run_model_mueller();
  } else if (strcmp(parsed[1], "Kd_490_Mueller_Secchi") == 0) {
    run_model_mueller_secchi();
  } else if (strcmp(parsed[1], "Chl_Morel") == 0) {
    run_model_morel();
  } else if (strcmp(parsed[1], "Lee_Kd_LS8") == 0) {
    run_model_lee_kd();
  } else if (strcmp(parsed[1], "Lee_Secchi_LS8") == 0) {
    run_model_lee_zsd();
  } else if (strcmp(parsed[1], "Lee_Semi_Analytical") == 0) {
    run_model_lee_semi_analytical();
  } else {
    printf("\nERROR: unknown model name '%s'.\n", parsed[1]);
    return ;    
  }
}



//  MODEL Lee_Semi_Analytical IN /scene 1/ /scene 2/ ... OUT /bathymetry grid/ 
//    /bathymetry one sigma error estimate/ /model error grid/ /bottom albedo grid/ 
//    /bottom sand percentage/ /bottom seagrass percentage/ /bottom coral percentage/
//    /Kd min grid/ /bottom type grid/ /index optical depth/ [optional]
//    DEPTHS /depth estimates/ BOTTOMONLY

void run_model_lee_semi_analytical() {

  int k, n, m, nrows, ncols, nscenes, nbands, ps_in[MAX_GRIDS], ps_bathy, ps_model_error, ps_empirical_depths,
    ps_bottom_albedo, ps_bottom_sand, ps_bottom_seagrass, ps_bottom_coral, ps_kd_min, ps_bottom_type, 
    ps_index_optical_depth, ps_bathy_sigma;
  bool success, empirical_depth_present, bottom_only = false;

  // Can we allocate the 8 new grids? 

  if (n_grids_in_use + 8 >= MAX_GRIDS) {
    printf("\nERROR: unable to store grids.\n");
    return ;
  }    

  // Input scenes. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: MODEL Lee_Semi_Analytical IN /scene 1/ /scene 2/ ... OUT /bathymetry grid/ \
/bathymetry one sigma error estimate/ \
/model error grid/ /bottom albedo grid/ \
/bottom sand percentage/ /bottom seagrass percentage/ /bottom coral percentage/ \
/Kd min grid/ /bottom type grid/ /index optical depth/ [optional] \
DEPTHS /depth estimates/ BOTTOMONLY \n");
    return ;
  }  

  nscenes = 0;
  k = 3;

  while (k < MAX_SCENES && strcmp(parsed[k], "OUT") != 0) {

    success = parse_old_scene(&ps_in[nscenes++], k++);
    if (! success) return ;
  }

  // Check input scenes are spatially commensurate. 

  for (n = 0; n < nscenes; n++) {
    for (m = 0; m < nscenes; m++) {
      if (! commensurate_grids(
              gridded_data[scene_data[ps_in[n]].band_indexes[0]], 
              gridded_data[scene_data[ps_in[m]].band_indexes[0]])) {
        printf("\nERROR: input grids are not commensurate.\n");
        return ; 
      }
    }
  }

  // Output grids. 

  if (strcmp(parsed[k], "OUT") != 0) {
    printf("\nERROR: MODEL Lee_Semi_Analytical IN /scene 1/ /scene 2/ ... OUT /bathymetry grid/ \
/bathymetry one sigma error estimate/ \
/model error grid/ /bottom albedo grid/ \
/bottom sand percentage/ /bottom seagrass percentage/ /bottom coral percentage/ \
/Kd min grid/ /bottom type grid/ /index optical depth/ [optional] \
DEPTHS /depth estimates/  BOTTOMONLY \n");
    return ;
  }

  success = parse_new_grid(&ps_bathy, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_bathy_sigma, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_model_error, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_bottom_albedo, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_bottom_sand, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_bottom_seagrass, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_bottom_coral, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_kd_min, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_bottom_type, ++k);
  if (! success) return ;

  success = parse_new_grid(&ps_index_optical_depth, ++k);
  if (! success) return ;

  k++;

  // Read Empirically-derived depths. 

  if (strcmp(parsed[k], "DEPTHS") == 0) {
    empirical_depth_present = true; 
    success = parse_old_grid(&ps_empirical_depths, ++k);
    k++;
    if (! success) return ;

    if (! commensurate_grids(
              gridded_data[scene_data[ps_in[0]].band_indexes[0]], 
              gridded_data[ps_empirical_depths])) {
        printf("\nERROR: input grids and approximate depth grid are not commensurate.\n");
        return ; 
    }
  } else {
    printf("\n\nWARNING!! Empirical depths not present!! This greatly effects the model optimisation!!\n\n");
    empirical_depth_present = false;
    ps_empirical_depths = 0;
  }

// Allocate memory for output grids. 

  float **bathy, **bathy_sigma, **model_error, **bottom_albedo, **bottom_sand, **bottom_seagrass, 
    **bottom_coral, **kd_min, **bottom_type, **index_optical_depth;

  nrows = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nrows;
  ncols = gridded_data[scene_data[ps_in[0]].band_indexes[0]].ncols;

  allocate_float_array_2d(&bathy, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&bathy_sigma, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&model_error, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&bottom_albedo, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&bottom_sand, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&bottom_seagrass, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&bottom_coral, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&kd_min, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&bottom_type, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  allocate_float_array_2d(&index_optical_depth, nrows, ncols);
  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  gridded_data[ps_bathy].array = bathy;
  gridded_data[ps_bathy].ncols = ncols;
  gridded_data[ps_bathy].nrows = nrows;
  gridded_data[ps_bathy].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bathy].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bathy].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bathy].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bathy].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bathy].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value;  

  gridded_data[ps_bathy_sigma].array = bathy_sigma;
  gridded_data[ps_bathy_sigma].ncols = ncols;
  gridded_data[ps_bathy_sigma].nrows = nrows;
  gridded_data[ps_bathy_sigma].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bathy_sigma].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bathy_sigma].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bathy_sigma].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bathy_sigma].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bathy_sigma].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value;  

  gridded_data[ps_model_error].array = model_error;
  gridded_data[ps_model_error].ncols = ncols;
  gridded_data[ps_model_error].nrows = nrows;
  gridded_data[ps_model_error].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_model_error].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_model_error].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_model_error].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_model_error].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_model_error].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value; 

  gridded_data[ps_bottom_albedo].array = bottom_albedo;
  gridded_data[ps_bottom_albedo].ncols = ncols;
  gridded_data[ps_bottom_albedo].nrows = nrows;
  gridded_data[ps_bottom_albedo].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bottom_albedo].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bottom_albedo].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bottom_albedo].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bottom_albedo].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bottom_albedo].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value; 

  gridded_data[ps_bottom_sand].array = bottom_sand;
  gridded_data[ps_bottom_sand].ncols = ncols;
  gridded_data[ps_bottom_sand].nrows = nrows;
  gridded_data[ps_bottom_sand].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bottom_sand].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bottom_sand].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bottom_sand].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bottom_sand].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bottom_sand].nodata_value = -9999.0; 

  gridded_data[ps_bottom_seagrass].array = bottom_seagrass;
  gridded_data[ps_bottom_seagrass].ncols = ncols;
  gridded_data[ps_bottom_seagrass].nrows = nrows;
  gridded_data[ps_bottom_seagrass].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bottom_seagrass].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bottom_seagrass].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bottom_seagrass].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bottom_seagrass].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bottom_seagrass].nodata_value = -9999.0; 

  gridded_data[ps_bottom_coral].array = bottom_coral;
  gridded_data[ps_bottom_coral].ncols = ncols;
  gridded_data[ps_bottom_coral].nrows = nrows;
  gridded_data[ps_bottom_coral].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bottom_coral].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bottom_coral].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bottom_coral].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bottom_coral].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bottom_coral].nodata_value = -9999.0; 

  gridded_data[ps_kd_min].array = kd_min;
  gridded_data[ps_kd_min].ncols = ncols;
  gridded_data[ps_kd_min].nrows = nrows;
  gridded_data[ps_kd_min].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_kd_min].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_kd_min].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_kd_min].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_kd_min].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_kd_min].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value; 

  gridded_data[ps_bottom_type].array = bottom_type;
  gridded_data[ps_bottom_type].ncols = ncols;
  gridded_data[ps_bottom_type].nrows = nrows;
  gridded_data[ps_bottom_type].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_bottom_type].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_bottom_type].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_bottom_type].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_bottom_type].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_bottom_type].nodata_value = -9999.0; 

  gridded_data[ps_index_optical_depth].array = index_optical_depth;
  gridded_data[ps_index_optical_depth].ncols = ncols;
  gridded_data[ps_index_optical_depth].nrows = nrows;
  gridded_data[ps_index_optical_depth].wlon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].wlon;
  gridded_data[ps_index_optical_depth].elon = gridded_data[scene_data[ps_in[0]].band_indexes[0]].elon;
  gridded_data[ps_index_optical_depth].slat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].slat;
  gridded_data[ps_index_optical_depth].nlat = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nlat;
  gridded_data[ps_index_optical_depth].cellsize = gridded_data[scene_data[ps_in[0]].band_indexes[0]].cellsize;
  gridded_data[ps_index_optical_depth].nodata_value = gridded_data[scene_data[ps_in[0]].band_indexes[0]].nodata_value; 


// Call the SA model. 

  samodel(scene_data, gridded_data, ps_in, nscenes, // inputs
    empirical_depth_present, gridded_data[ps_empirical_depths],
    samodel_n_smoothing_radius, samodel_n_spatial, samodel_n_bottoms,
    bathy, bathy_sigma, model_error, bottom_albedo, bottom_sand, bottom_seagrass, bottom_coral, kd_min, bottom_type, index_optical_depth, // outputs
    pagesize, background, line_width // graphics 
  );

}


// MODEL Lee_Secchi_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ ...
//  /red spectral band/ OUT /Secchi disk depth/ [optional] THETA /solar elevation in degrees/


void run_model_lee_zsd() {

  int n, ps_coastal, ps_blue, ps_green, ps_red, ps_out;
  float theta, **array;
  bool success;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

  // Input grids. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: MODEL Lee_Secchi_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ \
/red spectral band/ OUT /Secchi disk depth/ [optional] THETA /solar elevation in degrees/\n");
    return ;
  }

  // Read coastal spectral band. 

  success = parse_old_grid(&ps_coastal, 3);
  if (! success) return ;

  // Read blue spectral band. 

  success = parse_old_grid(&ps_blue, 4);
  if (! success) return ;

  // Read green spectral band. 

  success = parse_old_grid(&ps_green, 5);
  if (! success) return ;

  // Read red spectral band. 

  success = parse_old_grid(&ps_red, 6);
  if (! success) return ;

  // Check input grids are commensurate. 

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_blue].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_blue].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_blue].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_blue].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_blue].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_blue].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_blue].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_red].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_red].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_red].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_red].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_red].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_red].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_red].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Output grid. 

  if (strcmp(parsed[7], "OUT") != 0) {
    printf("\nERROR: MODEL Lee_Secchi_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ \
/red spectral band/ OUT /Secchi disk depth/ [optional] THETA /solar elevation in degrees/\n");
    return ;
  }

// Read THETA. 

  if (strcmp(parsed[9], "THETA") != 0) {
    if (thetas[ps_coastal] != 0.0) {
      theta = thetas[ps_coastal];
    } else {
      printf("\nERROR: MODEL Lee_Secchi_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ \
/red spectral band/ OUT /Secchi disk depth/ [optional] THETA /solar elevation in degrees/\n");
      return ;
    }
  } else {
    theta = atof(parsed[10]);
  }


// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (! allocated_grids[n]) {
      ps_out = n;
      break;
    }
  }

// Do we already have a grid with this name? 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[8]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[8]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[8]);

  int nrows = gridded_data[ps_blue].nrows;
  int ncols = gridded_data[ps_blue].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute the Secchi disk depth. 

  secchi_disk_depth(
    gridded_data[ps_coastal].array, gridded_data[ps_blue].array, 
    gridded_data[ps_green].array, gridded_data[ps_red].array, 
    array, nrows, ncols, 
    gridded_data[ps_coastal].nodata_value, gridded_data[ps_blue].nodata_value, 
    gridded_data[ps_green].nodata_value, gridded_data[ps_red].nodata_value,
    theta); 

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].wlon = gridded_data[ps_blue].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_blue].elon;
  gridded_data[ps_out].slat = gridded_data[ps_blue].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_blue].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_blue].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_blue].nodata_value;

}


// MODEL Lee_Kd_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ ...
//  /red spectral band/ OUT /Kd grid/ 

void run_model_lee_kd() {


  int n, ps_coastal, ps_blue, ps_green, ps_red, ps_out;
  float theta, **array;
  bool success;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

  // Input grids. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: MODEL Lee_Kd_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ \
/red spectral band/ OUT /Kd grid/\n");
    return ;
  }

  // Read coastal spectral band. 

  success = parse_old_grid(&ps_coastal, 3);
  if (! success) return ;

  // Read blue spectral band. 

  success = parse_old_grid(&ps_blue, 4);
  if (! success) return ;

  // Read green spectral band. 

  success = parse_old_grid(&ps_green, 5);
  if (! success) return ;

  // Read red spectral band. 

  success = parse_old_grid(&ps_red, 6);
  if (! success) return ;

  // Check input grids are commensurate. 

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_blue].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_blue].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_blue].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_blue].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_blue].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_blue].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_blue].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (gridded_data[ps_coastal].nrows != gridded_data[ps_red].nrows || 
      gridded_data[ps_coastal].ncols != gridded_data[ps_red].ncols || 
      fabs(gridded_data[ps_coastal].slat - gridded_data[ps_red].slat) > extent_eps || 
      fabs(gridded_data[ps_coastal].nlat - gridded_data[ps_red].nlat) > extent_eps || 
      fabs(gridded_data[ps_coastal].wlon - gridded_data[ps_red].wlon) > extent_eps || 
      fabs(gridded_data[ps_coastal].elon - gridded_data[ps_red].elon) > extent_eps || 
      fabs(gridded_data[ps_coastal].cellsize - gridded_data[ps_red].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Output grid. 

  if (strcmp(parsed[7], "OUT") != 0) {
    printf("\nERROR: MODEL Lee_Kd_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ \
/red spectral band/ OUT /Kd grid/\n");
    return ;
  }

// Read THETA. 

  if (strcmp(parsed[9], "THETA") != 0) {
    if (thetas[ps_coastal] != 0.0) {
      theta = thetas[ps_coastal];
    } else {
    printf("\nERROR: MODEL Lee_Kd_LS8 IN /coastal spectral band/ /blue spectral band/ /green spectral band/ \
/red spectral band/ OUT /Kd grid/\n");
      return ;
    }
  } else {
    theta = atof(parsed[10]);
  }


// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (! allocated_grids[n]) {
      ps_out = n;
      break;
    }
  }

// Do we already have a grid with this name? 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[8]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[8]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[8]);

  int nrows = gridded_data[ps_blue].nrows;
  int ncols = gridded_data[ps_blue].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute Kd(490). 

  Kd_LS8(
    gridded_data[ps_coastal].array, gridded_data[ps_blue].array, 
    gridded_data[ps_green].array, gridded_data[ps_red].array, 
    array, nrows, ncols, 
    gridded_data[ps_coastal].nodata_value, gridded_data[ps_blue].nodata_value, 
    gridded_data[ps_green].nodata_value, gridded_data[ps_red].nodata_value,
    theta); 

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].wlon = gridded_data[ps_blue].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_blue].elon;
  gridded_data[ps_out].slat = gridded_data[ps_blue].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_blue].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_blue].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_blue].nodata_value;


}


// MODEL Chl_Morel IN /blue spectral band/ /green spectral band/ OUT /Chl grid/ 

void run_model_morel() {

  int i, j, n, ps_blue, ps_green, ps_out, nrows, ncols;
  bool success;
  float **array, blue, green, x, p;

  // Input grids. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: MODEL Chl_Morel IN /blue spectral band/ /green spectral band/ OUT /Chl grid/ \n");
    return ;
  }

  // Read blue spectral band. 

  success = parse_old_grid(&ps_blue, 3);
  if (! success) return ;

  // Read green spectral band. 

  success = parse_old_grid(&ps_green, 4);
  if (! success) return ;

  // Check input grids are commensurate. 

  if (gridded_data[ps_blue].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_blue].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_blue].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_blue].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_blue].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_blue].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_blue].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Output grid. 

  if (strcmp(parsed[5], "OUT") != 0) {
    printf("\nERROR: MODEL Chl_Morel IN /blue spectral band/ /green spectral band/ OUT /Chl grid/ \n");
    return ;
  }


// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (! allocated_grids[n]) {
      ps_out = n;
      break;
    }
  }

// Do we already have a grid with this name? 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[6]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[6]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[6]);
  allocated_grids[ps_out] = true;

  nrows = gridded_data[ps_blue].nrows;
  ncols = gridded_data[ps_blue].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].wlon = gridded_data[ps_blue].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_blue].elon;
  gridded_data[ps_out].slat = gridded_data[ps_blue].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_blue].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_blue].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_blue].nodata_value;

// Compute Kd at 490nm

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (gridded_data[ps_blue].array[i][j] != gridded_data[ps_blue].nodata_value && 
          gridded_data[ps_green].array[i][j] != gridded_data[ps_green].nodata_value) {
        blue = gridded_data[ps_blue].array[i][j];
        green = gridded_data[ps_green].array[i][j];
        // per https://oceancolor.gsfc.nasa.gov/atbd/chlor_a/
        x = log10(blue/green);
        p = 0.1977 - 1.8117*x + 1.9743*pow(x, 2) - 2.5635*pow(x, 3) - 0.7218*pow(x, 4);
        array[i][j] = pow(10.0, p);
      } else {
        array[i][j] = gridded_data[ps_blue].nodata_value;
      }
    }
  }

}


// MODEL Kd_490_Mueller_Secchi IN /blue spectral band/ /green spectral band/ OUT /Secchi disk depth grid/

void run_model_mueller_secchi() {

  int i, j, n, ps_blue, ps_green, ps_out, nrows, ncols;
  bool success;
  float **array, blue, green, x, p;

  // Input grid. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: MODEL Kd_490_Mueller_Secchi IN /blue spectral band/ /green spectral band/ OUT /Kd grid/\n");
    return ;
  }

  // Read blue spectral band. 

  success = parse_old_grid(&ps_blue, 3);
  if (! success) return ;

  // Read green spectral band. 

  success = parse_old_grid(&ps_green, 4);
  if (! success) return ;

  // Check input grids are commensurate. 

  if (gridded_data[ps_blue].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_blue].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_blue].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_blue].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_blue].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_blue].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_blue].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Output grid. 

  if (strcmp(parsed[5], "OUT") != 0) {
    printf("\nERROR: MODEL Kd_490_Mueller_Secchi IN /blue spectral band/ /green spectral band/ OUT /Kd grid/\n");
    return ;
  }


// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (! allocated_grids[n]) {
      ps_out = n;
      break;
    }
  }

// Do we already have a grid with this name? 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[6]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[6]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[6]);
  allocated_grids[ps_out] = true;

  nrows = gridded_data[ps_blue].nrows;
  ncols = gridded_data[ps_blue].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].wlon = gridded_data[ps_blue].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_blue].elon;
  gridded_data[ps_out].slat = gridded_data[ps_blue].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_blue].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_blue].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_blue].nodata_value;

// Compute Secchi disk depth (crude approximation) at 490nm

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (gridded_data[ps_blue].array[i][j] != gridded_data[ps_blue].nodata_value && 
          gridded_data[ps_green].array[i][j] != gridded_data[ps_green].nodata_value) {
        blue = gridded_data[ps_blue].array[i][j];
        green = gridded_data[ps_green].array[i][j];
        x = log10(blue/green);
        p = -0.8813 - 2.0584*x + 2.5878*pow(x, 2) - 3.4885*pow(x, 3) - 1.5061*pow(x, 4);
        array[i][j] = 1.0/(0.0166 + pow(10, p));
      } else {
        array[i][j] = gridded_data[ps_blue].nodata_value;
      }
    }
  }

}


// MODEL Kd_490_Mueller IN /blue spectral band/ /green spectral band/ OUT /Kd grid/

void run_model_mueller() {

  int i, j, n, ps_blue, ps_green, ps_out, nrows, ncols;
  bool success;
  float **array, blue, green, x, p;

  // Input grid. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: MODEL Kd_490_Mueller IN /blue spectral band/ /green spectral band/ OUT /Kd grid/\n");
    return ;
  }

  // Read blue spectral band. 

  success = parse_old_grid(&ps_blue, 3);
  if (! success) return ;

  // Read green spectral band. 

  success = parse_old_grid(&ps_green, 4);
  if (! success) return ;

  // Check input grids are commensurate. 

  if (gridded_data[ps_blue].nrows != gridded_data[ps_green].nrows || 
      gridded_data[ps_blue].ncols != gridded_data[ps_green].ncols || 
      fabs(gridded_data[ps_blue].slat - gridded_data[ps_green].slat) > extent_eps || 
      fabs(gridded_data[ps_blue].nlat - gridded_data[ps_green].nlat) > extent_eps || 
      fabs(gridded_data[ps_blue].wlon - gridded_data[ps_green].wlon) > extent_eps || 
      fabs(gridded_data[ps_blue].elon - gridded_data[ps_green].elon) > extent_eps || 
      fabs(gridded_data[ps_blue].cellsize - gridded_data[ps_green].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Output grid. 

  if (strcmp(parsed[5], "OUT") != 0) {
    printf("\nERROR: MODEL Kd_490_Mueller IN /blue spectral band/ /green spectral band/ OUT /Kd grid/\n");
    return ;
  }


// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (! allocated_grids[n]) {
      ps_out = n;
      break;
    }
  }

// Do we already have a grid with this name? 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[6]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[6]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[6]);
  allocated_grids[ps_out] = true;

  nrows = gridded_data[ps_blue].nrows;
  ncols = gridded_data[ps_blue].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].wlon = gridded_data[ps_blue].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_blue].elon;
  gridded_data[ps_out].slat = gridded_data[ps_blue].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_blue].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_blue].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_blue].nodata_value;

// Compute Kd at 490nm.

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (gridded_data[ps_blue].array[i][j] != gridded_data[ps_blue].nodata_value && 
          gridded_data[ps_green].array[i][j] != gridded_data[ps_green].nodata_value) {
        blue = gridded_data[ps_blue].array[i][j];
        green = gridded_data[ps_green].array[i][j];
        x = log10(blue/green);
        p = -0.8813 - 2.0584*x + 2.5878*pow(x, 2) - 3.4885*pow(x, 3) - 1.5061*pow(x, 4);
        array[i][j] = 0.0166 + pow(10, p);
      } else {
        array[i][j] = gridded_data[ps_blue].nodata_value;
      }
    }
  }

}


void run_empirical_model() {
	
	int i, j, n, m, k, ps_out, ps_bottom, ps_error, ps_boundary, nrows, ncols, n_valid_soundings, 
    shallow_i, shallow_j, boundary_i, boundary_j, ns, nsamples; 
	float **bathymetry, **bottom, **error_estimate, depth, xmin, ymin, xmax, ymax, 
    err, L_sk, L_k, bk, si, lon, lat, mean_rel_error = 0.0, mean_abs_error = 0.0, nobs = 0.0, 
    *valid_model_depths, *valid_sounding_depths, a, b, r, error, Xmin, Xmax, Ymin, Ymax, sigma,
    shallow_lat, shallow_lon;
	bool present, land_present = false, shallow_present = false, 
    return_bottom_type = false, return_error_estimate = false, success;
  FILE *fp;
  char model_file_name[] = "model.bam", error_file_name[] = "output.csv";

	line_xa = 0.0; line_ya = 0.0; line_xb = 0.0; line_yb = 0.0; // Global.
	xmin = 0.0; ymin = 0.0; xmax = 0.0; ymax = 0.0;

// Has the model been previously run?

	if (zdepth != NULL) {
		free(mparams);
		zdepth = NULL;
	}

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
  	printf("\nERROR: unable to store grid.\n");
  	return ;
  }

// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
  	if (! allocated_grids[n]) {
  		ps_out = n;
  		break;
  	}
  }

// Read model name. 

	if (! (strcmp(parsed[1], "lyzenga_log_linear") == 0 || 
			strcmp(parsed[1], "polcyn_log_ratio") == 0 || 
			strcmp(parsed[1], "stumpf_log_ratio") == 0 || 
			strcmp(parsed[1], "blake_local") == 0) ) {
		printf("\nERROR: unknown model name '%s'.\n", parsed[1]);
		return ; 
	} else {
		strcpy(model_name, parsed[1]);
	}

// Read spectral banded grid names.

	if (strcmp(parsed[2], "GRIDS") != 0) {
    printf("\nERROR: MODEL /model name/ GRIDS /band A/ /band B/ /.../ LAND /land grid name/ ... \
SHALLOW /shallow water grid name/ BATHYMETRY /output bathymetry grid name/ [optional] LINE  /xa/ /ya/ /xb/ /yb/ \
WEIGHTS /sounding weight/ /line weight/ BOTTOM /output bottom type grid name/ ERROR /depth error estimate/ LSM \
/deep water radiances/ LSB /radiances at zero depth/ Alpha /attenuation coefficients/\n");
		return ;
	}

	k = 3; 
	n_model_bands = 0; // Global. 

	while (true) {

		if (n_model_bands >= MAX_GRIDS || strcmp(parsed[k], "BATHYMETRY") == 0)
			break ;

		present = false;

		for (n = 0; n < MAX_GRIDS; n++) {
  		if (strcmp(grid_names[n], parsed[k]) == 0) {
  			bands_indexes[n_model_bands++] = n;
  			present = true;
  			break ;
  		}
  	}

  	if (! present) {
  		printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
  		return ;
  	}

  	k++;
	}

	if (n_model_bands == 0) {
		printf("\nERROR: grids incorrectly specified.\n");
		return ; 
	}

// Read bathymetry grid name. 

  if (strcmp(parsed[k], "BATHYMETRY") != 0) {
    printf("\nERROR: MODEL /model name/ GRIDS /band A/ /band B/ /.../ LAND /land grid name/ ... \
SHALLOW /shallow water grid name/ BATHYMETRY /output bathymetry grid name/ [optional] LINE  /xa/ /ya/ /xb/ /yb/ \
WEIGHTS /sounding weight/ /line weight/ BOTTOM /output bottom type grid name/ ERROR /depth error estimate/ LSM \
/deep water radiances/ LSB /radiances at zero depth/ Alpha /attenuation coefficients/\n");
    return ;
  }

// Do we already have a grid with this name? 

  k++;
  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[k]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[k]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[k]);
  allocated_grids[ps_out] = true;
  n_grids_in_use++;
  k++;

// [optional]

// Read LAND cutoff grid name. 

  if (strcmp(parsed[k], "LAND") == 0) {
	  k++;
	  land_present = false;
	  for (n = 0; n < MAX_GRIDS; n++) {
    	if (strcmp(grid_names[n], parsed[k]) == 0) {
    		land_index = n; // Global.
    		land_present = true;
     		break ;
    	}
    }

    if (! land_present) {
    	printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
    	return ;
    }
    k++;
  }

// Read shallow water grid name. 

  if (strcmp(parsed[k], "SHALLOW") == 0) {
    k++;
    shallow_present = false;
    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        shallow_index = n; // Global.
        shallow_present = true;
        break ;
      }
    }

    if (! shallow_present) {
      printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
      return ;
    }
    k++;
  }

// Weights. 

  if (strcmp(parsed[k], "WEIGHTS") == 0) {
  	weight_soundings = atof(parsed[++k]);
  	weight_constraints = atof(parsed[++k]);
    k++;
  }

// Bottom-type grid name. 

  if (strcmp(parsed[k], "BOTTOM") == 0) {

    return_bottom_type = true;

    // Can we allocate another grid?

    if (n_grids_in_use == MAX_GRIDS) {
      printf("\nERROR: unable to store grid.\n");
      return ;
    }

    // Allocate grid into the first empty slot. 

    for (n = 0; n < MAX_GRIDS; n++) {
      if (! allocated_grids[n]) {
        ps_bottom = n;
        break;
      }
    }

    allocated_grids[ps_bottom] = true;

    // Do we already have a grid with this name? 

    k++;
    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        printf("\nERROR: grid '%s' already in use.\n", parsed[k]);
        return ;
      }
    }

    strcpy(grid_names[ps_bottom], parsed[k]);
    k++;
  }

// ERROR (dz) estimate a la Lyzenga 1985. 

  if (strcmp(parsed[k], "ERROR") == 0) {

    return_error_estimate = true;

    // Can we allocate another grid?

    if (n_grids_in_use == MAX_GRIDS) {
      printf("\nERROR: unable to store grid.\n");
      return ;
    }

    // Allocate grid into the first empty slot. 

    for (n = 0; n < MAX_GRIDS; n++) {
      if (! allocated_grids[n]) {
        ps_error = n;
        break;
      }
    }

    allocated_grids[ps_error] = true;

    // Do we already have a grid with this name? 

    k++;
    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        printf("\nERROR: grid '%s' already in use.\n", parsed[k]);
        return ;
      }
    }

    strcpy(grid_names[ps_error], parsed[k]);

    k++;
  }

// Lsm. 

  printf("\n%s\n",parsed[k]);

  if (strcmp(parsed[k], "LSM") == 0) {
    k++;
    for (n = 0; n < n_model_bands; n++) {
      Lsm[bands_indexes[n]] = atof(parsed[k++]);
      printf("\nLsm[%d] = %f", bands_indexes[n], Lsm[bands_indexes[n]]);
    }
  } else {
  
    // Compute Lsm for each grid. 

    for (n = 0; n < n_model_bands; n++) {
      if (Lsm[bands_indexes[n]] != 0.0)
        continue ;

      printf("\nApproximating Lsm for each grid...\n");
      Lsm[bands_indexes[n]] = Lsm_shallow(
          gridded_data[bands_indexes[n]].array,
          gridded_data[bands_indexes[n]].nrows, 
          gridded_data[bands_indexes[n]].ncols, 
          gridded_data[bands_indexes[n]].nodata_value,
          gridded_data[shallow_index].array,
          gridded_data[shallow_index].nodata_value,
          gridded_data[land_index].array,
          gridded_data[land_index].nodata_value);
      printf("\nLsm[%d] = %f", bands_indexes[n], Lsm[bands_indexes[n]]);
    }
    printf("\n");
  }

// Lsb. Radiance at zero depth - the sea floor radiance/reflectance.

  if (strcmp(parsed[k], "LSB") == 0) {
    k++;
    for (n = 0; n < n_model_bands; n++) {
      Lsb[bands_indexes[n]] = atof(parsed[k++]);
      // printf("\nLsb[%d] = %f", bands_indexes[n], Lsb[bands_indexes[n]]);
    }
    printf("\n");
  } else if (strcmp(model_name, "polcyn_log_ratio") == 0) {
  
    // Compute Lsb for each grid. 
    
    for (n = 0; n < n_model_bands; n++) {
      if (Lsb[bands_indexes[n]] != 0.0)
        continue ;      
      
      printf("\nApproximating Lsb for each grid...\n");

      Lsb[bands_indexes[n]] = Lsb_shallow(
          gridded_data[bands_indexes[n]].array,
          gridded_data[bands_indexes[n]].nrows, 
          gridded_data[bands_indexes[n]].ncols, 
          gridded_data[bands_indexes[n]].nodata_value,
          gridded_data[shallow_index].array,
          gridded_data[shallow_index].nodata_value,
          gridded_data[land_index].array,
          gridded_data[land_index].nodata_value, 
          lsb_spread, &sigma);
      // printf("\nLsb[%d] = %f", bands_indexes[n], Lsb[bands_indexes[n]]);
    }
    printf("\n");
  }

// Alpha. Attenuation coefficients. 

  printf("\n%s\n",parsed[k]);

  if (strcmp(parsed[k], "ALPHA") == 0) {
    k++;
    for (n = 0; n < n_model_bands; n++) {
      Alpha[bands_indexes[n]] = atof(parsed[k++]);
    }
  }

// Read REGION. 

  if (strcmp(parsed[k], "REGION") == 0) {
  	xmin = atof(parsed[++k]);
  	ymin = atof(parsed[++k]);
  	xmax = atof(parsed[++k]);
  	ymax = atof(parsed[++k]);
    k++;
  }

// Read BOUNDARY grid. 

  if (strcmp(parsed[k], "BOUNDARY") == 0) {

    k++;
    present = false;

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        ps_boundary = n;
        present = true;
        break ;
      }
    }

    if (! present) {
      printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
      return ;
    }

    k++;

    // From the boundary grid, create NBOUNDARY artificial soundings. 

    if (read_soundings) {
      free_float_array_2d(soundings, nsoundings);
      free(model_depths);
      free(sounding_errors);
      free(sounding_depths);
      meminuse -= (1.0e-6*(float) 5*nsoundings)*((float) sizeof(float));
    }

    nsoundings = nboundary; 
    allocate_float_array_2d(&soundings, nsoundings, 3);
    meminuse += (1.0e-6*(float) 3*nsoundings)*((float) sizeof(float));

    sounding_errors = (float*) malloc(nsoundings*sizeof(float));
    model_depths = (float*) malloc(nsoundings*sizeof(float));
    sounding_depths = (float*) malloc(nsoundings*sizeof(float));

    ns = 0;
    nsamples = 0;

//    printf("\nnsoundings = %d\n", nsoundings);

    while (true) {
      nsamples++;

      // Check if a randomly generated point on the boundary grid and the shallow grid is not nodata. 

      shallow_i = random_in_range(0, gridded_data[shallow_index].nrows);
      shallow_j = random_in_range(0, gridded_data[shallow_index].ncols);

//      printf("\nshallow_i,j = %d, %d", shallow_i, shallow_j);

      shallow_lat = gridded_data[shallow_index].slat + ((float) shallow_i)*gridded_data[shallow_index].cellsize;
      shallow_lon = gridded_data[shallow_index].wlon + ((float) shallow_j)*gridded_data[shallow_index].cellsize;

      // Check shallow grid point is within the extents of the boundary grid. 

      if (shallow_lat < gridded_data[ps_boundary].slat || shallow_lat > gridded_data[ps_boundary].nlat || 
        shallow_lon < gridded_data[ps_boundary].wlon || shallow_lon > gridded_data[ps_boundary].elon) {
        continue ;
      }

      boundary_i = (int) round( ((float) gridded_data[ps_boundary].nrows)*(shallow_lat - gridded_data[ps_boundary].slat)/(gridded_data[ps_boundary].nlat - gridded_data[ps_boundary].slat) );
      boundary_j = (int) round( ((float) gridded_data[ps_boundary].ncols)*(shallow_lon - gridded_data[ps_boundary].wlon)/(gridded_data[ps_boundary].elon - gridded_data[ps_boundary].wlon) );

//      printf("\nboundary_i,j = %d, %d", boundary_i, boundary_j);

      if (boundary_i < 0 || boundary_i >= gridded_data[ps_boundary].nrows || 
          boundary_j < 0 || boundary_j >= gridded_data[ps_boundary].ncols)
        continue ;

      if (gridded_data[shallow_index].array[shallow_i][shallow_j] != gridded_data[shallow_index].nodata_value && 
        gridded_data[ps_boundary].array[boundary_i][boundary_j] != gridded_data[ps_boundary].nodata_value) {
        soundings[ns][0] = shallow_lon;
        soundings[ns][1] = shallow_lat;
        soundings[ns][2] = -1.0*gridded_data[ps_boundary].array[boundary_i][boundary_j];
        ns++;
//        printf("\nboundary depth = %.3f", soundings[ns-1][2]);
      }

      if (ns == nsoundings) 
        break ;
    }

    printf("\nsamples on boundary grid = %d\n\n", nsamples);

  }

// Check soundings are present for optimisation driven model. 

  if (strcmp(model_name, "polcyn_log_ratio") != 0 && nsoundings == 0) {
    printf("\nERROR: depth soundings are not present.\n");
    return ;
  }

// Check uniform grid extents. 

  for (n = 1; n < n_model_bands; n++) {
  	if (gridded_data[bands_indexes[0]].nrows != gridded_data[bands_indexes[n]].nrows || 
  		gridded_data[bands_indexes[0]].ncols != gridded_data[bands_indexes[n]].ncols || 
  		fabs(gridded_data[bands_indexes[0]].slat - gridded_data[bands_indexes[n]].slat) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].nlat - gridded_data[bands_indexes[n]].nlat) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].wlon - gridded_data[bands_indexes[n]].wlon) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].elon - gridded_data[bands_indexes[n]].elon) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].cellsize - gridded_data[bands_indexes[n]].cellsize) > extent_eps) {
  		printf("\nERROR: input grids are not commensurate.\n");
  		return ;
  	}
  }

  if (gridded_data[bands_indexes[0]].nrows != gridded_data[land_index].nrows || 
  		gridded_data[bands_indexes[0]].ncols != gridded_data[land_index].ncols || 
  		fabs(gridded_data[bands_indexes[0]].slat - gridded_data[land_index].slat) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].nlat - gridded_data[land_index].nlat) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].wlon - gridded_data[land_index].wlon) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].elon - gridded_data[land_index].elon) > extent_eps || 
  		fabs(gridded_data[bands_indexes[0]].cellsize - gridded_data[land_index].cellsize) > extent_eps) {
  	printf("\nERROR: input grids are not commensurate.\n");
  	return ;
  }

  if (gridded_data[bands_indexes[0]].nrows != gridded_data[shallow_index].nrows || 
      gridded_data[bands_indexes[0]].ncols != gridded_data[shallow_index].ncols || 
      fabs(gridded_data[bands_indexes[0]].slat - gridded_data[shallow_index].slat) > extent_eps || 
      fabs(gridded_data[bands_indexes[0]].nlat - gridded_data[shallow_index].nlat) > extent_eps || 
      fabs(gridded_data[bands_indexes[0]].wlon - gridded_data[shallow_index].wlon) > extent_eps || 
      fabs(gridded_data[bands_indexes[0]].elon - gridded_data[shallow_index].elon) > extent_eps || 
      fabs(gridded_data[bands_indexes[0]].cellsize - gridded_data[shallow_index].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

  if (nsoundings != 0) {

    // Allocate memory for storing the radiance of each spectral 
    // band at each of the soundings and the sounding weights. 

    allocate_float_array_2d(&radiance_at_soundings, nsoundings, n_model_bands);
    sounding_weights = malloc(nsoundings*sizeof(float));

    // Compute the sounding weights. 

    compute_weights(&sounding_weights, soundings, nsoundings);
  }

// Run the model. 

  if (strcmp(model_name, "lyzenga_log_linear") == 0) {
  	lyzenga_log_linear();
  } else if (strcmp(model_name, "polcyn_log_ratio") == 0) {
		polcyn_log_ratio();
  } else if (strcmp(model_name, "stumpf_log_ratio") == 0) {
  	// stumpf_log_ratio();
  } else if (strcmp(model_name, "blake_local") == 0) {
  	blake_local();
  } else {
  	printf("\nERROR: unknown model '%s'.\n", model_name);
  	return ;
  }

// Create bathymetric grid. 

	nrows = gridded_data[bands_indexes[0]].nrows;
  ncols = gridded_data[bands_indexes[0]].ncols;

  allocate_float_array_2d(&bathymetry, nrows, ncols);

  float bathy_spv = gridded_data[bands_indexes[0]].nodata_value;
  float land_spv = gridded_data[land_index].nodata_value;
  float shallow_spv = gridded_data[shallow_index].nodata_value;

  for (i = 0; i < nrows; i++) {
  	for (j = 0; j < ncols; j++) {
  		if ((land_present == false || gridded_data[land_index].array[i][j] == land_spv) && 
        (shallow_present == false || gridded_data[shallow_index].array[i][j] != shallow_spv)) {

  			depth = zdepth(i, j);
#if 0
  			if (depth < depthmin)
					depth = depthmin;
#endif
        if (depth > depthmax)
					depth = depthmax;

				bathymetry[i][j] = -1.0*depth;
#if 0

// The following modification was per Greg Ashcroft (WKC) wanting depthmax at optically deep water pixels. I 
prefer to keep them as no data values. 

			} else if ((shallow_present == true && gridded_data[shallow_index].array[i][j] == shallow_spv) && 
          (land_present == false || gridded_data[land_index].array[i][j] == land_spv)) {
        // optically deep water
        bathymetry[i][j] = -1.0*depthmax;
#endif
      } else {
        // land 
				bathymetry[i][j] = bathy_spv;
			}
  	}
  }

  gridded_data[ps_out].array = bathymetry;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].cellsize = gridded_data[bands_indexes[0]].cellsize;
  gridded_data[ps_out].slat = gridded_data[bands_indexes[0]].slat;
  gridded_data[ps_out].nlat = gridded_data[bands_indexes[0]].nlat;
  gridded_data[ps_out].wlon = gridded_data[bands_indexes[0]].wlon;
  gridded_data[ps_out].elon = gridded_data[bands_indexes[0]].elon;
  gridded_data[ps_out].nodata_value = bathy_spv;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));  // bathymetry grid

// Create bottom type grid. 

  if (return_bottom_type) {

    if (nsoundings != 0)
      model_error();

    // Allocate memory for the bottom type grid. 

    allocate_float_array_2d(&bottom, nrows, ncols);
    n_grids_in_use++;

    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        if (bathymetry[i][j] != bathy_spv) {
          bottom[i][j] = bottom_type(i, j);
        } else {
          // land or optically deep water
          bottom[i][j] = bathy_spv;
        }
      }
    }

    gridded_data[ps_bottom].array = bottom;
    gridded_data[ps_bottom].ncols = ncols;
    gridded_data[ps_bottom].nrows = nrows;
    gridded_data[ps_bottom].cellsize = gridded_data[bands_indexes[0]].cellsize;
    gridded_data[ps_bottom].slat = gridded_data[bands_indexes[0]].slat;
    gridded_data[ps_bottom].nlat = gridded_data[bands_indexes[0]].nlat;
    gridded_data[ps_bottom].wlon = gridded_data[bands_indexes[0]].wlon;
    gridded_data[ps_bottom].elon = gridded_data[bands_indexes[0]].elon;
    gridded_data[ps_bottom].nodata_value = gridded_data[bands_indexes[0]].nodata_value;

    meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float)); 

  }

// Create error estimate grid. 

  if (return_error_estimate) {

    if (nsoundings != 0)
      model_error();

    // Check that the standard deviation of the deep water (Lsm) signal has been computed. 

    for (k = 0; k < n_model_bands; k++) {
      if ( approx_equal(Lsm_sigma[bands_indexes[k]], 0.0, 1.0e-4) ) 
        printf("\nWARNING: the standard deviation of the deep water (Lsm) signal for band %d is zero!\n", bands_indexes[k]);
    }

    // Allocate memory for the error estimate (dz (m)) array. 

    allocate_float_array_2d(&error_estimate, nrows, ncols);
    allocated_grids[ps_error] = true;
    n_grids_in_use++;

    float max_error = sqrt(fabs(depthmax));

    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        if ((land_present == false || gridded_data[land_index].array[i][j] == land_spv) && 
          (shallow_present == false || gridded_data[shallow_index].array[i][j] != shallow_spv) && 
          bathymetry[i][j] != bathy_spv) {

          err = 0.0;

          for (k = 0; k < n_model_bands; k++) {
            L_k = gridded_data[bands_indexes[k]].array[i][j];
            L_sk = Lsm[bands_indexes[k]];
            if (L_k - L_sk < epsilon) {
              continue ;
            }
            bk = Hj[bands_indexes[k]];
            si = Lsm_sigma[bands_indexes[k]];
            // err += pow(bk*si/(L_k - L_sk + log_epsilon), 2);
            err += fabs( bk*si/(L_k - L_sk + log_epsilon) );
          }

          // err = sqrt( fabs( err/((float) n_model_bands) ) );

          err = sqrt( fabs(err) );

          if (err > max_error) {
            error_estimate[i][j] = max_error;
          } else {
            error_estimate[i][j] = err;
          }

        } else {

          // Land or optically deep water.

          error_estimate[i][j] = bathy_spv;
        }
      }
    }

    gridded_data[ps_error].array = error_estimate;
    gridded_data[ps_error].ncols = ncols;
    gridded_data[ps_error].nrows = nrows;
    gridded_data[ps_error].cellsize = gridded_data[bands_indexes[0]].cellsize;
    gridded_data[ps_error].slat = gridded_data[bands_indexes[0]].slat;
    gridded_data[ps_error].nlat = gridded_data[bands_indexes[0]].nlat;
    gridded_data[ps_error].wlon = gridded_data[bands_indexes[0]].wlon;
    gridded_data[ps_error].elon = gridded_data[bands_indexes[0]].elon;
    gridded_data[ps_error].nodata_value = gridded_data[bands_indexes[0]].nodata_value;

    meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));  

    // Compute dz at the soundings. 

    if (nsoundings != 0) {

      dz_at_soundings = malloc(nsoundings*sizeof(float)); // Global. 

//      printf("\nncols,nrows = %d, %d", ncols, nrows);
//      printf("\nwlon,slat   = %f, %f", gridded_data[ps_error].wlon, gridded_data[ps_error].slat);
//      printf("\ncellsize    = %f", gridded_data[ps_error].cellsize);

      for (n = 0; n < nsoundings; n++) {
        lon = soundings[n][0];
        lat = soundings[n][1];
        if (lon > gridded_data[ps_error].wlon && lon < gridded_data[ps_error].elon && 
            lat > gridded_data[ps_error].slat && lat < gridded_data[ps_error].nlat) {
          dz_at_soundings[n] = interp_bicubic(
              error_estimate, 
              ncols,
              nrows,  
              gridded_data[ps_error].wlon, 
              gridded_data[ps_error].slat, 
              gridded_data[ps_error].cellsize, 
              gridded_data[ps_error].nodata_value, 
              lon, 
              lat);
        } else {
          dz_at_soundings[n] = 0.0;
        }
      }
    }
  }

// Display model.

  if (nsoundings != 0) {
	  model_error();

	  printf("\n");
    if (return_error_estimate) {
      printf("\n         LON         |          LAT        |         WEIGHT      |  SOUNDING DEPTH  |    MODEL DEPTH   |    ABS ERROR (m)    |    REL %% ERROR      |      DZ (m)         |        RADIANCE ");
      printf("\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
      for (n = 0; n < nsoundings; n++) {
          if (approx_equal(model_depths[n], 0.0, 1e-3)) {
            continue ;
          }
        printf("\n %14.6f      | %14.6f      | %14.6f      |    %8.2f      |    %8.2f      |    %8.2f         |    %8.2f         |    %8.2f         |   ", 
          soundings[n][0], soundings[n][1], sounding_weights[n], soundings[n][2], model_depths[n], 
          fabs(model_depths[n] - soundings[n][2]), 100.0*fabs(model_depths[n] - soundings[n][2])/fabs(soundings[n][2]), dz_at_soundings[n]);
        for (m = 0; m < n_model_bands; m++)
          printf(" %8.2f,", radiance_at_soundings[n][m]);

        if (model_depths[n] != bathy_spv) {
          if (! approx_equal(fabs(soundings[n][2]), 0.0, 1e-3)) {
            mean_rel_error += 100.0*fabs(model_depths[n] - soundings[n][2])/fabs(soundings[n][2]);
          }
          mean_abs_error += fabs(model_depths[n] - soundings[n][2]);
          nobs += 1.0;
        }
      }
    } else {
      printf("\n         LON         |          LAT        |         WEIGHT      |  SOUNDING DEPTH  |    MODEL DEPTH   |    ABS ERROR (m)    |    REL %% ERROR      |        RADIANCE ");
      printf("\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
      for (n = 0; n < nsoundings; n++) {
        if (approx_equal(model_depths[n], 0.0, 1e-3)) 
          continue ;
        printf("\n %14.6f      | %14.6f      | %14.6f      |    %8.2f      |    %8.2f      |    %8.2f         |    %8.2f         |   ", 
          soundings[n][0], soundings[n][1], sounding_weights[n], soundings[n][2], model_depths[n], 
          model_depths[n] - soundings[n][2], 100.0*fabs(model_depths[n] - soundings[n][2])/fabs(soundings[n][2]));
        for (m = 0; m < n_model_bands; m++)
          printf(" %8.2f,", radiance_at_soundings[n][m]);

        if (model_depths[n] != bathy_spv) {
          mean_rel_error += 100.0*fabs(model_depths[n] - soundings[n][2])/fabs(soundings[n][2]);
          mean_abs_error += fabs(model_depths[n] - soundings[n][2]);
          nobs += 1.0;
        }
      }
    }

    mean_rel_error /= nobs;
    mean_abs_error /= nobs;

    printf("\n\nMEAN ABSOLUTE ERROR = %.2f  (m)", mean_abs_error);
    printf("\nMEAN RELATIVE ERROR = %.2f (%%)\n\n", mean_rel_error);

    fp = fopen(error_file_name, "w");

    if (return_error_estimate) {
      fprintf(fp, "LON,LAT,WEIGHT,SOUNDING DEPTH,MODEL DEPTH,ABS ERROR (m),REL %% ERROR,DZ (m),RADIANCE");
      for (n = 0; n < nsoundings; n++) {
        if (approx_equal(model_depths[n], 0.0, 1e-3)) {
          continue ;
        }
        fprintf(fp, "\n%14.6f,%14.6f,%14.6f,%8.2f,%8.2f,%8.2f,%8.2f,%8.2f,", 
          soundings[n][0], 
          soundings[n][1],
          sounding_weights[n],
          soundings[n][2],
          model_depths[n], 
          fabs(model_depths[n] - soundings[n][2]), 
          100.0*fabs(model_depths[n] - soundings[n][2])/fabs(soundings[n][2]), 
          dz_at_soundings[n]);
        for (m = 0; m < n_model_bands - 1; m++) {
          fprintf(fp, "%8.2f,", radiance_at_soundings[n][m]);
        }
        fprintf(fp, "%8.2f", radiance_at_soundings[n][n_model_bands - 1]);
      }
    } else {
      fprintf(fp, "LON,LAT,WEIGHT,SOUNDING DEPTH,MODEL DEPTH,ABS ERROR (m),REL %% ERROR,RADIANCE");
      for (n = 0; n < nsoundings; n++) {
        if (approx_equal(model_depths[n], 0.0, 1e-3)) {
          continue ;
        }
        fprintf(fp, "\n%14.6f,%14.6f,%14.6f,%8.2f,%8.2f,%8.2f,%8.2f,", 
          soundings[n][0], 
          soundings[n][1], 
          sounding_weights[n], 
          soundings[n][2], 
          model_depths[n], 
          fabs(model_depths[n] - soundings[n][2]), 
          100.0*fabs(model_depths[n] - soundings[n][2])/fabs(soundings[n][2]));
        for (m = 0; m < n_model_bands - 1; m++) {
          fprintf(fp, "%8.2f,", radiance_at_soundings[n][m]);
          }
        fprintf(fp, "%8.2f", radiance_at_soundings[n][n_model_bands - 1]);
      }
    }  

    fprintf(fp, "\n\nMEAN ABSOLUTE ERROR,%.2f,(m)", mean_abs_error);
    fprintf(fp, "\nMEAN RELATIVE ERROR,%.2f,(%%)", mean_rel_error);

    fclose(fp);
  }

  // Write model to file. 

  fp = fopen(model_file_name, "w");

  if (strcmp(model_name, "lyzenga_log_linear") == 0)
    fprintf(fp, "MODEL=LYZENGA2005\n");
  else if (strcmp(model_name, "polcyn_log_ratio") == 0)
    fprintf(fp, "MODEL=POLCYN1972\n");
  else 
    fprintf(fp, "MODEL=BAM2017\n");

  fprintf(fp, "NBANDS=%d\n", n_model_bands);

  fprintf(fp, "BANDS=");
  for (m = 0; m < n_model_bands - 1; m++)
    fprintf(fp, "%s,", grid_names[bands_indexes[m]]);
  fprintf(fp, "%s\n", grid_names[bands_indexes[n_model_bands - 1]]);

  fprintf(fp, "LSM=");
  for (m = 0; m < n_model_bands - 1; m++)
    fprintf(fp, "%.3f,", Lsm[bands_indexes[m]]);
  fprintf(fp, "%.3f\n", Lsm[bands_indexes[n_model_bands - 1]]);

  fprintf(fp, "LSB=");
  for (m = 0; m < n_model_bands - 1; m++)
    fprintf(fp, "%.3f,", Lsb[bands_indexes[m]]);
  fprintf(fp, "%.3f\n", Lsb[bands_indexes[n_model_bands - 1]]);

  fprintf(fp, "ALPHA=");
  for (m = 0; m < n_model_bands - 1; m++)
    fprintf(fp, "%.6f,", Alpha[bands_indexes[m]]);
  fprintf(fp, "%.6f\n", Alpha[bands_indexes[n_model_bands - 1]]);

  if (strcmp(model_name, "lyzenga_log_linear") == 0) {
    fprintf(fp, "H=%.3f,", mparams[0]);
    for (m = 0; m < n_model_bands - 1; m++)
      fprintf(fp, "%.3f,", Hj[bands_indexes[m]]);
    fprintf(fp, "%.3f\n\n", Hj[bands_indexes[n_model_bands - 1]]);
  } else {
    fprintf(fp, "H=%.3f,%.3f,%.3f\n", H0, Hj[bands_indexes[0]], Hj[bands_indexes[1]]);
  }

  fclose(fp);

  // Correlation plot of model depth vs sounding depth. 

#if PGPLOT
  if (nsoundings != 0) {

    char plot_title[] = "BAM Depth Correlation Plot", xlabel[] = "sounding depths (m)", ylabel[] = "model depths (m)";

    valid_model_depths = malloc(nsoundings*sizeof(float));
    valid_sounding_depths = malloc(nsoundings*sizeof(float));

    n_valid_soundings = 0;

    for (n = 0; n < nsoundings; n++) {
      if (model_depths[n] != bathy_spv && soundings[n][2] < depthmax) {
        valid_model_depths[n_valid_soundings]    = model_depths[n];
        valid_sounding_depths[n_valid_soundings] = soundings[n][2];
        n_valid_soundings++;
      }
    }

    error = cpgopen("?");
    if (error < 1) {
      printf("\nERROR: cannot load graphics library.\n");
      return ;
    }

    // Set the page size and aspect ratio to 1.0. 

    cpgpap(pagesize, 1.0);

    // Background and text colours. 

    if (background == BLACK) {
      cpgscr(0, 0.0, 0.0, 0.0);
      cpgscr(1, 1.0, 1.0, 1.0);
    } else {
      cpgscr(0, 1.0, 1.0, 1.0);
      cpgscr(1, 0.0, 0.0, 0.0);
    }

    // Plot the line of best fit.

    Xmin = vec_min(valid_sounding_depths, n_valid_soundings);
    Xmax = vec_max(valid_sounding_depths, n_valid_soundings);
    // Ymin = vec_min(valid_model_depths, n_valid_soundings);
    // Ymax = vec_max(valid_model_depths, n_valid_soundings);
    Ymin = Xmin; 
    Ymax = Xmax;

    plot_scatter(plot_title, xlabel, ylabel, 2, 16, 
      n_valid_soundings, valid_sounding_depths, valid_model_depths, line_width,
      Xmin, Xmax, Xmin, Xmax, true);

    // Compute the lines of best fit and correlation coefficient.

    success = linear_fit(valid_sounding_depths, valid_model_depths, n_valid_soundings, &a, &b, &r);

    printf("\nm, b, r, r^2 = %.3f, %.3f, %.3f, %.3f\n", a, b, r, pow(r, 2));

    cpgsls(2);
    cpgsci(8);
    cpgmove(Xmin, a*Xmin + b);
    cpgdraw(Xmax, a*Xmax + b);
    cpgsls(1);

    // Plot the line y == x.

    cpgsls(1);
    cpgsci(15);
    cpgmove(Xmin, Xmin);
    cpgdraw(Xmax, Xmax);
        
    cpgsci(1);

    // Write r^2, slope, mean and relative error to screen. 

    char r_label[64], slope_label[64], abs_err_label[64], rel_err_label[64];
    sprintf(r_label,       "R\\u2\\d              = %.2f", pow(r, 2));
    sprintf(slope_label,   "slope           = %.2f", a);
    sprintf(abs_err_label, "mean abs error = %.2f (m)",  mean_abs_error);
    sprintf(rel_err_label, "mean rel error  = %.2f (%%)", mean_rel_error);

    cpgslw(line_width);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.025*(Ymax - Ymin), r_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.075*(Ymax - Ymin), slope_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.125*(Ymax - Ymin), abs_err_label);
    cpgtext(Xmin + 0.025*(Xmax - Xmin), Ymax - 0.175*(Ymax - Ymin), rel_err_label);

    free(valid_model_depths);
    free(valid_sounding_depths);

    // Close the graphics device. 

    cpgclos();
  }
#endif

  // Free memory. 

  if (return_error_estimate && nsoundings != 0)
    free(dz_at_soundings);

  if (nsoundings != 0)
    free_float_array_2d(radiance_at_soundings, nsoundings);
}


/* Compute a weight for each sounding. This prevents bias towards commonly charted depths. (Usually shoal depths.) */

/* The weighting scheme used is a decreasing exponential. The weight, W, associated with a depth, d, from N depths, D is 
given by: 

  w = sum(exp((d - D_i)^2), i=0,N-1)

  W = 1 - w/(ws*w_max)

  where the weight scale, ws, used is usually around 5% (1.05). 
 */

void compute_weights(float **weights, float **soundings, int nsoundings) {

  int n, k; 
  float weight_scale = 1.05, sum, max;

  for (n = 0; n < nsoundings; n++) {
    sum = 0.0; 
    for (k = 0; k < nsoundings; k++) {
      sum += exp(-pow(soundings[n][2] - soundings[k][2], 2));
    }
    (*weights)[n] = sum;
  }

  max = vec_max(*weights, nsoundings);

  for (n = 0; n < nsoundings; n++) {
    (*weights)[n] = 1.0 - ((*weights)[n])/(weight_scale*max); 
  }

}


/* Notes on the analytical Lyzenga/Polcyn log ratio algorithm with spectral 
attenuation coefficients approximated from Jerlov's table. */

void polcyn_log_ratio() {

/* Presently this is just a placeholder as all the parameters required for this 
  model are user-defined. Eventually this will be automatically computed by the model. */


// Set model depth function.

  zdepth = &zdepth_polcyn_log_ratio;
}


float zdepth_polcyn_log_ratio(int i, int j) {
  float Li, Lj, Lsmi, Lsmj, Lsbi, Lsbj, Ki, Kj, z, a, b, c, xi, xj, xbi, xbj;

  // This should be an average over depth_radius cells, as in the other zdepth routines. 

  Li = gridded_data[bands_indexes[0]].array[i][j];
  Lj = gridded_data[bands_indexes[1]].array[i][j];

  Lsmi = Lsm[bands_indexes[0]];
  Lsmj = Lsm[bands_indexes[1]];
  
  Lsbi = Lsb[bands_indexes[0]];
  Lsbj = Lsb[bands_indexes[1]];
  
  Ki = Alpha[bands_indexes[0]];
  Kj = Alpha[bands_indexes[1]];

  if (Li < Lsmi || Lj < Lsmj)
    return 0.;

  xi = log(Li - Lsmi + log_epsilon);
  xj = log(Lj - Lsmj + log_epsilon);

  xbi = log(Lsbi - Lsmi + log_epsilon);
  xbj = log(Lsbj - Lsmj + log_epsilon);

  a = (xbj - xbi)/(2.0*Kj - 2.0*Ki);
  b = 1.0/(2.0*Kj - 2.0*Ki);
  c = -b;

  H0 = a;
  Hj[bands_indexes[0]] = b;
  Hj[bands_indexes[1]] = c;

  z = a + b*xi + c*xj;
//  printf("\na,b,c,z = %8.2f, %8.2f, %8.2f, %8.2f", a, b, c, z);
  return z;
}


/* The bottom type as defined by Lyzenga "Remote sensing of bottom reflectance and water 
attenuation parameters in shallow water using aircraft and LANDSAT data", International Journal 
of Remote Sensing, January 1981]. The equation for the index of bottom reflectance is given by: 

    Y_i = (K_j log(L_i - Lsm_i) - K_i log(L_j - Lsm_j))/(sqrt(K_i^2 + K_j^2))

where K_i, K_j are the attenuation coefficients for bands i and j. This relies on the 
observation that bottom-reflected radiance is approximately a linear function of 
the bottom reflectance and an exponential function of the water depth.  */

float bottom_type(int i, int j) {

  float Li, Lj, Lsmi, Lsmj, Ki, Kj, Yij;

  Li = gridded_data[bands_indexes[0]].array[i][j];
  Lj = gridded_data[bands_indexes[1]].array[i][j];
  
  Lsmi = Lsm[bands_indexes[0]];
  Lsmj = Lsm[bands_indexes[1]];

  if (Li < Lsmi || Lj < Lsmj) {
    return 0.;
  }

  Ki = Alpha[bands_indexes[0]];
  Kj = Alpha[bands_indexes[1]]; 

  Yij = Kj*log(Li - Lsmi + log_epsilon) - Ki*log(Lj - Lsmj + log_epsilon);
  Yij /= sqrt(Ki*Ki + Kj*Kj);

  return Yij;
}


/* Notes on my regionally-adaptive, sounding-directed model. We construct a global 
model, Rm, which is composed on Nm local models (Nm is often set to the number of 
soundings.) such that 

		Rm(i,j) = sum(w_k(i,j)*m_k(i,j), k=1,Nm)/sum(w_k(i,j), k=1,Nm)

where the weight function is given by 

		w_k(i,j) = exp(-d^2/R)

where d is the distance from (i,j) to the origin of the k-th local model. The local 
model, m_k(i,j), can be any of the known global models. 

*/

void blake_local() {

	int n, nparams, randindex;
	float scalex, scaley;

// Set model depth function.

	zdepth = &zdepth_blake;

// Set the number of local models. 

	n_local_models = nsoundings;  // TODO: this should be adaptive

// Set the number of model parameters.

	nparams = n_local_models*(1 + n_model_bands);
	mparams = (float*) malloc(nparams*sizeof(float));
	
// Set the origin of the model soundings. 

	allocate_int_array_2d(&model_origins, n_local_models, 2);

	for (n = 0; n < n_local_models; n++) {
		randindex = random_in_range(0,nsoundings);
		model_origins[n][0] = round( ((float) gridded_data[bands_indexes[0]].ncols)*
				(soundings[randindex][0] - gridded_data[bands_indexes[0]].wlon)/
				(gridded_data[bands_indexes[0]].elon - gridded_data[bands_indexes[0]].wlon) );

		model_origins[n][1] = round( ((float) gridded_data[bands_indexes[0]].nrows)*
				(soundings[randindex][1] - gridded_data[bands_indexes[0]].slat)/
				(gridded_data[bands_indexes[0]].nlat - gridded_data[bands_indexes[0]].slat) );
	}

	scalex = geodistance(
		gridded_data[bands_indexes[0]].wlon, gridded_data[bands_indexes[0]].slat, 
		gridded_data[bands_indexes[0]].elon, gridded_data[bands_indexes[0]].slat);
	scalex += geodistance(
		gridded_data[bands_indexes[0]].wlon, gridded_data[bands_indexes[0]].nlat, 
		gridded_data[bands_indexes[0]].elon, gridded_data[bands_indexes[0]].nlat);
	scalex /= 2.0*gridded_data[bands_indexes[0]].ncols;

	scaley = geodistance(
		gridded_data[bands_indexes[0]].wlon, gridded_data[bands_indexes[0]].slat, 
		gridded_data[bands_indexes[0]].wlon, gridded_data[bands_indexes[0]].nlat);
	scaley += geodistance(
		gridded_data[bands_indexes[0]].elon, gridded_data[bands_indexes[0]].slat, 
		gridded_data[bands_indexes[0]].elon, gridded_data[bands_indexes[0]].nlat);
	scaley /= 2.0*gridded_data[bands_indexes[0]].nrows;

	map_scale = 0.5*(scalex + scaley);

	// printf("\nmap resolution = %f\n", map_scale);

// Run the optimisation strategy. 

	simulated_annealing(nparams);
}


float zdepth_blake(int i, int j) {
	int pin, n, m, ii, jj, iii, jjj, nobs, nrows, ncols;
	float rad, Lm, depth = 0.0, weight, weightsum,
		distx, disty, dist, local_depth;

	nrows = gridded_data[bands_indexes[m]].nrows;
	ncols = gridded_data[bands_indexes[m]].ncols;

	weightsum = 0.0; // Total of all the weights.
	pin = 0; // Model parameter index.

	for (n = 0; n < n_local_models; n++) {
	
		// Compute the weight.

		distx = j - model_origins[n][0];
		disty = i - model_origins[n][1];
		dist = 1.0e-3*map_scale*sqrt(pow(distx, 2) + pow(disty, 2));

		weight = exp(-pow(dist, 2)/pow(falloff_param, 2));

		// printf("\ndist = %f", dist);
		// printf("\nweight = %f", weight);

//		if (weight < model_weight_cutoff)
//			continue ;

		weightsum += weight;

		// Compute the depth for this model component. 

		local_depth = mparams[pin++];

		for (m = 0; m < n_model_bands; m++) {

			nobs = 0;
			Lm = 0.0;

			for (ii = 1 - depth_radius; ii < depth_radius; ii++) {

				if (i + ii < 0)
					iii = 0;
				else if (i + ii > nrows - 1)
					iii = nrows - 1;
				else 
					iii = i + ii;

				for (jj = 1 - depth_radius; jj < depth_radius; jj++) {

					if (j + jj < 0)	
						jjj = 0;
					else if (j + jj > ncols - 1)
						jjj = ncols - 1;
					else 
						jjj = j + jj;

					rad = gridded_data[bands_indexes[m]].array[iii][jjj];

					if (rad != gridded_data[bands_indexes[m]].nodata_value) {
						Lm += rad;
						nobs++;
					}
				}
			}
		
			if (nobs != 0)
				Lm /= (float) nobs;

			if (Lm > Lsm[bands_indexes[m]])
				local_depth -= mparams[pin++]*log(Lm - Lsm[bands_indexes[m]] + log_epsilon); 
		}

		depth += weight*local_depth;
	}

	if (weightsum == 0.0)
		return depth;
	else
		return depth/weightsum;
}

/* 

Notes on my sounding directed implementation of Lyzenga's 1985 method: We construct 
a depth function, z, which is a function of radiance, L, for the m-th 
spectral band, at position i,j:

		z_m_{i,j} = a_{m,0} + a_{m,1} log( L_m_{i,j} - L_sm )

where L_sm is the radiance of deep water for band m. Then for a M-banded model: 

		z_{i,j} = z_{m_0}_{i,j} + z_{m_1}_{i,j} + ... + z_{m_{M-1}}_{i,j}, 

where the coefficients a_{m,0}, a_{m,1}, for each m, are found using a global 
optimisation scheme (simulated annealing or random search with backtracking). The 
error model for this scheme follows a weighted rms of the relative error: 

		error(z) = sqrt( sum( (w_k (z_{i,j} - S_k)/S_k)^2, k=1,N)/N )

where we have N soundings (obtained from nautical charts, lidar, side scan sonar, etc). 
Generally the weights, w_k, are inversly proportional to the depth and can be set 
and modified by the user. 

An estimate, due to Lyzenga (1985 pp. 123), of the error due to combined system and 
environmental noise is given by 

  dz_{i,j} = sqrt(sum( (a_{k,1}*sqrt(L_sk)/(L_k_{i,j} - L_sk))^2, k=1,N ))


*/

void lyzenga_log_linear() {

	int i, k = 1, n, nparams, lsb_sigma_present, linear_present;

	zdepth = &zdepth_lyzenga_log_linear;
	n_local_models = 1;

#if LINEAR
	linear_present = 2;
#else
  linear_present = 0;
#endif

#if LSB_SIGMA
  lsb_sigma_present = 2;
#else
  lsb_sigma_present = 1;
#endif

  nparams = 1 + lsb_sigma_present*n_model_bands + linear_present;

	mparams = (float*) malloc(nparams*sizeof(float));

	simulated_annealing(nparams);

  for (i = 0; i < n_model_bands; i++)
    Hj[bands_indexes[i]] = mparams[i + 1];

  printf("\n Z = %8.4f", mparams[0]);

	for (n = 0; n < n_model_bands; n++) {
		if (mparams[k] > 0.0) {	
      // Signs correctly reversed below. See eqn. 9, pp. 2253 Lyzenga 2005. 
			printf(" - %8.4f log(L_%d - %8.4f)",  mparams[k], n, Lsm[bands_indexes[n]]);
		} else { 
			printf(" + %8.4f log(L_%d - %8.4f)", -mparams[k], n, Lsm[bands_indexes[n]]);
		}
		k++;
	}

#if LINEAR
  printf("\n a*z^b: a,b = %.3f, %.3f\n", mparams[k], mparams[k + 1]);
#endif
}

// zdepth_lyzenga_log_linear computes the bathy model depth at the grid position i, j, 
// for a n-banded model using Lyzenga's log-linear aproximation. 

float zdepth_lyzenga_log_linear(int i, int j) {

	int k, m, ii, jj, iii, jjj, nobs, nrows, ncols;
	float rad, Lm, depth = mparams[0], a, b, X, h0, h1, h2; 

	nrows = gridded_data[bands_indexes[0]].nrows;
	ncols = gridded_data[bands_indexes[0]].ncols;

	k = 1;

	for (m = 0; m < n_model_bands; m++) {

		nobs = 0;
		Lm = 0.0;
		for (ii = 1 - depth_radius; ii < depth_radius; ii++) {

			if (i + ii < 0)
				iii = 0;
			else if (i + ii > nrows - 1)
				iii = nrows - 1;
			else 
				iii = i + ii;

			for (jj = 1 - depth_radius; jj < depth_radius; jj++) {

				if (j + jj < 0)
					jjj = 0;
				else if (j + jj > ncols - 1)
					jjj = ncols - 1;
				else 
					jjj = j + jj;

				rad = gridded_data[bands_indexes[m]].array[iii][jjj];

				if (rad != gridded_data[bands_indexes[m]].nodata_value) {
					Lm += rad;
					nobs++;
				}
			}
		}
		
		if (nobs != 0)
			Lm /= (float) nobs;

    if (n_model_bands == 2) {
      if (m == 1) {
        h0 = 1.0 + fabs(mparams[k++]);
      } else {
        h0 = -1.0 - fabs(mparams[k++]);
      }
    } else {
      h0 = mparams[k++];
    }
    
		if (Lm > Lsm[bands_indexes[m]]) {
      X = log(Lm - Lsm[bands_indexes[m]] + log_epsilon);
			depth -= h0*X;
		}
	}

#if LINEAR
  a = mparams[k];
  b = mparams[++k];
  return a*depth + b;
#else
  return depth;
#endif
}


// RMS relative error of the model at the depth soundings and constraints from Lyzenga, 2005. 

float model_error() {

	int n, m, i, j, nobs = 0, nsoundings_inside_shallow = 0, n_shallow_soundings;
	float c, m_fit, c_fit, r_fit, slope_error, del, xpts, ypts, slope_depth,
		lon, lat, dep, z, sum_of_weights = 0.0, sum, lsb, cost, ws;
	bool zero_radiance, success;

	int ncols = gridded_data[bands_indexes[0]].ncols;
	int nrows = gridded_data[bands_indexes[0]].nrows;

// Check model parameter bounds. 

  bound_error = 0.0;

  for (i = 0; i <= n_model_bands; i++) {
    if (mparams[i] > 1.5*depthmax) {
      bound_error += fabs(mparams[i] - 1.5*depthmax); 
      mparams[i] = 1.5*depthmax;
    } else if (mparams[i] < -1.5*depthmax) {
      bound_error += fabs(1.5*depthmax - mparams[i]);
      mparams[i] = -1.5*depthmax;
    }
  }

  bound_error /= (float) 2*n_model_bands;

  if (bound_error > 2.0*depthmax) {
    return 6.0*bound_error;
  }

// Larger coefficients are preferable. 

  coeff_size_preference = 0.0;

  for (i = 0; i < n_model_bands; i++) {
    coeff_size_preference += (fabs(depthmax) - fabs(mparams[i]))/fabs(depthmax);  
  }

// Store the model coefficients. 

  for (i = 0; i < n_model_bands; i++)
    Hj[bands_indexes[i]] = mparams[i + 1];

// Compute the corrected value of h_0 = mparams[0] per Lyzenga 2005: sum(h_j*log(Lsb[j]), j=1,N)

#if 0
  mparams[0] = 0.0; 

  for (i = 0; i < n_model_bands; i++) {
//    printf("\n %d, %d, %.2f, %.2f", n_model_bands, bands_indexes[i], Lsb[bands_indexes[i]], Hj[bands_indexes[i]]);

#if LSB_SIGMA
    if (mparams[1 + n_model_bands + i] > 2.0)
      mparams[1 + n_model_bands + i] = 2.0;
    else if (mparams[1 + n_model_bands + i] < -2.0)
      mparams[1 + n_model_bands + i] = -2.0;

    c = mparams[1 + n_model_bands + i];
#else
    c = 0.0;
#endif
    if (Lsb[bands_indexes[i]] > 0.0) {
      lsb = Lsb[bands_indexes[i]] + c*Lsb_sigma[bands_indexes[i]];
      if (lsb > log_epsilon) {
        mparams[0] += Hj[bands_indexes[i]]*log(lsb);
      }
    }
  }
#else 
  mparams[0] = fabs(mparams[0]);
#endif

// Sounding errors. 

  sounding_error = 0.0; // Global.

	if (weight_soundings != 0.0) {

		for (n = 0; n < nsoundings; n++) {

			model_depths[n] = 0.0;
			sounding_errors[n] = 0.0;
      sounding_depths[n] = soundings[n][2];

			for (m = 0; m < n_model_bands; m++)
				radiance_at_soundings[n][m] = 0.0;

			// Position of sounding. 

			lon = soundings[n][0];
			lat = soundings[n][1];
			dep = soundings[n][2];

			if (dep < depthmin || dep > depthmax)
				continue ;

			// Grid coordinates of sounding. 

			i = round( ((float) nrows)*(lat - gridded_data[bands_indexes[0]].slat)/
				(gridded_data[bands_indexes[0]].nlat - gridded_data[bands_indexes[0]].slat) );

			j = round( ((float) ncols)*(lon - gridded_data[bands_indexes[0]].wlon)/
				(gridded_data[bands_indexes[0]].elon - gridded_data[bands_indexes[0]].wlon) );

			// Check sounding is within the domain of the grid. 

			if (i < 0 || i >= gridded_data[bands_indexes[0]].nrows || j < 0 || j >= gridded_data[bands_indexes[0]].ncols)
				continue ;

      // Store the radiance at the sounding. 

      zero_radiance = true;

      for (m = 0; m < n_model_bands; m++) {
        radiance_at_soundings[n][m] += gridded_data[bands_indexes[m]].array[i][j];
        if (radiance_at_soundings[n][m] > epsilon) {
          zero_radiance = false;
        }
      }

      if (zero_radiance)
        continue ;

			// Compute model depth. 

			z = zdepth(i, j);

      if (approx_equal(z, mparams[0], 1e-3)) // All radiances are approx Lsm up to log_epsilon. 
        continue ;

			model_depths[n] = z;

      nsoundings_inside_shallow++; // Sounding is in optically shallow water and sufficiently shallow. 

			// Compute error at depth sounding. 

			if (dep == 0.0) {
				sounding_errors[n] = z;
			} else {
        if (relative_error) {
				  sounding_errors[n] = (z - dep)/dep;
        } else {
          sounding_errors[n] = (z - dep);
        }
      }

      // EXPERIMENTAL!! 17/12/2018
      if (z < 0.0) {
        sounding_errors[n] = pow(sounding_errors[n], 3.0);
      }

			// Weighted RMS error.

			sounding_error += sounding_weights[n]*pow(sounding_errors[n], 2);
			nobs++;

      sum_of_weights += sounding_weights[n];
		}

		sounding_error = sqrt( sounding_error/sum_of_weights );

    if (first_model_error_soundings) {
      first_model_error_soundings = false; // Global.
      n_shallow_soundings = (int) 100.0*((float) nsoundings_inside_shallow)/((float) nsoundings);
      printf("\nPercentage of shallow soundings in optically shallow water = %d%%\n", n_shallow_soundings);
    }
	}

// Deviation from constraints from Lyzenga, 2005.

	constraint_error = 0.0; // Global
  local_constraint_error = 0.0; // Global

	if (weight_constraints != 0.0) {

    // Constraint 2, Lyzenga 2005, pp. 2253, eqn 12: sum(h_j*alpha_j, j=1,N) == 1.0

    for (m = 0; m < n_model_bands; m++) {
//      printf("\nAlpha, bi = %f, %f", Alpha[bands_indexes[m]], Hj[bands_indexes[m]]);
      if (Alpha[bands_indexes[m]] > 0.0) {
        constraint_error += Hj[bands_indexes[m]]*Alpha[bands_indexes[m]];
        local_constraint_error += fabs(Hj[bands_indexes[m]]*Alpha[bands_indexes[m]] - 1.0);
      }
    }

    constraint_error = fabs(constraint_error - 1.0); 
    // constraint_error = fabs(constraint_error); 
	}

// Linear regression slope, m x + c. We use the deviation from m = 1 as an indicator of the correct model slope. 

regression_error = 12.0;

#if 1
  if (weight_regression != 0.0 && ! all_approx_zero(model_depths, nsoundings, epsilon)) {
    success = linear_fit2(sounding_depths, model_depths, nsoundings, &m_fit, &c_fit, &r_fit, 0.0);
    if (success) {
      regression_error = 12.0*fabs(m_fit - 1.0) + 0.5*fabs(c_fit);
    }
  }
#endif


//  printf("\nsounding_error         = %.2f", sounding_error);
//  printf("\nconstraint_error       = %.2f", constraint_error);
//  printf("\nlocal_constraint_error = %.2f", local_constraint_error);
//  printf("\ncoeff_size_preference  = %.2f", coeff_size_preference);
//  printf("\nbound_error            = %.2f", bound_error);
//  printf("\nregression_error           = %.2f", regression_error);

  weight_local_constraints = weight_local_constraints/((float) n_model_bands);
  // weight_local_constraints = 0.0;

  cost = weight_soundings*sounding_error;
  cost += weight_constraints*constraint_error;
  cost += weight_local_constraints*local_constraint_error;
  cost += weight_coeff_size_preference*coeff_size_preference;
  cost += weight_bound_error*bound_error;
  cost += weight_regression*regression_error;

  ws = weight_soundings + weight_constraints + weight_local_constraints + weight_coeff_size_preference + weight_bound_error + weight_regression;
  cost /= ws;

//  printf("\ncost                  = %.2f\n", cost);

  return cost;

//  float partial_cost = (weight_soundings*sounding_error + weight_constraints*constraint_error)/(weight_soundings + weight_constraints);
// 	return (1.0 - main_cost_weight)*(coeff_size_preference + bound_error) + main_cost_weight*partial_cost;
}


void random_search() { }

#define SA_INITIAL_TEMP 1000.0
#define SA_COOLING_RATE 0.00001
#define SA_PARAMETER_RANGE 1.0
#define SA_MAX_ITERATIONS 16384

float initial_temperature = SA_INITIAL_TEMP;
float cooling_rate = SA_COOLING_RATE;
float parameter_range = SA_PARAMETER_RANGE;
int sa_max_iterations = SA_MAX_ITERATIONS;


void simulated_annealing(int nparams) {

	int i, n, nsteps = 0, nexplore = 0, niter = 0, nrestarts = 0, n_low = 0, max_low_repeats = 8;
	float *mparams_best, *mparams_solution, temp, best_cost, neighbour_cost, solution_cost, 
    sounding_error_best, constraint_error_best, coeff_size_preference_best, bound_error_best, regression_error_best, cast;

	srand(time(NULL));

	mparams_best = (float*) malloc(nparams*sizeof(float));
	mparams_solution = (float*) malloc(nparams*sizeof(float));

	temp = initial_temperature;

	// Random initial solution. 

	for (i = 0; i < nparams; i++) {
      mparams_solution[i] = frand2(parameter_range);
  }

	copy_vec(mparams_solution, mparams_best, nparams);
	copy_vec(mparams_solution, mparams, nparams);

	best_cost = model_error();
	solution_cost = best_cost;
	sounding_error_best = sounding_error;
  constraint_error_best = constraint_error; 
  coeff_size_preference_best = coeff_size_preference;
  bound_error_best = bound_error;
  regression_error_best = regression_error;

	printf("INITIAL COST = %6.2f (%%)\n", 100*best_cost);

	while (temp > 1.0) {

		copy_vec(mparams_solution, mparams, nparams);

		// 	Pertubate candidate. 
    if (niter++ > sa_max_iterations) {
      nrestarts++;
      niter = 0;
//      for (i = 0; i < nparams; i++)
//        mparams[i] = frand2(parameter_range);
      copy_vec(mparams, mparams_solution, nparams);
//      copy_vec(mparams_best, mparams, nparams);
//      copy_vec(mparams_best, mparams_solution, nparams);
      neighbour_cost = best_cost;
      solution_cost = neighbour_cost;
    } else {
      cast = parameter_range*(1.0 - ((float) niter)/((float) sa_max_iterations));
		  pertubate_vec(mparams, nparams, cast, parameter_range);
      neighbour_cost = model_error();
    }		

		// Display. 
		if (nsteps++%1000 == 0) {
			printf("\nITERATION   = %d\n", nsteps);
			printf("TEMPERATURE = %10.2f\n", temp);
			printf("COST        = %10.2f\n", best_cost);
      printf("cast        = %10.2f\n", cast);
      printf("meanders    = %d\n", nexplore);
      printf("restarts    = %d\n", nrestarts);
      if (relative_error)
			  printf("RMS relative error = %6.2f (%%)\n", 100.0*sounding_error_best);
      else
        printf("RMS absolute error = %6.2f (m)\n", sounding_error_best);
      printf("regression error           = %6.2f\n", regression_error_best);
      printf("parameter constraint error = %6.2f\n", constraint_error_best); 
      printf("coeff size preference      = %6.2f\n", coeff_size_preference_best);
      printf("bound error                = %6.2f\n", bound_error_best);
#if 1
      for (i = 0; i < nparams; i++) 
        printf("%10.3f ", mparams_best[i]);
      printf("\n");
      for (i = 0; i < nparams; i++) 
        printf("%10.3f ", mparams_solution[i]);
      printf("\n");
#endif
		}

		if (neighbour_cost < solution_cost) {
			solution_cost = neighbour_cost;
			copy_vec(mparams, mparams_solution, nparams);
			if (neighbour_cost < best_cost) {
				best_cost = neighbour_cost;
				copy_vec(mparams, mparams_best, nparams);
				sounding_error_best = sounding_error;
        constraint_error_best = constraint_error;
        coeff_size_preference_best = coeff_size_preference;
        bound_error_best = bound_error;
        regression_error_best = regression_error;
			}
		} else if (exp((solution_cost - neighbour_cost)/temp) > frand()) {
      nexplore++;
			solution_cost = neighbour_cost; 
			copy_vec(mparams, mparams_solution, nparams);
		}

		// Decrease temperature. 
		temp *= 1.0 - cooling_rate;

    if (temp < 1.0) {
      if (++n_low < max_low_repeats) {
        temp = 3.0;
      }
    }
	}

	rms_error = best_cost;
	copy_vec(mparams_best, mparams, nparams);

	free(mparams_best);
	free(mparams_solution);

}

void pertubate_vec(float *vec, int len, float scale, float scale0) {
  int n, discretized_depth;
#if 0
  n = random_in_range(0, len);
	float rndnrm, x1, x2;
	x1 = frand();
	x2 = frand();
	rndnrm = sqrt(-2.0*log(x1))*cos(2.0*PI*x2);
	vec[n] += rndnrm*scale;
#else
  n = random_in_range(0, len);
  vec[n] += scale*(2.0*frand() - 1.0);
  
  // if (n == 1) {
  // discretized_depth = random_in_range(-40*((int) depthmax), 0); 
  // } else {
  // discretized_depth = random_in_range(0, 40*((int) depthmax));
  // }
  
  // discretized_depth = random_in_range(-20*((int) depthmax), 20*((int) depthmax));
  // vec[n] = 0.05*((float) discretized_depth);
#endif
}

void pertubate_vec2(float *vec, int len, float scale) {

	int n = random_in_range(0, len);
	float rndnrm, x1, x2;
	x1 = frand();
	x2 = frand();
	rndnrm = sqrt(-2.0*log(x1))*cos(2.0*PI*x2);
	vec[n] += rndnrm*scale;
	if (vec[n] <= 0.0) {
		vec[n] = fabs(vec[n]) + 0.0001;
	}
}


//
//    DELETE
//

// DELETE /grid name/

void run_delete() {

	int n, np, nrows, ncols; 
	bool present = false;

// Check grid name is present. 

	for (n = 0; n < MAX_GRIDS; n++) {
  	if (strcmp(grid_names[n], parsed[1]) == 0) {
  		np = n;
  		present = true;
  		break ;
  	}
  }

	if (! present) {
	  printf("\nERROR: unknown grid '%s'.\n", parsed[1]);
	  return ;
	}

	nrows = gridded_data[np].nrows;
	ncols = gridded_data[np].ncols;

// Clear any spectral grid stats. 


  thetas[np] = 0.0;
  Lsm[np] = 0.0;
  Lsm_sigma[np] = 0.0;
  Lsb[np] = 0.0;
  Lsb_sigma[np] = 0.0;
  Alpha[np] = 0.0;
  Hj[np]  = 0.0;
  NAlpha[np] = 0;
  model_wavelengths[np] = 0.0;
  radiance_mult_band[np]    = BAM_NODATA;
  radiance_add_band[np]     = BAM_NODATA;
  quantize_cal_min_band[np] = BAM_NODATA;
  quantize_cal_max_band[np] = BAM_NODATA;
  radiance_minimum_band[np] = BAM_NODATA;
  radiance_maximum_band[np] = BAM_NODATA;

// Delete grid name. 

	grid_names[np][0] = '\0';

// Mark as unallocated. 

	allocated_grids[np] = false;

// Update memory usage (in units of megabytes). 

  meminuse -= 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));	

// Free memory. 

	free_float_array_2d(gridded_data[np].array, nrows);

	n_grids_in_use--;
}

//
//    HELP
//
void run_help() {

}


//
//    READ
//

// READ COMMANDS /command file/
// READ GRID /grid name/ /file/ SCALE /scale/ OFFSET /offset/
// READ SOUNDINGS /file/ CHARTDATUM /chart datum/
// READ LANDSAT METADATA /file/


void run_read() {

	if (strcmp(parsed[1], "COMMANDS") == 0) {
		run_read_commands();
	} else if (strcmp(parsed[1], "GRID") == 0) {
		run_read_grid();
	} else if (strcmp(parsed[1], "SOUNDINGS") == 0) {
		run_read_soundings();
  } else if (strcmp(parsed[1], "LANDSAT") == 0) {
    run_read_landsat();
  } else if (strcmp(parsed[1], "BOTTOM") == 0) {
    run_read_bottom();
	} else {
		printf("\nERROR: unknown command '%s'\n", parsed[1]);
	}
}


// READ BOTTOM /bottom spectra file/

void run_read_bottom() {

  int nrows, ncols; 

// Check bottom spectra has not been read in already. 

  if (n_bottom_spectra > 0) {
    printf("\n\nERROR: bottom spectra already specified.\n\n");
    return ;
  }

// Read file name. 

  if (! file_exists(parsed[2])) {
    printf("\nERROR: missing file '%s'\n", parsed[2]);
    return ;
  }

// Check csv file. 

  if (strstr(parsed[2], ".csv") == NULL) { // Crude check.
    printf("\nERROR: output file should have the extension .csv\n");
    return ;
  }

  read_csv(parsed[2], &bottom_spectra, &nrows, &ncols);

  n_bottom_spectra = nrows; 
}




// READ LANDSAT METADATA /file/

void run_read_landsat() {

  int i;
  bool success;

// Read METADATA.

  if (strcmp(parsed[2], "METADATA") != 0) {
    printf("\nERROR: READ LANDSAT METADATA /file/\n");
    return ;
  }

// Read file name. 

  if (! file_exists(parsed[3])) {
    printf("\nERROR: missing file '%s'\n", parsed[3]);
    return ;
  }

// Read landsat metadata. 

#if 0
  success = read_landsat_metadata(parsed[3], &sun_azimuth, &sun_elevation, &earth_sun_distance, 
      radiance_mult_band, radiance_add_band, quantize_cal_min_band, quantize_cal_max_band,
      radiance_minimum_band, radiance_maximum_band, MAX_GRIDS);

  if (! success) {
    printf("\nERROR: cannot read LANDSAT metadata from %s.\n", parsed[3]);

    // Clear metadata variables. 

    sun_azimuth = BAM_NODATA;
    sun_elevation = BAM_NODATA;
    earth_sun_distance = BAM_NODATA;

    for (i = 0; i < MAX_GRIDS; i++) {
      radiance_mult_band[i] = BAM_NODATA;
      radiance_add_band[i] = BAM_NODATA;
      quantize_cal_min_band[i] = BAM_NODATA;
      quantize_cal_max_band[i] = BAM_NODATA;
      radiance_minimum_band[i] = BAM_NODATA;
      radiance_maximum_band[i] = BAM_NODATA;
    }
  }

#endif
}


void run_read_commands() {
	
	FILE *fp;
  char line[MAX_LINE];
  int code;

	// Read command file name. 

	if (parsed[2][0] == '\0') {
		printf("\nERROR: READ COMMANDS /command file/\n");
		return ;
	}

	if (! file_exists(parsed[2])) {
  	printf("\nERROR: missing file '%s'\n", parsed[2]);
  	return ;
  }

  // Run commands. 

  fp = fopen(parsed[2], "r");

  while (!feof(fp)) {

    // Read line.

    fgets(line, sizeof(line), fp);

    trim(line);

    // Check for %, %n.

    if (line[0] == '%') {
      if (line[1] == '\0' || line[1] == '1')
        strcpy(line, prev1);
      else if (line[1] == '2')
        strcpy(line, prev2);
      else if (line[1] == '3')
        strcpy(line, prev3);
      else if (line[1] == '4')
        strcpy(line, prev4);
      else if (line[1] == '5')
        strcpy(line, prev5);
      else if (line[1] == '6')
        strcpy(line, prev6);
      else if (line[1] == '7')
        strcpy(line, prev7);
      else if (line[1] == '8')
        strcpy(line, prev8);
      else if (line[1] == '9')
        strcpy(line, prev9);
    }

    strcpy(prev9, prev8);
    strcpy(prev8, prev7);
    strcpy(prev7, prev6);
    strcpy(prev6, prev5);
    strcpy(prev5, prev4);
    strcpy(prev4, prev3);
    strcpy(prev3, prev2);
    strcpy(prev2, prev1);
    strcpy(prev1, line);

    // Display input. 

  	printf(">> %s\n", line);
    fflush(stdout);

    // Add to history. 

    linenoiseHistoryAdd(line);
    linenoiseHistorySave("history.bam");

    // Parse input.

  	code = parse_input(line);

  	switch (code) {
      case STOP:
        goto finish;
        break ;
  		case COMMENT:
  			break ;
  		case EXIT:
        // goto finish;
        exit(0);
#if 0
  			printf("\nAre you sure you want to exit? (yes/no) ");
  			fgets(line, sizeof(line), stdin);
 		 		if (strncmp(line, "yes", 3) == 0) {
  				goto finish; 
        }
#endif
  			break ;
  		default:
 				run_command();
 		}

    fflush(stdout);
 	}

finish:

  fclose(fp);
}



void run_read_soundings() {

  int i, k;
  float *x, *y, *z, scale = 1.0, chart_datum = 0.0, offsetX = 0.0, offsetY = 0.0;

// Check if we have previously read in sounding data. 

  if (read_soundings) {
  	free_float_array_2d(soundings, nsoundings);
  	meminuse -= (1.0e-6*(float) 3*nsoundings)*((float) sizeof(float));
  }

// Check file exists.  

  if (! file_exists(parsed[2])) {
  	printf("\nERROR: missing file '%s'\n", parsed[2]);
  	return ;
  }

// Read SCALE. 

  k = 3;

  if (strcmp(parsed[k], "SCALE") == 0) {
    scale = atof(parsed[++k]);
    k++;
  }  


// Check for chart datum.

  if (strcmp(parsed[k], "CHARTDATUM") == 0) {
  	chart_datum = atof(parsed[++k]);
    k++;
  }

// Check for offsets. 

  if (strcmp(parsed[k], "OFFSET") == 0) {
    offsetX = atof(parsed[++k]);
    offsetY = atof(parsed[++k]);
    k++;
  }

// Read sounding data. 

  read_xyz(parsed[2], &x, &y, &z, &nsoundings);

  allocate_float_array_2d(&soundings, nsoundings, 3);

  sounding_errors = (float*) malloc(nsoundings*sizeof(float));
  model_depths = (float*) malloc(nsoundings*sizeof(float));
  sounding_depths = (float*) malloc(nsoundings*sizeof(float));

  for (i = 0; i < nsoundings; i++) {
  	soundings[i][0] = x[i] + offsetX;
  	soundings[i][1] = y[i] + offsetY;
  	z[i] = (z[i] + chart_datum)*scale;
  	soundings[i][2] = z[i];
  	sounding_errors[i] = 0.0;
  }

  read_soundings = true;

// Display simple stats. 

  float soundmin = vec_min(z, nsoundings);
  float soundmax = vec_max(z, nsoundings);
  float soundmean = vec_mean(z, nsoundings);

  printf("\n  READ %d SOUNDINGS.\n", nsoundings);
  printf("  MINIMUM DEPTH = %f\n", soundmin);
  printf("  MAXIMUM DEPTH = %f\n", soundmax);
  printf("  MEAN DEPTH    = %f\n", soundmean);

// Update memory usage (in units of megabytes). 

  meminuse += (1.0e-6*(float) 3*nsoundings)*((float) sizeof(float));

// Free memory. 

  free(x);
  free(y);
  free(z);
}


void run_read_grid() {

  int i, j, k, n, ps, ncols, nrows;
  double wlon, slat, cellsize, spval;
  float **array, scale, offset;
  geogrid grid;
  bool convert2reflectance = false;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
  	printf("\nERROR: unable to store grid.\n");
  	return ;
  }

// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
  	if (! allocated_grids[n]) {
  		ps = n;
  		break;
  	}
  }

// Set name. 

  if (parsed[2][0] == '\0' || parsed[3][0] == '\0') {
  	printf("\nERROR: READ GRID /name/ /file/ [optional] SCALE /scale factor/ \
OFFSET /add offset/ THETA /sun elevation angle/\n");
  	return ;
  } else {

	// Do we already have a grid with this name? 

  	for (n = 0; n < MAX_GRIDS; n++) {
  		if (strcmp(grid_names[n], parsed[2]) == 0) {
  			printf("\nERROR: grid '%s' already in use.\n", parsed[2]);
  			return ;
  		}
  	}

  	strcpy(grid_names[ps], parsed[2]);
  }

// Check file exists.  

  if (! file_exists(parsed[3])) {
  	printf("\nERROR: missing file '%s'\n", parsed[3]);
  	return ;
  }

// Read in gridded data. 

  read_grid(parsed[3], &array, &ncols, &nrows, &wlon, 
              &slat, &cellsize, &spval);

  grid.ncols = ncols;
  grid.nrows = nrows;
  grid.cellsize = (float) cellsize;
  grid.wlon = (float) wlon;
  grid.slat = (float) slat;
  grid.elon = wlon + cellsize*((float) ncols - 1);
  grid.nlat = slat + cellsize*((float) nrows - 1);
  grid.nodata_value = (float) spval;
  grid.array = array;

// Check for SCALE. 

  k = 4; 

  if (strcmp(parsed[k], "SCALE") == 0) {
    scale = atof(parsed[++k]);
    k++;
  } else {
    scale = 1.0;
  }

// Check for OFFSET. 

  if (strcmp(parsed[k], "OFFSET") == 0) {
    offset = atof(parsed[++k]);
    k++;
  } else {
    offset = 0.0;
  }

// Set allocated. 

  allocated_grids[ps] = true;

// Transform array. 

  for (i = 0; i < nrows; i++) {
  	for (j = 0; j < ncols; j++) {
  		if (grid.array[i][j] == grid.nodata_value) {
  		  grid.array[i][j] = grid.nodata_value;
  		} else {
  		  grid.array[i][j] = grid.array[i][j]*scale + offset;  // radiance scale and offset. 
      }
  	}
  }

#if 0
  printf("\nmin, max = %f, %f\n", 
  	array_min(grid.array, nrows, ncols), 
  	array_max(grid.array, nrows, ncols));
#endif

  gridded_data[ps] = grid;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;
}



//
//    WRITE
//

// WRITE GRID /grid name/ /output file name/ REGION /xmin/ /ymin/ /xmax/ /ymax/

// WRITE LINE /grid name/ /output file name/ COORDINATES /xa/ /ya/ /xb/ /yb/

// WRITE POINTS /grid name/ /output file name/

void run_write() {

  if (strcmp(parsed[1], "GRID") == 0) {
    run_write_grid();
  } else if (strcmp(parsed[1], "LINE") == 0) {
    run_write_line();
  } else if (strcmp(parsed[1], "POINTS") == 0) {
    run_write_points();
  } else {
    printf("\nERROR: unknown command '%s'\n", parsed[1]);
  }
}

// WRITE POINTS /grid name/ /output file name/

void run_write_points() {

  int n, ps;
  float *xpts, *ypts, *zpts;
  bool present = false; 

  // Get position of grid name.
      
  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[2]) == 0) {
      present = true;
      ps = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
    return ;
  }

  // Check soundings are already in memory. 

  if (! read_soundings) {
    printf("\nERROR: soundings have not been read into memory.\n");
    return ;
  }

  // Interpolate grid at soundings. 

  xpts = malloc(nsoundings*sizeof(float));
  ypts = malloc(nsoundings*sizeof(float));
  zpts = malloc(nsoundings*sizeof(float));

  for (n = 0; n < nsoundings; n++) {
    xpts[n] = soundings[n][0];
    ypts[n] = soundings[n][1];

    if (xpts[n] > gridded_data[ps].wlon && 
    xpts[n] < gridded_data[ps].elon && 
    ypts[n] > gridded_data[ps].slat && 
    ypts[n] < gridded_data[ps].nlat) {
      zpts[n] = interp_bicubic(gridded_data[ps].array, 
              gridded_data[ps].ncols,
              gridded_data[ps].nrows, 
              gridded_data[ps].wlon, 
              gridded_data[ps].slat, 
              gridded_data[ps].cellsize, 
              gridded_data[ps].nodata_value, 
              xpts[n], 
              ypts[n]);
    } else {
      zpts[n] = gridded_data[ps].nodata_value;
    }
  }

  // Write points to file. 

  write_xyz(parsed[3], xpts, ypts, zpts, nsoundings);

  free(xpts);
  free(ypts);
  free(zpts);
}

// WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/ FILE //

void run_write_line() {

  int n, ps, npoints, i, j, ii, jj, np;
  float xa, ya, xb, yb, rad, grad, del, *xpts, *ypts, *radiance;
  bool present = false;

  if (parsed[2][0] == '\0'){
    printf("\nERROR: WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

  // Get position of grid name.
      
  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[2]) == 0) {
      present = true;
      ps = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
    return ;
  }

  // Read output file name. 

  if (parsed[3][0] == '\0'){
    printf("\nERROR: WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

  // Check output grid name is a netCDF (.nc or .nc4 extension)

  if (strstr(parsed[3], ".csv") == NULL && strstr(parsed[3], ".xyz") == NULL) { // Crude check.
    printf("\nERROR: output file should have the extension .csv or .xyz\n");
    return ;
  }

  // Read LINE.

  if (strcmp(parsed[4], "COORDINATES") != 0) {
    printf("\nERROR: WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

  if (strcmp(parsed[5], "SAND") == 0) {
    xa = sand_lon;
    ya = sand_lat;
    if (strcmp(parsed[6], "DEEP") == 0) {
      xb = deep_lon;
      yb = deep_lat;
      n = 7;
    } else {
      xb = atof(parsed[6]);
      yb = atof(parsed[7]);
      n = 8;
    }
  } else if (strcmp(parsed[5], "DEEP") == 0) {
    xa = deep_lon;
    ya = deep_lat;
    if (strcmp(parsed[6], "SAND") == 0) {
      xb = sand_lon;
      yb = sand_lat;
      n = 7;
    } else {
      xb = atof(parsed[6]);
      yb = atof(parsed[7]);
      n = 8;
    }
  } else {
    xa = atof(parsed[5]);
    ya = atof(parsed[6]);
    xb = atof(parsed[7]);
    yb = atof(parsed[8]);
    n = 9;
  }

// Extract points on the line from the grid. 

  float dist = sqrt(pow(xb - xa, 2) + pow(yb - ya, 2));
  npoints = 1 + floor(dist/gridded_data[ps].cellsize);

  xpts = (float*) malloc(npoints*sizeof(float));
  ypts = (float*) malloc(npoints*sizeof(float));
  radiance = (float*) malloc(npoints*sizeof(float));

// Gradient of line. 

  grad = (yb - ya)/(xb - xa); // TODO: check for vertical line.

  for (n = 0; n < npoints; n++) {
    // Compute the point on the line. 
    del = ((float)n)/((float) npoints);
    xpts[n] = xa + del*(xb - xa);
    ypts[n] = ya + grad*del*(xb - xa);
    // Compute grid coordinates of point. 
    // i = (ypts[n] - gridded_data[ps].slat)/gridded_data[ps].cellsize;
    // j = (xpts[n] - gridded_data[ps].wlon)/gridded_data[ps].cellsize;
    i = round( ((float) gridded_data[ps].nrows)*(ypts[n] - gridded_data[ps].slat)/
          (gridded_data[ps].nlat - gridded_data[ps].slat) );
    j = round( ((float) gridded_data[ps].ncols)*(xpts[n] - gridded_data[ps].wlon)/
          (gridded_data[ps].elon - gridded_data[ps].wlon) );
    // Radiance at i,j.
    rad = 0.0; 
    np = 0;
    for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
      for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
        if (gridded_data[ps].array[i + ii][j + jj] != gridded_data[ps].nodata_value) {
          rad += gridded_data[ps].array[i + ii][j + jj];
          np++;
        }
      }
    }
    if (np == 0)
      radiance[n] = 0.0;
    else
      radiance[n] = rad/((float) np);
  }

//  Write line to xyz file. 

  write_xyz(parsed[3], xpts, ypts, radiance, npoints);

// Free memory.

  free(xpts);
  free(ypts);
  free(radiance);
}

// WRITE GRID /grid name/ /output file name/ REGION /xmin/ /ymin/ /xmax/ /ymax/

void run_write_grid() {

  int n, ps;
  float *lats, *lons, xmin, ymin, xmax, ymax;
  bool present = false, region = false;

  if (parsed[2][0] == '\0'){
    printf("\nERROR: WRITE GRID /grid name/ /output file name/\n");
    return ;
  }
  	
  // Get position of grid name.
  	  
	for (n = 0; n < MAX_GRIDS; n++) {
	 	if (strcmp(grid_names[n], parsed[2]) == 0) {
    	present = true;
			ps = n;
			break;
	 	}
	}

	if (! present) {
	  printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
	  return ;
	}

  // Read output file name. 

  if (parsed[3][0] == '\0'){
    printf("\nERROR: WRITE GRID /grid name/ /output file name/\n");
    return ;
  }

  // Check output grid name is a netCDF (.nc or .nc4 extension)

  if (strstr(parsed[3], ".nc") == NULL && strstr(parsed[3], ".asc") == NULL) { // Crude check.
    printf("\nERROR: output file should have the extension .nc, .nc4 or .asc\n");
    return ;
  }

  if (strcmp(parsed[4], "REGION") == 0) {
  	region = true;
  	xmin = atof(parsed[5]);
  	ymin = atof(parsed[6]);
  	xmax = atof(parsed[7]);
  	ymax = atof(parsed[8]);
  }

  // Write output grid. 

  if (region) {
  	int ncols = (xmax - xmin)/gridded_data[ps].cellsize;
  	int nrows = (ymax - ymin)/gridded_data[ps].cellsize;
  	lons = (float*) malloc(ncols*sizeof(float));
  	lats = (float*) malloc(nrows*sizeof(float));

  	for (n = 0; n < ncols; n++)
  		lons[n] = xmin + ((float) n)*gridded_data[ps].cellsize;

  	for (n = 0; n < nrows; n++)
  		lats[n] = ymin + ((float) n)*gridded_data[ps].cellsize;
  } else {
  	lons = (float*) malloc(gridded_data[ps].ncols*sizeof(float));
  	lats = (float*) malloc(gridded_data[ps].nrows*sizeof(float));

  	for (n = 0; n < gridded_data[ps].ncols; n++)
  		lons[n] = gridded_data[ps].wlon + ((float) n)*gridded_data[ps].cellsize;

  	for (n = 0; n < gridded_data[ps].nrows; n++)
  		lats[n] = gridded_data[ps].slat + ((float) n)*gridded_data[ps].cellsize;

  	write_nc(parsed[3], 
  		gridded_data[ps].array, 
  		gridded_data[ps].ncols,
  		gridded_data[ps].nrows,
  		lons, lats, 
  		gridded_data[ps].nodata_value);
  }

  free(lons);
  free(lats);
}


//
//    SET
//

void run_set() {

  int n; 
  bool grid_present; 

// NSMOOTH (default 1)

  if (strcmp(parsed[1], "NSMOOTH") == 0) {
    samodel_n_smoothing_radius = atoi(parsed[2]);
    printf("\n    NSMOOTH = %d\n", samodel_n_smoothing_radius);
    return ;
  }  

// NSPATIAL (default 2)

  if (strcmp(parsed[1], "NSPATIAL") == 0) {
    samodel_n_spatial = atoi(parsed[2]);
    printf("\n    NSPATIAL = %d\n", samodel_n_spatial);
    return ;
  }  

// NBOTTOMS (default 3)

  if (strcmp(parsed[1], "NBOTTOMS") == 0) {
    samodel_n_bottoms = atoi(parsed[2]);
    printf("\n    NBOTTOMS = %d\n", samodel_n_bottoms);
    return ;
  }  

// K (ie. SET K coastal 0.09)

  if (strcmp(parsed[1], "K") == 0) {
    grid_present = false; 
    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[2]) == 0) {
        Alpha[n] = atof(parsed[3]);
        grid_present = true; 
        break;
      }
    }
    if (! grid_present) {
      printf("\nERROR: grid '%s' not found.", parsed[2]);
    }
    return ;
  }

// LSM (ie. SET K coastal 0.09)

  if (strcmp(parsed[1], "LSM") == 0) {
    grid_present = false; 
    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[2]) == 0) {
        Lsm[n] = atof(parsed[3]);
        grid_present = true; 
        break;
      }
    }
    if (! grid_present) {
      printf("\nERROR: grid '%s' not found.", parsed[2]);
    }
    return ;
  }

// REGRESSIONWEIGHT

  if (strcmp(parsed[1], "REGRESSIONWEIGHT") == 0) {
    weight_regression = atof(parsed[2]);
    return ;
  }

// PARAMBOUNDWEIGHT

  if (strcmp(parsed[1], "PARAMBOUNDWEIGHT") == 0) {
    weight_bound_error = atof(parsed[2]);
    return ;
  }

// COEFFSIZEWEIGHT

  if (strcmp(parsed[1], "COEFFSIZEWEIGHT") == 0) {
    weight_coeff_size_preference = atof(parsed[2]);
    return ;
  }

// NBOUNDARY 

  if (strcmp(parsed[1], "NBOUNDARY") == 0) {
    nboundary = atof(parsed[2]);
    return ;
  }

// CONSTRAINTSWEIGHT

  if (strcmp(parsed[1], "CONSTRAINTSWEIGHT") == 0) {
    weight_constraints = atof(parsed[2]);
    return ;
  }

// SOUNDINGSWEIGHT

  if (strcmp(parsed[1], "SOUNDINGSWEIGHT") == 0) {
    weight_soundings = atof(parsed[2]);
    return ;
  }

// RELATIVEERROR

  if (strcmp(parsed[1], "RELATIVEERROR") == 0) {
    if (strcmp(parsed[2], "TRUE") == 0) {
      relative_error = true;
      printf("\n    RELATIVEERROR = true\n");
    } else if (strcmp(parsed[2], "FALSE") == 0) {
      relative_error = false;
      printf("\n    RELATIVEERROR = false\n");
    } else {
      printf("\nERROR: expected true/false input.\n");
    }
    return ;
  }

// SAMAXITERATIONS

  if (strcmp(parsed[1], "SAMAXITERATIONS") == 0) {
    sa_max_iterations = atof(parsed[2]);
    printf("\n    SAMAXITERATIONS = %d\n", sa_max_iterations);
    return ;
  }

// LSBSPREAD

  if (strcmp(parsed[1], "LSBSPREAD") == 0) {
    lsb_spread = atof(parsed[2]);
    printf("\n    LSBSPREAD = %6.3f\n", lsb_spread);
    return ;
  }

// PANNORMALISE

  if (strcmp(parsed[1], "PANNORMALISE") == 0) {
    if (strcmp(parsed[2], "true") == 0) {
      pan_normalise = true;
      printf("\n    PAN NORMALISE = true\n");
    } else if (strcmp(parsed[2], "false") == 0)  {
      pan_normalise = false;
      printf("\n    PAN NORMALISE = false\n");
    } else {
      printf("\nERROR: expected true/false input.\n");
    }
    return ;
  }

// LINEWIDTH

  if (strcmp(parsed[1], "LINEWIDTH") == 0) {
    line_width = atoi(parsed[2]);
    printf("\n    LINE WIDTH = %d\n", line_width);
    return ;
  }

// TEXTSIZE

  if (strcmp(parsed[1], "TEXTSIZE") == 0) {
    textsize = atof(parsed[2]);
    printf("\n     TEXT SIZE = %.2f\n", textsize);
    return ;
  }

// DEEP

  if (strcmp(parsed[1], "DEEP") == 0) {
    deep_lon = atof(parsed[2]);
    deep_lat = atof(parsed[3]);
    printf("\n    DEEP LON/LAT = %f, %f\n", deep_lon, deep_lat);
    return ;
  }

// SAND

  if (strcmp(parsed[1], "SAND") == 0) {
    sand_lon = atof(parsed[2]);
    sand_lat = atof(parsed[3]);
    printf("\n    SAND LON/LAT = %f, %f\n", sand_lon, sand_lat);
    return ;
  }

// FALLOFF

	if (strcmp(parsed[1], "FALLOFF") == 0) {
		falloff_param = atof(parsed[2]);
		printf("\n    FALLOFF = %f\n", falloff_param);
		return ;
	}

// MODELRADIUS 

	if (strcmp(parsed[1], "MODELRADIUS") == 0) {
		depth_radius = atoi(parsed[2]);
		printf("\n    MODEL RADIUS = %d\n", depth_radius);
		return ;
	}

// SATEMP

	if (strcmp(parsed[1], "SAINITIALTEMPERATURE") == 0) {
		initial_temperature = atof(parsed[2]);
		printf("\n    SET SA INITIAL TEMPERATURE = %f\n", initial_temperature);
		return ;
	}

// COOLINGRATE

	if (strcmp(parsed[1], "SACOOLINGRATE") == 0) {
		cooling_rate = atof(parsed[2]);
		printf("\n    SET SA COOLING RATE = %f\n", cooling_rate);
		return ;
	}

// SAPARAMRANGE

	if (strcmp(parsed[1], "SAPARAMRANGE") == 0) {
		parameter_range = atof(parsed[2]);
		printf("\n    SET SA PARAMETER RANGE = %f\n", parameter_range);
		return ;
	}

// BACKGROUND

	if (strcmp(parsed[1], "BACKGROUND") == 0) {
		if (strcmp(parsed[2], "WHITE") == 0) {
			background = WHITE;
			printf("\n    SET BACKGROUND WHITE\n");
		} else if (strcmp(parsed[2], "BLACK") == 0){
			background = BLACK;
			printf("\n    SET BACKGROUND BLACK\n");
		} else {
			printf("\nERROR: unknown command '%s'\n", parsed[2]);
		}
		return ;
	}

// NOISETHRESHOLD

	if (strcmp(parsed[1], "NOISETHRESHOLD") == 0) {
		noise_threshold = atof(parsed[2]);
		printf("\n    NOISETHRESHOLD = %f\n", noise_threshold);
		return ;
	}

// MINRADIANCETHRESHOLD

  if (strcmp(parsed[1], "MINRADIANCETHRESHOLD") == 0) {
  	min_radiance_threshold = atof(parsed[2]);
  	printf("\n    MINRADIANCETHRESHOLD = %f\n", min_radiance_threshold);
  	return ;
  }

// MAXRADIANCETHRESHOLD

  if (strcmp(parsed[1], "MAXRADIANCETHRESHOLD") == 0) {
  	max_radiance_threshold = atof(parsed[2]);
  	printf("\n    MAXRADIANCETHRESHOLD = %f\n", max_radiance_threshold);
  	return ;
  }

// PAGESIZE

		if (strcmp(parsed[1], "PAGESIZE") == 0) {
			pagesize = atof(parsed[2]);
			printf("\n  PAGESIZE = %f\n", pagesize);
			return ;
		}

// DEPTHMIN

  if (strcmp(parsed[1], "DEPTHMIN") == 0) {
  	depthmin = atof(parsed[2]);
  	printf("\n  DEPTHMIN = %f\n", depthmin);
  	return ;
  }

// DEPTHMAX

  if (strcmp(parsed[1], "DEPTHMAX") == 0) {
  	depthmax = atof(parsed[2]);
  	printf("\n  DEPTHMAX = %f\n", depthmax);
  	return ;
  }

// MODEL OPTIMISATION

  if (strcmp(parsed[1], "MODEL") == 0 && 
  		(strcmp(parsed[2], "OPTIMISATION") == 0 || 
  		 strcmp(parsed[2], "OPTIMIZATION") == 0)) {
  	if (strcmp(parsed[3], "SimulatedAnnealing") == 0 || strcmp(parsed[3], "SA") == 0) {
  		model_optimisation = SIMULATED_ANNEALING;
  		printf("\nMODEL OPTIMISATION = SimulatedAnnealing\n");
  	} else if (strcmp(parsed[3], "RandomSearch") == 0 || strcmp(parsed[3], "RS") == 0) {
  		model_optimisation = RANDOM_SEARCH;
  		printf("\nMODEL OPTIMISATION = RandomSearch\n");
  	} else {
  		printf("\nERROR: unknown MODEL OPTIMISATION '%s'\n", parsed[3]);
  		return ;
  	}
  }

  printf("\nERROR: unknown command '%s'\n", parsed[1]);
}


//
//    PRINT
//

void run_print() {

  int n, ps;
  bool present = false;

// REGRESSIONWEIGHT

  if (strcmp(parsed[1], "REGRESSIONWEIGHT") == 0) {
    printf("\n    REGRESSION WEIGHT = %f\n", weight_regression);
    return ;
  }

// PARAMBOUNDWEIGHT

  if (strcmp(parsed[1], "PARAMBOUNDWEIGHT") == 0) {
    printf("\n    PARAMETER BOUND WEIGHT = %f\n", weight_bound_error);
    return ;
  }

// COEFFSIZEWEIGHT

  if (strcmp(parsed[1], "COEFFSIZEWEIGHT") == 0) {
    printf("\n    COEFFICIENT SIZE WEIGHT = %f\n", weight_coeff_size_preference);
    return ;
  }

// NBOUNDARY 

  if (strcmp(parsed[1], "NBOUNDARY") == 0) {
    printf("\n    N BOUNDARY = %d\n", nboundary);
    return ;
  }

// CONSTRAINTSWEIGHT

  if (strcmp(parsed[1], "CONSTRAINTSWEIGHT") == 0) {
    printf("\n    CONSTRAINTS WEIGHT = %f\n", weight_constraints);
    return ;
  }

// SOUNDINGSWEIGHT

  if (strcmp(parsed[1], "SOUNDINGSWEIGHT") == 0) {
    printf("\n    SOUNDING WEIGHT = %f\n", weight_soundings);
    return ;
  }

// K (Spectral attenuation coefficients.)

  if (strcmp(parsed[1], "K") == 0) {
    printf("\n");
    printf("    BAND NAME        K    \n");
    printf("--------------------------\n");
    for (n = 0; n < MAX_GRIDS; n++) {
      if (allocated_grids[n] && Alpha[n] != 0.0)
        printf("%12.12s        %.3f\n", grid_names[n], Alpha[n]);
    }
    printf("\n");
    return ;
  }

// LANDSAT8

  if (strcmp(parsed[1], "LANDSAT8") == 0) {
    printf("\n");
    printf(" BAND       NAME        WIDTH   MEAN\n");
    printf(" -----------------------------------\n");
    printf("    1    COASTAL    435 - 451    443\n");
    printf("    2       BLUE    452 - 512    482\n");
    printf("    3      GREEN    533 - 590    562\n");
    printf("    4        RED    636 - 673    655\n");
    printf("    5        NIR    851 - 879    865\n\n");
    return ;
  }

// LINEWIDTH

  if (strcmp(parsed[1], "LINEWIDTH") == 0) {
#if PGPLOT
    cpgqlw(&n);
#endif
    printf("\n    LINE WIDTH = %d\n", n);
    return ;
  }

// DEEP

  if (strcmp(parsed[1], "DEEP") == 0) {
    printf("\n    DEEP LON/LAT = %f, %f\n", deep_lon, deep_lat);
    return ;
  }

// SAND

  if (strcmp(parsed[1], "SAND") == 0) {
    printf("\n    SAND LON/LAT = %f, %f\n", sand_lon, sand_lat);
    return ;
  }

// FALLOFF

  if (strcmp(parsed[1], "FALLOFF") == 0) {
  	printf("\n    FALLOFF = %f\n", falloff_param);
  	return ;
  }

// MODELRADIUS

  if (strcmp(parsed[1], "MODELRADIUS") == 0) {
  	printf("\n    MODEL RADIUS = %d\n", depth_radius);
  	return ;
  }  

// BACKGROUND

  if (strcmp(parsed[1], "BACKGROUND") == 0) {
  	if (background == WHITE)
  		printf("\n    BACKGROUND = WHITE\n");
  	else
  		printf("\n    BACKGROUND = BLACK\n");
  	return ;
  }    

// NOISE THRESHOLD

  if (strcmp(parsed[1], "NOISETHRESHOLD") == 0) {
  	printf("\n    NOISETHRESHOLD = %f\n", noise_threshold);
  	return ;
  }    

// MINRADIANCETHRESHOLD

  if (strcmp(parsed[1], "MINRADIANCETHRESHOLD") == 0) {
  	printf("\n    MINRADIANCETHRESHOLD = %f\n", min_radiance_threshold);
  	return ;
  }  

// MAXRADIANCETHRESHOLD

  if (strcmp(parsed[1], "MAXRADIANCETHRESHOLD") == 0) {
  	printf("\n    MAXRADIANCETHRESHOLD = %f\n", max_radiance_threshold);
  	return ;
  }  

// DEPTHMIN

  if (strcmp(parsed[1], "DEPTHMIN") == 0) {
  	printf("\n    DEPTHMIN = %f\n", depthmin);
  	return ;
  }

// DEPTHMAX

  if (strcmp(parsed[1], "DEPTHMAX") == 0) {
  	printf("\n    DEPTHMAX = %f\n", depthmax);
  	return ;
  }

// MODEL OPTIMISATION

  if (strcmp(parsed[1], "MODEL") == 0 && 
  	strcmp(parsed[2], "OPTIMISATION") == 0) {
  	if (model_optimisation == SIMULATED_ANNEALING)
  		printf("\n    MODEL OPTIMISATION = SimulatedAnnealing\n");
  	else if (model_optimisation == RANDOM_SEARCH)
  		printf("\n    MODEL OPTIMISATION = RandomSearch\n");

  	return ;
  }

// MEMINUSE

  if (strcmp(parsed[1], "MEMINUSE") == 0) {
  	printf("\n    MEMORY IN USE = %f (Mb)\n", meminuse);
  	return ;
  } 

// GRIDSTATS

  if (strcmp(parsed[1], "GRIDSTATS") == 0) {

  	// Check grid is present. 

	for (n = 0; n < MAX_GRIDS; n++) {
		if (strcmp(grid_names[n], parsed[2]) == 0) {
			present = true;
			ps = n;
			break;
		}
	}

	if (! present) {
		printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
		return ;
	}

	printf("\n-------------------");
	printf("\n  GRID STATISTICS  ");
	printf("\n-------------------");
	printf("\nname         = %s", parsed[2]);
	printf("\nncols        = %d", gridded_data[ps].ncols);
	printf("\nnrows        = %d", gridded_data[ps].nrows);
	printf("\nW lon        = %f", gridded_data[ps].wlon);
	printf("\nS lat        = %f", gridded_data[ps].slat);
	printf("\nE lon        = %f", gridded_data[ps].elon);
	printf("\nN lat        = %f", gridded_data[ps].nlat);
	printf("\nresolution   = %f", gridded_data[ps].cellsize);
	printf("\nNODATA_value = %f\n", gridded_data[ps].nodata_value);

	return ;
  }

// GRIDSINUSE

  if (strcmp(parsed[1], "GRIDSINUSE") == 0) {

		for (n = 0; n < MAX_GRIDS; n++) {
			if (allocated_grids[n]) {
				printf("\n#%d : NAME =  %s, NCOLS = %d, NROWS = %d, RESOLUTION = %f", 
					n + 1,
					grid_names[n], 
					gridded_data[n].ncols,
					gridded_data[n].nrows,
					gridded_data[n].cellsize);
			} else {
				printf("\n#%d : EMPTY", n + 1);
			}
		}
		printf("\n");
  	return ;
  } 


  printf("\nERROR: unknown command '%s'\n", parsed[1]);
}

//
//    PLOT
//

// PLOT GRID /grid name/ [optional] SOUNDINGS ...
// 		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/

// PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ ... 
//    [optional] THRESHOLD /min radiance/ /max radiance/ ... 
//		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/

// PLOT LOGSCATTER /band name 1/ /band name 2/ TYPE /water type/ ... 
//    COORDINATES /x1/ /y1/ /x2/ /y2/ ...

// PLOT QA /qa grid name/ [optional] REGION /xmin/ /ymin/ /xmax/ /ymax/

void run_plot() { 

	if (strcmp(parsed[1], "GRID") == 0) {
		run_plot_grid();
		return ;
	} else if (strcmp(parsed[1], "LINE") == 0) {
		run_plot_line();
		return ;
  } else if (strcmp(parsed[1], "LOGSCATTER") == 0) {
    run_plot_logscatter();
	} else {
		printf("\nERROR: unknown command '%s'\n", parsed[1]);
		return ;
	}
}



// PLOT LOGSCATTER /band name 1/ /band name 2/ ... 
//    COORDINATES /x1/ /y1/ /x2/ /y2/ ... [optional] Lsm /Lsm band 1/
//    /Lsm band 2/

void run_plot_logscatter() {

  int n, m, ii, jj, i, j, k, ps1, ps2, segpts, np, npoints, nlines, 
    ncoords, segsizes[LINE_MAX_POINTS], 
    ncolours = 9, colours[] = {2,6,8,4,3,11,12,13,5};
  float coords[LINE_MAX_POINTS], dist, x, y, *radiance1, *radiance2, del, 
    grad, rad, xa, xb, ya, yb, Lsm1, Lsm2, error, *seg1, *seg2,
    xmin, xmax, ymin, ymax, a, b, r;
  bool present, success;

  for (n=0; n < LINE_MAX_POINTS; n++)
    segsizes[n] = 0;

  if (parsed[2][0] == '\0') {
    printf("\nERROR: PLOT LOGSCATTER /band name 1/ /band name 2/ \
COORDINATES /x1/ /y1/ /x2/ /y2/ ... [optional] LSM /Lsm band 1/ /Lsm band 2/\n");
    return ;
  }

// Check first grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[2]) == 0) {
      present = true;
      ps1 = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
   return ;
  }

// Check second grid is present. 

  present = false;
  
  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps2 = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
   return ;
  }

// Read coordinates of lines. 

  if (strcmp(parsed[4], "COORDINATES") != 0) {
    printf("\nERROR: PLOT LOGSCATTER /band name 1/ /band name 2/ \
COORDINATES /x1/ /y1/ /x2/ /y2/ ... [optional] LSM /Lsm band 1/ /Lsm band 2/\n");
    return ;
  }

  n = 5;
  nlines = 0;
  ncoords = 0;

  while(parsed[n][0] != '\0' && strcmp(parsed[n], "LSM") != 0) {

    nlines++;

    if ((ncoords + 4) >= LINE_MAX_POINTS) {
      printf("\nERROR: too many lines.\n");
      return ;
    }

    if (parsed[n][0] == '\0' || strcmp(parsed[n], "LSM") == 0) { 
      break ;
    } else if (strcmp(parsed[n], "SAND") == 0) {
      coords[ncoords++] = sand_lon;
      coords[ncoords++] = sand_lat;
      if (strcmp(parsed[n + 1], "DEEP") == 0) {
        coords[ncoords++] = deep_lon;
        coords[ncoords++] = deep_lat;
        n += 2;
      } else {
        coords[ncoords++] = atof(parsed[n + 1]);
        coords[ncoords++] = atof(parsed[n + 2]);
        n += 3;
      }
    } else if (strcmp(parsed[n], "DEEP") == 0) {
      coords[ncoords++] = deep_lon;
      coords[ncoords++] = deep_lat;
      if (strcmp(parsed[n + 1], "SAND") == 0) {
        coords[ncoords++] = sand_lon;
        coords[ncoords++] = sand_lat;
        n += 2;
      } else {
        coords[ncoords++] = atof(parsed[n + 1]);
        coords[ncoords++] = atof(parsed[n + 2]);
        n += 3;  
      }
    } else {
      coords[ncoords++] = atof(parsed[n]);
      coords[ncoords++] = atof(parsed[n + 1]);
      coords[ncoords++] = atof(parsed[n + 2]);
      coords[ncoords++] = atof(parsed[n + 3]);
      n += 4;
    }
  }

// Check grids are on the same extents. 

  if (gridded_data[ps1].nrows != gridded_data[ps2].nrows || 
    gridded_data[ps1].ncols != gridded_data[ps2].ncols || 
    fabs(gridded_data[ps1].slat - gridded_data[ps2].slat) > extent_eps || 
    fabs(gridded_data[ps1].nlat - gridded_data[ps2].nlat) > extent_eps || 
    fabs(gridded_data[ps1].wlon - gridded_data[ps2].wlon) > extent_eps || 
    fabs(gridded_data[ps1].elon - gridded_data[ps2].elon) > extent_eps || 
    fabs(gridded_data[ps1].cellsize - gridded_data[ps2].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Compute Lsm if not specified by the user.

  if (strcmp(parsed[n], "LSM") == 0) {
    Lsm1 = atof(parsed[++n]);
    Lsm2 = atof(parsed[++n]);
  } else if (Lsm[ps1] != 0.0 && Lsm[ps2] != 0.0) {
    Lsm1 = Lsm[ps1];
    Lsm2 = Lsm[ps2];
  } else {
    Lsm1 = array_min2(
        gridded_data[ps1].array, 
        gridded_data[ps1].nrows, 
        gridded_data[ps1].ncols, 
        gridded_data[ps1].nodata_value);
    Lsm2 = array_min2(
        gridded_data[ps2].array, 
        gridded_data[ps2].nrows, 
        gridded_data[ps2].ncols, 
        gridded_data[ps2].nodata_value);
  }

//  printf("\nLsm1,Lsm2 = %f, %f\n", Lsm1, Lsm2);

// Extract points on each line on each grid. 

  ii = 0;
  npoints = 0;

  for (k = 0; k < nlines; k++) {
    xa = coords[ii++];
    ya = coords[ii++];
    xb = coords[ii++];
    yb = coords[ii++];
    dist = sqrt(pow(xb - xa,2) + pow(yb - ya,2));
    npoints += floor(dist/gridded_data[ps1].cellsize);
  }

  radiance1 = (float*) malloc(npoints*sizeof(float));
  radiance2 = (float*) malloc(npoints*sizeof(float));

  ii = 0;
  m = 0;

  for (k = 0; k < nlines; k++) {
    xa = coords[ii++];
    ya = coords[ii++];
    xb = coords[ii++];
    yb = coords[ii++];

    dist = sqrt(pow(xb - xa,2) + pow(yb - ya,2));
    segpts = floor(dist/gridded_data[ps1].cellsize);
    segsizes[k] = segpts;
    grad = (yb - ya)/(xb - xa);

    for (n = 0; n < segpts; n++) {

      if (m > npoints - 1)
        break ;

      // Compute the point on the line. 

      del = ((float)n)/((float) segpts);
      x = xa + del*(xb - xa);
      y = ya + grad*del*(xb - xa);

      // Compute grid coordinates of point. 

      i = round( ((float) gridded_data[ps1].nrows)*(y - gridded_data[ps1].slat)/
          (gridded_data[ps1].nlat - gridded_data[ps1].slat) );
      j = round( ((float) gridded_data[ps1].ncols)*(x - gridded_data[ps1].wlon)/
          (gridded_data[ps1].elon - gridded_data[ps1].wlon) );

      // Radiance at i,j for band 1. 

      rad = 0.0; 
      np = 0;

      for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
        for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
          if (gridded_data[ps1].array[i + ii][j + jj] != gridded_data[ps1].nodata_value) {
            rad += gridded_data[ps1].array[i + ii][j + jj];
            np++;
          }
        }
      }

      if (np == 0)
        radiance1[m] = 0.0;
      else
        radiance1[m] = rad/((float) np);

      // Radiance at i,j for band 2. 
      
      rad = 0.0; 
      np = 0;

      for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
        for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
          if (gridded_data[ps2].array[i + ii][j + jj] != gridded_data[ps2].nodata_value) {
            rad += gridded_data[ps2].array[i + ii][j + jj];
            np++;
          }
        }
      }

      if (np == 0)
        radiance2[m] = 0.0;
      else
        radiance2[m] = rad/((float) np);

      m++;
    }
  }

// Load PGPLOT. 
#if PGPLOT
  error = cpgopen("?");
  if (error < 1) {
    printf("\nERROR: cannot load graphics library.\n");
    return ;
  }

// Set the page size and aspect ratio to 1.0. 

  cpgpap(pagesize, 1.0);

// Background and text colours. 

  if (background == BLACK) {
    cpgscr(0, 0.0, 0.0, 0.0);
    cpgscr(1, 1.0, 1.0, 1.0);
  } else {
    cpgscr(0, 1.0, 1.0, 1.0);
    cpgscr(1, 0.0, 0.0, 0.0);
  }

// Plot the log-transformed scatter plot of radiances. 

  for (n = 0; n < npoints; n++) {
    if (radiance1[n] < Lsm1)
      radiance1[n] = 0.0; 
    else
      radiance1[n] = log(radiance1[n] - Lsm1);

    if (radiance2[n] < Lsm2)
      radiance2[n] = 0.0; 
    else
      radiance2[n] = log(radiance2[n] - Lsm2);
  }

  int offset = 0;

  xmin = vec_min(radiance1, npoints);
  xmax = vec_max(radiance1, npoints);
  ymin = vec_min(radiance2, npoints);
  ymax = vec_max(radiance2, npoints);


  seg1 = (float*) malloc(npoints*sizeof(float));
  seg2 = (float*) malloc(npoints*sizeof(float));

  bool first = true;
  char xlabel[128], ylabel[128];

  strcpy(xlabel, "log(Ls\\d");
  strcat(xlabel, parsed[2]);
  strcat(xlabel, "\\u - Lsm\\d");
  strcat(xlabel, parsed[2]);
  strcat(xlabel, "\\u)");

  strcpy(ylabel, "log(Ls\\d");
  strcat(ylabel, parsed[3]);
  strcat(ylabel, "\\u - Lsm\\d");
  strcat(ylabel, parsed[3]);
  strcat(ylabel, "\\u)");

  for (k = 0; k < nlines; k++) {

    for (n = 0; n < segsizes[k]; n++) {
      seg1[n] = radiance1[offset + n];
      seg2[n] = radiance2[offset + n];
    }

    plot_scatter("LOG-LOG SCATTER PLOT", xlabel, ylabel, 
      colours[k%ncolours], 16, 
      segsizes[k], seg1, seg2, line_width,
      xmin, xmax, ymin, ymax, first);

    offset += segsizes[k];
    first = false;
  }

  if (nlines == 1) {

    // Compute the lines of best fit.

    success = linear_fit(radiance1, radiance2, segsizes[0], &a, &b, &r);

    printf("\nm, b, r = %.3f, %.3f, %.3f\n", a, b, r);

    // Plot the lines of best fit.
    cpgsls(2);
    cpgsci(3);
    cpgmove(xmin, a*xmin + b);
    cpgdraw(xmax, a*xmax + b);
    cpgsls(1);
  }

// Close the graphics device. 

  cpgclos();

  free(seg1);
  free(seg2);
#endif

  free(radiance1);
  free(radiance2);
}

// PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ ... 
//    [optional] THRESHOLD /min radiance/ /max radiance/ ... 
//		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/

void run_plot_line() {

	int n, ps, error, palette, i, j, np, npoints, ii, jj;
	bool present = false;
	float x[2], y[2], xa, ya, xb, yb, xmin = 0.0, xmax = 0.0, 
	  ymin = 0.0, ymax = 0.0, *xpts, *ypts, *radiance, grad, 
	  rad, del, min_threshold = 0.0, max_threshold = 0.0, 
	  max_slope, slope, radmin, radmax, noise;
	char title[256], snum[32];

// Read the grid. 

  if (parsed[2][0] == '\0') {
  	printf("\nERROR: PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ \
[optional] THRESHOLD /min radiance/ /max radiance/ REGION /xmin/ /ymin/ \
/xmax/ /ymax/ PALETTE /palette name/\n");
  	return ;
  } else {

  	// Check grid is present. 

		for (n = 0; n < MAX_GRIDS; n++) {
			if (strcmp(grid_names[n], parsed[2]) == 0) {
				present = true;
				ps = n;
				break;
			}
		}

		if (! present) {
			printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
			return ;
		}
  }

// Read the line coordinates. 

  if (strcmp(parsed[3], "COORDINATES") != 0) {
  	printf("\nERROR: PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ \
[optional] THRESHOLD /min radiance/ /max radiance/ REGION /xmin/ /ymin/ \
/xmax/ /ymax/ PALETTE /palette name/\n");
  	return ;
  } else {
    if (strcmp(parsed[4], "SAND") == 0) {
      xa = sand_lon;
      ya = sand_lat;
      if (strcmp(parsed[5], "DEEP") == 0) {
        xb = deep_lon;
        yb = deep_lat;
        n = 6;
      } else {
        xb = atof(parsed[5]);
        yb = atof(parsed[6]);
        n = 7;       
      }
    } else if (strcmp(parsed[4], "DEEP") == 0) {
      xa = deep_lon;
      ya = deep_lat;
      if (strcmp(parsed[5], "SAND") == 0) {
        xb = sand_lon;
        yb = sand_lat;
        n = 6;
      } else {
        xb = atof(parsed[5]);
        yb = atof(parsed[6]);
        n = 7;       
      }
    } else {
  	  xa = atof(parsed[4]);
  	  ya = atof(parsed[5]);
  	  xb = atof(parsed[6]);
  	  yb = atof(parsed[7]);
      n = 8;
    }
  }

  // printf("\nxa,ya = %f, %f", xa, ya);
  // printf("\nxb,yb = %f, %f\n", xb, yb);

// Read optional arguments. 

  if (strcmp(parsed[n], "THRESHOLD") == 0) {
  	min_threshold = atof(parsed[++n]);
  	max_threshold = atof(parsed[++n]);
  	n++;
  }

  if (strcmp(parsed[n], "REGION") == 0) {
  	xmin = atof(parsed[++n]);
  	ymin = atof(parsed[++n]);
  	xmax = atof(parsed[++n]);
  	ymax = atof(parsed[++n]);
  	n++;
  }

  palette = GRAYSCALE; // Default.

  if (strcmp(parsed[n], "PALETTE") == 0) {
    if (strcmp(parsed[++n], "jet") == 0)
      palette = JET;
    else if (strcmp(parsed[n], "inverse_jet") == 0)
      palette = INVERSEJET;
    else if (strcmp(parsed[n], "grayscale") == 0)
      palette = GRAYSCALE;
    else if (strcmp(parsed[n], "ocean") == 0)
      palette = OCEAN;
    else if (strcmp(parsed[n], "gmt") == 0)
      palette = GMT;
    else if (strcmp(parsed[n], "rygb") == 0)
      palette = RYGB;
    else if (strcmp(parsed[n], "wor") == 0)
      palette = WOR;
    else if (strcmp(parsed[n], "viridis") == 0)
      palette = VIRIDIS;
    else if (strcmp(parsed[n], "plasma") == 0)
      palette = PLASMA;
    else {
      printf("\nERROR: unknown colour palette '%s'\n", parsed[n]);
    }
  }

// Extract points on the line from the grid. 

  npoints = 1 + floor(sqrt(pow(xb - xa, 2) + pow(yb - ya, 2))/gridded_data[ps].cellsize);

  xpts = (float*) malloc(npoints*sizeof(float));
  ypts = (float*) malloc(npoints*sizeof(float));
  radiance = (float*) malloc(npoints*sizeof(float));

// Gradient of line. 

  grad = (yb - ya)/(xb - xa); // TODO: check for vertical line.

  for (n = 0; n < npoints; n++) {
  	// Compute the point on the line. 
  	del = ((float)n)/((float) npoints);
  	xpts[n] = xa + del*(xb - xa);
  	ypts[n] = ya + grad*del*(xb - xa);
  	// Compute grid coordinates of point. 
  	// i = (ypts[n] - gridded_data[ps].slat)/gridded_data[ps].cellsize;
  	// j = (xpts[n] - gridded_data[ps].wlon)/gridded_data[ps].cellsize;
  	i = round( ((float) gridded_data[ps].nrows)*(ypts[n] - gridded_data[ps].slat)/
          (gridded_data[ps].nlat - gridded_data[ps].slat) );
  	j = round( ((float) gridded_data[ps].ncols)*(xpts[n] - gridded_data[ps].wlon)/
          (gridded_data[ps].elon - gridded_data[ps].wlon) );
  	// Radiance at i,j.
  	rad = 0.0; 
  	np = 0;
  	for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
  		for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
  			if (gridded_data[ps].array[i + ii][j + jj] != gridded_data[ps].nodata_value) {
  				rad += gridded_data[ps].array[i + ii][j + jj];
  				np++;
  			}
  		}
		}
		if (np == 0)
			radiance[n] = 0.0;
		else
			radiance[n] = rad/((float) np);
  }

// Estimate deep water radiance. 

  radmin = vec_min(radiance, npoints);
  radmax = vec_max(radiance, npoints);

  if (min_threshold == 0.0)
	  min_threshold = radmin + noise_threshold*(radmax - radmin);
  
  min_radiance_threshold = min_threshold; // Global.

// Estimate land-sea radiance. 
  
  max_slope = 0.0;

  if (max_threshold == 0.0) {
  	for (n = 1; n < npoints - 1; n++) {
   		slope = radiance[n + 1] - radiance[n]; 
			if (fabs(slope) > max_slope) {
				max_slope = fabs(slope);
				if (radiance[n + 1] > radiance[n])
					max_threshold = radiance[n];
				else 
					max_threshold = radiance[n + 1];
			}
  	}
  }

  max_radiance_threshold = max_threshold; // Global.

#if PGPLOT
// Load PGPLOT. 

  error = cpgopen("?");
  if (error < 1) {
  	printf("\nERROR: cannot load graphics library.\n");
  	return ;
  }

// Background and text colours. 

	if (background == BLACK) {
  	cpgscr(0, 0.0, 0.0, 0.0);
  	cpgscr(1, 1.0, 1.0, 1.0);
  } else {
  	cpgscr(0, 1.0, 1.0, 1.0);
  	cpgscr(1, 0.0, 0.0, 0.0);
  }

// Create two (horizontally-alligned) graphics windows. 

  cpgsubp(2,1);

// Create pixel-array plot.

  x[0] = xa;
  x[1] = xb; 
  y[0] = ya;
  y[1] = yb;

  plot_array(gridded_data[ps], parsed[2], pagesize, textsize, pointsize,
  	xmin, ymin, xmax, ymax, palette, 2, x, y, 0.5, line_width, 
    false, gridded_data[ps], 0.0, 0.0);

// Plot land/sea intersection on pixel-array plot. 

  cpgsci(8);

  for (n = 0; n < npoints - 1; n++) {
  	if ((radiance[n] <= max_threshold && radiance[n + 1] >= max_threshold) || 
  		(radiance[n] >= max_threshold && radiance[n + 1] <= max_threshold) ) {
  		cpgpt1(xpts[n], ypts[n], 20);
  	}
  }

// Plot deep water threshold on pixel-array plot. 

  cpgsci(4);

  for (n = 0; n < npoints - 1; n++) {
  	if ((radiance[n] <= min_threshold && radiance[n + 1] >= min_threshold) || 
  		(radiance[n] >= min_threshold && radiance[n + 1] <= min_threshold) ) {
  		cpgpt1(xpts[n], ypts[n], 20);
  	}
  }


// Reset colour. 

  cpgsci(1);

// Create scatter plot. 

  strcpy(title, "RADIANCE MIN/MAX THRESHOLDS = ");
  sprintf(snum, "%5.2f, %5.2f", min_threshold, max_threshold);
  strcat(title, snum);

  plot_scatter(title, "", "", 3, 20, npoints, xpts, radiance, line_width, 0.0, 0.0, 0.0, 0.0, true);

// Plot the radiance thresholds on the scatter plot.

  cpgsci(4);
  cpgsls(2);
  cpgmove(xpts[0], min_threshold);
  cpgdraw(xpts[npoints - 1], min_threshold);
  cpgsls(1);

  cpgsci(8);
  cpgsls(2);
  cpgmove(xpts[0], max_threshold);
  cpgdraw(xpts[npoints - 1], max_threshold);
  cpgsls(1);

// Close the graphics device. 

  cpgclos();
#endif

// Free memory.

  free(xpts);
  free(ypts);
  free(radiance);

}

// PLOT GRID /grid name/ [optional] SOUNDINGS ...
// 		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/ ...
//    LAND /land grid name/ MIN /min/ MAX /max/

void run_plot_grid() {
	
	int i, k, n, ps, ps_land = 0, error;
	bool present = false, plot_soundings = false, plot_land = false;
	float xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, res = 0.0,
		*x, *y, min = 0.0, max = 0.0;

// Read the grid name. 

  if (parsed[2][0] == '\0') {
  	printf("\nERROR: PLOT /grid name/ [optional] REGION /xmin/ /ymin/ /xmax/ /ymax/ \
PALETTE /palette name/ LAND /land grid name/\n");
  	return ;
  } else {

  	// Check grid is present. 

		for (n = 0; n < MAX_GRIDS; n++) {
			if (strcmp(grid_names[n], parsed[2]) == 0) {
				present = true;
				ps = n;
				break;
			}
		}

		if (! present) {
			printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
			return ;
		}
  }

  n = 3;

// Read SOUNDINGS

  if (strcmp(parsed[n], "SOUNDINGS") == 0) {
  	plot_soundings = true;
  	n++;
  	x = (float*) malloc(nsoundings*sizeof(float));
  	y = (float*) malloc(nsoundings*sizeof(float));
  	for (i = 0; i < nsoundings; i++) {
  		x[i] = soundings[i][0];
  		y[i] = soundings[i][1];
  	}
  }

// Read REGION

  if (strcmp(parsed[n], "REGION") == 0) {
  	xmin = atof(parsed[++n]);
  	ymin = atof(parsed[++n]);
  	xmax = atof(parsed[++n]);
  	ymax = atof(parsed[++n]);
  	n++;
  } 

// Read PALETTE 

  palette = JET; // Default.

  if (strcmp(parsed[n], "PALETTE") == 0) {
    if (strcmp(parsed[++n], "jet") == 0)
      palette = JET;
    else if (strcmp(parsed[n], "inverse_jet") == 0)
      palette = INVERSEJET;
    else if (strcmp(parsed[n], "grayscale") == 0)
      palette = GRAYSCALE;
    else if (strcmp(parsed[n], "ocean") == 0)
      palette = OCEAN;
    else if (strcmp(parsed[n], "gmt") == 0)
      palette = GMT;
    else if (strcmp(parsed[n], "rygb") == 0)
      palette = RYGB;
    else if (strcmp(parsed[n], "wor") == 0)
      palette = WOR;
    else if (strcmp(parsed[n], "viridis") == 0)
      palette = VIRIDIS;
    else if (strcmp(parsed[n], "plasma") == 0)
      palette = PLASMA;
    else {
      printf("\nERROR: unknown colour palette '%s'\n", parsed[n]);
    }
    n++;
  }

// Read LAND

  if (strcmp(parsed[n], "LAND") == 0) {
    plot_land = true;

    // Check grid is present. 

    present = false;
    n++;
    for (k = 0; k < MAX_GRIDS; k++) {
      if (strcmp(grid_names[k], parsed[n]) == 0) {
        present = true;
        ps_land = k;
        break;
      }
    }

    if (! present) {
      printf("\nERROR: unknown grid '%s'.\n", parsed[n]);
      return ;
    }
    n++;
  }

// Read MIN and MAX.

  if (strcmp(parsed[n], "MIN") == 0) {
    min = atof(parsed[++n]);
    n++;
  }

  if (strcmp(parsed[n], "MAX") == 0) {
    max = atof(parsed[++n]);
    n++;
  }

#if PGPLOT
// Load PGPLOT. 

  error = cpgopen("?");
  if (error < 1) {
  	printf("\nERROR: cannot load graphics library.\n");
  	return ;
  }

// Background and text colours. 

	if (background == BLACK) {
  		cpgscr(0, 0.0, 0.0, 0.0);
  		cpgscr(1, 1.0, 1.0, 1.0);
  	} else {
  		cpgscr(0, 1.0, 1.0, 1.0);
  		cpgscr(1, 0.0, 0.0, 0.0);
  	}

// Create plot.

  if (plot_soundings) {
  	plot_array(gridded_data[ps], parsed[2], pagesize, textsize, pointsize,
  		xmin, ymin, xmax, ymax, palette, nsoundings, x, y, 1.0, line_width, 
      plot_land, gridded_data[ps_land], min, max);
  	free(x);
  	free(y);
  } else {
  	plot_array(gridded_data[ps], parsed[2], pagesize, textsize, pointsize,
  		xmin, ymin, xmax, ymax, palette, 0, NULL, NULL, 1.0, line_width,
      plot_land, gridded_data[ps_land], min, max);
  }

// Close the graphics device. 

  cpgclos();
#endif

}


//
//    GEODETIC COMPUTATIONS
//

float geodistance(float lon1, float lat1, float lon2, float lat2) {

	float d2r, radius, dlat, dlon, a, c;

	d2r = PI/180.0;
	radius = 6.371e6;
	dlat = lat2 - lat1;
	dlon = lon2 - lon1;
	a = pow(sin(0.5*dlat*d2r), 2) + cos(lat1*d2r)*cos(lat2*d2r)*pow(sin(0.5*dlon*d2r), 2);
	c = 2.0*asin(sqrt(a));
	return c*radius;
}

//
//    INPUT PARSER
//

int parse_input(char line[MAX_LINE]) {

	int n, m;

// Reset parsed. 

	for (n = 0; n < MAX_ARGS; n++)
		for (m = 0; m < MAX_ARG_SIZE; m++)
			parsed[n][m] = '\0';

// Tokenize. 

	tokenize(line);

// Check for high priority commands. 

	if (strncmp(parsed[0], "//", 2) == 0)
		return COMMENT;
	else if (strncmp(parsed[0], "EXIT", 4) == 0)
		return EXIT;
  else if (strncmp(parsed[0], "STOP", 4) == 0)
    return STOP;
	else
		return OK;
}

void tokenize(char line[MAX_LINE]) {

	int n, m = 0, tok = 0; 

	for (n = 0; n < strlen(line); n++) {
//		printf("\n'%c' ", line[n]);
		if (line[n] == ' ' || line[n] == '\0' || line[n] == '\n') {
			parsed[tok][m++] = '\0';
			tok++;
			m = 0;
//			printf("%s\n", parsed[tok - 1]);
		} else {
			parsed[tok][m++] = line[n];
		}
	}
}
