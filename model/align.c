//
//  ALIGN
//

// ALIGN IN /grid 1/ /grid 2/ ... OUT /output grid/ SOUNDINGS /soundings file or grid/ [optional] CHARTDATUM 
//  /chart datum for soundings/ MAXSOUNDINGS /maximum number of soundins/ NEGATEGRID NEGATESOUNDINGS

// Finds c[], a[] and b[] to align the grids with the soundings file. Initially we 
// find a[] and b[], such that each of 
//
//  a[0]*Z[0]^b[0], ..., a[n-1]*Z[n-1]^b[n-1]
//
//  minimises the error between each grid and the soundings file. Then we find 
//  an optimal weighted linear combination of all the grids
//
//  (a[0]*grid[0] + ... + a[0]*grid[0])/(a[0] + a[1] + ... + a[n-1])
//

#include "align.h"

void run_align() {

  int i, j, k, l, n, kk, n_grids, ps_in[MAX_GRIDS], ps_out, ps_grid2, nrows, ncols, 
    n_soundings, max_soundings, r_index, best_index;
  float *x, *y, *z_soundings_float, chart_datum, spval;
  double *z_soundings, error, lowest, z_depth;
  bool present, negate_grid, negate_soundings, nodata;

  srand(time(NULL));

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  } 

// Read input grid names.   

  if (strcmp(parsed[1], "IN") != 0) {
    printf("\nERROR: ALIGN IN /grid 1/ /grid 2/ ... OUT /output grid/ SOUNDINGS \
/soundings file or grid/ [optional] CHARTDATUM /chart datum for soundings/ MAXSOUNDINGS \
/maximum number of soundins/ NEGATEGRID NEGATESOUNDINGS\n");
    return ;
  }

  k = 2;
  n_grids = 0;

  while (true) {

    if (k == MAX_ARGS || n_grids >= MAX_GRIDS || strcmp(parsed[k], "OUT") == 0)
      break ;

    present = false;

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        ps_in[n_grids++] = n;
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

  printf("\nn_grids = %d\n", n_grids);

// Check grids are on the same extents. 

  for (n = 1; n < n_grids; n++) {
    if (gridded_data[ps_in[0]].nrows != gridded_data[ps_in[n]].nrows || 
      gridded_data[ps_in[0]].ncols != gridded_data[ps_in[n]].ncols || 
      fabs(gridded_data[ps_in[0]].slat - gridded_data[ps_in[n]].slat) > extent_eps || 
      fabs(gridded_data[ps_in[0]].nlat - gridded_data[ps_in[n]].nlat) > extent_eps || 
      fabs(gridded_data[ps_in[0]].wlon - gridded_data[ps_in[n]].wlon) > extent_eps || 
      fabs(gridded_data[ps_in[0]].elon - gridded_data[ps_in[n]].elon) > extent_eps || 
      fabs(gridded_data[ps_in[0]].cellsize - gridded_data[ps_in[n]].cellsize) > extent_eps) {
      printf("\nERROR: input grids are not commensurate.\n");
      return ;
    }
  }

// Read output grid name.

  if (strcmp(parsed[k], "OUT") != 0) {
    printf("\nERROR: ALIGN IN /grid 1/ /grid 2/ ... OUT /output grid/ SOUNDINGS \
/soundings file or grid/ [optional] CHARTDATUM /chart datum for soundings/ MAXSOUNDINGS \
/maximum number of soundins/ NEGATEGRID NEGATESOUNDINGS\n");
    return ;
  }

// Check output grid name is not in use. 

  k++;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[k]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[k]);
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

  strcpy(grid_names[ps_out], parsed[k]);

  k++;


// Read soundings files.

  if (strcmp(parsed[k], "SOUNDINGS") != 0) {
    printf("\nERROR: ALIGN IN /grid 1/ /grid 2/ ... OUT /output grid/ SOUNDINGS \
/soundings file or grid/ [optional] CHARTDATUM /chart datum for soundings/ \
SPVAL /no data value/ MAXSOUNDINGS /maximum number of soundins/ NEGATEGRID NEGATESOUNDINGS\n");
    return ;
  } 

  k++;

// Read sounding data. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[k]) == 0) {
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
    z_soundings_float = malloc(nrows*ncols*sizeof(float));

    l = 0; 
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        if (! approx_equal(gridded_data[ps_grid2].array[i][j], gridded_data[ps_grid2].nodata_value, 1.0e-4)) {
          x[l] = gridded_data[ps_grid2].wlon + ((float) j)*gridded_data[ps_grid2].cellsize;
          y[l] = gridded_data[ps_grid2].slat + ((float) i)*gridded_data[ps_grid2].cellsize;
          z_soundings_float[l] = gridded_data[ps_grid2].array[i][j];
          l++;
        }
      }
    }
    n_soundings = l;
    printf("\nn_soundings = %d\n", n_soundings);
  } else if (! file_exists(parsed[k])) {
    printf("\nERROR: File %s does not exist!\n\n", parsed[k]);
    return ;
  } else {

    // Read soundings from CSV file. 

    read_xyz(parsed[k], &x, &y, &z_soundings_float, &n_soundings); 

    // Select soundings within model domain. 

    spval = gridded_data[ps_in[0]].nodata_value;
    kk = 0;
    for (k = 0; k < n_soundings; k++) {

      nodata = false; 

      for (n = 0; n < n_grids; n++) {
        if (x[k] < gridded_data[ps_in[n]].wlon || x[k] > gridded_data[ps_in[n]].elon || 
          y[k] < gridded_data[ps_in[n]].slat || y[k] > gridded_data[ps_in[n]].nlat) {
          nodata = true; 
          break ;
        } else {
          z_depth = (double) interp_bilinear(
                gridded_data[ps_in[n]].array, 
                gridded_data[ps_in[n]].ncols, 
                gridded_data[ps_in[n]].nrows, 
                gridded_data[ps_in[n]].wlon, 
                gridded_data[ps_in[n]].slat, 
                gridded_data[ps_in[n]].cellsize, 
                gridded_data[ps_in[n]].nodata_value, 
                x[k], y[k]); 
          // printf("\n\t\t%.2f\n", z_depth);
          if (approx_equal_double(z_depth, (double) spval, 1.0e-4)) {
            nodata = true; 
            break ;
          }
        }
      }

      if (! nodata) {
        // printf("\n%.2f  %.2f", z_depth, z_soundings_float[k]);
        x[kk] = x[k];
        y[kk] = y[k];
        z_soundings_float[kk] = z_soundings_float[k];
        kk++;
      }
    }

    n_soundings = kk;
  }

  z_soundings = malloc(n_soundings*sizeof(double));
  for (i = 0; i < n_soundings; i++) {
    z_soundings[i] = (double) z_soundings_float[i];
  }

  k++;

  // [optional] Read CHARTDATUM

  if (strcmp(parsed[k], "CHARTDATUM") == 0) {
    chart_datum = atof(parsed[++k]);
    k++; 
  } else {
    chart_datum = 0.0;
  }

  // Read SPVAL. 

  if (strcmp(parsed[k], "SPVAL") == 0) {
    spval = atof(parsed[++k]);
    k++; 
  } else {
    spval = 0.0;
  }

  // Read MAXSOUNDINGS. 

  if (strcmp(parsed[k], "MAXSOUNDINGS") == 0) {
    max_soundings = atof(parsed[++k]);
    k++; 
  } else {
    max_soundings = n_soundings;
  }

  printf("\nmax_soundings = %d\n", max_soundings);

  if (max_soundings != 0) {
    for (n = 0; n < max_soundings; n++) {
      r_index = random_in_range(0, n_soundings);
      x[n] = x[r_index];
      y[n] = y[r_index];
      z_soundings[n] = z_soundings[r_index];
    }

    n_soundings = max_soundings;
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

  // Create grid data at each sounding. 

  double **z_model;

  spval = gridded_data[ps_in[0]].nodata_value;

  allocate_double_array_2d(&z_model, n_grids, n_soundings);

  for (k = 0; k < n_soundings; k++) {
    for (n = 0; n < n_grids; n++) {
      
      if (x[k] < gridded_data[ps_in[n]].wlon || x[k] > gridded_data[ps_in[n]].elon || 
          y[k] < gridded_data[ps_in[n]].slat || y[k] > gridded_data[ps_in[n]].nlat) {
        z_model[n][k] = (double) spval;  
      } else {
      z_model[n][k] = (double) interp_bilinear(
                gridded_data[ps_in[n]].array, 
                gridded_data[ps_in[n]].ncols, 
                gridded_data[ps_in[n]].nrows, 
                gridded_data[ps_in[n]].wlon, 
                gridded_data[ps_in[n]].slat, 
                gridded_data[ps_in[n]].cellsize, 
                gridded_data[ps_in[n]].nodata_value, 
                x[k], y[k]); 
      }

      if (approx_equal_double(z_model[n][k], (double) spval, 1.0e-4)) {
        continue ;
      }

      if (negate_soundings) {
        z_soundings[k] *= -1.0;
      }
      if (negate_grid) {
        z_model[n][k] *= -1.0; 
      }

      // printf("\nmodel/sounding depths = %.2f, %.2f\n", z_model[n][k], z_soundings[k]);
    }
  }

  // Compute weights for soundings. 

  double *sounding_weights; 
  sounding_weights = malloc(n_soundings*sizeof(double)); 
  align_compute_weights(&sounding_weights, z_soundings, n_soundings);

#if 0
  printf("\nsoundings, weights = \n");
  for (i = 0; i < n_soundings; i++) {
    printf("%.2f\t%.4f \n", z_soundings[i], sounding_weights[i]);
  }
  printf("\n");
#endif

  // Simplex algorithm parameters. 

  double reqmin = 1.0e-4; // Terminating limit for the variance of function values.
  int konvge = 250;  // Number of iterations for each convergence check. 
  int kcount = 5000; // Maximum number of function evaluations. 
  int icount = 0;   // Total number of function evaluations. 
  int numres = 0;   // Number of restarts. 
  int ifault = 0;   // 0 - no errors, 1 - bad inputs, 2 - failed to converge. 

  // Optimise each grid individually. 

  double **transformed_params, *start, *min, *step;

  start = malloc(2*sizeof(double));
  min = malloc(2*sizeof(double));
  step = malloc(2*sizeof(double));

  allocate_double_array_2d(&transformed_params, n_grids, 2);

  printf("\nAligning each grid...");

  best_index = 0;
  lowest = 1.0e4;

  for (i = 0; i < n_grids; i++) {

    // error_aligned(float *params, float *z_model, float *z_soundings, int n_soundings, float spval)

    start[0] = 1.0;
    start[1] = 1.0;
    step[0] = 0.5;
    step[1] = 0.5;

    error = error_aligned(start, z_model[i], z_soundings, sounding_weights, n_soundings, (double) spval);
    printf("\ndataset %d: initial error = %.3f\n", i, error);

    nelmin2(&error_aligned, 
          z_model[i], z_soundings, sounding_weights, n_soundings, (double) spval, 
          2, start, min, 
          &error, reqmin, step, konvge, kcount, 
          &icount, &numres, &ifault);

    printf("\nfinal error = %.3f\n", error);
    printf("\na*Z^b: a = %.3f, b = %.3f\n", fabs(min[0]), fabs(min[1]));

    if (error < lowest) {
      lowest = error;
      best_index = i;
    }

    transformed_params[i][0] = fabs(min[0]);
    transformed_params[i][1] = fabs(min[1]);

    // Update modelled depths. 

    for (j = 0; j < n_soundings; j++) {
      z_model[i][j] = transformed_params[i][0]*pow(fabs(z_model[i][j]), transformed_params[i][1]);
    }
  }

  free(start);
  free(min);
  free(step);

  // Compute optimal weighted sum. 

  float *weighted_params;

  weighted_params = malloc(n_grids*sizeof(double));

  start = malloc(n_grids*sizeof(double));
  min = malloc(n_grids*sizeof(double));
  step = malloc(n_grids*sizeof(double));

  int k_restart, n_restarts = 1;
  
  lowest = 1.0e4;

  printf("\nComputed weighted sum of all grids...");

  for (k_restart = 0; k_restart < n_restarts; k_restart++) {

    for (i = 0; i < n_grids; i++) {
      start[i] = 1.0; 
      step[i] = 0.1;
    }

    // double error_weighted_sum(double *params, double **z_model, double *z_soundings, int n_grids, int n_soundings, double spval)

    error = error_weighted_sum(start, z_model, z_soundings, sounding_weights, n_grids, n_soundings, spval);
    printf("\ndataset initial error = %.3f\n", error);

    icount = 0;   // Total number of function evaluations. 
    numres = 0;   // Number of restarts. 
    ifault = 0;   // 0 - no errors, 1 - bad inputs, 2 - failed to converge. 

    nelmin3(&error_weighted_sum, 
          z_model, z_soundings, sounding_weights, n_grids, n_soundings, (double) spval, 
          n_grids, start, min, 
          &error, reqmin, step, konvge, kcount, 
          &icount, &numres, &ifault);

    if (error < lowest) {
      lowest = error;
      for (i = 0; i < n_grids; i++) {
        weighted_params[i] = fabs(min[i]); 
      }
    }

    printf("\npost_optimisation error = %.3f\n", error);
    printf("\nc = ");
    for (i = 0; i < n_grids; i++) {
      printf("%.3f\t", fabs(min[i]));
    }
    printf("\n");
  }

  printf("\nBEST error = %.3f\n", lowest);
  printf("\nc = ");
  for (i = 0; i < n_grids; i++) {
    printf("%.3f\t", weighted_params[i]);
  }
  printf("\n");

  free(start);
  free(min);
  free(step);

  // Create gridded dataset using the optimised dataset. 

  float **array, z, depth, sum;

  nrows = gridded_data[ps_in[0]].nrows;
  ncols = gridded_data[ps_in[0]].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

  allocated_grids[ps_out] = true;
  n_grids_in_use++;

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {

      depth = 0.0; 
      sum = 0.0;

      nodata = false;
      for (k = 0; k < n_grids; k++) {
        z = gridded_data[ps_in[k]].array[i][j];
        if (approx_equal(z, spval, 1.0e-4)) {
          nodata = true;
          break ;
        } else {
          depth += fabs(weighted_params[k])*fabs(transformed_params[k][0])*pow(fabs(z), fabs(transformed_params[k][1]));
          sum += fabs(weighted_params[k]);
        }
      }

      if (nodata) {
        z = gridded_data[ps_in[best_index]].array[i][j];
        if (! approx_equal(z, spval, 1.0e-4)) {
          nodata = false;
          depth = fabs(transformed_params[best_index][0])*pow(fabs(z), fabs(transformed_params[best_index][1]));
          array[i][j] = depth;
          if (z < 0.0) {
            array[i][j] *= -1.0;
          }
        } else {
          array[i][j] = spval;
        }
      } else if (! nodata && sum > epsilon) {
        array[i][j] = depth/sum;
        if (z < 0.0) {
          array[i][j] *= -1.0;
        }
      } else {
        array[i][j] = spval;
      }
    }
  }


  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in[0]].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in[0]].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in[0]].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in[0]].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in[0]].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in[0]].nodata_value;

  // Deallocate memory. 

  free(x);
  free(y);
  free(z_soundings);
  free(z_soundings_float);
  free(weighted_params);
  free(sounding_weights);
  free_double_array_2d(z_model, n_grids);
  free_double_array_2d(transformed_params, n_grids);

}



void align_compute_weights(double **weights, double *soundings, int nsoundings) {

  int n, k; 
  double weight_scale = 1.0, sum, max;

  for (n = 0; n < nsoundings; n++) {
    sum = 0.0; 
    for (k = 0; k < nsoundings; k++) {
      sum += exp(-pow(soundings[k] - soundings[n], 2));
    }
    (*weights)[n] = sum;
  }

  max = vec_max_double(*weights, nsoundings);

  for (n = 0; n < nsoundings; n++) {
    (*weights)[n] = 1.0 - ((*weights)[n])/(weight_scale*max); 
  }

}



// RMS error of the scaled depths a*Z^b
//

double error_aligned(double *params, double *z_model, double *z_soundings, double *sounding_weights, int n_soundings, double spval) {

  int i; 
  double error, n_obs, transformed_depth, soundings_sum;

  error = 0.0;
  n_obs = 0.0;
  soundings_sum = 0.0; 

  for (i = 0; i < n_soundings; i++) {

    if (approx_equal_double(z_soundings[i], spval, 1.0e-4) || approx_equal_double(z_model[i], spval, 1.0e-4)) {
      continue ;
    }

    transformed_depth = fabs(params[0])*pow(fabs(z_model[i]), fabs(params[1]));
    error += sounding_weights[i]*fabs((transformed_depth - fabs(z_soundings[i]))/fabs(z_soundings[i]));
    // soundings_sum += fabs(z_soundings[i]);
    n_obs += sounding_weights[i];
  }

  // soundings_sum /= n_obs;
  // error = 100.0*sqrt(error/n_obs)/soundings_sum;
  error = 100.0*error/n_obs;
  return error;
}


// RMS error of the optimal weighted linear combination. 
//
//  a[0]*Z[0]^b[0], ..., a[n-1]*Z[n-1]^b[n-1]
//

double error_weighted_sum(double *params, double **z_model, double *z_soundings, double *sounding_weights, int n_grids, int n_soundings, double spval) {

  int i, j; 
  double error, weighted_depth, weight_sum, soundings_sum, n_obs;
  bool nodata;

  n_obs = 0.0;
  error = 0.0;
  soundings_sum = 0.0;

  for (i = 0; i < n_soundings; i++) {

    if (approx_equal_double(z_soundings[i], spval, 1.0e-4)) {
      continue ;
    }

    nodata = false;
    for (j = 0; j < n_grids; j++) {
      if (approx_equal_double(z_model[j][i], spval, 0.0)) {
        nodata = true;
        break ;
      }
    }

    if (nodata) {
      continue ;
    }

    weighted_depth = 0.0;
    weight_sum = 0.0;

    for (j = 0; j < n_grids; j++) {
      weighted_depth += fabs(params[j])*fabs(z_model[j][i]);
      weight_sum += fabs(params[j]);
    }

    if (weight_sum > epsilon) {
      weighted_depth /= weight_sum;
      error += sounding_weights[i]*fabs((weighted_depth - fabs(z_soundings[i]))/fabs(z_soundings[i]));
      // error += (weighted_depth - fabs(z_soundings[i]))/fabs(z_soundings[i]);
      // printf("\n model, sounding, error = %.2f, %.2f, %.4f", weighted_depth, z_soundings[i], pow((weighted_depth - fabs(z_soundings[i]))/fabs(z_soundings[i]), 2));
      // soundings_sum += fabs(z_soundings[i]);
      n_obs += sounding_weights[i];
    }
  }

  // soundings_sum /= n_obs;
  // error = 100.0*sqrt(error/n_obs)/soundings_sum;
  error = 100.0*error/n_obs;
  //printf("\n error = %.2f", error);
  return error;
}







