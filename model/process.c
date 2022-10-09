
#include "process.h"

//
//    PROCESS
//

// PROCESS INTERPOLATE IN /input grid/ OUT /output grid/ ... 
//    [optional] REGION /xmin/ /ymin/ /xmax/ /ymax/ RESOLUTION /resolution/

// PROCESS SMOOTH IN /input grid/ OUT /output grid/

// PROCESS SHALLOW IN /input grid name/ OUT /output grid name/ ...
//    THRESHOLD /min radiance level/ /max radiance level/

// PROCESS PANSHARPEN IN /input grid name/ OUT /output grid name/ ...
//    SPECTRAL /spectral band name 1/ /spectral band name 2/ ... 
//    PAN /panchromatic band grid name/

// PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ SPECTRAL ...
//    /spectral band name 1/ /spectral band name 2/ ... NIR ... 
//    /nir band grid name/ LAND /land grid name/

// PROCESS NOISE IN /input grid name/ OUT /output grid name/ ...
//    MAXPOINTS /maximum number of noise points/ FILTER NODATA

// PROCESS UNITE IN /input grid name 1/ /input grid name 2/ ... OUT /output grid name/ ...
//    [optional] REGION /wlon/ /slat/ /elon/ /nlat/ RESOLUTION /resolution/ ...
//    METHOD MIN || MAX || MEAN

// PROCESS OVERLAY IN /input grid name/ OUT /output grid name/ LAYERS /layer 1/ /layer 2/ ...

// PROCESS CLOUD IN /input grid name/ OUT /output grid name/ CLOUD /cloud mask grid name/ ...
//     [optional] RADIUS /integer radius/ 

// PROCESS MASK IN /input grid name/ OUT /output grid name/ MASK /mask grid/

// PROCESS POLYGONS IN /input grid name/ OUT /output grid name/ [optional] POLYGONS /file name 1/ /file name 2/ ... 
//    FILL /fill value 1/ /fill value 2/ ... 

// PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ ERROR /error grid name/ 
//    THRESHOLD /error threshold/

// PROCESS RELATIVEERROR IN /input bathymetry grid/ OUT /output relative error grid/ ERROR /absolute error grid/ 

// PROCESS MEDIAN IN /input grid name/ OUT /output grid name/ [optional] WINDOWSIZE /window size/ NSIGMA /n sigma threshold/

void run_process() {
  
  if (strcmp(parsed[1], "INTERPOLATE") == 0) {
    run_process_interpolate();
  } else if (strcmp(parsed[1], "SMOOTH") == 0) {
    run_process_smooth();
  } else if (strcmp(parsed[1], "SHALLOW") == 0) {
    run_process_shallow();
  } else if (strcmp(parsed[1], "PANSHARPEN") == 0) {
    run_process_pansharpen();
  } else if (strcmp(parsed[1], "NOISE") == 0) {
    run_process_noise();
  } else if (strcmp(parsed[1], "DEGLINT") == 0) {
    run_process_deglint();
  } else if (strcmp(parsed[1], "UNITE") == 0) {
    run_process_unite();
  } else if (strcmp(parsed[1], "CLOUD") == 0) {
    run_process_cloud();
  } else if (strcmp(parsed[1], "MASK") == 0) {
    run_process_mask();
  } else if (strcmp(parsed[1], "LAND") == 0) {
    run_process_land();
  } else if (strcmp(parsed[1], "OVERLAY") == 0) {
    run_process_overlay();
  } else if (strcmp(parsed[1], "POLYGONS") == 0) {
    run_process_polygons();
  } else if (strcmp(parsed[1], "ERROR") == 0) {
    run_process_error();
  } else if (strcmp(parsed[1], "RELATIVEERROR") == 0) {
    run_process_relativeerror();
  } else if (strcmp(parsed[1], "MNDWI") == 0) {
    run_process_mndwi();
  } else if (strcmp(parsed[1], "FILL") == 0) {
    run_process_fill();
  } else if (strcmp(parsed[1], "ADD") == 0) {
    run_process_add();
  } else if (strcmp(parsed[1], "MEDIAN") == 0) {
    run_process_median();
  } else {
    printf("\nERROR: unknown command '%s'\n", parsed[1]);
    return ;
  }
}



//	GRID PROCESSING ROUTINES


#include "process.h"

float land_cutoff(float **array, int nrows, int ncols, float spval, 
    int i_min, int j_min, int i_max, int j_max);

float lyzenga_extinction_cutoff(float **array, int nrows, int ncols, float spval, int n_sigma);


// PROCESS MEDIAN IN /input grid name/ OUT /output grid name/ [optional] WINDOWSIZE /window size/ NSIGMA /n sigma threshold/

void run_process_median() {

  int i, j, ii, jj, len, n, ps_in, ps_out, nrows, ncols, window_size = 1;
  float **array, *vec, n_sigma = 0.0, center, mean, median, stddev;
  bool present = false;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Check input grid.  

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS MEDIAN IN /input grid name/ OUT /output grid name/ [optional] \
WINDOWSIZE /window size/ NSIGMA /n sigma threshold/ \n");
    return ;
  }

// Check input grid is present. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS MEDIAN IN /input grid name/ OUT /output grid name/ [optional] \
WINDOWSIZE /window size/ NSIGMA /n sigma threshold/ \n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Read WINDOWSIZE. 

  n = 6;
  if (strcmp(parsed[n], "WINDOWSIZE") == 0) {
    window_size = atoi(parsed[++n]);
    n++;
  }

// Read NSIGMA. 

  if (strcmp(parsed[n], "NSIGMA") == 0) {
    n_sigma = atof(parsed[++n]);
    n++;
  }

// Apply adaptive median filter on the array. 

  copy_float_array_2d(gridded_data[ps_in].array, array, nrows, ncols);

  vec = malloc(pow(2*window_size + 1, 2)*sizeof(float)); 

  for (i = window_size; i < nrows - window_size - 1; i++) {
    for (j = window_size; j < ncols - window_size - 1; j++) {

      if (approx_equal(gridded_data[ps_in].array[i][j], gridded_data[ps_in].nodata_value, 1.0e-6)) {
        continue ;
      } else {
        center = gridded_data[ps_in].array[i][j];
      }

      len = 0;
      for (ii = -window_size; ii < window_size + 1; ii++) {
        for (jj = -window_size; jj < window_size + 1; jj++) {
          if (ii == 0 && jj == 0) {
            continue ;
          }
          if (! approx_equal(gridded_data[ps_in].array[i + ii][j + jj], gridded_data[ps_in].nodata_value, 1.0e-6)) {
            vec[len++] = gridded_data[ps_in].array[i + ii][j + jj];
          }
        }
      }

      // Sort, median, stddev.

      qsort(vec, len, sizeof(float), compare_floats);

      if (len%2 == 0) {
        median = (vec[(len-1)/2] + vec[len/2])/2.0;
      } else {
        median = vec[len/2];
      }

      mean = vec_mean(vec, len);
      stddev = vec_stddev(vec, len, 0.0, mean);

      // Filter. 

      if (center < median - n_sigma*stddev || center > median + n_sigma*stddev) {
        array[i][j] = median;
      } 
    }
  }


  free(vec);

  gridded_data[ps_out].array = array; 
  gridded_data[ps_out].ncols = gridded_data[ps_in].ncols;
  gridded_data[ps_out].nrows = gridded_data[ps_in].nrows;
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;


}



// PROCESS FILL IN /input grid name/ OUT /output grid name/ DENSITY /required surrounding data density/ 

void run_process_fill() {

  int i, j, ps_in, ps_out, density, nrows, ncols, ndata;
  float **array, spval, tot;
  bool success;

  // Read input grid. 

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS FILL IN /input grid name/ OUT /output grid name/ DENSITY \
/required surrounding data density/ \n");
    return ;
  }

  success = parse_old_grid(&ps_in, 3);
  if (! success) return ;

  // Read in output grid.

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS FILL IN /input grid name/ OUT /output grid name/ DENSITY \
/required surrounding data density/ \n");
    return ;
  }

  success = parse_new_grid(&ps_out, 5);
  if (! success) return ;

  // Read data density. 

  if (strcmp(parsed[6], "DENSITY") != 0) {
    printf("\nERROR: PROCESS FILL IN /input grid name/ OUT /output grid name/ DENSITY \
/required surrounding data density/ \n");
    return ;
  }

  density = atoi(parsed[7]);

  if (density == 0) {
    printf("\nERROR: DENSITY cannot be zero. Setting DENSITY to 1.\n");
    density = 1;
  }

  // Allocate output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Copy input grid to output grid. 

  copy_float_array_2d(gridded_data[ps_in].array, array, nrows, ncols);

// Fill no-data points. 

  spval = gridded_data[ps_in].nodata_value;

  for (i = 1; i < nrows - 1; i++) {
    for (j = 1; j < ncols - 1; j++) {
      if (gridded_data[ps_in].array[i][j] == spval) {
        ndata = 0;
        tot = 0.0;
        if (gridded_data[ps_in].array[i - 1][j - 1] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i - 1][j - 1];
        }
        if (gridded_data[ps_in].array[i - 1][j] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i - 1][j];
        }
        if (gridded_data[ps_in].array[i - 1][j + 1] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i - 1][j + 1];
        }
        if (gridded_data[ps_in].array[i][j - 1] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i][j - 1];
        }
        if (gridded_data[ps_in].array[i][j + 1] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i][j + 1];
        }
        if (gridded_data[ps_in].array[i + 1][j - 1] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i + 1][j - 1];
        }
        if (gridded_data[ps_in].array[i + 1][j] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i + 1][j];
        }
        if (gridded_data[ps_in].array[i + 1][j + 1] != spval) {
          ndata++;
          tot += gridded_data[ps_in].array[i + 1][j + 1];
        }

        if (ndata >= density) {
          array[i][j] = tot/((float) ndata);
        }
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}



// PROCESS ADD IN /input grid 1/ /input grid 2/ OUT /output grid name/ [optional] COEFFICIENTS /c1/ /c2/ /c3/ ABS

// We compute c1*/grid 1/ + c2*/grid 2/ + c3, where by default c1 = c2 = 1.0, c3 = 0.0, and ABS == false. 

void run_process_add() {

  int i, j, k, n, ps_in_1, ps_in_2, ps_out, nrows, ncols;
  float **array, c1, c2, c3;
  bool present, abs_p = false;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS ADD IN /input grid 1/ /input grid 2/ OUT /output grid name/ [optional] COEFFICIENTS /c1/ /c2/ /c3/ ABS \n");
    return ;
  }

// Check grid 1 is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in_1 = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[4]) == 0) {
      present = true;
      ps_in_2 = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[4]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[5], "OUT") != 0) {
    printf("\nERROR: PROCESS ADD IN /input grid 1/ /input grid 2/ OUT /output grid name/ [optional] COEFFICIENTS /c1/ /c2/ /c3/ ABS \n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[6]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[6]);
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

  strcpy(grid_names[ps_out], parsed[6]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// [optional] Read WEIGHTS

  k = 7;

  if (strcmp(parsed[k], "COEFFICIENTS") != 0) {
    c1 = atof(parsed[++k]);
    c2 = atof(parsed[++k]);
    c3 = atof(parsed[++k]);
  } else {
    c1 = 1.0; 
    c2 = 1.0;
    c3 = 0.0;
  }

// [optional] Read ABS. 

  if (strcmp(parsed[k], "ABS") != 0) {
    abs_p = true; 
  }

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

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in_1].nrows;
  ncols = gridded_data[ps_in_1].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute w1*grid1 + w2*grid2

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {

      if (approx_equal(gridded_data[ps_in_1].array[i][j], gridded_data[ps_in_1].nodata_value, 1.0e-4) || 
          approx_equal(gridded_data[ps_in_2].array[i][j], gridded_data[ps_in_2].nodata_value, 1.0e-4)) {
        array[i][j] = gridded_data[ps_in_1].nodata_value;
      } else {
        array[i][j] = c1*gridded_data[ps_in_1].array[i][j] + c2*gridded_data[ps_in_2].array[i][j] + c3;
        if (abs_p) {
          array[i][j] = fabs(array[i][j]);
        }
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in_1].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in_1].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in_1].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in_1].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in_1].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in_1].nodata_value;
}




// PROCESS RELATIVEERROR IN /input bathymetry grid/ OUT /output relative error grid/ ERROR /absolute error grid/ 

void run_process_relativeerror() {

  int i, j, n, ps_in, ps_error, ps_out, nrows, ncols;
  float **array;
  bool present;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS RELATIVEERROR IN /input bathymetry grid/ OUT /output relative error grid/ ERROR /absolute error grid/ \n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ \
ERROR /error grid name/ THRESHOLD /error thereshold/\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Read ERROR grid. 

  if (strcmp(parsed[6], "ERROR") != 0) {
    printf("\nERROR: PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ \
ERROR /error grid name/ THRESHOLD /error thereshold/\n");
    return ;
  }

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[7]) == 0) {
      present = true;
      ps_error = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[7]);
    return ;
  }

// Check grids are commensurate. 

  if (gridded_data[ps_in].nrows != gridded_data[ps_error].nrows || 
    gridded_data[ps_in].ncols != gridded_data[ps_error].ncols || 
    fabs(gridded_data[ps_in].slat - gridded_data[ps_error].slat) > extent_eps || 
    fabs(gridded_data[ps_in].nlat - gridded_data[ps_error].nlat) > extent_eps || 
    fabs(gridded_data[ps_in].wlon - gridded_data[ps_error].wlon) > extent_eps || 
    fabs(gridded_data[ps_in].elon - gridded_data[ps_error].elon) > extent_eps || 
    fabs(gridded_data[ps_in].cellsize - gridded_data[ps_error].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute relative percentage error. 

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {

      if (gridded_data[ps_in].array[i][j] == gridded_data[ps_in].nodata_value || 
          gridded_data[ps_error].array[i][j] == gridded_data[ps_error].nodata_value) {
        array[i][j] = gridded_data[ps_in].nodata_value;
      } else if (fabs(gridded_data[ps_in].array[i][j]) < 0.1) {
        array[i][j] = 0.0;
      } else if (fabs(gridded_data[ps_in].array[i][j]) < 1.0) {
        array[i][j] = 100.0*gridded_data[ps_error].array[i][j]*fabs(gridded_data[ps_in].array[i][j]);
      } else {
        array[i][j] = 100.0*gridded_data[ps_error].array[i][j]/fabs(gridded_data[ps_in].array[i][j]);
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;
}



// PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ ERROR /error grid name/ 
//    THRESHOLD /error threshold/ [optional] REGION /wlon/ /slat/ /elon/ /nlat/ 

void run_process_error() {

  int n, i, j, ps_in, ps_out, ps_error, ncols, nrows;
  bool present, region = false;
  float threshold, **array, wlon, slat, elon, nlat, lon, lat;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ \
ERROR /error grid name/ THRESHOLD /error thereshold/\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ \
ERROR /error grid name/ THRESHOLD /error thereshold/\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Read ERROR grid. 

  if (strcmp(parsed[6], "ERROR") != 0) {
    printf("\nERROR: PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ \
ERROR /error grid name/ THRESHOLD /error thereshold/\n");
    return ;
  }

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[7]) == 0) {
      present = true;
      ps_error = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[7]);
    return ;
  }

// Read THRESHOLD. 

  if (strcmp(parsed[8], "THRESHOLD") != 0) {
    printf("\nERROR: PROCESS ERROR IN /bathymetry grid name/ OUT /output grid name/ \
ERROR /error grid name/ THRESHOLD /error thereshold/\n");
    return ;
  }

  threshold = atof(parsed[9]);

//  printf("\nthreshold = %.2f\n\n", threshold);

// [optional] REGION.

  if (strcmp(parsed[10], "REGION") == 0) {
    region = true;
    wlon = atof(parsed[11]);
    slat = atof(parsed[12]);
    elon = atof(parsed[13]);
    nlat = atof(parsed[14]);
  }  

// Check grids are commensurate. 

  if (gridded_data[ps_in].nrows != gridded_data[ps_error].nrows || 
    gridded_data[ps_in].ncols != gridded_data[ps_error].ncols || 
    fabs(gridded_data[ps_in].slat - gridded_data[ps_error].slat) > extent_eps || 
    fabs(gridded_data[ps_in].nlat - gridded_data[ps_error].nlat) > extent_eps || 
    fabs(gridded_data[ps_in].wlon - gridded_data[ps_error].wlon) > extent_eps || 
    fabs(gridded_data[ps_in].elon - gridded_data[ps_error].elon) > extent_eps || 
    fabs(gridded_data[ps_in].cellsize - gridded_data[ps_error].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Copy input grid to output grid. 

  copy_float_array_2d(gridded_data[ps_in].array, array, nrows, ncols);

// Propagate error estimates on the ERROR grid which are greater than the user-specified threshold 
// to nodata values on the bathymetry (input) grid. 

  int nupdated = 0;

  for (i = 0; i < nrows; i++) {
    lat = gridded_data[ps_in].slat + gridded_data[ps_in].cellsize*i;
    for (j = 0; j < ncols; j++) {
      lon = gridded_data[ps_in].wlon + gridded_data[ps_in].cellsize*j;
      if (gridded_data[ps_error].array[i][j] > threshold) {
        if (region) {
          if (lon > wlon && lon < elon && lat > slat && lat < nlat) {
            array[i][j] = gridded_data[ps_in].nodata_value;
            nupdated++;
          }
        } else {
          array[i][j] = gridded_data[ps_in].nodata_value;
          nupdated++;
        }
      }
    }
  }

  printf("\nUpdated %d grid points.\n\n", nupdated);

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}



// PROCESS POLYGONS IN /input grid name/ OUT /output grid name/ [optional] POLYGONS /file name 1/ /file name 2/ ... 
//    FILL /fill value 1/ /fill value 2/ ... MODE INTERIOR|EXTERIOR RANGE min max

// Default fill value is the nodata value of the input grid. 

// Default mode is INTERIOR. 

void run_process_polygons() {

  int n, i, j, k, ps_in, ps_out, nrows, ncols, n_polygons, n_vertices;
  bool present, interior = true, range_present = false;
  float **array, fill[MAX_ARGS], slat, wlon, nlat, elon, lon, lat, *polygon_X, *polygon_Y, 
    range_min, range_max, spval;
  char polygon_files[MAX_ARGS][MAX_LINE];

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS POLYGONS IN /input grid name/ OUT /output grid name/ [optional] \
POLYGONS /file name 1/ /file name 2/ ... FILL /fill value 1/ /fill value 2/ ... \n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS POLYGONS IN /input grid name/ OUT /output grid name/ [optional] \
POLYGONS /file name 1/ /file name 2/ ... FILL /fill value 1/ /fill value 2/ ... \n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

// Read POLYGONS. 

  k = 6;
  n_polygons = 0;

  if (strcmp(parsed[k], "POLYGONS") == 0) {

    k++;

    while (true) {

      if (parsed[k][0] == '\0' || strcmp(parsed[k], "FILL") == 0 || strcmp(parsed[k], "MODE") == 0 || strcmp(parsed[k], "RANGE") == 0)
        break ;

      // Check polygon file exists. 

      if (! file_exists(parsed[k]) ) {
        printf("\nERROR: missing file '%s'.\n", parsed[k]);
        continue ;
      }

      // Read polygon file. 

      strcpy(polygon_files[n_polygons++], parsed[k]);
      k++;
    }
  }  

// Read FILL. 

  if (strcmp(parsed[k], "FILL") == 0) {

    k++;
    int nfills = 0;

    while ( parsed[k][0] != '\0' && strcmp(parsed[k], "MODE") != 0 && strcmp(parsed[k], "RANGE") != 0 ) {
      fill[nfills++] = atof(parsed[k++]);
    }

    if (nfills != n_polygons) {
      printf("\nERROR: the number of fills should equal the number of polygons.\n");
    }

  } else {
    for (i = 0; i < MAX_ARGS; i++) {
      fill[i] = gridded_data[ps_in].nodata_value;
    }
  }

// Read MODE. 

  if (strcmp(parsed[k], "MODE") == 0) {
    k++;
    if (strcmp(parsed[k], "INTERIOR") == 0) {
      interior = true;
    } else if (strcmp(parsed[k], "EXTERIOR") == 0) {
      interior = false;
    } else {
      printf("\nERROR: MODE should be either INTERIOR or EXTERIOR.\n");
      return ;
    }
    k++;
  }

// Read RANGE. 

  if (strcmp(parsed[k], "RANGE") == 0) {
    range_present = true;
    range_min = atof(parsed[++k]);
    range_max = atof(parsed[++k]);
    printf("\nrange_min, range_max = %f, %f", range_min, range_max);
    k++;
  }

// Copy grid name. 

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
//  printf("\nnrows,ncols = %d, %d\n", nrows, ncols);
  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Copy input grid to output grid. 

  copy_float_array_2d(gridded_data[ps_in].array, array, nrows, ncols);

// Process each polygon. 

  spval = gridded_data[ps_in].nodata_value;

  for (n = 0; n < n_polygons; n++) {

    // Read polygon data. 

    read_xy(polygon_files[n], &polygon_X, &polygon_Y, &n_vertices);

    printf("\npolygon file = %s, number of vertices = %d\n", polygon_files[n], n_vertices);

    // Determine bounding box of vertices. 

    wlon = vec_min(polygon_X, n_vertices);
    elon = vec_max(polygon_X, n_vertices);
    slat = vec_min(polygon_Y, n_vertices);
    nlat = vec_max(polygon_Y, n_vertices);

    printf("\nwlon,elon,slat,nlat = %.3f,%.3f,%.3f,%.3f\n", wlon, elon, slat, nlat);

    // Overlay polgon. 

    int nupdated = 0;

    if (interior) {
      for (i = 0; i < nrows; i++) {
        lat = gridded_data[ps_in].slat + i*gridded_data[ps_in].cellsize;
        for (j = 0; j < ncols; j++) {
          lon = gridded_data[ps_in].wlon + j*gridded_data[ps_in].cellsize;
          if (lon >= wlon && lon <= elon && lat >= slat && lat <= nlat) { // Trivial optimisation. 
            if ( point_in_polygon(polygon_X, polygon_Y, n_vertices, lon, lat) ) {
              if ((range_present && array[i][j] > range_min && array[i][j] < range_max) || ! range_present) {
                array[i][j] = fill[n];
                nupdated++;
              }
            }
          }
        }
      }
    } else {
      // exterior
      for (i = 0; i < nrows; i++) {
        lat = gridded_data[ps_in].slat + i*gridded_data[ps_in].cellsize;
        for (j = 0; j < ncols; j++) {
          lon = gridded_data[ps_in].wlon + j*gridded_data[ps_in].cellsize;
          if (lon >= wlon && lon <= elon && lat >= slat && lat <= nlat) { // Trivial optimisation. 
            if ( ! point_in_polygon(polygon_X, polygon_Y, n_vertices, lon, lat) ) {
              if ((range_present && array[i][j] > range_min && array[i][j] < range_max) || ! range_present) {
                array[i][j] = fill[n];
                nupdated++;
              }
            }
          } else {
            array[i][j] = fill[n];
            nupdated++;            
          }
        }
      }
    }

    printf("\nUpdated %d grid points.\n\n", nupdated);

    free(polygon_X);
    free(polygon_Y);
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}


// PROCESS OVERLAY IN /input grid name/ OUT /output grid name/ LAYERS /layer 1/ /layer 2/ ...

void run_process_overlay() {

  int i, j, k, n, ps_in, ps_out, n_overlays, overlay_indexes[MAX_GRIDS], nrows, ncols;
  bool present;
  float **array, slat, wlon, nlat, elon, cellsize, lon, lat, z;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS OVERLAY IN /input grid name/ OUT /output grid name/ LAYERS /layer 1/ /layer 2/ ...\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }


// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS OVERLAY IN /input grid name/ OUT /output grid name/ LAYERS /layer 1/ /layer 2/ ...\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Read LAYERS

  if (strcmp(parsed[6], "LAYERS") != 0) {
    printf("\nERROR: PROCESS OVERLAY IN /input grid name/ OUT /output grid name/ LAYERS /layer 1/ /layer 2/ ...\n");
    return ;
  }

// Read OVERLAY grids. 

  k = 7; 
  n_overlays = 0; // Global. 

  while (true) {

    if (n_overlays >= MAX_GRIDS || parsed[k][0] == '\0')
      break ;

    present = false;

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        overlay_indexes[n_overlays++] = n;
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

  if (n_overlays == 0) {
    printf("\nERROR: overlay grids incorrectly specified.\n");
    return ; 
  }

// Create output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute the overlays. 

  slat = gridded_data[ps_in].slat; 
  wlon = gridded_data[ps_in].wlon;
  nlat = gridded_data[ps_in].nlat;
  elon = gridded_data[ps_in].elon;
  cellsize = gridded_data[ps_in].cellsize;

  for (i = 0; i < nrows; i++) {
    lat = slat + i*cellsize;
    for (j = 0; j < ncols; j++) {
      lon = wlon + j*cellsize;
      array[i][j] = gridded_data[ps_in].array[i][j];
      for (k = 0; k < n_overlays; k++) {
        if (commensurate_grids(gridded_data[ps_in], gridded_data[overlay_indexes[k]])) {
          z = gridded_data[overlay_indexes[k]].array[i][j];
        } else {
        z = interp_bicubic(gridded_data[overlay_indexes[k]].array, 
              gridded_data[overlay_indexes[k]].ncols,
              gridded_data[overlay_indexes[k]].nrows, 
              gridded_data[overlay_indexes[k]].wlon, 
              gridded_data[overlay_indexes[k]].slat, 
              gridded_data[overlay_indexes[k]].cellsize, 
              gridded_data[overlay_indexes[k]].nodata_value, 
              lon, 
              lat);
        }
        if (! approx_equal(z, gridded_data[overlay_indexes[k]].nodata_value, 1.0e-6)) { 
          array[i][j] = z; // overlay
        }
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = wlon; 
  gridded_data[ps_out].slat = slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}



// PROCESS LAND IN /input grid name/ OUT /output grid name/ [optional] THRESHOLD /radiance threshold/ COORDINATES /xa/ /ya/ /xb/ /yb/

void run_process_land() {

  int i, j, k, n, ps_in, ps_out, nrows, ncols;
  float **array, xa = 0.0, ya = 0.0, xb = 0.0, yb = 0.0, threshold;
  bool present, coords_present = false, threshold_present = false;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS LAND IN /input grid name/ OUT /output grid name/ [optional] \
THRESHOLD /radiance threshold/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS LAND IN /input grid name/ OUT /output grid name/ [optional] \
THRESHOLD /radiance threshold/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// [optional] Read threshold. 

  k = 6;

  if (strcmp(parsed[k], "THRESHOLD") == 0) {
    threshold_present = true;
    threshold = atof(parsed[++k]);
    k++;
  }

// [optional] Read COORDINATES. 

  if (strcmp(parsed[k], "COORDINATES") == 0) {
    coords_present = true;
    xa = atof(parsed[++k]);
    ya = atof(parsed[++k]);
    xb = atof(parsed[++k]);
    yb = atof(parsed[++k]);
    k++;
  } else {
    xa = sand_lon;
    ya = sand_lat;
    xb = deep_lon;
    yb = deep_lat;
  }

// Do we have sufficient information to compute land/sea threshold? 

  if (! threshold_present && coords_present && 
      (deep_lon == 0.0 || deep_lat == 0.0 || sand_lon == 0.0 || sand_lat == 0.0)) {
    printf("\nERROR: THRESHOLD, COORDINATES, and SAND/DEEP location not present.\n");
    allocated_grids[ps_out] = false;
    grid_names[ps_out][0] = '\0';
    return ;
  }

// Compute threshold. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;

  if (coords_present || !threshold_present) {

    int j_min = ((float) ncols)*(xa - gridded_data[ps_in].wlon)/(gridded_data[ps_in].elon - gridded_data[ps_in].wlon);
    int i_min = ((float) nrows)*(ya - gridded_data[ps_in].slat)/(gridded_data[ps_in].nlat - gridded_data[ps_in].slat);

    int j_max = ((float) ncols)*(xb - gridded_data[ps_in].wlon)/(gridded_data[ps_in].elon - gridded_data[ps_in].wlon);
    int i_max = ((float) nrows)*(yb - gridded_data[ps_in].slat)/(gridded_data[ps_in].nlat - gridded_data[ps_in].slat);

//    printf("\n%.5f,%.5f,%.5f,%.5f\n", xa, ya, xb, yb);
//    printf("\n%d,%d,%d,%d\n", j_min, i_min, j_max, i_max);

    if (j_min < 0 || j_min >= ncols || i_min < 0 || i_min >= nrows || 
        j_max < 0 || j_max >= ncols || i_max < 0 || i_max >= nrows) {
      printf("\nERROR: coordinates outside grid limits.\n");
      allocated_grids[ps_out] = false;
      grid_names[ps_out][0] = '\0';
      return ;
    }

    threshold = land_cutoff(
      gridded_data[ps_in].array, 
      gridded_data[ps_in].nrows, 
      gridded_data[ps_in].ncols, 
      gridded_data[ps_in].nodata_value, 
      i_min, j_min, i_max, j_max);
  }

  printf("\nland/sea threshold = %.2f\n\n", threshold);

// Allocate memory for the output grid. 

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute threshold.

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (gridded_data[ps_in].array[i][j] == gridded_data[ps_in].nodata_value || 
        gridded_data[ps_in].array[i][j] < threshold) {
        array[i][j] = gridded_data[ps_in].nodata_value;
      } else {
        array[i][j] = gridded_data[ps_in].array[i][j];
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;
}




// PROCESS MNDWI IN /green spectral grid name/ /nir spectral grid name/ OUT /output grid name/ ... 
//     [optional] THRESHOLD /radiance threshold/ COORDINATES /xa/ /ya/ /xb/ /yb/

void run_process_mndwi() {

  int i, j, k, n, ps_spec, ps_nir, ps_out, nrows, ncols;
  float **array, xa = 0.0, ya = 0.0, xb = 0.0, yb = 0.0, threshold, spec, nir;
  bool present, coords_present = false, threshold_present = false;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS MNDWI IN /green spectral grid name/ /nir spectral grid name/ OUT /output grid name/ [optional] \
THRESHOLD /radiance threshold/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_spec = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[4]) == 0) {
      present = true;
      ps_nir = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[4]);
    return ;
  }

// Check uniform grid extents. 

  if (gridded_data[ps_spec].nrows != gridded_data[ps_nir].nrows || 
    gridded_data[ps_spec].ncols != gridded_data[ps_nir].ncols || 
    fabs(gridded_data[ps_spec].slat - gridded_data[ps_nir].slat) > extent_eps || 
    fabs(gridded_data[ps_spec].nlat - gridded_data[ps_nir].nlat) > extent_eps || 
    fabs(gridded_data[ps_spec].wlon - gridded_data[ps_nir].wlon) > extent_eps || 
    fabs(gridded_data[ps_spec].elon - gridded_data[ps_nir].elon) > extent_eps || 
    fabs(gridded_data[ps_spec].cellsize - gridded_data[ps_nir].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[5], "OUT") != 0) {
    printf("\nERROR: PROCESS MNDWI IN /green spectral grid name/ /nir spectral grid name/ OUT /output grid name/ [optional] \
THRESHOLD /radiance threshold/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[6]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[6]);
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

  strcpy(grid_names[ps_out], parsed[6]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// [optional] Read threshold. 

  k = 7;

  if (strcmp(parsed[k], "THRESHOLD") == 0) {
    threshold_present = true;
    threshold = atof(parsed[++k]);
    k++;
  }

// [optional] Read COORDINATES. 

  if (strcmp(parsed[k], "COORDINATES") == 0) {
    coords_present = true;
    xa = atof(parsed[++k]);
    ya = atof(parsed[++k]);
    xb = atof(parsed[++k]);
    yb = atof(parsed[++k]);
    k++;
  } else {
    xa = sand_lon;
    ya = sand_lat;
    xb = deep_lon;
    yb = deep_lat;
  }

// Do we have sufficient information to compute land/sea threshold? 

  if (! threshold_present && coords_present && 
      (deep_lon == 0.0 || deep_lat == 0.0 || sand_lon == 0.0 || sand_lat == 0.0)) {
    printf("\nERROR: THRESHOLD, COORDINATES, and SAND/DEEP location not present.\n");
    allocated_grids[ps_out] = false;
    grid_names[ps_out][0] = '\0';
    return ;
  }

// Compute threshold. 

  nrows = gridded_data[ps_spec].nrows;
  ncols = gridded_data[ps_spec].ncols;

  if (coords_present || !threshold_present) {

    int j_min = ((float) ncols)*(xa - gridded_data[ps_spec].wlon)/(gridded_data[ps_spec].elon - gridded_data[ps_spec].wlon);
    int i_min = ((float) nrows)*(ya - gridded_data[ps_spec].slat)/(gridded_data[ps_spec].nlat - gridded_data[ps_spec].slat);

    int j_max = ((float) ncols)*(xb - gridded_data[ps_spec].wlon)/(gridded_data[ps_spec].elon - gridded_data[ps_spec].wlon);
    int i_max = ((float) nrows)*(yb - gridded_data[ps_spec].slat)/(gridded_data[ps_spec].nlat - gridded_data[ps_spec].slat);

//    printf("\n%.5f,%.5f,%.5f,%.5f\n", xa, ya, xb, yb);
//    printf("\n%d,%d,%d,%d\n", j_min, i_min, j_max, i_max);

    if (j_min < 0 || j_min >= ncols || i_min < 0 || i_min >= nrows || 
        j_max < 0 || j_max >= ncols || i_max < 0 || i_max >= nrows) {
      printf("\nERROR: coordinates outside grid limits.\n");
      allocated_grids[ps_out] = false;
      grid_names[ps_out][0] = '\0';
      return ;
    }

    threshold = land_cutoff(
      gridded_data[ps_spec].array, 
      gridded_data[ps_spec].nrows, 
      gridded_data[ps_spec].ncols, 
      gridded_data[ps_spec].nodata_value, 
      i_min, j_min, i_max, j_max);
  }

  printf("\nland/sea threshold = %.2f\n\n", threshold);

// Allocate memory for the output grid. 

  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute threshold.

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (gridded_data[ps_spec].array[i][j] == gridded_data[ps_spec].nodata_value) {
        array[i][j] = gridded_data[ps_spec].nodata_value;
      } else {
        spec = gridded_data[ps_spec].array[i][j];
        nir  = gridded_data[ps_nir].array[i][j];

        array[i][j] = (spec - nir)/(spec + nir); // MNDWI

        if (array[i][j] > threshold) {
          array[i][j] = gridded_data[ps_spec].nodata_value;
        }
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_spec].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_spec].slat;
  gridded_data[ps_out].elon = gridded_data[ps_spec].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_spec].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_spec].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_spec].nodata_value;
}


// PROCESS MASK IN /input grid name/ OUT /output grid name/ MASK /mask grid/ [optional] REVERSE

void run_process_mask() {

// Basically this routine just copies the input grid to the output grid and propogates 
// no data values on the mask grid to the output grid. 

  int n, i, j, ii, jj, iii, jjj, nrows, ncols, ps_in, ps_out, ps_mask, radius = 1; 
  float **array;
  bool present, reverse, no_data;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS MASK IN /input grid name/ OUT /output grid name/ MASK /mask grid/ [optional] REVERSE\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS MASK IN /input grid name/ OUT /output grid name/ MASK /mask grid/ [optional] REVERSE\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Read mask grid name. 

  if (strcmp(parsed[6], "MASK") != 0) {
    printf("\nERROR: PROCESS MASK IN /input grid name/ OUT /output grid name/ MASK /mask grid/ [optional] REVERSE\n");
    return ;
  }

// Check mask grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[7]) == 0) {
      present = true;
      ps_mask = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[7]);
    return ;
  }

// Check uniform grid extents. 

  if (gridded_data[ps_in].nrows != gridded_data[ps_mask].nrows || 
    gridded_data[ps_in].ncols != gridded_data[ps_mask].ncols || 
    fabs(gridded_data[ps_in].slat - gridded_data[ps_mask].slat) > extent_eps || 
    fabs(gridded_data[ps_in].nlat - gridded_data[ps_mask].nlat) > extent_eps || 
    fabs(gridded_data[ps_in].wlon - gridded_data[ps_mask].wlon) > extent_eps || 
    fabs(gridded_data[ps_in].elon - gridded_data[ps_mask].elon) > extent_eps || 
    fabs(gridded_data[ps_in].cellsize - gridded_data[ps_mask].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// [optional] REVERSE (Mask in the opposite direction.) 

  int k = 8;

  if (strcmp(parsed[k], "REVERSE") == 0) {
    reverse = true;
    k++;
  } else {
    reverse = false;
  }

// [optional] RADIUS. 

  if (strcmp(parsed[k], "RADIUS") == 0) {
    radius = atoi(parsed[++k]);
    // printf("\nRADIUS = %d\n", radius);
    k++;
  }

// Allocate memory for the output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Mask. 

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {

      no_data = false; 

      for (ii = 1 - radius; ii < radius; ii++) {

        if (i + ii < 0)
          iii = 0;
        else if (i + ii > nrows - 1)
          iii = nrows - 1;
        else 
          iii = i + ii;

        for (jj = 1 - radius; jj < radius; jj++) {

          if (j + jj < 0) 
            jjj = 0;
          else if (j + jj > ncols - 1)
            jjj = ncols - 1;
          else 
            jjj = j + jj;

          if (approx_equal(gridded_data[ps_mask].array[iii][jjj], gridded_data[ps_mask].nodata_value, 1.0e-4)) {
            no_data = true;
          }
        }
      }

      if (reverse) {
        if (! no_data) {
          array[i][j] = gridded_data[ps_in].nodata_value;
        } else {
          array[i][j] = gridded_data[ps_in].array[i][j];
        }
      } else {
        if (! no_data) {
          array[i][j] = gridded_data[ps_in].array[i][j];
        } else {
          array[i][j] = gridded_data[ps_in].nodata_value;
        }
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

 }


// PROCESS CLOUD IN /input grid name/ OUT /output grid name/ CLOUD /cloud mask grid name/ ...
//     [optional] RADIUS /integer radius/ 

void run_process_cloud() {

  int i, j, n, nrows, ncols, ps_in, ps_out, ps_cloud, radius = 0;
  float **array;
  bool present;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Read input grid name.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS CLOUD IN /input grid name/ OUT /output grid name/ \
CLOUD /cloud mask grid name/ [optional] RADIUS /integer radius/ \n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS CLOUD IN /input grid name/ OUT /output grid name/ \
CLOUD /cloud mask grid name/ [optional] RADIUS /integer radius/ \n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// Read cloud mask grid name. 

  if (strcmp(parsed[6], "CLOUD") != 0) {
    printf("\nERROR: PROCESS CLOUD IN /input grid name/ OUT /output grid name/ \
CLOUD /cloud mask grid name/ [optional] RADIUS /integer radius/ \n");
    return ;
  }

// Check cloud mask grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[7]) == 0) {
      present = true;
      ps_cloud = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[7]);
    return ;
  }

// [optional] Read radius. 

  if (strcmp(parsed[8], "RADIUS") == 0) {
    radius = atoi(parsed[9]);
  }

// Check uniform grid extents. 

  if (gridded_data[ps_in].nrows != gridded_data[ps_cloud].nrows || 
    gridded_data[ps_in].ncols != gridded_data[ps_cloud].ncols || 
    fabs(gridded_data[ps_in].slat - gridded_data[ps_cloud].slat) > extent_eps || 
    fabs(gridded_data[ps_in].nlat - gridded_data[ps_cloud].nlat) > extent_eps || 
    fabs(gridded_data[ps_in].wlon - gridded_data[ps_cloud].wlon) > extent_eps || 
    fabs(gridded_data[ps_in].elon - gridded_data[ps_cloud].elon) > extent_eps || 
    fabs(gridded_data[ps_in].cellsize - gridded_data[ps_cloud].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Allocate memory for the output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Mask out the clouds. 

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      
      if (gridded_data[ps_in].array[i][j] == gridded_data[ps_in].nodata_value || 
        gridded_data[ps_cloud].array[i][j] == gridded_data[ps_cloud].nodata_value) {
        array[i][j] = gridded_data[ps_in].nodata_value;
      } else {
        array[i][j] = gridded_data[ps_in].array[i][j];
      }
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;
}


// PROCESS UNITE IN /input grid name 1/ /input grid name 2/ ... OUT /output grid name/ ...
//    [optional] REGION /wlon/ /slat/ /elon/ /nlat/ RESOLUTION /resolution/ ...
//    METHOD MIN || MAX || MEAN || EXP || MEDIAN || MEDIAN_NO_SPVAL COMBINE AND||OR WEIGHTS /weight for grid 1/ ... 
//    /weight for grid 2/ ... 

void run_process_unite() {
  int i, j, k, n, ps_in[MAX_GRIDS], n_in, ps_out, ncols, nrows, nobs,
    method, method_min = 1, method_mean = 2, method_max = 3, method_exp = 4, 
    method_median = 5, method_median_no_spval = 6;
  float **array, *scene_list, res, wlon = 360.0, elon = 0.0, slat = 90.0, nlat = -90.0, 
    min, max, mean, lon, lat, z, z0, spval, weights[MAX_GRIDS], ws, w, stddev, avg, median;
  bool present, region = false, and, first, no_data;

// Can we allocate another grid?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }  

// Read input grid names.   

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS UNITE IN /input grid name 1/ /input grid name 2/ ... \
OUT /output grid name/ [optional] REGION /wlon/ /slat/ /elon/ /nlat/ \
RESOLUTION /resolution/ METHOD MIN || MAX || MEAN || EXP || MEDIAN || MEDIAN_NO_SPVAL \
COMBINE AND||OR WEIGHTS /weight for grid 1/ ... /weight for grid 2/\n");
    return ;
  }

  k = 3;
  n_in = 0;

  while (true) {

    if (k == MAX_ARGS || n_in >= MAX_GRIDS || strcmp(parsed[k], "OUT") == 0)
      break ;

    present = false;

    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        ps_in[n_in++] = n;
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

  // printf("\nn_in = %d\n", n_in);

// Read output grid name.

  if (strcmp(parsed[k], "OUT") != 0) {
    printf("\nERROR: PROCESS UNITE IN /input grid name 1/ /input grid name 2/ ... \
OUT /output grid name/ [optional] REGION /wlon/ /slat/ /elon/ /nlat/ \
RESOLUTION /resolution/ METHOD MIN || MAX || MEAN || EXP || MEDIAN || MEDIAN_NO_SPVAL \
COMBINE AND||OR WEIGHTS /weight for grid 1/ ... /weight for grid 2/\n");
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

// [optional] REGION.

  if (strcmp(parsed[k], "REGION") == 0) {
    region = true;
    wlon = atof(parsed[++k]);
    slat = atof(parsed[++k]);
    elon = atof(parsed[++k]);
    nlat = atof(parsed[++k]);
    k++;
  }  

// [optional] RESOLUTION.

  if (strcmp(parsed[k], "RESOLUTION") == 0) {
    res = atof(parsed[++k]);
    k++;
  } else {
    res = gridded_data[ps_in[0]].cellsize;
  }

  // printf("\nresolution = %f\n", res);

// [optional] METHOD. 

  if (strcmp(parsed[k], "METHOD") == 0) {
    k++;
    if (strcmp(parsed[k], "MIN") == 0) {
      method = method_min;
    } else if (strcmp(parsed[k], "MEAN") == 0) {
      method = method_mean;
    } else if (strcmp(parsed[k], "MAX") == 0) {
      method = method_max;
    } else if (strcmp(parsed[k], "EXP") == 0) {
      method = method_exp;
    } else if (strcmp(parsed[k], "MEDIAN") == 0) {
      method = method_median;
    } else if (strcmp(parsed[k], "MEDIAN_NO_SPVAL") == 0) {
      method = method_median_no_spval;
    } else {
      method = method_mean;
    }
    k++;
  } else {
    method = method_mean;
  }

// [optional] COMBINE. AND or OR

  and = false; // OR by default.
  first = false;

  if (strcmp(parsed[k], "COMBINE") == 0) {
    k++;
    if (strcmp(parsed[k], "AND") == 0) {
      and = true;
    }

    if (strcmp(parsed[k], "FIRST") == 0) {
      first = true;
    }

    k++;
  }

// [optional] WEIGHTS.

  if (strcmp(parsed[k], "WEIGHTS") == 0) {
    for (n = 0; n < n_in; n++) {
      weights[n] = atof(parsed[++k]);
    }
    k++;
  } else {
    for (n = 0; n < n_in; n++) {
      weights[n] = 1.0;
    }
  }

// Determine extent of output grid.

  if (! region) {
    for (n = 0; n < n_in; n++) {
      if (gridded_data[ps_in[n]].wlon < wlon)
        wlon = gridded_data[ps_in[n]].wlon;
      if (gridded_data[ps_in[n]].elon > elon)
        elon = gridded_data[ps_in[n]].elon;
      if (gridded_data[ps_in[n]].slat < slat)
        slat = gridded_data[ps_in[n]].slat;
      if (gridded_data[ps_in[n]].nlat > nlat)
        nlat = gridded_data[ps_in[n]].nlat;
    }
  }

  // printf("\nwlon = %f, slat = %f, elon = %f, nlat = %f\n", wlon, slat, elon, nlat);

// Allocate memory for output grid. 

  ncols = 1 + round((elon - wlon)/res);
  nrows = 1 + round((nlat - slat)/res);

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// "Unite" the grids. 

  spval = gridded_data[ps_in[0]].nodata_value; // TODO: this should be documented. 

  scene_list = malloc(n_in*sizeof(float)); 

  for (i = 0; i < nrows; i++) {
    lat = slat + i*res;
    for (j = 0; j < ncols; j++) {
      lon = wlon + j*res;
      min = BIG;
      max = -BIG;
      mean = 0.0;
      ws = 0.0; 
      nobs = 0;
      no_data = false;
      for (k = 0; k < n_in; k++) {
        if (lon > gridded_data[ps_in[k]].wlon && 
            lon < gridded_data[ps_in[k]].elon && 
            lat > gridded_data[ps_in[k]].slat && 
            lat < gridded_data[ps_in[k]].nlat) {

          z = interp_bilinear(gridded_data[ps_in[k]].array, 
              gridded_data[ps_in[k]].ncols,
              gridded_data[ps_in[k]].nrows, 
              gridded_data[ps_in[k]].wlon, 
              gridded_data[ps_in[k]].slat, 
              gridded_data[ps_in[k]].cellsize, 
              gridded_data[ps_in[k]].nodata_value, 
              lon, 
              lat);

        } else {
          z = gridded_data[ps_in[k]].nodata_value;
        }

#if 0
        if (z < -depthmax || z > -depthmin) {
          continue;
        }
#endif
        if (first && k == 0 && approx_equal(z, gridded_data[ps_in[k]].nodata_value, 1.0e-6)) {
          no_data = true;
        }

        if (and && approx_equal(z, gridded_data[ps_in[k]].nodata_value, 1.0e-6)) {
          no_data = true;
        }

        if (! approx_equal(z, gridded_data[ps_in[k]].nodata_value, 1.0e-6) || method_median_no_spval) {

          if (! approx_equal(z, gridded_data[ps_in[k]].nodata_value, 1.0e-6)) {
            nobs++;
          }

          if (method == method_min && z < min) {
            min = z;
          }

          if (method == method_max && z > max) {
            max = z;
          }

          if (method == method_mean) {
            mean += weights[k]*z;
            ws += weights[k];
          }

          if (method == method_median || method == method_median_no_spval) {
            scene_list[nobs - 1] = z;
          }

          if (method == method_exp) {
            w = weights[k]*exp(fabs(z)/2.0);
            mean += w*z;
            ws += w;
          }
        }
      }

      if (method == method_mean || method == method_exp) {
        if (nobs == 0 || no_data) {
          array[i][j] = spval;
        } else { 
          array[i][j] = mean/ws;
        }
      } else if (method == method_median || method == method_median_no_spval) {
        if (((nobs == 0 || no_data) && method != method_median_no_spval) || (method == method_median_no_spval && ((float) nobs) < (0.667*(float) n_in))) {
          array[i][j] = spval;
        } else { 
          qsort(scene_list, nobs, sizeof(float), compare_floats);

          avg = vec_mean(scene_list, nobs);
          stddev = vec_stddev(scene_list, nobs, 0.0, avg);
          median = scene_list[(nobs-1)/2];
          while (nobs > 0 && scene_list[nobs-1] > median + 1.645*stddev) {
            nobs--;
          }

          if (method == method_median) {
            array[i][j] = 0.5*(scene_list[nobs/2] + scene_list[(nobs-1)/2]);
          } else {
            array[i][j] = scene_list[(nobs-1)/2];
          }
        }   
      } else if (method == method_min) {
        if (nobs == 0 || no_data) {
          array[i][j] = spval;
        } else { 
          array[i][j] = min;
        }
      } else {
        if (nobs == 0 || no_data) {
          array[i][j] = spval;
        } else {
          array[i][j] = max;
        }
      }
    }
  }

  free(scene_list);

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = wlon; 
  gridded_data[ps_out].slat = slat;
  gridded_data[ps_out].elon = elon;
  gridded_data[ps_out].nlat = nlat;
  gridded_data[ps_out].cellsize = res;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in[0]].nodata_value;

}



// PROCESS NOISE IN /input grid name/ OUT /output grid name/ ...
//    MAXPOINTS /maximum number of noise points/ FILTER NODATA

void run_process_noise() {

  int n, ps_in, ps_out, maxpoints = 1, nrows, ncols;
  bool filter_nodata = false, present;
  float **array;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Check input grid.  

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS NOISE IN /input grid name/ OUT /output grid name/ \
MAXPOINTS /maximum number of noise points/ FILTER NODATA\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS NOISE IN /input grid name/ OUT /output grid name/ \
MAXPOINTS /maximum number of noise points/ FILTER NODATA\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// Read MAXPOINTS. 

  if (strcmp(parsed[6], "MAXPOINTS") != 0) {
    printf("\nERROR: PROCESS NOISE IN /input grid name/ OUT /output grid name/ \
MAXPOINTS /maximum number of noise points/ FILTER NODATA\n");
    return ;
  }

  maxpoints = atoi(parsed[7]);

  if (strcmp(parsed[8], "FILTER") == 0 && strcmp(parsed[9], "NODATA") == 0) {
    filter_nodata = true;
  }

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;

// Allocate memory for the filtered grid. 

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Copy input array to output array. 

  copy_float_array_2d(gridded_data[ps_in].array, array, nrows, ncols);

// Filter output array. 

  fill_lakes(array, ncols, nrows, maxpoints, true, 
    filter_nodata, gridded_data[ps_in].nodata_value);

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}



// PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ NIR ... 
//    /nir band grid name/ REGION /wlon/ /slat/ /elon/ /nlat/ [optional] SCALE /scale/

void run_process_deglint() {

  int i, j, k, n, ps_in, ps_out, ps_nir, nrows, ncols;
  float **array, rho_ij, rho_jj, nir_mean, spectral_mean, lon, lat, nobs, 
    xi, xj, r_ij, deglinted, scale = 1.0;
  bool present;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }  

// Check input grid.  

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ \
NIR /nir band grid name/ REGION /wlon/ /slat/ /elon/ /nlat/ [optional] SCALE /scale/\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ \
NIR /nir band grid name/ REGION /wlon/ /slat/ /elon/ /nlat/ [optional] SCALE /scale/\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// Read NIR grid. 

  if (strcmp(parsed[6], "NIR") != 0) {
    printf("\nERROR: PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ \
NIR /nir band grid name/ REGION /wlon/ /slat/ /elon/ /nlat/ [optional] SCALE /scale/\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[7]) == 0) {
      present = true;
      ps_nir = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[7]);
    return ;
  }

// Read REGION. 

  if (strcmp(parsed[8], "REGION") != 0) {
    printf("\nERROR: PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ \
NIR /nir band grid name/ REGION /wlon/ /slat/ /elon/ /nlat/ [optional] SCALE /scale/\n");
    return ;
  }

  float wlon = atof(parsed[9]);
  float slat = atof(parsed[10]);
  float elon = atof(parsed[11]);
  float nlat = atof(parsed[12]);

  if (slat >= nlat || wlon >= elon) {
    printf("\nERROR: PROCESS DEGLINT IN /input grid name/ OUT /output grid name/ \
NIR /nir band grid name/ REGION /wlon/ /slat/ /elon/ /nlat/ [optional] SCALE /scale/\n");
    return ;
  }

// Read SCALE. 

  if (strcmp(parsed[13], "SCALE") == 0) {
    scale = atof(parsed[14]);
  }

// Check uniform grid extents. 

  if (gridded_data[ps_in].nrows != gridded_data[ps_nir].nrows || 
    gridded_data[ps_in].ncols != gridded_data[ps_nir].ncols || 
    fabs(gridded_data[ps_in].slat - gridded_data[ps_nir].slat) > extent_eps || 
    fabs(gridded_data[ps_in].nlat - gridded_data[ps_nir].nlat) > extent_eps || 
    fabs(gridded_data[ps_in].wlon - gridded_data[ps_nir].wlon) > extent_eps || 
    fabs(gridded_data[ps_in].elon - gridded_data[ps_nir].elon) > extent_eps || 
    fabs(gridded_data[ps_in].cellsize - gridded_data[ps_nir].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Allocate memory for the output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute NIR and spectral means.

  nir_mean = 0.0;
  spectral_mean = 0.0;
  nobs = 0.0;

  for (i = 0; i < nrows; i++) {
    lat = gridded_data[ps_in].slat + i*gridded_data[ps_in].cellsize;
    for (j = 0; j < ncols; j++) {
      lon = gridded_data[ps_in].wlon + j*gridded_data[ps_in].cellsize;
      if (lon > wlon && lon < elon && lat > slat && lat < nlat && 
          gridded_data[ps_in].array[i][j] != gridded_data[ps_in].nodata_value && 
          gridded_data[ps_nir].array[i][j] != gridded_data[ps_nir].nodata_value) {
        spectral_mean += gridded_data[ps_in].array[i][j];
        nir_mean += gridded_data[ps_nir].array[i][j];
        nobs++;
      }
    }
  }

  nir_mean /= nobs;
  spectral_mean /= nobs;

  printf("\nSamples       = %d", (int) nobs);
  printf("\nNIR mean      = %f", nir_mean);
  printf("\nSpectral mean = %f", spectral_mean);

// Compute covariances. 

  rho_ij = 0.0;
  rho_jj = 0.0;

  for (i = 0; i < nrows; i++) {
    lat = gridded_data[ps_in].slat + i*gridded_data[ps_in].cellsize;
    for (j = 0; j < ncols; j++) {
      lon = gridded_data[ps_in].wlon + j*gridded_data[ps_in].cellsize;
      if (lon > wlon && lon < elon && lat > slat && lat < nlat  && 
          gridded_data[ps_in].array[i][j] != gridded_data[ps_in].nodata_value && 
          gridded_data[ps_nir].array[i][j] != gridded_data[ps_nir].nodata_value) {
        xi = gridded_data[ps_in].array[i][j];
        xj = gridded_data[ps_nir].array[i][j];
        rho_ij += (xi - spectral_mean)*(xj - nir_mean); // covariance
        rho_jj += pow(xj - nir_mean, 2);                // variance
      }
    }
  }

  rho_ij /= nobs;
  rho_jj /= nobs;
  r_ij = rho_ij/rho_jj;

  printf("\nrho_ij, rho_jj, r_ij = %f, %f, %f\n\n", rho_ij, rho_jj, r_ij);

// Deglint.

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      
      // Handle nodata value. 

      if (gridded_data[ps_in].array[i][j] == gridded_data[ps_in].nodata_value) {
        array[i][j] = gridded_data[ps_in].nodata_value;
        continue ;
      }

      // Compute deglinting correction.
#if 0
      if (j == ncols/2) {
        printf("\n %f, %f, %f", 
          gridded_data[ps_in].array[i][j], 
          gridded_data[ps_nir].array[i][j], 
          r_ij*(gridded_data[ps_nir].array[i][j] - nir_mean));
      }
#endif
      array[i][j] = gridded_data[ps_in].array[i][j] - scale*r_ij*(gridded_data[ps_nir].array[i][j] - nir_mean);
    }
  }

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}




// PROCESS PANSHARPEN IN /input grid name/ OUT /output grid name/ ...
//    SPECTRAL /spectral band name 1/ /spectral band name 2/ ...
//    PAN /panchromatic band grid name/

void run_process_pansharpen() {

  int i, j, k, ii, jj, n, m, ps_in, ps_out, ps_pan, nrows, ncols, nspec;
  float **array, res, lon, lat, smin, smax, pmin, pmax, sharp, pan, spec, spectot;
  bool present;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Check input grid.  

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS PANSHARPEN IN /input grid/ OUT /output grid/ \
SPECTRAL /spectral band name 1/ /spectral band name 2/ ... \
PAN /panchromatic band grid name/\n");
    return ;
  }

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS PANSHARPEN IN /input grid/ OUT /output grid/ \
SPECTRAL /spectral band name 1/ /spectral band name 2/ ... \
PAN /panchromatic band grid name/\n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// Read the spectral band grid names. 

  if (strcmp(parsed[6], "SPECTRAL") != 0) {
    printf("\nERROR: PROCESS PANSHARPEN IN /input grid/ OUT /output grid/ \
SPECTRAL /spectral band name 1/ /spectral band name 2/ ... \
PAN /panchromatic band grid name/\n");
    return ;    
  }

  nspec = 0;
  k = 7;

  while (true) {

    if (k == MAX_ARGS || nspec >= MAX_GRIDS || strcmp(parsed[k], "PAN") == 0)
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

// Read the panchromatic band grid name. 

  if (strcmp(parsed[k], "PAN") != 0) {
    printf("\nERROR: PROCESS PANSHARPEN IN /input grid/ OUT /output grid/ \
SPECTRAL /spectral band name 1/ /spectral band name 2/ ... \
PAN /panchromatic band grid name/\n");
    return ;
  }

  k++;

// Check grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[k]) == 0) {
      present = true;
      ps_pan = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
    return ;
  }

// Allocate memory for the pan-sharpened grid. 

  res = gridded_data[ps_pan].cellsize;
  ncols = 1 + round((gridded_data[ps_in].elon - gridded_data[ps_in].wlon)/res);
  nrows = 1 + round((gridded_data[ps_in].nlat - gridded_data[ps_in].slat)/res);

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Set grid parameters.

  gridded_data[ps_out].array = array; 
  gridded_data[ps_out].ncols = ncols; 
  gridded_data[ps_out].nrows = nrows; 
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon; 
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = res;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

// Compute min/max of input grid. 

  smin = array_min2(
    gridded_data[ps_in].array, 
    gridded_data[ps_in].nrows, 
    gridded_data[ps_in].ncols, 
    gridded_data[ps_in].nodata_value);

  smax = array_max2(
    gridded_data[ps_in].array, 
    gridded_data[ps_in].nrows, 
    gridded_data[ps_in].ncols, 
    gridded_data[ps_in].nodata_value);

  // printf("\nspectral min,max = %f, %f", smin,smax);

// Interpolate to pan resolution. 

  grid2gridinterp(
    // Input grid. 
    gridded_data[ps_in].array, 
    gridded_data[ps_in].ncols, 
    gridded_data[ps_in].nrows, 
    (double) gridded_data[ps_in].wlon, 
    (double) gridded_data[ps_in].slat, 
    (double) gridded_data[ps_in].cellsize,
    (double) gridded_data[ps_in].nodata_value,
    // Output grid. 
    array, 
    gridded_data[ps_out].ncols, 
    gridded_data[ps_out].nrows, 
    (double) gridded_data[ps_out].wlon, 
    (double) gridded_data[ps_out].slat, 
    (double) gridded_data[ps_out].cellsize, 
    (double) gridded_data[ps_out].nodata_value);

// Pan-sharpen the grid. 

  for (i = 0; i < nrows; i++) {
    // Latitude on the output grid. 
    lat = gridded_data[ps_out].slat + i*gridded_data[ps_out].cellsize;

    for (j = 0; j < ncols; j++) {

      if (approx_equal(array[i][j], gridded_data[ps_in].nodata_value, 1.0e-4)) {
        continue ;
      }

      // Longitude on the output grid. 
      lon = gridded_data[ps_out].wlon + j*gridded_data[ps_out].cellsize;

      pan = interp_bicubic(gridded_data[ps_pan].array, 
                    gridded_data[ps_pan].ncols,
                    gridded_data[ps_pan].nrows, 
                    gridded_data[ps_pan].wlon, 
                    gridded_data[ps_pan].slat, 
                    gridded_data[ps_pan].cellsize, 
                    gridded_data[ps_pan].nodata_value, 
                    lon, 
                    lat);

      if (approx_equal(pan, gridded_data[ps_pan].nodata_value, 1.0e-4)) {
        continue ;
      }

      sharp = nspec*fabs(array[i][j])*fabs(pan);

      if (pan_normalise) {
        spectot = 0.0;
        for (k = 0; k < nspec; k++) {
          spec = interp_bicubic(gridded_data[spectral_indexes[k]].array, 
                    gridded_data[spectral_indexes[k]].ncols,
                    gridded_data[spectral_indexes[k]].nrows, 
                    gridded_data[spectral_indexes[k]].wlon, 
                    gridded_data[spectral_indexes[k]].slat, 
                    gridded_data[spectral_indexes[k]].cellsize, 
                    gridded_data[spectral_indexes[k]].nodata_value, 
                    lon, 
                    lat);
          if (! approx_equal(spec, gridded_data[spectral_indexes[k]].nodata_value, 1.0e-4)) {
            spectot += fabs(spec);
          }
        }

        array[i][j] = sharp/spectot;
#if 0        
        if (array[i][j] > smax) 
          array[i][j] = smax;
#endif
      } else {
        array[i][j] = sharp/((float) nspec);
      }

    }
  }

  pmin = array_min2(
    array, 
    gridded_data[ps_out].nrows, 
    gridded_data[ps_out].ncols, 
    gridded_data[ps_out].nodata_value);

  pmax = array_max2(
    array, 
    gridded_data[ps_out].nrows, 
    gridded_data[ps_out].ncols, 
    gridded_data[ps_out].nodata_value);


  // printf("\npan min,max = %f, %f\n", pmin,pmax);

  // Rescale min radiance levels to spectral levels. 

#if 0
  if (smin < 0.0) smin = 0.0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (! approx_equal(array[i][j], gridded_data[ps_in].nodata_value, 1.0e-4)) {
        array[i][j] = array[i][j] + smin - pmin;
      }

      if (array[i][j] > smax) {
        array[i][j] = smax;
      }
    }
  }
#endif

#if 0

// Rescale grid to original radiance levels. 

  float alpha = (pmax*smin - pmin*smax)/(pmax - pmin);
  float beta = (smax - smin)/(pmax - pmin);

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (array[i][j] != gridded_data[ps_in].nodata_value)
        array[i][j] = alpha + beta*array[i][j];
    }
  }
#endif
}





// PROCESS INTERPOLATE IN /input grid/ OUT /output grid/ ... 
//		[optional] REGION /xmin/ /ymin/ /xmax/ /ymax/ RESOLUTION /resolution/

void run_process_interpolate() {

  int k, n, ps_in, ps_out, ps_res, ncols, nrows;
  bool present, grid_region, gridres, ps_region;
  float xmin, xmax, ymin, ymax, res, **array;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
  	printf("\nERROR: unable to store grid.\n");
  	return ;
  }

// Check input grid.  

  if (strcmp(parsed[2], "IN") != 0) {
  	printf("\nERROR: PROCESS INTERPOLATE IN /input grid/ OUT /output grid/ \
[optional] REGION /xmin/ /ymin/ /xmax/ /ymax/ RESOLUTION /resolution/\n");
  	return ;
  }

// Check grid is present. 

  present = false;

	for (n = 0; n < MAX_GRIDS; n++) {
	 	if (strcmp(grid_names[n], parsed[3]) == 0) {
	   	present = true;
	   	ps_in = n;
	   	break;
	 	}
	}

	if (! present) {
	  printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
  	return ;
	}

  xmin = gridded_data[ps_in].wlon;
  xmax = gridded_data[ps_in].elon;
  ymin = gridded_data[ps_in].slat;
  ymax = gridded_data[ps_in].nlat;
  res =  gridded_data[ps_in].cellsize;

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS INTERPOLATE IN /input grid/ OUT /output grid/ \
[optional] REGION /xmin/ /ymin/ /xmax/ /ymax/ RESOLUTION /resolution/\n");
    return ;
  }
  	
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
  	if (strcmp(grid_names[n], parsed[5]) == 0) {
  		printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// Read region. 

  n = 6;
  if (strcmp(parsed[n], "REGION") == 0) {

    // Check for a grid name. 
    n++;
    grid_region = false;
    for (k = 0; k < MAX_GRIDS; k++) {
      if (strcmp(grid_names[k], parsed[n]) == 0) {
        grid_region = true;
        ps_region = k;
        break;
      }
    }

    if (grid_region) {
      xmin = gridded_data[ps_region].wlon;
      xmax = gridded_data[ps_region].elon;
      ymin = gridded_data[ps_region].slat;
      ymax = gridded_data[ps_region].nlat;
    } else {
  	  xmin = atof(parsed[n]);
  	  ymin = atof(parsed[++n]);
  	  xmax = atof(parsed[++n]);
  	  ymax = atof(parsed[++n]);
    }
  	n++;
  }

// TODO: Sense check the region. 

// Read resolution. 

  if (strcmp(parsed[n], "RESOLUTION") == 0) {
  	// Check for a grid name.
  	n++;
    gridres = false;
		for (k = 0; k < MAX_GRIDS; k++) {
	 		if (strcmp(grid_names[k], parsed[n]) == 0) {
	   		gridres = true;
	   		ps_res = k;
	   		break;
	 		}
		}

  	if (gridres) {
  		res = gridded_data[ps_res].cellsize;
  	} else {
  		res = atof(parsed[n]);
    }
  } else {
    res = gridded_data[ps_in].cellsize;
  }

  // printf("\nxmin,ymin,xmax,ymax,res = %f,%f,%f,%f,%f", xmin,ymin,xmax,ymax,res);

// Allocate memory for the interpolated grid. 

	ncols = 1 + round((xmax - xmin)/res);
	nrows = 1 + round((ymax - ymin)/res);

  // printf("\nncols, nrows = %d, %d\n\n", ncols, nrows);

  allocate_float_array_2d(&array, nrows, ncols);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Interpolate. 

  grid2gridinterp(
  	// Input grid. 
  	gridded_data[ps_in].array, gridded_data[ps_in].ncols, 
  	gridded_data[ps_in].nrows, gridded_data[ps_in].wlon, 
  	gridded_data[ps_in].slat, gridded_data[ps_in].cellsize,
  	gridded_data[ps_in].nodata_value,
  	// Output grid. 
  	array, ncols, nrows, xmin, ymin, res, gridded_data[ps_in].nodata_value);

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].wlon = xmin;
  gridded_data[ps_out].elon = xmax;
  gridded_data[ps_out].slat = ymin;
  gridded_data[ps_out].nlat = ymax;
  gridded_data[ps_out].cellsize = res;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}


// PROCESS SMOOTH IN /input grid name/ OUT /output grid name/ [optional]
//    NPASSES /number of passes/ WEIGHT /weight/ 

void run_process_smooth() {

  int n, ps_in, ps_out, nrows, ncols, npasses = 1;
  float **array, central_weight = 1.0;
  bool present = false;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
    printf("\nERROR: unable to store grid.\n");
    return ;
  }

// Check input grid.  

  if (strcmp(parsed[2], "IN") != 0) {
    printf("\nERROR: PROCESS SMOOTH IN /input grid name/ OUT /output grid name/ [optional] \
NPASSES /number of passes/ WEIGHT /weight/ \n");
    return ;
  }

// Check input grid is present. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps_in = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
    return ;
  }

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
    printf("\nERROR: PROCESS SMOOTH IN /input grid name/ OUT /output grid name/ [optional] \
NPASSES /number of passes/ WEIGHT /weight/ \n");
    return ;
  }
    
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[5]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Read NPASSES. 

  n = 6;
  if (strcmp(parsed[n], "NPASSES") == 0) {
    npasses = atoi(parsed[++n]);
    n++;
  }

  if (strcmp(parsed[n], "WEIGHT") == 0) {
    central_weight = atoi(parsed[++n]);
    n++;
  }

// Smooth the array. 

  smooth(
    gridded_data[ps_in].array,
    array,
    ncols,
    nrows,
    central_weight,
    npasses,
    gridded_data[ps_in].nodata_value);

  gridded_data[ps_out].array = array; 
  gridded_data[ps_out].ncols = gridded_data[ps_in].ncols;
  gridded_data[ps_out].nrows = gridded_data[ps_in].nrows;
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}


// PROCESS SHALLOW IN /input grid name/ OUT /output grid name/ ...
//		[optional] THRESHOLD /min radiance level/ /max radiance level/

void run_process_shallow() {

	int i, j, k, n, n_grids, ps, ps_in, ps_grids[MAX_GRIDS], ps_out, nrows, ncols, n_sigma;
	float **array;
	bool present = false;

  double gr, ggr, nobs, deep_radiance_total, mean_deep_radiance, stddev_a, stddev_b, stddev;

  n_sigma = 3;

// Can we allocate another band?

  if (n_grids_in_use == MAX_GRIDS) {
  	printf("\nERROR: unable to store grid.\n");
  	return ;
  }

// Check input grid.  

	if (strcmp(parsed[2], "IN") != 0) {
		printf("\nERROR: PROCESS SHALLOW IN /input grid name/ OUT /output grid name/ \
THRESHOLD /min radiance level/ /max radiance level/\n");
		return ;
	}

// Check input grid is present. 

	for (n = 0; n < MAX_GRIDS; n++) {
	 	if (strcmp(grid_names[n], parsed[3]) == 0) {
	   	present = true;
	   	ps_in = n;
	   	break;
	 	}
	}

	if (! present) {
	  printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
  	return ;
	}

// Read output grid name. 

  if (strcmp(parsed[4], "OUT") != 0) {
		printf("\nERROR: PROCESS SHALLOW IN /input grid name/ OUT /output grid name/ \
THRESHOLD /min radiance level/ /max radiance level/\n");
    return ;
  }
  	
// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
  	if (strcmp(grid_names[n], parsed[5]) == 0) {
  		printf("\nERROR: grid '%s' already in use.\n", parsed[5]);
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

  strcpy(grid_names[ps_out], parsed[5]);

// This grid is now in use. 

  allocated_grids[ps_out] = true;

// Read min/max radiance threshold.

  if (strcmp(parsed[6], "THRESHOLD") == 0) {
    if (strcmp(parsed[7], "MIN") == 0) {
      min_radiance_threshold = array_min2(
        gridded_data[ps_in].array,
        gridded_data[ps_in].nrows,
        gridded_data[ps_in].ncols,
        gridded_data[ps_in].nodata_value);
    } else if (strcmp(parsed[7], "AUTOMATIC") == 0) {

      min_radiance_threshold = lyzenga_extinction_cutoff(
        gridded_data[ps_in].array,
        gridded_data[ps_in].nrows,
        gridded_data[ps_in].ncols,
        gridded_data[ps_in].nodata_value, 
        n_sigma);

    } else {
      min_radiance_threshold = atof(parsed[7]); // Global.       
    }

    if (strcmp(parsed[8], "MAX") == 0 || strcmp(parsed[8], "AUTOMATIC") == 0) {
      max_radiance_threshold = array_max2(
        gridded_data[ps_in].array,
        gridded_data[ps_in].nrows,
        gridded_data[ps_in].ncols,
        gridded_data[ps_in].nodata_value);
    } else {
      max_radiance_threshold = atof(parsed[8]); // Global.       
    }
  }

// Read additional grids for which to compute optically deep water signal. 

  k = 10;
  n_grids = 0;

  if (strcmp(parsed[9], "GRIDS") == 0) {

    while (true) {

      if (k == MAX_ARGS || n_grids >= MAX_GRIDS || parsed[k][0] == '\0') {
        break ;
      }

      present = false;

      for (n = 0; n < MAX_GRIDS; n++) {
        if (strcmp(grid_names[n], parsed[k]) == 0) {
          ps_grids[n_grids++] = n;
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
  }

// Allocate memory for output grid. 

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;
  allocate_float_array_2d(&array, nrows, ncols);

// Update memory usage (in units of megabytes). 

  meminuse += 1.0e-6*((float) nrows*ncols)*((float) sizeof(float));

// Update number of bands in use. 

  n_grids_in_use++;

// Compute the shallow water grid, mean deep water signal and standard deviation. 

  printf("\nmin_radiance_threshold = %.9f", min_radiance_threshold);
  printf("\nmax_radiance_threshold = %.9f\n", max_radiance_threshold);

  deep_radiance_total = 0.0;
  nobs = 0.0;

  for (i = 0; i < nrows; i++) {
  	for (j = 0; j < ncols; j++) {
  		gr = gridded_data[ps_in].array[i][j];
  		if (gr >= min_radiance_threshold && gr <= max_radiance_threshold){
  			array[i][j] = gr;
  		} else if (gr > 0.0 && ! approx_equal(gr, gridded_data[ps_in].nodata_value, 1.0e-4)) {
        deep_radiance_total += gr;
        nobs += 1.0;
  			array[i][j] = gridded_data[ps_in].nodata_value;
  		} else {
        array[i][j] = gridded_data[ps_in].nodata_value;
      }
  	}
  }

  mean_deep_radiance = deep_radiance_total/nobs;

  nobs = 0.0;
  stddev_a = 0.0;
  stddev_b = 0.0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      gr = gridded_data[ps_in].array[i][j];
      if (gr >= min_radiance_threshold && gr <= max_radiance_threshold){
        continue ;
      } else if (gr > 0.0 && ! approx_equal(gr, gridded_data[ps_in].nodata_value, 1.0e-4)) {
        stddev_a += gr - mean_deep_radiance; 
        stddev_b += pow(gr - mean_deep_radiance, 2);
        nobs += 1.0;
      }
    }
  }

  stddev = sqrt((stddev_b - pow(stddev_a, 2)/nobs)/nobs);

  Lsm[ps_in] = mean_deep_radiance;
  Lsm_sigma[ps_in] = stddev;

  printf("\nMean deep water radiance of %s = %f +/- %f\n", grid_names[ps_in], mean_deep_radiance, stddev);

// Compute Lsm and Lsm_sigma in any additional grids. 

  if (n_grids > 0) {

    for (n = 0; n < n_grids; n++) {
      ps = ps_grids[n];
      nobs = 0.0;
      deep_radiance_total = 0.0;

      for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
          gr = gridded_data[ps_in].array[i][j];
          if (gr >= min_radiance_threshold && gr <= max_radiance_threshold) {
            continue ;
          } else if (gr > 0.0 && ! approx_equal(gr, gridded_data[ps_in].nodata_value, 1.0e-4)) {
            deep_radiance_total += gridded_data[ps].array[i][j];
            nobs += 1.0;
          }
        }
      }

      mean_deep_radiance = deep_radiance_total/nobs;

      stddev_a = 0.0;
      stddev_b = 0.0;
      nobs = 0.0;

      for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
          gr = gridded_data[ps_in].array[i][j];
          if (gr >= min_radiance_threshold && gr <= max_radiance_threshold) {
            continue ;
          } else if (gr > 0.0 && ! approx_equal(gr, gridded_data[ps_in].nodata_value, 1.0e-4)) {
            ggr = gridded_data[ps].array[i][j];
            stddev_a += ggr - mean_deep_radiance; 
            stddev_b += pow(ggr - mean_deep_radiance, 2);
            nobs += 1.0;
          }
        }
      }

      stddev = sqrt((stddev_b - pow(stddev_a, 2)/nobs)/nobs);
      Lsm[ps] = mean_deep_radiance;
      Lsm_sigma[ps] = stddev;
      printf("Mean deep water radiance of %s = %f +/- %f\n", grid_names[ps], mean_deep_radiance, stddev);
    }
  }

  printf("\n\n");

  gridded_data[ps_out].array = array;
  gridded_data[ps_out].nrows = nrows;
  gridded_data[ps_out].ncols = ncols;
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;

}


float lyzenga_extinction_cutoff(float **array, int nrows, int ncols, float spval, int n_sigma) {

  int i, j, ii, jj, k = 0, indx, len, nbelow, partial_indx;
  float min, max, *vec, threshold, pc, partial_mean, partial_stddev, base_stddev, partial_threshold;

  double deep_water_mean = 0.0, stddev, stddev_a = 0.0, stddev_b = 0.0, nobs = 0.0;

  float percentile = 0.1; // 10%

  // Compute the nth-percentile radiance. 

  len = nrows*ncols;

  vec = malloc(len*sizeof(float));

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (array[i][j] > 0.0 && ! approx_equal(array[i][j], spval, 1.0e-6)) {
        vec[k++] = array[i][j];
      }
    }
  }

  qsort(vec, k, sizeof(float), compare_floats);

  indx = 0.01*((float) k);
  partial_mean = vec_mean(vec, indx);
  base_stddev = vec_stddev(vec, indx, 0.0, partial_mean);

  for (pc = 0.01; pc < 1.01; pc += 0.01) {
    partial_indx = pc*((float) k);
    partial_threshold = vec[partial_indx];
    partial_mean = vec_mean(vec, partial_indx);
    partial_stddev = vec_stddev(vec, partial_indx, 0.0, partial_mean);
    printf("\n%.0fth-percentile radiance = %.4f, stddev = %.4f\n", 100.0*pc, partial_threshold, partial_stddev);

#if 0
    if (partial_stddev > 3.0*base_stddev) {
      percentile = pc;
      threshold = partial_threshold;
      break ;
    }
#endif
  }

  indx = percentile*((float) k);
  threshold = vec[indx];

  free(vec);

  printf("\n%.0fth-percentile radiance = %.9f\n", 100.0*percentile, threshold);

  int window_size = 6; 

  for (i = 0; i < nrows - window_size; i++) {
    for (j = 0; j < ncols - window_size; j++) {

      nbelow = 0;

      for (ii = i; ii < i + window_size; ii++) {
        for (jj = j; jj < j + window_size; jj++) {
          if (array[ii][jj] > 0.0 && ! approx_equal(array[ii][jj], spval, 1.0e-6) && array[ii][jj] <= threshold)
            nbelow++;
        }
      }

      // Deep water region? 

      if (nbelow > window_size*window_size/2) {

        for (ii = i; ii < i + window_size; ii++) {
          for (jj = j; jj < j + window_size; jj++) {
            if (array[ii][jj] > 0.0 && ! approx_equal(array[ii][jj], spval, 1.0e-6)) {
              deep_water_mean += array[ii][jj];
              nobs += 1.0;
            }
          }
        }
      }
    }
  }

  deep_water_mean /= nobs;
  nobs = 0.0;

  for (i = 0; i < nrows - window_size; i++) {
    for (j = 0; j < ncols - window_size; j++) {

      nbelow = 0;

      for (ii = i; ii < i + window_size; ii++) {
        for (jj = j; jj < j + window_size; jj++) {
          if (array[ii][jj] > 0.0 && ! approx_equal(array[ii][jj], spval, 1.0e-6) && array[ii][jj] <= threshold)
            nbelow++;
        }
      }

      // Deep water region? 

      if (nbelow > window_size*window_size/2) {

        for (ii = i; ii < i + window_size; ii++) {
          for (jj = j; jj < j + window_size; jj++) {
            if (array[ii][jj] > 0.0 && ! approx_equal(array[ii][jj], spval, 1.0e-6)) {
              stddev_a += array[ii][jj] - deep_water_mean; 
              stddev_b += pow(array[ii][jj] - deep_water_mean, 2);
              nobs += 1.0;
            }
          }
        }
      }
    }
  }

  stddev = sqrt((stddev_b - pow(stddev_a, 2)/nobs)/nobs);

  printf("\ndeep_water_mean = %.9f\n", deep_water_mean);
  printf("\nstddev          = %.9f\n", stddev);  

  return threshold;
}



float land_cutoff(float **array, int nrows, int ncols, float spval, 
    int i_min, int j_min, int i_max, int j_max) {

  int i, j, k, nx, ny, npoints;
  float alpha, radA, radB, radC, slope, cutoff = 0.0, maxrad = 0.0;

// Straight line from min radiance to max radiance. 

  nx = abs(j_max - j_min);
  ny = abs(i_max - i_min);
  npoints = MAX(nx, ny);

  radA = array[i_min][j_min];
  radB = radA;

  for (k = 1; k < npoints; k++) {
    // Grid coordinates on the line. 
    alpha = ((float) k)/((float) (npoints - 1));
    i = ((float) i_min) + alpha*((float) i_max - i_min);
    j = ((float) j_min) + alpha*((float) j_max - j_min);
    // Difference between successive points. 
    radC = array[i][j];
    if (radA < 0.0 || radB < 0.0 || radC < 0.0) 
      continue ;
    slope = fabs(radA - 2.0*radB + radC);
    // Max slope. 
    if (slope > cutoff) {
      cutoff = slope;
      maxrad = radB;
    }
    //
    radA = radB; 
    radB = radC;
  }

  return maxrad;
}


