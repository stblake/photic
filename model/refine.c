
# include "refine.h"

//
//    REFINE
//

// REFINE GRID IN /input grid name/ OUT /output grid name/ LAND /land grid name/ ... 
// SHALLOW /shallow water grid name/ [optional] SCRAP /scrap min/ /scrap max/ ...
// CLIP /clipped min/ /clipped max/ SCALE /rescaled min/ /rescaled max/ POWER /power scale/

void run_refine() {

	int i, j, k, n, ps_in, ps_out, nrows, ncols, shallow_index, land_index;
	bool present = false, shallow_present = false, land_present = false, 
    clip = false, rescale = false, linear = false, scrap = false, power = false;
	float **array, v, depth, dmin, dmax, oldmin, oldmax, scale = 1.0, 
		smin, smax, sca, scb, linear_m, linear_c, alpha, beta, scrapmin, 
    scrapmax, power_a, power_b;

#if 0
	if (zdepth == NULL) {
		printf("\nERROR: Uninitialised MODEL.\n");
		return ;
	}
#endif

	if (strcmp(parsed[1], "IN") != 0) {
    printf("\nERROR: REFINE IN /input grid name/ OUT /output grid name/ [optional] LAND \
/land grid name/ SHALLOW /shallow water grid name/ SCRAP /scrap min/ /scrap max/\
CLIP /clipped min/ /clipped max/ SCALE /rescaled min/ /rescaled max/ POWER /exponent/\n");
		return ;
	}

// Check input grid name is present. 

	for (n = 0; n < MAX_GRIDS; n++) {
  	if (strcmp(grid_names[n], parsed[2]) == 0) {
  		ps_in = n;
  		present = true;
  		break ;
  	}
  }

	if (! present) {
	  printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
	  return ;
	}

  if (strcmp(parsed[3], "OUT") != 0) {
    printf("\nERROR: REFINE IN /input grid name/ OUT /output grid name/ [optional] LAND \
/land grid name/ SHALLOW /shallow water grid name/ SCRAP /scrap min/ /scrap max/\
CLIP /clipped min/ /clipped max/ SCALE /rescaled min/ /rescaled max/ POWER /exponent/\n");
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
    if (strcmp(grid_names[n], parsed[4]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[4]);
      return ;
    }
  }

  strcpy(grid_names[ps_out], parsed[4]);
  allocated_grids[ps_out] = true;

  nrows = gridded_data[ps_in].nrows;
  ncols = gridded_data[ps_in].ncols;

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
  gridded_data[ps_out].wlon = gridded_data[ps_in].wlon;
  gridded_data[ps_out].elon = gridded_data[ps_in].elon;
  gridded_data[ps_out].slat = gridded_data[ps_in].slat;
  gridded_data[ps_out].nlat = gridded_data[ps_in].nlat;
  gridded_data[ps_out].cellsize = gridded_data[ps_in].cellsize;
  gridded_data[ps_out].nodata_value = gridded_data[ps_in].nodata_value;


// Read LAND grid name.

  k = 5;

	if (strcmp(parsed[k], "LAND") == 0) {

    // Check grid name is present. 
    k++;
    for (n = 0; n < MAX_GRIDS; n++) {
      if (strcmp(grid_names[n], parsed[k]) == 0) {
        land_index = n;
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

// Read SHALLOW grid name.

  if (strcmp(parsed[k], "SHALLOW") == 0) {

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

// Read SCRAP. 

  if (strcmp(parsed[k], "SCRAP") == 0) {
    scrap = true;
    scrapmin = atof(parsed[++k]);
    scrapmax = atof(parsed[++k]);
    k++;
//    printf("\nscrapmin,scrapmax = %f, %f\n", scrapmin, scrapmax);
  }

// Read CLIP. 

	if (strcmp(parsed[k], "CLIP") == 0) {
		clip = true;
		oldmin = atof(parsed[++k]);
		oldmax = atof(parsed[++k]);
    k++;
	}

// Read SCALE. 

	if (strcmp(parsed[k], "SCALE") == 0) {
		rescale = true;
		dmin = atof(parsed[++k]);
		dmax = atof(parsed[++k]);
    k++;
    // printf("\nscalemin, scalemax = %.3f, %.3f\n", dmin, dmax);
	}

// Read SHAPE. 

	if (strcmp(parsed[k], "SHAPE") == 0) {
		scale = atof(parsed[++k]);
    k++;
	}

// Read LINEAR.

  if (strcmp(parsed[k], "LINEAR") == 0) {
    // printf("\nlinear = true\n");
    linear = true;
    linear_m = atof(parsed[++k]);
    linear_c = atof(parsed[++k]);
    // printf("\nm,c = %f,%f\n", linear_m, linear_c);
    k++;
  }

// Read POWER

  if (strcmp(parsed[k], "POWER") == 0) {
    power = true; 
    power_a = atof(parsed[++k]);
    power_b = atof(parsed[++k]);
    printf("\npower_a,power_b = %f, %f\n\n", power_a, power_b);
    k++;
  }

	if (! clip) {
		oldmin = array_min2(gridded_data[ps_in].array, 
			gridded_data[ps_in].nrows,
			gridded_data[ps_in].ncols,
			gridded_data[ps_in].nodata_value);

		oldmax = array_max2(gridded_data[ps_in].array, 
			gridded_data[ps_in].nrows,
			gridded_data[ps_in].ncols,
			gridded_data[ps_in].nodata_value);
	}

  // printf("\ninput min,max = %10.2f, %10.2f\n", oldmin, oldmax);

// Refine model. 

	if (rescale) {
		smin = pow(fabs(oldmin), scale);
		if (oldmin < 0.0) smin *= -1.0;
		smax = pow(fabs(oldmax), scale);
		if (oldmax < 0.0) smax *= -1.0;
		sca = (oldmax - oldmin)/(smax - smin);
		scb = (oldmax*smin - oldmin*smax)/(smin - smax);
	}

  for (i = 0; i < nrows; i++) {
  	for (j = 0; j < ncols; j++) {
  		if ((land_present == false && shallow_present == false) || 
        ((land_present == true && gridded_data[land_index].array[i][j] != gridded_data[land_index].nodata_value) && 
        (shallow_present == true && gridded_data[shallow_index].array[i][j] != gridded_data[shallow_index].nodata_value))) {

        depth = gridded_data[ps_in].array[i][j];
  			
        if (depth == gridded_data[ps_in].nodata_value) {
          array[i][j] = gridded_data[ps_in].nodata_value;
          continue;
        }

        if (clip) {
          if (depth < oldmin) 
            depth = oldmin;
          else if (depth > oldmax)
            depth = oldmax;
        }

        if (linear) {
          depth = linear_m*depth + linear_c;
        }

  			if (rescale) {
          if (scale != 1.0) {
            v = (dmin*oldmax - dmax*oldmin + dmax*depth - dmin*depth)/(oldmax - oldmin);
            depth = sca*pow(fabs(v), scale) + scb;
            if (v < 0.0) depth *= -1.0;
          } else {
            alpha = (oldmax*dmin - oldmin*dmax)/(oldmax - oldmin);
            beta = (dmax - dmin)/(oldmax - oldmin);
            depth = alpha + beta*depth;
          }
  			}

        if (scrap) {
          if (depth < scrapmin) {
            array[i][j] = gridded_data[ps_in].nodata_value;
            continue;
          }
          if (depth > scrapmax) {
            array[i][j] = gridded_data[ps_in].nodata_value;
            continue;
          }
        }
        
        if (power) {
          if (depth < 0.0) {
            depth = -1.0*power_a*pow(fabs(depth), power_b);
          } else {
            depth = power_a*pow(fabs(depth), power_b);
          }
        }

				array[i][j] = depth;

			} else {
				array[i][j] = gridded_data[ps_in].nodata_value;
			}
  	}
  }
}

