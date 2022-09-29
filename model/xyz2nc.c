/*
 *
 *    xyz2nc
 *
 */

/* Written by Sam Blake. */


#include "xyz2nc.h"

#define MXFILE 512

int main(int argc, char **argv) {

  int n, np, i, j, ii, jj, len, nrows, ncols, nremoved;
  float spval, *x, *y, *z, *lons, *lats, wlon, elon, slat, nlat, 
  	xa, ya, dx, dy, dxy, **grid, **smoothed, interp, s[8], eps, 
    scale = 1.0, offset = 0.0, resolution;
  char inputfile[MXFILE], outputfile[MXFILE];
  bool smoothspval = false, user_resolution = false;

// Set defaults. 

  spval = -99999.0; 
  eps = 1.0e-4;

// Read command line inputs. 

  for (n = 1; n < argc; n++) {
  	if (strcmp(argv[n], "-input") == 0) {
      strcpy(inputfile, argv[++n]);
      printf("-input %s\n", inputfile);
    } else if (strcmp(argv[n], "-output") == 0) {
      strcpy(outputfile, argv[++n]);
      printf("-output %s\n", outputfile);
    } else if (strcmp(argv[n], "-spval") == 0 || strcmp(argv[n], "-NODATA_value") == 0) {
      spval = atof(argv[++n]);
      printf("-spval %f\n", spval);
    } else if (strcmp(argv[n], "-smoothspval") == 0) {
      smoothspval = true;
      printf("-smoothspval \n");
    } else if (strcmp(argv[n], "-scale") == 0) {
      scale = atof(argv[++n]);;
      printf("-scale %.7f \n", scale);
    } else if (strcmp(argv[n], "-offset") == 0) {
      offset = atof(argv[++n]);;
      printf("-offset %.7f \n", offset);
    } else if (strcmp(argv[n], "-resolution") == 0) {
      user_resolution = true; 
      resolution = atof(argv[++n]);;
      printf("-resolution %.7f \n", resolution);
    } else {
      printf("\nERROR: unknown option:  %s\n\n", argv[n]);
      exit(1);
    }
  }

// Check input file exists. 

  if (! file_exists(inputfile)) {
    printf("\nERROR: missing file: %s\n", inputfile);
    exit(1);
  }

// Read xyz file.

  read_xyz(inputfile, &x, &y, &z, &len);

  printf("\nnumber of points = %d", len);

// Determine slat, nlat, wlon, elon. 

  wlon = vec_min(x, len);
  elon = vec_max(x, len);

  slat = vec_min(y, len);
  nlat = vec_max(y, len);

  printf("\nwlon, elon = %f, %f", wlon, elon);
  printf("\nslat, nlat = %f, %f", slat, nlat);

// Determine dx, dy.  

  ya = y[0];
  n = 1;
  while (ya == y[n]){
  	n++;
  }
  dy = fabs(y[n] - ya);

  xa = x[0];
  n = 1;
  while (xa == x[n]) {
  	n++;
  }
  dx = fabs(x[n] - xa);

  if (user_resolution) {
    dx = resolution;
    dy = resolution; 
    printf("\ndx,dy = %f, %f", dx, dy);
  }

// Determine output grid size. 

  dxy = (dx < dy) ? dx : dy;

  nrows = 1 + floor((nlat - slat)/dxy);
  ncols = 1 + floor((elon - wlon)/dxy);

  printf("\nncols, nrows = %d, %d", ncols, nrows);

// Create output grid and fill with spval. 

  allocate_float_array_2d(&grid, nrows, ncols);

  for (i = 0; i < nrows; i++) {
  	for (j = 0; j < ncols; j++) {
  		grid[i][j] = spval;
  	}
  }

// Convert (gridded) xyz data to grid.

  for (n = 0; n < len; n++) {
  	// Compute index on the grid.
  	i = floor((y[n] - slat)/dxy);
  	j = floor((x[n] - wlon)/dxy);
    if (z[n] == spval)
      grid[i][j] = spval;
    else
      grid[i][j] = scale*z[n] + offset;
  }

// Free xyz memory. 

  free(x);
  free(y);
  free(z);

// Smooth out landsat artifacts. 

  if (smoothspval) {

  	nremoved = 0;
  	allocate_float_array_2d(&smoothed, nrows, ncols);

  	for(i = 0; i < nrows; i++) {
  	  for(j = 0; j < ncols; j++) {

        if (fabs(grid[i][j] - spval) > eps) {
          smoothed[i][j] = grid[i][j];
          continue ;
        }

        if (i == 0)
          ii = i + 1;
        else if (i == nrows - 1)
          ii = i - 1;
        else 
          ii = i;

        if (j == 0)
          jj = j + 1;
        else if (j == ncols - 1)
          jj = j - 1;
        else 
          jj = j;

        s[0] = grid[ii + 1][jj - 1];
        s[1] = grid[ii + 1][jj];
        s[2] = grid[ii + 1][jj + 1];
        s[3] = grid[ii][jj - 1];
        s[4] = grid[ii][jj + 1];
        s[5] = grid[ii - 1][jj - 1];
        s[6] = grid[ii - 1][jj];
        s[7] = grid[ii - 1][jj + 1];

        interp = 0.0;
        np = 0;

        for (n = 0; n < 8; n++) {
          if (fabs(s[n] - spval) > eps) {
            interp += s[n];
            np++;
          }
        }

    	if (np != 0) {
    	  interp /= ((float) np);

    	  if (fabs(interp - spval) < eps)
    	  	interp = spval;
    	  else
    	    nremoved++;

          smoothed[i][j] = interp;
        } else {
          smoothed[i][j] = spval;
        }
   	  }
    }

    free_float_array_2d(grid, nrows);
    printf("\nRemoved %d special values.\n\n", nremoved);
  }

// Export gridded data. 

  lons = (float*) malloc(ncols*sizeof(float));

  for (n = 0; n < ncols; n++) {
    lons[n] = wlon + dxy*n;
  }

  lats = (float*) malloc(nrows*sizeof(float));

  for (n = 0; n < nrows; n++) {
    lats[n] = slat + dxy*n;
  }

  if (smoothspval)
  	write_nc(outputfile, smoothed, ncols, nrows, lons, lats, (double) spval);
  else
  	write_nc(outputfile, grid, ncols, nrows, lons, lats, (double) spval);

// Free memory. 

  free(lons);
  free(lats);
  if (smoothspval) 
  	free_float_array_2d(smoothed, nrows);
  else
  	free_float_array_2d(grid, nrows);

  return 0;
}
