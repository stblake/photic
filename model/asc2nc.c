/*
 *
 *    asc2nc
 *
 */

/* Written by Sam Blake. */

/* Started on 25 Feb 2016. */


#include "asc2nc.h"

#define MXFILE 512

int main(int argc, char **argv) {
  
  int i, j, n, nrows, ncols;
  double wlon, slat, cellsize, spval;
  float **grid, *lons, *lats;
  char inputfile[MXFILE], outputfile[MXFILE];


// Read command line inputs. 

  for (n = 1; n < argc; n++) {
  	if (strcmp(argv[n], "-input") == 0) {
      strcpy(inputfile, argv[++n]);
      printf("-input %s\n", inputfile);
    } else if (strcmp(argv[n], "-output") == 0) {
      strcpy(outputfile, argv[++n]);
      printf("-output %s\n", outputfile);
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

// Read input asc grid file. 

  read_asc(inputfile, &grid, &ncols, &nrows, &wlon, &slat, &cellsize, &spval);

// Create vector of longitudes. 

  lons = (float*) malloc(ncols*sizeof(float));

  for (n = 0; n < ncols; n++) {
    lons[n] = wlon + cellsize*n;
  }

// Create vector of latitudes. 

  lats = (float*) malloc(nrows*sizeof(float));

  for (n = 0; n < nrows; n++) {
    lats[n] = slat + cellsize*n;
  }

// Write out netCDF file. 

  write_nc(outputfile, grid, ncols, nrows, lons, lats, spval);

// Free memory. 

  free(lons);
  free(lats);
  free_float_array_2d(grid, nrows);

  return 0;
}

