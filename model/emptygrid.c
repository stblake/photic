/*
 *
 *    emptygrid
 *
 */

/* Written by Sam Blake. */

/* Started on 23 Feb 2016. */


#include "common.h"
#include "io.h"


#define MXFILE 512

int main(int argc, char **argv) {
  
  int n, i, j, nrows, ncols;
  double wlon, slat, cellsize, spval;
  float **grid;
  char inputfile[MXFILE];

  // Read command line inputs. 

  for (n = 1; n < argc; n++) {
  	if (strcmp(argv[n], "-input") == 0) {
      strcpy(inputfile, argv[++n]);
      printf("-input %s\n", inputfile);
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

  // Read header data from output file. 
  read_asc_header(inputfile, &ncols, &nrows, &wlon, &slat, &cellsize, &spval);

  // Allocate memory. 
  allocate_float_array_2d(&grid, nrows, ncols);

  // Fill grid with spval. 
  for(i = 0; i < nrows; i++) {
  	for(j = 0; j < ncols; j++) {
  		grid[i][j] = spval;
  	}
  }

  write_asc(inputfile, grid, ncols, nrows, wlon, slat, cellsize, spval);

  // Free memory. 
  free_float_array_2d(grid, nrows);

  return 0;
}