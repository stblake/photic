/*
 *
 *    smooth
 *
 */

/* Written by Sam Blake. */

/* Started on 28 Feb 2016. */

#include "smooth.h"
#include "smoothutils.h"

#define MXFILE 512

int main(int argc, char **argv) {
  
  int n, i, j, nrows, ncols, npasses;
  double wlon, slat, cellsize, spval;
  float **unsmoothed, **smoothed, weight;
  char inputfile[MXFILE], outputfile[MXFILE];

// Defaults. 

  npasses = 2; 
  weight  = 3;

// Read command line inputs. 

  for (n = 1; n < argc; n++) {
  	if (strcmp(argv[n], "-input") == 0) {
      strcpy(inputfile, argv[++n]);
      printf("-input %s\n", inputfile);
    } else if (strcmp(argv[n], "-output") == 0) {
      strcpy(outputfile, argv[++n]);
      printf("-output %s\n", outputfile);
    } else if (strcmp(argv[n], "-npasses") == 0) {
      npasses = atoi(argv[++n]);
      printf("\n-npasses %d", npasses);
    } else if (strcmp(argv[n], "-weight") == 0) {
      weight = atof(argv[++n]);
      printf("\n-weight %f", weight);
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

// Read in asc grid. 

  read_asc(inputfile, &unsmoothed, &ncols, &nrows, &wlon, &slat, &cellsize, &spval); 

// Allocate memory for smoothed array.   

  allocate_float_array_2d(&smoothed, nrows, ncols);

// Smooth the array. 

  smooth(unsmoothed, smoothed, ncols, nrows, weight, npasses, (float) spval);

// Write out smoothed array. 

  write_asc(outputfile, smoothed, ncols, nrows, wlon, slat, cellsize, spval);

// Free memory. 

  free_float_array_2d(unsmoothed, nrows);
  free_float_array_2d(smoothed, nrows);

  return 0;
}

