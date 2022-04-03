
//
//    filllakes.c
//

/* Written by Sam Blake. */

/* Started on 28 Feb 2016. */

#include "filllakes.h"
#include "filllakesutils.h"

#define MXFILE 512

int main(int argc, char **argv) {
  
  int n, i, j, nrows, ncols, maxwet;
  double wlon, slat, cellsize, spval;
  float **grid;
  char inputfile[MXFILE], outputfile[MXFILE];
  bool spvalpos, fillspval;

// Default. 

  spvalpos = true;  // Assumes spval is positive. (For bounding purposes.)
  fillspval = false;
  maxwet = 1;

// Read command line inputs. 

  for (n = 1; n < argc; n++) {
  	if (strcmp(argv[n], "-input") == 0) {
      strcpy(inputfile, argv[++n]);
      printf("\n-input %s", inputfile);
    } else if (strcmp(argv[n], "-output") == 0) {
      strcpy(outputfile, argv[++n]);
      printf("\n-output %s", outputfile);
    } else if (strcmp(argv[n], "-maxpoints") == 0) {
      maxwet = atoi(argv[++n]);
      printf("\n-maxpoints %d", maxwet);
    } else if (strcmp(argv[n], "-spvalpositive") == 0) {
      spvalpos = true;
      printf("\n-spvalpositive");
    } else if (strcmp(argv[n], "-spvalnegative") == 0) {
      spvalpos = false;
      printf("\n-spvalnegative");
    } else if (strcmp(argv[n], "-fillspval") == 0) {
      fillspval = true;
      printf("\n-fillspval");
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

  read_asc(inputfile, &grid, &ncols, &nrows, &wlon, &slat, &cellsize, &spval); 

// Fill in small lakes. 

  fill_lakes(grid, ncols, nrows, maxwet, spvalpos, fillspval, (float) spval);

// Write out filled array. 

  write_asc(outputfile, grid, ncols, nrows, wlon, slat, cellsize, spval);

// Free memory. 

  free_float_array_2d(grid, nrows);

  return 0;
}





