/*
 *
 *    smoothspval
 *
 */

/* Written by Sam Blake. */

/* Started on 24 Feb 2016. */


// This program only smoothes-out special values. 

#include "common.h"
#include "io.h"


#define MXFILE 512

int main(int argc, char **argv) {
  
  int n, np, i, j, ii, jj, nrows, ncols, nremoved;
  double wlon, slat, cellsize, spval;
  float **grid, s[8], interp;
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

  // Read header data from output file. 
  read_asc(inputfile, &grid, &ncols, &nrows, &wlon, &slat, &cellsize, &spval);

  nremoved = 0;

  // Smooth-out spvals. 

  for(i = 0; i < nrows; i++) {
  	for(j = 0; j < ncols; j++) {
      if (grid[i][j] != spval) 
        continue ;

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
        if (s[n] != spval) {
          interp += s[n];
          np++;
        }
      }

    	if (np != 0) {
        grid[i][j] = interp/((float) np);
        nremoved++;
      } else {
        printf("\nWARNING: cannot remove spval at i,j = %d,%d", i, j);
      }
  	}
  }

  printf("\nSmoothed out %d special values.\n", nremoved);

  write_asc(outputfile, grid, ncols, nrows, wlon, slat, cellsize, spval);

  // Free memory. 
  free_float_array_2d(grid, nrows);

  return 0;
}
