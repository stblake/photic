/*
 *
 *    gridinterp
 *
 */

/* Written by Sam Blake. */

/* Started on 20 Feb 2016. */

 /* This program interpolates from one ESRI ascii grid to another. The input file
 must be an ESRI ascii formatted file and the output file must contain the first 
 6 header lines with the appropriate output parameters. */

/*
    Command line inputs: 

      -input /name of input .asc formatted file/

      -output /name of the output .asc formatted file/

*/

#include "common.h"
#include "io.h"
#include "interp.h"



#define MXFILE 512

int main(int argc, char **argv) {

  char inputfile[MXFILE], outputfile[MXFILE];
  int n, nx_cg, ny_cg, nx_fg, ny_fg;
  double wlon_cg, elon_cg, slat_cg, nlat_cg, spval_cg, cellsize_cg,
         wlon_fg, elon_fg, slat_fg, nlat_fg, spval_fg, cellsize_fg,
         resolution = 0.0;
  float **grid_cg, **grid_fg;

// Read command line inputs. 

  for (n = 1; n < argc; n++) {
  	if (strcmp(argv[n], "-input") == 0) {
      strcpy(inputfile, argv[++n]);
      printf("-input %s\n", inputfile);
    } else if (strcmp(argv[n], "-output") == 0) {
      strcpy(outputfile, argv[++n]);
      printf("-output %s\n", outputfile);
    } else if (strcmp(argv[n], "-resolution") == 0) {
      resolution = atof(argv[++n]);
      printf("-resolution %lf", resolution);
	  } else {
      printf("\nERROR: unknown option:  %s\n\n", argv[n]);
      exit(1);
    }
  }

// Check input and output files exist. 

  if (file_exists(inputfile)) {
    // Read coarse grid file. 
    read_asc(inputfile, &grid_cg, &nx_cg, &ny_cg, &wlon_cg, &slat_cg, &cellsize_cg, &spval_cg);
    nlat_cg = slat_cg + (ny_cg - 1)*cellsize_cg;
    elon_cg = wlon_cg + (nx_cg - 1)*cellsize_cg;
  } else {
    // Input file missing. 
    printf("\nERROR: missing file: %s\n", inputfile);
    exit(1);
  }

  if (file_exists(outputfile)) {
    // Read fine grid file. 
    read_asc(outputfile, &grid_fg, &nx_fg, &ny_fg, &wlon_fg, &slat_fg, &cellsize_fg, &spval_fg); 
    nlat_fg = slat_fg + (ny_fg - 1)*cellsize_fg; 
    elon_fg = wlon_fg + (nx_fg - 1)*cellsize_fg;   
  } else if (resolution == 0.0) {
    printf("\nERROR: missing output resolution.\n");
    exit(1);
  } else {
    wlon_fg = wlon_cg;
    slat_fg = slat_cg;
    nlat_fg = nlat_cg; 
    elon_fg = elon_cg;  
    cellsize_fg = resolution; 
    spval_fg = spval_cg;  
    nx_fg = floor((elon_fg - wlon_fg)/cellsize_fg) + 1;
    ny_fg = floor((nlat_fg - slat_fg)/cellsize_fg) + 1;
    allocate_float_array_2d(&grid_fg, ny_fg, nx_fg);
  }

// Check fine grid overlaps the coarse grid. 

  if (! overlap(wlon_cg, elon_cg, slat_cg, nlat_cg,
  				wlon_fg, elon_fg, slat_fg, nlat_fg)) {
  	printf("\nERROR: output grid is not a subset of the input grid.\n");
  	exit(1);
  }

// Interpolate from coarse grid to fine grid. 

  grid2gridinterp(grid_cg, nx_cg, ny_cg, wlon_cg, slat_cg, cellsize_cg, spval_cg,
  				  grid_fg, nx_fg, ny_fg, wlon_fg, slat_fg, cellsize_fg, spval_fg);

// Write out fine grid. 

  write_asc(outputfile, grid_fg, nx_fg, ny_fg, wlon_fg, slat_fg, cellsize_fg, spval_fg);

// Free memory.

  free_float_array_2d(grid_cg, ny_cg);
  free_float_array_2d(grid_fg, ny_fg);

  return 0;
}

