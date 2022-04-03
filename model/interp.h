
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

void grid2gridinterp(
  float **grid_cg, int nx_cg, int ny_cg, 
  double wlon_cg, double slat_cg, double cellsize_cg, double spval_cg,
  float **grid_fg, int nx_fg, int ny_fg, 
  double wlon_fg, double slat_fg, double cellsize_fg, double spval_fg);


float interp_bilinear(float **grid, int nx, int ny, double wlon, double slat, double cellsize, 
      double spval, double lon, double lat);

float interp_bilinear_raw(float **grid, int nx, int ny, float xf, float yf, float spval);

float interp_bicubic(float **grid, int nx, int ny, double wlon, double slat, double cellsize, 
      double spval, double lon, double lat);

float interp_bicubic_raw(float **grid, int nx, int ny, float x, float y, float spval);