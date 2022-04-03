#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

void read_asc_header(char *file, int *ncols, int *nrows, double *wlon, double *slat, 
        double *cellsize, double *spval);
void read_asc(char *file, float ***grid, int *ncols, int *nrows, double *wlon, double *slat, 
	      double *cellsize, double *spval);
void write_asc(char *file, float **grid, int ncols, int nrows, double xllcorner,
	       double yllcorner, double cellsize, double spval);
void read_xyz(char *file, float **x, float **y, float **z, int *len);
void read_xy(char *file, float **x, float **y, int *len);
void write_xyz(char *file, float *x, float *y, float *z, int len);
void read_grid(char *file, float ***array, int *ncols, int *nrows, double *wlon, 
              double *slat, double *cellsize, double *spval);
void read_csv(char *file, float ***array, int *nrows, int *ncols);

#define MEMCHECK(p) \
  if (p == NULL)    \
    {\
      printf("\n\nERROR: Out of memory.\n\n");\
      exit(1);\
    }

