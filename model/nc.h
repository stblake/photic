
#include "common.h"
#include <netcdf.h>
#include <float.h>
#include <limits.h>

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}

void write_nc(char *file, float **grid, int ncols, int nrows, 
	float *lons, float *lats, double spval);

void read_nc(char *file, float ***grid, int *ncols, int *nrows, 
	float **lons, float **lats, double *spval);

