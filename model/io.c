
#include "io.h"
#include "common.h"
#include "nc.h"

#define NETCDF_GRID 19
#define ESRI_GRID   23
#define GEOTIFF     25
#define UNKNOWN_GRID 27

int grid_type(char *file); // Not externally visible. 

// Generic grid reader. 

void read_grid(char *file, float ***array, int *ncols, int *nrows, double *wlon, 
              double *slat, double *cellsize, double *spval) {
  int type;
  float cellsizeX, cellsizeY, *lons, *lats, elon, nlat; 

// Inspect grid type. 

  type = grid_type(file);

  switch (type) {
    case NETCDF_GRID: 
      // Read netCDF grid.
      read_nc(file, array, ncols, nrows, &lons, &lats, spval);
      *wlon = vec_min(lons, *ncols);
      *slat = vec_min(lats, *nrows);
      elon = vec_max(lons, *ncols);
      nlat = vec_max(lats, *nrows);
      cellsizeX = fabs( ((double) elon) - ((double) *wlon) )/((double) *ncols - 1);
      cellsizeY = fabs( ((double) nlat) - ((double) *slat) )/((double) *nrows - 1);
      *cellsize = 0.5*(cellsizeX + cellsizeY);
      free(lons);
      free(lats);
      break;
    case ESRI_GRID:
      // Read ESRI ASCII grid. 
      read_asc(file, array, ncols, nrows, wlon, slat, cellsize, spval);
      break;
    default:
      printf("\nERROR: unknown grid type: %s", trim(file));
      exit(1);
  }

}

int grid_type(char *file) {
  if (strstr(file, ".nc") != NULL) {
    return NETCDF_GRID;
  } else if (strstr(file, ".asc") != NULL) {
    return ESRI_GRID;
  } else if (strstr(file, ".tif") != NULL) {
    return GEOTIFF;
  } else {
    return UNKNOWN_GRID;
  }
}

// ESRI grid spec: https://en.wikipedia.org/wiki/Esri_grid

void read_asc_header(char *file, int *ncols, int *nrows, double *wlon, double *slat, 
        double *cellsize, double *spval) {

  char str[32];
  FILE *fp;  

  fp = fopen(file, "r");

// ncols

  fscanf(fp, "%s %d", str, ncols);

  if ( strcmp(str, "ncols") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n");
    exit(1);
  }

  // printf("\nncols %d", *ncols);
  
// nrows

  fscanf(fp, "%s %d", str, nrows);  

  if ( strcmp(str, "nrows") != 0 ) { 
    printf("\nERROR: malformed ascii grid.\n\n");
    exit(1);
  }

  // printf("\nnrows %d", *nrows); 

// xllcorner

  fscanf(fp, "%s %lf", str, wlon); 

  if ( strcmp(str, "xllcorner") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n"); 
    exit(1);
  }

  // printf("\nxllcorner %lf", *wlon);

// yllcorner

  fscanf(fp, "%s %lf", str, slat);

  if ( strcmp(str, "yllcorner") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n"); 
    exit(1);
  }

  // printf("\nyllcorner %lf", *slat);

// cellsize

  fscanf(fp, "%s %lf", str, cellsize);

  if ( strcmp(str, "cellsize") != 0 ) { 
    printf("\nERROR: malformed ascii grid.\n\n");
    exit(1);
  }

  // printf("\ncellsize %lf", *cellsize);

// spval

  fscanf(fp, "%s %lf", str, spval);

  if ( strcmp(str, "NODATA_value") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n"); 
    exit(1);
  }
 
  // printf("\nNODATA_value %lf", *spval);

  // printf("\n\n");
  fclose(fp);
  return ;

}

void read_asc(char *file, float ***grid, int *ncols, int *nrows, double *wlon, double *slat, 
	      double *cellsize, double *spval) {

  int m, row, col;
  char str[32];
  FILE *fp;  
  float height;


  fp = fopen(file, "r");

// ncols

  fscanf(fp, "%s %d", str, ncols);

  if ( strcmp(str, "ncols") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n");
    exit(1);
  }

  // printf("\nncols %d", *ncols);
  
// nrows

  fscanf(fp, "%s %d", str, nrows);  

  if ( strcmp(str, "nrows") != 0 ) { 
    printf("\nERROR: malformed ascii grid.\n\n");
    exit(1);
  }

  // printf("\nnrows %d", *nrows); 

// xllcorner

  fscanf(fp, "%s %lf", str, wlon); 

  if ( strcmp(str, "xllcorner") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n"); 
    exit(1);
  }

  // printf("\nxllcorner %lf", *wlon);

// yllcorner

  fscanf(fp, "%s %lf", str, slat);

  if ( strcmp(str, "yllcorner") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n"); 
    exit(1);
  }

  // printf("\nyllcorner %lf", *slat);

// cellsize

  fscanf(fp, "%s %lf", str, cellsize);

  if ( strcmp(str, "cellsize") != 0 ) { 
    printf("\nERROR: malformed ascii grid.\n\n");
    exit(1);
  }

  // printf("\ncellsize %lf", *cellsize);

// spval

  fscanf(fp, "%s %lf", str, spval);

  if ( strcmp(str, "NODATA_value") != 0 ) {
    printf("\nERROR: malformed ascii grid.\n\n"); 
    exit(1);
  }
 
  // printf("\nNODATA_value %lf", *spval);

// Allocate memory for the array. 

  allocate_float_array_2d(grid, *nrows, *ncols);

  // gridded data. We flip the grid in the vertical direction
  // so the [0][0] entry corresponds to the south-west corner
  // of the grid. 

  row = 0; col = 0;
  while ( ! feof(fp) ) {
    m = fscanf(fp, "%f ", &height);
    if (m == 1) {
        (*grid)[*nrows - row - 1][col++] = height;
      if (col == (*ncols)) {
        col = 0;
        row++;
      }
    }
  }

  // printf("\n\n");
  fclose(fp);
  return ;
}


void write_asc(char *file, float **grid, int ncols, int nrows, double xllcorner, 
	       double yllcorner, double cellsize, double spval) {
  FILE *fp;
  int i, j;

  fp = fopen(file, "w"); 

  fprintf(fp, "ncols        %d\n", ncols);
  fprintf(fp, "nrows        %d\n", nrows);
  fprintf(fp, "xllcorner    %lf\n", xllcorner);
  fprintf(fp, "yllcorner    %lf\n", yllcorner);
  fprintf(fp, "cellsize     %lf\n", cellsize);
  fprintf(fp, "NODATA_value %lf\n", spval);
 
  for (i = nrows - 1; i >= 0; i--) {
    for (j = 0; j < ncols; j++) {
      fprintf(fp, "%.2f ", grid[i][j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return ;
}


void read_csv(char *file, float ***array, int *nrows, int *ncols) {

  int nr, nc, row, col, m;
  float v;
  char c;

  FILE *fp;

// Read nrows. 

  fp = fopen(file, "r");

  nr = 0;
  while ( (c = fgetc(fp)) != EOF ) {
    if ( c == '\n' )
      nr++;
  }

  *nrows = nr + 1;

  fclose(fp);

// Read ncols.

  fp = fopen(file, "r");

  nc = 0;
  while ( (c = fgetc(fp)) != '\n' && c != EOF) {
    if ( c == ',' )
      nc++;
  }

  *ncols = nc + 1;

  fclose(fp);

// Allocate memory. 

  allocate_float_array_2d(array, *nrows, *ncols);

// Read file data into array. 

  fp = fopen(file, "r");

  row = 0;
  col = 0;

  while ( ! feof(fp) ) {

    if (col == (*ncols)) {
      m = fscanf(fp, "%f", &v);
    } else {
      m = fscanf(fp, "%f,", &v);
    }

    if (m == 1) {
      (*array)[row][col++] = v;
      if (col == (*ncols)) {
        col = 0;
        row++;
      }
    }
  }

  fclose(fp);

}


void read_xyz(char *file, float **x, float **y, float **z, int *len) {

  int i, nlines = 0, count, matched;
  bool csv;
  FILE *fp;
  char c, line[256];

// Count the number of points.  

  fp = fopen(file, "r");

  while ( (c = fgetc(fp)) != EOF ) {
    if ( c == '\n' )
      nlines++;
  }

  *len = nlines;

  fclose(fp);

// Allocate memory. 

  (*x) = (float*) malloc(nlines*sizeof(float));
  MEMCHECK((*x));

  (*y) = (float*) malloc(nlines*sizeof(float));
  MEMCHECK((*y));

  (*z) = (float*) malloc(nlines*sizeof(float));
  MEMCHECK((*z));

// Check for space or comma separated lines. 

  fp = fopen(file, "r");

  fgets(line, sizeof(line), fp);

  if(strstr(line, ",") != NULL) {
    csv = true;
  } else {
    csv = false;
  }

  fclose(fp);

  fp = fopen(file, "r");  

// Read in xyz data. 

  if (csv) {
    for (i = 0; i < nlines; i++) {
      fgets(line, sizeof(line), fp);
      matched = sscanf(line, "%f, %f, %f\n", &((*x)[i]), &((*y)[i]), &((*z)[i]));
      if (matched != 3) {
        printf("\nERROR: xyz file reader error.\n");
        exit(1);
      }    }
  } else {
    for (i = 0; i < nlines; i++) {
      fgets(line, sizeof(line), fp);
      matched = sscanf(line, "%f %f %f\n", &((*x)[i]), &((*y)[i]), &((*z)[i]));
      if (matched != 3) {
        printf("\nERROR: xyz file reader error.\n");
        exit(1);
      }
    }
  }

  fclose(fp);
}


void write_xyz(char *file, float *x, float *y, float *z, int len) {

  FILE *fp;
  int i;

  fp = fopen(file, "a");

  for (i = 0; i < len; i++) {
    fprintf(fp, "%.6f, %.6f, %.6f \n", x[i], y[i], z[i]);
  }

  fclose(fp);
  return ;
}


void read_xy(char *file, float **x, float **y, int *len) {

  int i, nlines = 0, count, matched;
  bool csv;
  FILE *fp;
  char c, line[256];

// Count the number of points.  

  fp = fopen(file, "r");

  while ( (c = fgetc(fp)) != EOF ) {
    if ( c == '\n' )
      nlines++;
  }

  *len = nlines;

  fclose(fp);

// Allocate memory. 

  (*x) = (float*) malloc(nlines*sizeof(float));
  MEMCHECK((*x));

  (*y) = (float*) malloc(nlines*sizeof(float));
  MEMCHECK((*y));

// Check for space or comma separated lines. 

  fp = fopen(file, "r");

  fgets(line, sizeof(line), fp);

  if(strstr(line, ",") != NULL) {
    csv = true;
  } else {
    csv = false;
  }

  fclose(fp);

  fp = fopen(file, "r");  

// Read in xy data. 

  if (csv) {
    for (i = 0; i < nlines; i++) {
      fgets(line, sizeof(line), fp);
      matched = sscanf(line, "%f, %f\n", &((*x)[i]), &((*y)[i]));
      if (matched != 2) {
        printf("\nERROR: xy file reader error.\n");
        exit(1);
      }    }
  } else {
    for (i = 0; i < nlines; i++) {
      fgets(line, sizeof(line), fp);
      matched = sscanf(line, "%f %f\n", &((*x)[i]), &((*y)[i]));
      if (matched != 2) {
        printf("\nERROR: xy file reader error.\n");
        exit(1);
      }
    }
  }

  fclose(fp);
}





