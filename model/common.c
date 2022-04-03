
#include "common.h"

#define BIG 1.0e10



#if PLOTTING
int n_fn_evals = 0;
double Error_path[6000];
double P_path[6000]; 
double G_path[6000]; 
double X_path[6000];
double B_path[6000];
double H_path[6000];
double Delta_path[6000];
#endif


// Globals. 
int n_grids_in_use = 0;
int n_images_in_use = 0;
float meminuse = 0.0; // Mb

char command_file[MAX_STRING_LEN];
bool allocated_grids[MAX_GRIDS];
bool allocated_images[MAX_IMAGES];
char grid_names[MAX_GRIDS][MAX_STRING_LEN];
char image_names[MAX_IMAGES][MAX_STRING_LEN];
geogrid gridded_data[MAX_GRIDS];
image image_data[MAX_IMAGES];
line_segment line_data[MAX_GRIDS];
char parsed[MAX_ARGS][MAX_ARG_SIZE];


scene scene_data[MAX_SCENES];
bool allocated_scenes[MAX_SCENES];
char scene_names[MAX_SCENES][MAX_STRING_LEN];
int n_scenes_in_use = 0;

char prev1[MAX_STRING_LEN];
char prev2[MAX_STRING_LEN];
char prev3[MAX_STRING_LEN];
char prev4[MAX_STRING_LEN];
char prev5[MAX_STRING_LEN];
char prev6[MAX_STRING_LEN];
char prev7[MAX_STRING_LEN];
char prev8[MAX_STRING_LEN];
char prev9[MAX_STRING_LEN];

// Depth soundings. 

bool read_soundings = false;
bool first_model_error_soundings = true;
int nsoundings = 0;
float **soundings;
float *model_depths, *sounding_depths, *sounding_errors, *sounding_weights, **radiance_at_soundings, 
  *dz_at_soundings, *lsb_at_soundings;

// Graphics defaults. 

int background = BLACK;
int palette = GRAYSCALE;
float pagesize = 8, textsize = 1.0, pointsize = 1.0;
int line_width = 2;

// Satellite metadata. 

float sun_azimuth = BAM_NODATA;
float sun_elevation = BAM_NODATA;
float earth_sun_distance = BAM_NODATA;
float radiance_mult_band[MAX_GRIDS];
float radiance_add_band[MAX_GRIDS];
float quantize_cal_min_band[MAX_GRIDS];
float quantize_cal_max_band[MAX_GRIDS];
float radiance_minimum_band[MAX_GRIDS];
float radiance_maximum_band[MAX_GRIDS];

// Model defaults.

float H0;

bool relative_error = false;

float (*zdepth)(int, int);

char model_name[MAX_STRING_LEN];

bool pan_normalise = true;

float extent_eps = 0.005;

float epsilon = 1.0e-4;
// float log_epsilon = 0.00001; 
float log_epsilon = 1.0;
float falloff_param = 1.0; // units km.
float model_weight_cutoff = 0.0001;
float map_scale = 1.0;

int nboundary = 256;

float deep_lon = 0.0;
float deep_lat = 0.0;
float sand_lon = 0.0;
float sand_lat = 0.0;

int n_local_models = 1;
int depth_radius = 2;
int n_depth_bins = 6;
float lsb_spread = 0.0;

float sounding_error = 0.0;
float constraint_error = 0.0;
float local_constraint_error = 0.0;

int n_bottom_spectra = 0;

float noise_threshold = 0.025;
float min_radiance_threshold = 0.0;
float max_radiance_threshold = 0.0;

int bands_indexes[MAX_GRIDS];
int spectral_indexes[MAX_GRIDS];
int land_index = 0;
int shallow_index = 0;

int n_model_bands = 0, **model_origins, NAlpha[MAX_GRIDS]; 
float thetas[MAX_GRIDS], model_wavelengths[MAX_GRIDS], Lsm[MAX_GRIDS], Lsb[MAX_GRIDS], Lsb_sigma[MAX_GRIDS], 
  Alpha[MAX_GRIDS], Hj[MAX_GRIDS], Lsm_sigma[MAX_GRIDS], *mparams;

float K_ratios[MAX_GRIDS][4];

float jerlov_water_type = 2.0; // Default is OIA + 0.00
int n_water_type_estimates = 0;

float slope_line[LINE_MAX_POINTS];

float depthmin = 0.0;
float depthmax = 30.0;

float line_xa, line_ya, line_xb, line_yb;

float bound_error = 0.0, coeff_size_preference = 0.0, regression_error = 0.0;

float weight_soundings = 0.775, weight_constraints = 0.3, weight_local_constraints = 0.15, 
  weight_coeff_size_preference = 0.05, weight_bound_error = 0.05, weight_regression = 0.7;

float rms_error;

int model_optimisation = SIMULATED_ANNEALING;

int samodel_n_smoothing_radius = 1, samodel_n_spatial = 2, samodel_n_bottoms = 3;

bool commensurate_grids(geogrid a, geogrid b) {

  return a.nrows == b.nrows && 
         a.ncols == b.ncols && 
    fabs(a.slat - b.slat) < extent_eps && 
    fabs(a.nlat - b.nlat) < extent_eps && 
    fabs(a.wlon - b.wlon) < extent_eps && 
    fabs(a.elon - b.elon) < extent_eps && 
    fabs(a.cellsize - b.cellsize) < extent_eps;
}

int cmp_floats(const void *a, const void *b) {
    
  struct sort_index_float_struct *a1 = (struct sort_index_float_struct *)a;
  struct sort_index_float_struct *a2 = (struct sort_index_float_struct *)b;

  if ((*a1).value < (*a2).value) {
    return -1;
  } else if ((*a1).value > (*a2).value) {
    return 1;
  } else {
    return 0;
  }
}


void sort_indexes_2d(float **array, int nrows, int ncols, int **indexes) {

  int i, j, k, nelems;
  struct sort_index_float_struct *temp;

  nelems = nrows*ncols;
  temp = malloc(nelems*sizeof(struct sort_index_float_struct));

  k = 0;
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      temp[k].value = array[i][j];
      temp[k].indx = k;
      k++;
    }
  }

  qsort(temp, nelems, sizeof(struct sort_index_float_struct), cmp_floats);

  k = 0;
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      indexes[i][j] = temp[k++].indx;
    }
  }

  free(temp);

  return ;
}



// RANDOM FLOATS

float frand() {
  return ((float) rand())/((float) RAND_MAX);
}


float frand2(float max) {
  if (frand() < 0.5)
    return ((float) rand())/((float) RAND_MAX/max);
  else
    return -1.0*((float) rand())/((float) RAND_MAX/max);
}


//
//  SMOOTHED ARRAY ELEMENT
//

float smoothed_index_2d(float **array, int i, int j, int nrows, int ncols, int smoothing_radius, float spval) {

  int ii, jj, iii, jjj;
  float v, smoothed, nobs; 

  smoothed = 0.0;
  nobs = 0.0;
        
  for (ii = 1 - smoothing_radius; ii < smoothing_radius; ii++) {

    if (i + ii < 0) {
      iii = 0;
    } else if (i + ii > nrows - 1) {
      iii = nrows - 1;
    } else { 
      iii = i + ii;
    }

    for (jj = 1 - smoothing_radius; jj < smoothing_radius; jj++) {

      if (j + jj < 0) {
        jjj = 0;
      } else if (j + jj > ncols - 1) {
        jjj = ncols - 1;
      } else {
        jjj = j + jj;
      }

      v = array[iii][jjj];

      if (! approx_equal(v, spval, 1.0e-6)) {
        smoothed += v;
        nobs += 1.0;
      }
    }
  }

  if (nobs < 0.5) {
    smoothed = spval;
  } else {
    smoothed /= nobs;
  }

  return smoothed;
}


//
//  RANDOM NORMAL 
//

double random_norm(double mean, double stddev) {

  double x1, x2, rn;   
  
  x1 = ((double) rand())/((double) RAND_MAX);
  x2 = ((double) rand())/((double) RAND_MAX);
  rn = sqrt(-2.0*log(x1))*cos(2.0*PI*x2);

  return mean + stddev*rn; 
}

//
//  1D INTERPOLATION
//

double interp_1d(double *X, double *Y, int n, double x) {

  int i;
  double x0, x1, y0, y1, alpha;

  if (approx_equal(x, X[0], 1.0e-5)) {
    return Y[0];
  } else if (approx_equal(x, X[n - 1], 1.0e-5)) {
    return Y[n - 1];
  } else if (X[0] < X[n - 1] && x < X[0]) {
    // printf("\nWARNING: interp_1d is using extrapolation!\n");
    x0 = X[0];
    x1 = X[1];
    y0 = Y[0];
    y1 = Y[1];
  } else if (X[n - 1] > X[0] && x > X[n - 1]) {
    // printf("\nWARNING: interp_1d is using extrapolation!\n");
    x0 = X[n - 2];
    x1 = X[n - 1];
    y0 = Y[n - 2];
    y1 = Y[n - 1];
  } else {
    for (i = 0; i < n - 1; i++) {
      if ((approx_less_equal(X[i], x, 1.0e-5) && approx_greater_equal(X[i + 1], x, 1.0e-5)) || 
        (approx_greater_equal(X[i], x, 1.0e-5) && approx_less_equal(X[i + 1], x, 1.0e-5))) {
        x0 = X[i];
        x1 = X[i + 1];
        y0 = Y[i];
        y1 = Y[i + 1];
        break ;
      }
    }
  }

  alpha = (x - x0)/(x1 - x0);
  
  return y0*(1.0 - alpha) + y1*alpha;
}


float interp_float_1d(float *X, float *Y, int n, float x) {

  int i;
  float x0, x1, y0, y1, alpha;

  if (approx_equal(x, X[0], 1.0e-5)) {
    return Y[0];
  } else if (approx_equal(x, X[n - 1], 1.0e-5)) {
    return Y[n - 1];
  } else if (X[0] < X[n - 1] && x < X[0]) {
    // printf("\nWARNING: interp_1d is using extrapolation!\n");
    x0 = X[0];
    x1 = X[1];
    y0 = Y[0];
    y1 = Y[1];
  } else if (X[n - 1] > X[0] && x > X[n - 1]) {
    // printf("\nWARNING: interp_1d is using extrapolation!\n");
    x0 = X[n - 2];
    x1 = X[n - 1];
    y0 = Y[n - 2];
    y1 = Y[n - 1];
  } else {
    for (i = 0; i < n - 1; i++) {
      if ((approx_less_equal(X[i], x, 1.0e-5) && approx_greater_equal(X[i + 1], x, 1.0e-5)) || 
        (approx_greater_equal(X[i], x, 1.0e-5) && approx_less_equal(X[i + 1], x, 1.0e-5))) {
        x0 = X[i];
        x1 = X[i + 1];
        y0 = Y[i];
        y1 = Y[i + 1];
        break ;
      }
    }
  }

  alpha = (x - x0)/(x1 - x0);
  
  return y0*(1.0 - alpha) + y1*alpha;
}


bool all_approx_zero(float *x, int n, float epsilon) {

  int i;

  for (i = 0; i < n; i++) {
    if (! approx_equal(x[i], 0.0, epsilon))
      return false; 
  }

  return true;
}

// Ref: Knuth, Semi-Numerical Algorithms, 4.2.2, eqn. 22. 

bool approx_equal(float a, float b, float epsilon) {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool approx_equal_double(double a, double b, double epsilon) {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool approx_less_equal(float a, float b, float epsilon) {
    return a < b || fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool approx_greater_equal(float a, float b, float epsilon) {
    return a > b || fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

// Ref: Knuth, Semi-Numerical Algorithms, 4.2.2, eqn. 24. 

bool essentially_equal(float a, float b, float epsilon) {
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

// Ref: Numerical Recipies in FORTRAN, ch 15, p. 656

// The linear fit is y = m*x + b

bool linear_fit(float *x, float *y, int n, float *m, float *b, float *r) {

  int i;
  double Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0, Syy = 0.0, delta; 

  for (i = 0; i < n; i++) {
    Sx += x[i];
    Sy += y[i];
    Sxx += x[i]*x[i];
    Sxy += x[i]*y[i];
    Syy += y[i]*y[i];
  }

  delta = ((double) n)*Sxx - Sx*Sx;

  if (approx_equal(delta, 0.0, 1.0e-6)) {
    printf("\nERROR: singular matrix - cannot compute a linear fit.\n");
    *m = 0.0; 
    *b = 0.0; 
    if (r) *r = 0.0;
    return false;
  }

  *m = (float) (((double) n)*Sxy - Sx*Sy)/delta;
  *b = (float) (Sy*Sxx - Sx*Sxy)/delta;

  if (r != NULL) {
    *r = (float) (Sxy - Sx*Sy/((double) n)) /
          sqrt((Sxx - Sx*Sx/((double) n))*(Syy - Sy*Sy/((double) n)));
    }

  return true;
}



bool linear_fit2(float *x, float *y, int n, float *m, float *b, float *r, float spval) {

  int i, nv = 0;
  double Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0, Syy = 0.0, delta; 

  for (i = 0; i < n; i++) {
    if (! approx_equal(x[i], spval, 1.0e-6) && ! approx_equal(y[i], spval, 1.0e-6)) {
      Sx += x[i];
      Sy += y[i];
      Sxx += x[i]*x[i];
      Sxy += x[i]*y[i];
      Syy += y[i]*y[i];
      nv++;
    }
  }

  delta = ((double) nv)*Sxx - Sx*Sx;

  if (approx_equal(delta, 0.0, 1.0e-6)) {
    printf("\nERROR: singular matrix - cannot compute a linear fit.\n");
//    printf("\n");
//    for (i = 0; i < n; i++) 
//      printf("%.3f ", x[i]);
//    printf("\n");
//    for (i = 0; i < n; i++) 
//      printf("%.3f ", y[i]);
//    printf("\n");
    *m = 0.0; 
    *b = 0.0; 
    if (r) *r = 0.0;
    return false;
  }

  *m = (float) (((double) nv)*Sxy - Sx*Sy)/delta;
  *b = (float) (Sy*Sxx - Sx*Sxy)/delta;

  if (r != NULL) {
    *r = (double) (Sxy - Sx*Sy/((double) nv)) /
          sqrt((Sxx - Sx*Sx/((double) nv))*(Syy - Sy*Sy/((double) nv)));
    }

  return true;
}

// Check if a point, (pt_x, pt_y), is inside a polygon with n vertices specified by x and y. 

bool point_in_polygon(float *x, float *y, int n, float pt_x, float pt_y) {

  int i, j;
  bool inside = false;

  for (i = 0, j = n - 1; i < n; j = i++) {
    if (( (y[i] >= pt_y) != (y[j] >= pt_y) ) && (pt_x <= (x[j] - x[i]) * (pt_y - y[i]) / (y[j] - y[i]) + x[i])) {
      inside = !inside;
    }
  }

  return inside;
}


int compare_ints(const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}


int compare_floats(const void *a, const void *b) {
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}


int random_in_range(unsigned int min, unsigned int max)                                            
{                                                                                                  
  int base_random, range, remainder, bucket;                                                       
                                                                                                   
  base_random = rand(); /* in [0, RAND_MAX] */                                                     
  if (RAND_MAX == base_random) return random_in_range(min, max);                                   
  /* now guaranteed to be in [0, RAND_MAX) */                                                      
  range     = max - min;                                                                           
  remainder = RAND_MAX % range;                                                                    
  bucket    = RAND_MAX / range;                                                                    
  /* There are range buckets, plus one smaller interval                                            
     within remainder of RAND_MAX */                                                               
  if (base_random < RAND_MAX - remainder)                                                          
    return min + base_random/bucket;                                                               
  else                                                                                             
    return random_in_range(min, max);                                                              
}  

bool overlap(double a_x0, double a_x1, double a_y0, double a_y1,
       double b_x0, double b_x1, double b_y0, double b_y1) {
  return a_x0 < b_x1 && a_x1 > b_x0 && a_y0 < b_y1 && a_y1 > b_y0;
}

void allocate_int_array_2d(int ***grid, int nrows, int ncols) {
  int n;

  (*grid) = (int**) malloc(nrows*sizeof(int*));
  MEMCHECK((*grid));

  for (n = 0; n < nrows; n++) {
    (*grid)[n] = (int*) malloc(ncols*sizeof(int));
    MEMCHECK((*grid)[n]);
  }
}

void allocate_float_array_2d(float ***grid, int nrows, int ncols) {
  int n;

  (*grid) = (float**) malloc(nrows*sizeof(float*));
  MEMCHECK((*grid));

  for (n = 0; n < nrows; n++) {
    (*grid)[n] = (float*) malloc(ncols*sizeof(float));
    MEMCHECK((*grid)[n]);
  }
}

void allocate_double_array_2d(double ***grid, int nrows, int ncols) {
  int n;

  (*grid) = (double**) malloc(nrows*sizeof(double*));
  MEMCHECK((*grid));

  for (n = 0; n < nrows; n++) {
    (*grid)[n] = (double*) malloc(ncols*sizeof(double));
    MEMCHECK((*grid)[n]);
  }
}

void allocate_float_array_3d(float ****grid, int d0, int d1, int d2) {
  int n, m;

  (*grid) = (float***) malloc(d0*sizeof(float**));
  MEMCHECK((*grid));

  for (n = 0; n < d0; n++) {
    (*grid)[n] = (float**) malloc(d1*sizeof(float*));
    for (m = 0; m < d1; m++) {
      (*grid)[n][m] = (float*) malloc(d2*sizeof(float));
      MEMCHECK((*grid)[n][m]);
    }
  }
}

void allocate_double_array_3d(double ****grid, int d0, int d1, int d2) {
  int n, m;

  (*grid) = (double***) malloc(d0*sizeof(double**));
  MEMCHECK((*grid));

  for (n = 0; n < d0; n++) {
    (*grid)[n] = (double**) malloc(d1*sizeof(double*));
    for (m = 0; m < d1; m++) {
      (*grid)[n][m] = (double*) malloc(d2*sizeof(double));
      MEMCHECK((*grid)[n][m]);
    }
  }
}

void free_double_array_3d(double ***grid, int d0, int d1) {
  int n, m;

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      free(grid[n][m]);
    }
    free(grid[n]);
  }
  free(grid);
}

void free_float_array_3d(float ***grid, int d0, int d1) {
  int n, m;

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      free(grid[n][m]);
    }
    free(grid[n]);
  }
  free(grid);
}

void free_double_array_2d(double **grid, int nrows) {
  int n;
  for (n = 0; n < nrows; n++) 
    free(grid[n]);
  free(grid);
}

void free_int_array_2d(int **grid, int nrows) {
  int n;
  for (n = 0; n < nrows; n++) 
    free(grid[n]);
  free(grid);
}

void free_float_array_2d(float **grid, int nrows) {
  int n;
  for (n = 0; n < nrows; n++) 
    free(grid[n]);
  free(grid);
}

void copy_vec(float *source, float *dest, int len) {
  int n;

  for (n = 0; n < len; n++)
      dest[n] = source[n];
}

void copy_double_vec(double *source, double *dest, int len) {
  int n;

  for (n = 0; n < len; n++)
      dest[n] = source[n];
}

void copy_float_array_2d(float **src, float **dest, int nrows, int ncols) {
  int i, j; 

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      dest[i][j] = src[i][j];
    }
  }
}

void copy_float_array_3d(float ***src, float ***dest, int d0, int d1, int d2) {
  int i, j, k; 

  for (i = 0; i < d0; i++) {
    for (j = 0; j < d1; j++) {
      for (k = 0; k < d2; k++) {
        dest[i][j][k] = src[i][j][k];
      }
    }
  }
}

void copy_double_array_3d(double ***src, double ***dest, int d0, int d1, int d2) {
  int i, j, k; 

  for (i = 0; i < d0; i++) {
    for (j = 0; j < d1; j++) {
      for (k = 0; k < d2; k++) {
        dest[i][j][k] = src[i][j][k];
      }
    }
  }
}

int file_exists(const char *filename) {
  struct stat st;
  int result = stat(filename, &st);
  return result == 0;
}  

char* trim(char *str) {
  size_t len = 0;
  char *frontp = str;
  char *endp = NULL;

  if( str == NULL ) { return NULL; }
  if( str[0] == '\0' ) { return str; }

  len = strlen(str);
  endp = str + len;

  /* Move the front and back pointers to address the first non-whitespace
   * characters from each end.
   */
  while( isspace(*frontp) ) { ++frontp; }
  if( endp != frontp )
    {
      while( isspace(*(--endp)) && endp != frontp ) {}
    }

  if( str + len - 1 != endp )
    *(endp + 1) = '\0';
  else if( frontp != str &&  endp == frontp )
    *str = '\0';

  /* Shift the string so that it starts at str so that if it's dynamically
   * allocated, we can still free it on the returned pointer.  Note the reuse
   * of endp to mean the front of the string buffer now.
   */
  endp = str;
  if( frontp != str )
    {
      while( *frontp ) { *endp++ = *frontp++; }
      *endp = '\0';
    }
  
  return str;
}


float vec_min(float *vec, int len) {
  int n; 
  float min;

  min = vec[0];

  for (n = 0; n < len; n++) {
    if (vec[n] < min)
      min = vec[n];
  }

  return min;
}

double vec_min_double(double *vec, int len) {
  int n; 
  double min;

  min = vec[0];

  for (n = 0; n < len; n++) {
    if (vec[n] < min)
      min = vec[n];
  }

  return min;
}

float vec_min2(float *vec, int len, float spval) {
  int n; 
  float min;

  min = BIG;

  for (n = 0; n < len; n++) {
    if (! approx_equal(vec[n], spval, 1.0e-4) && vec[n] < min)
      min = vec[n];
  }

  return min;
}

float vec_max(float *vec, int len) {
  int n; 
  float max;

  max = vec[0];

  for (n = 0; n < len; n++) {
    if (vec[n] > max)
      max = vec[n];
  }

  return max;
}

double vec_max_double(double *vec, int len) {
  int n; 
  double max;

  max = vec[0];

  for (n = 0; n < len; n++) {
    if (vec[n] > max)
      max = vec[n];
  }

  return max;
}

float vec_max2(float *vec, int len, float spval) {
  int n; 
  float max;

  max = -BIG;

  for (n = 0; n < len; n++) {
    if (! approx_equal(vec[n], spval, 1.0e-4) && vec[n] > max)
      max = vec[n];
  }

  return max;
}

double vec_mean_double(double *vec, int len) {
  int n;
  double total = 0.0; 

  for (n = 0; n < len; n++)
    total += vec[n];

  return total/((double) len);
}

float vec_mean(float *vec, int len) {
  int n;
  float total = 0.0; 

  for (n = 0; n < len; n++)
    total += vec[n];

  return total/((float) len);
}

float vec_mean2(float *vec, int len, float spval) {
  int n, count = 0;
  float total = 0.0; 

  for (n = 0; n < len; n++) {
    if (! approx_equal(vec[n], spval, 1.0e-4)) {
      total += vec[n];
      count++;
    }
  }

  if (count == 0) {
    return spval;
  } else {
    return total/((float) count);
  }
}

float vec_mean2_double(double *vec, int len, double spval) {
  int n, count = 0;
  double total = 0.0; 

  for (n = 0; n < len; n++) {
    if (! approx_equal(vec[n], spval, 1.0e-4)) {
      total += vec[n];
      count++;
    }
  }

  if (count == 0) {
    return spval;
  } else {
    return total/((double) count);
  }
}

float vec_total(float *vec, int len) {
  int n;
  float total = 0.0;

  for (n = 0; n < len; n++)
    total += vec[n];

  return total;
}

float vec_stddev(float *vec, int len, float spval, float mean) {

  int n, count = 0;
  float sumdev = 0.0;

  for (n = 0; n < len; n++) {
    if (! approx_equal(vec[n], spval, 1.0e-4)) {
      sumdev += pow(vec[n] - mean, 2); 
      count++;
    }
  }  

  if (count == 0) {
    return spval; 
  } else {
    return sqrt(sumdev/((float) count));
  }
}

double vec_stddev_double(double *vec, int len, double spval, double mean) {

  int n, count = 0;
  double sumdev = 0.0;

  for (n = 0; n < len; n++) {
    if (! approx_equal(vec[n], spval, 1.0e-4)) {
      sumdev += pow(vec[n] - mean, 2); 
      count++;
    }
  }  

  if (count == 0) {
    return spval; 
  } else {
    return sqrt(sumdev/((double) count));
  }
}


float array_min(float **array, int nrows, int ncols) {
  int i, j;
  float min;

  min = array[0][0];

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (array[i][j] < min)
        min = array[i][j];
    }
  }

  return min;
}

double array_max_double2(double **array, int nrows, int ncols, double spval) {
  int i, j;
  double max;

  max = -BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (! approx_equal(array[i][j], spval, 1.0e-4) && array[i][j] > max)
        max = array[i][j];
    }
  }

  if (max == -BIG)
    return spval;
  else
    return max;
}

double array_min_double2(double **array, int nrows, int ncols, double spval) {
  int i, j;
  double min;

  min = BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (! approx_equal(array[i][j], spval, 1.0e-4) && array[i][j] < min)
        min = array[i][j];
    }
  }

  if (min == BIG)
    return spval;
  else
    return min;
}


float array_max(float **array, int nrows, int ncols) {
  int i, j;
  float max;

  max = array[0][0];

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (array[i][j] > max)
        max = array[i][j];
    }
  }

  return max;
}

float array_mean(float **array, int nrows, int ncols, float spval) {

  int i, j;
  float n = 0.0, tot = 0.0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4)) {
        continue ;
      }
      tot += array[i][j];
      n++;
    }
  }  

  if (n == 0.0)
    return 0.0;
  else
    return tot/n;
}

float array_mean_on_grid(float **array, int nrows, int ncols, float spval,
    float **mask, float mask_spval) {

  int i, j;
  float n = 0.0, tot = 0.0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4) || approx_equal(mask[i][j], mask_spval, 1.0e-4)) {
        continue ;
      }
      tot += array[i][j];
      n++;
    }
  }  

  if (n == 0.0)
    return 0.0;
  else
    return tot/n;
}


float Lsm_shallow(float **array, int nrows, int ncols, float arr_spval, 
    float **shallow, float shallow_spval, float **land, float land_spval) {
  int i, j;
  float min;

  min = BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (array[i][j] == arr_spval || shallow[i][j] == shallow_spval || 
          land[i][j] == land_spval)
        continue ;
      if (array[i][j] < min)
        min = array[i][j];
    }
  }

  return min;
}


float Lsb_shallow(float **array, int nrows, int ncols, float arr_spval, 
    float **shallow, float shallow_spval, float **land, float land_spval, float spread, float *sigma) {
  int i, j;
  float min, max, nobs, mean, stddev;

  min = BIG;
  max = -BIG;
  nobs = 0.0;
  mean = 0.0;
  stddev = 0.0;

// Estimate Lsb from near water points (hopefully sandy beach points). 

  // row-wise

  for (i = 0; i < nrows; i++) {
    // Water south of land.
    for (j = 0; j < ncols - 1; j++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i][j + 1] != land_spval) {
        // min
        if (array[i][j + 1] < min) {
          min = array[i][j + 1];
        }
        // max
        if (array[i][j + 1] > max) {
          max = array[i][j + 1];
        }
        // mean
        mean += array[i][j + 1];
        nobs++;
      }
    }
    // Water north of land. 
    for (j = 1; j < ncols; j++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i][j - 1] != land_spval) {
        // min
        if (array[i][j - 1] < min) {
          min = array[i][j - 1];
        }
        // max
        if (array[i][j - 1] > max) {
          max = array[i][j - 1];
        }
        // mean
        mean += array[i][j - 1];
        nobs++;
      }
    }
  }

  // column-wise

  for (j = 0; j < ncols; j++) {
    // Water west of land.
    for (i = 0; i < nrows - 1; i++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i + 1][j] != land_spval) {
        // min
        if (array[i + 1][j] < min) {
          min = array[i + 1][j];
        }
        // max
        if (array[i + 1][j] > max) {
          max = array[i + 1][j];
        }
        // mean
        mean += array[i + 1][j];
        nobs++;
      }
    }

    // Water east of land.
    for (i = 1; i < nrows; i++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i - 1][j] != land_spval) {
        // min
        if (array[i - 1][j] < min) {
          min = array[i - 1][j];
        }
        // max
        if (array[i - 1][j] > max) {
          max = array[i - 1][j];
        }
        // mean
        mean += array[i - 1][j];
        nobs++; 
      }
    }
  }

  mean /= nobs;

#if 0
  printf("\nLsb:");
  printf("\n      nobs   = %10.1f", nobs);
  printf("\n      min    = %7.2f", min);
  printf("\n      max    = %7.2f", max);
  printf("\n      mean   = %7.2f", mean);
#endif

  // Second pass - compute standard deviation. 

  // row-wise

  for (i = 0; i < nrows; i++) {
    // Water south of land.
    for (j = 0; j < ncols - 1; j++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i][j + 1] != land_spval) {
        stddev += pow(array[i][j + 1] - mean, 2);
      }
    }
    // Water north of land. 
    for (j = 1; j < ncols; j++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i][j - 1] != land_spval) {
        stddev += pow(array[i][j - 1] - mean, 2);
      }
    }
  }

  // column-wise

  for (j = 0; j < ncols; j++) {
    // Water west of land.
    for (i = 0; i < nrows - 1; i++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i + 1][j] != land_spval) {
        stddev += pow(array[i + 1][j] - mean, 2);
      }
    }

    // Water east of land.
    for (i = 1; i < nrows; i++) {
      if (shallow[i][j] != shallow_spval && land[i][j] == land_spval && land[i - 1][j] != land_spval) {
        stddev += pow(array[i - 1][j] - mean, 2);
      }
    }
  }

  stddev = sqrt(stddev/nobs);
  *sigma = stddev;

  printf("\n    sigma_i   = %7.2f", stddev);

  return mean + spread*stddev;
}


float array_min2(float **array, int nrows, int ncols, float spval) {
  int i, j;
  float min;

  min = BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4)) {
        continue ;
      }
      if (array[i][j] < min) {
        min = array[i][j];
      }
    }
  }

  return min;
}

float array_max2(float **array, int nrows, int ncols, float spval) {
  int i, j;
  float max;

  max = -BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4)) {
        continue ;
      }
      if (array[i][j] > max) {
        max = array[i][j];
      }
    }
  }

  return max;
}


void array_arg_min(float **array, int nrows, int ncols, float spval, int *i_min, int *j_min) {
  int i, j;
  float min;

  min = BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4)) {
        continue ;
      }
      if (array[i][j] < min) {
        min = array[i][j];
        *i_min = i;
        *j_min = j;
      }
    }
  }

}


void array_arg_max(float **array, int nrows, int ncols, float spval, int *i_max, int *j_max) {
  int i, j;
  float max;

  max = -BIG;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4)) {
        continue ;
      }
      if (array[i][j] > max) {
        max = array[i][j];
        *i_max = i;
        *j_max = j;
      }
    }
  }

}

float array_stddev(float **array, int nrows, int ncols, float spval, float mean) {

  int i, j;
  float n = 0.0, sumdev = 0.0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (approx_equal(array[i][j], spval, 1.0e-4)) {
        continue ;
      }
      sumdev += pow(array[i][j] - mean, 2); 
      n++;
    }
  }  

  return sqrt(sumdev/n);
}

float array_sample_min(float **array, int nrows, int ncols, float spval, int ntrials) {

  int i, j, n;
  float min = BIG;

  for (n = 0; n < ntrials; n++) {
    i = random_in_range(0, nrows);
    j = random_in_range(0, ncols);
    if (approx_equal(array[i][j], spval, 1.0e-4)) {
      continue ;
    }
    if (array[i][j] < min)
      min = array[i][j];
  }

  return min;
}


float array_sample_max(float **array, int nrows, int ncols, float spval, int ntrials) {
 
  int i, j, n;
  float max = -BIG;

  for (n = 0; n < ntrials; n++) {
    i = random_in_range(0, nrows);
    j = random_in_range(0, ncols);
    if (approx_equal(array[i][j], spval, 1.0e-4)) {
      continue ;
    }
    if (array[i][j] > max)
      max = array[i][j];
  }

  return max; 

}

void array_sample_arg_min(float **array, int nrows, int ncols, float spval, int *i_min, int *j_min, int ntrials) {

  int i, j, n;
  float min = BIG;

  for (n = 0; n < ntrials; n++) {
    i = random_in_range(0, nrows);
    j = random_in_range(0, ncols);
    if (approx_equal(array[i][j], spval, 1.0e-4)) {
      continue ;
    }
    if (array[i][j] < min) {
      min = array[i][j];
      *i_min = i;
      *j_min = j;
    }
  }

}

void array_sample_arg_max(float **array, int nrows, int ncols, float spval, int *i_max, int *j_max, int ntrials) {
 
  int i, j, n;
  float max = -BIG;

  for (n = 0; n < ntrials; n++) {
    i = random_in_range(0, nrows);
    j = random_in_range(0, ncols);
    if (approx_equal(array[i][j], spval, 1.0e-4)) {
      continue ;
    }
    if (array[i][j] > max) {
      max = array[i][j];
      *i_max = i;
      *j_max = j;
    }
  }

}




//
//    READ (NEW/OLD) GRID FROM COMMAND LINE
//

bool parse_old_grid(int *ps_in, int k) {

  int n;
  bool present = false; 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[k]) == 0) {
      *ps_in = n;
      present = true;
      break ;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[k]);
    return false;
  }

  return true;
}



bool parse_new_grid(int *ps_out, int k) {

  int n;

// Check output grid name is not in use. 

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[k]) == 0) {
      printf("\nERROR: grid '%s' already in use.\n", parsed[k]);
      return false;
    }
  }

// Allocate grid into the first empty slot. 

  for (n = 0; n < MAX_GRIDS; n++) {
   if (! allocated_grids[n]) {
    *ps_out = n;
    break;
    }
  }

// Copy grid name. 

  strcpy(grid_names[*ps_out], parsed[k]);

// This grid is now in use. 

  allocated_grids[*ps_out] = true;

// Update number of grids in use. 

  n_grids_in_use++;

  return true;
}

//
//    READ (NEW/OLD) SCENE FROM COMMAND LINE
//

bool parse_old_scene(int *ps_in, int k) {

  int n;
  bool present = false; 

  for (n = 0; n < MAX_SCENES; n++) {
    if (strcmp(scene_names[n], parsed[k]) == 0) {
      *ps_in = n;
      present = true;
      break ;
    }
  }

  if (! present) {
    printf("\nERROR: unknown scene '%s'.\n", parsed[k]);
    return false;
  }

  return true;
}

bool parse_new_scene(int *ps_out, int k) {

  int n;

// Check output scene name is not in use. 

  for (n = 0; n < MAX_SCENES; n++) {
    if (strcmp(scene_names[n], parsed[k]) == 0) {
      printf("\nERROR: scene '%s' already in use.\n", parsed[k]);
      return false;
    }
  }

// Allocate scene into the first empty slot. 

  for (n = 0; n < MAX_SCENES; n++) {
   if (! allocated_scenes[n]) {
    *ps_out = n;
    break;
    }
  }

// Copy scene name. 

  strcpy(scene_names[*ps_out], parsed[k]);

// This scene is now in use. 

  allocated_scenes[*ps_out] = true;

// Update number of scenes in use. 

  n_scenes_in_use++;

  return true;
}

