#ifndef _COMMONH_
#define _COMMONH_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <float.h>
#include <limits.h>
// #include <omp.h>

#define PGPLOT 1
#define PLOTTING 0

#define PI 3.141592653589793

#define BIG 1.0e10

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define JET 101
#define GRAYSCALE 102
#define OCEAN 103
#define GMT 104
#define INVERSEJET 105
#define RYGB 106
#define WOR 108
#define VIRIDIS 109
#define PLASMA 110

#define LYZENGALOGLINEAR1985 2001
#define LYZENGALOGRATIO1978 2002
#define STUMPF2003 2003
#define SPLINE 2004

#define SIMULATED_ANNEALING 2011
#define RANDOM_SEARCH 2012

#define BAM_NODATA -99999.0

#define MAX_SCENES 32
#define MAX_STRING_LEN 2048 
#define MAX_GRIDS 265
#define MAX_IMAGES 16
#define MAX_LINE 2048
#define MAX_ARGS 128
#define MAX_ARG_SIZE 128 
#define LINE_MAX_POINTS 256

#define BLACK 10
#define WHITE 11

typedef int bool;
#define true 1
#define false 0

#define MEMCHECK(p) \
  if (p == NULL)    \
    {\
      printf("\n\nERROR: Out of memory.\n\n");\
      exit(1);\
    }

struct grid_struct {
	int nrows; 
	int ncols;
	float cellsize; 
	float wlon; 
	float slat; 
	float elon; 
	float nlat;
	float nodata_value;
	float lambda;
	float theta_v;
	float theta_w;
	float **array;
};

typedef struct grid_struct geogrid;



struct line_segment_struct {
	char *name;
	int npoints;
	float *x;
	float *y;
};

typedef struct line_segment_struct line_segment;


struct image_struct {
	int ncols;
	int nrows;
	float cellsize; 
	float wlon; 
	float slat; 
	float elon; 
	float nlat;
	float nodata_value;
	unsigned short int **red;
	unsigned short int **green;
	unsigned short int **blue;
};

typedef struct image_struct image;



struct model_data_struct {
	int n_bottoms;
	int n_scenes;
	int n_regions;
	int max_n_bands;
	int max_n_regions;
	int *n_bands;
	int *n_raw_bands;
	int origin;
	int n_params;
	// inputs
	double ***Rrs_measured;
	double **wavelengths;
	double **raw_wavelengths;
	double ***bottom_reflectance;
	double **a_0; 
	double **a_1; 
	double **a_w; 
	double **b_bw;
	double a_w640;
	double *H_tide;
	double *theta_sun;
	double *theta_view;
	double *sec_theta_sun;
	double *sec_theta_view;
	double **Rrs440;
	double **Rrs490;
	double **Rrs550;
	double **Rrs640;
	double **Rrs750;
	double **rrs_noise;
	double h_empirical;
	int bottom_type_indexes[12];
	char B_type_name[12][MAX_STRING_LEN];
	bool empirical_depth_present;
	// restart
	int i;
	int j;
	int prev_i;
	int prev_j;
	bool start_at_previous;
	// outputs
	double *P;
	double *G;
	double *X;
	double *B;
	double *H;
	double *D;
	double *prev;
	double ***Rrs_modelled; 
	double ***rrs_modelled;
	double ***rrs_bottom;
	double ***rrs_dp;
	double **K;
	double ***rho;
	double depth;
	double depth_prev;
	double depth_sigma;
	double Rrs_error;
	double SAM_error;
	double K_error;
	double depth_error;
	double bottom_error;
	double model_error;
	double K_min;
	double K_min_prev;
	double bottom_albedo;
	double B_type_percent[12];
	double index_optical_depth;
	int bottom_type;
	bool converged; 
	int n_iterations;
};

typedef struct model_data_struct model_data;



struct scene_struct {
	char scene_name[MAX_STRING_LEN];
	int n_bands;
	int nrows;
	int ncols;
	int band_indexes[MAX_GRIDS];
	int wavelengths[MAX_GRIDS];
	double theta_w; 
	double theta_v;
	double H_tide;
	double R_inf[MAX_GRIDS];
	double R_sigma[MAX_GRIDS];
	double K[MAX_GRIDS];
	double K_sigma[MAX_GRIDS];
	double ratio_min;
	double ratio_max; 
	double slope_min;
	double slope_max;
  bool pgx_present;
  int ps_p;
  int ps_g; 
  int ps_x;
};

typedef struct scene_struct scene;


struct sort_index_float_struct
{
  float value;
  int indx;
};

void sort_indexes_2d(float **array, int nrows, int ncols, int **indexes);
float smoothed_index_2d(float **array, int i, int j, int nrows, int ncols, int smoothing_radius, float spval);
double interp_1d(double *X, double *Y, int n, double x);
float interp_float_1d(float *X, float *Y, int n, float x);
bool approx_equal(float a, float b, float epsilon);
bool approx_equal_double(double a, double b, double epsilon);
bool approx_less_equal(float a, float b, float epsilon);
bool approx_greater_equal(float a, float b, float epsilon);
bool linear_fit(float *x, float *y, int n, float *m, float *b, float *r);
bool linear_fit2(float *x, float *y, int n, float *m, float *b, float *r, float spval);
int random_in_range(unsigned int min, unsigned int max);
char* trim(char *str);
void allocate_int_array_2d(int ***grid, int nrows, int ncols);
void allocate_float_array_2d(float ***grid, int nrows, int ncols);
void free_int_array_2d(int **grid, int nrows);
void free_float_array_2d(float **grid, int nrows);
int file_exists(const char *filename);
bool overlap(double a_x0, double a_x1, double a_y0, double a_y1,
	     double b_x0, double b_x1, double b_y0, double b_y1);
void copy_float_array_2d(float **src, float **dest, int nrows, int ncols);
void copy_float_array_3d(float ***src, float ***dest, int d0, int d1, int d2);
void copy_double_array_3d(double ***src, double ***dest, int d0, int d1, int d2);
float vec_min(float *vec, int len);
double vec_min_double(double *vec, int len);
float vec_max(float *vec, int len);
double vec_max_double(double *vec, int len);
float vec_mean(float *vec, int len);
float array_min(float **array, int nrows, int ncols);
float array_max(float **array, int nrows, int ncols);
float array_min2(float **array, int nrows, int ncols, float spval);
float array_max2(float **array, int nrows, int ncols, float spval);
float vec_total(float *vec, int len);
void copy_vec(float *source, float *dest, int len);
void copy_double_vec(double *source, double *dest, int len);
void array_arg_min(float **array, int nrows, int ncols, float spval, int *i_min, int *j_min);
void array_arg_max(float **array, int nrows, int ncols, float spval, int *i_max, int *j_max);
float array_mean(float **array, int nrows, int ncols, float spval);
float array_stddev(float **array, int nrows, int ncols, float spval, float mean);
float array_sample_min(float **array, int nrows, int ncols, float spval, int ntrials);
float array_sample_max(float **array, int nrows, int ncols, float spval, int ntrials);
void array_sample_arg_min(float **array, int nrows, int ncols, float spval, int *i_min, int *j_min, int ntrials);
void array_sample_arg_max(float **array, int nrows, int ncols, float spval, int *i_max, int *j_max, int ntrials);
float Lsm_shallow(float **array, int nrows, int ncols, float arr_spval, 
    float **shallow, float shallow_spval, float **land, float land_spval);
float array_mean_on_grid(float **array, int nrows, int ncols, float spval,
    float **mask, float mask_spval);
float Lsb_shallow(float **array, int nrows, int ncols, float arr_spval, 
    float **shallow, float shallow_spval, float **land, float land_spval, float spread, float *sigma);
int compare_floats(const void *a, const void *b);
int compare_ints(const void * a, const void * b);
bool point_in_polygon(float *x, float *y, int n, float pt_x, float pt_y);
float vec_mean2(float *vec, int len, float spval);
float vec_stddev(float *vec, int len, float spval, float mean);
float vec_min2(float *vec, int len, float spval);
float vec_max2(float *vec, int len, float spval);
bool all_approx_zero(float *x, int n, float epsilon);
double random_norm(double mean, double stddev);
double vec_stddev_double(double *vec, int len, double spval, double mean);
double vec_mean_double(double *vec, int len);
double array_min_double2(double **array, int nrows, int ncols, double spval);
double array_max_double2(double **array, int nrows, int ncols, double spval);
void allocate_double_array_2d(double ***grid, int nrows, int ncols);
void free_double_array_2d(double **grid, int nrows);
bool parse_old_grid(int *ps, int k);
bool parse_new_grid(int *ps_out, int k);
bool parse_new_scene(int *ps_out, int k);
bool parse_old_scene(int *ps_in, int k);
void allocate_float_array_3d(float ****grid, int d0, int d1, int d2);
void free_float_array_3d(float ***grid, int d0, int d1);
void allocate_double_array_3d(double ****grid, int d0, int d1, int d2);
void free_double_array_3d(double ***grid, int d0, int d1);
float frand();
float frand2(float max);
float vec_mean2_double(double *vec, int len, double spval);
bool commensurate_grids(geogrid a, geogrid b);

#if PLOTTING
int n_fn_evals;
double Error_path[6000];
double P_path[6000]; 
double G_path[6000]; 
double X_path[6000];
double B_path[6000];
double H_path[6000];
double Delta_path[6000];
#endif


// Globals. 
int n_grids_in_use;
int n_images_in_use;
float meminuse; 

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
int n_scenes_in_use;

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

bool read_soundings;
bool first_model_error_soundings;
int nsoundings;
float **soundings;
float *model_depths, *sounding_depths, *sounding_errors, *sounding_weights, **radiance_at_soundings, 
  *dz_at_soundings, *lsb_at_soundings;

// Bottom spectra. 

float **bottom_spectra;
int n_bottom_spectra;

// Graphics defaults. 

int background, palette;
float pagesize, textsize, pointsize;
int line_width;

// Satellite metadata. 

float sun_azimuth;
float sun_elevation;
float earth_sun_distance;
float radiance_mult_band[MAX_GRIDS];
float radiance_add_band[MAX_GRIDS];
float quantize_cal_min_band[MAX_GRIDS];
float quantize_cal_max_band[MAX_GRIDS];
float radiance_minimum_band[MAX_GRIDS];
float radiance_maximum_band[MAX_GRIDS];

// Model defaults.

#define LINEAR 0
#define LSB_SIGMA 0

float H0;

bool relative_error;

float (*zdepth)(int, int);

char model_name[MAX_STRING_LEN];

bool pan_normalise;

float extent_eps;

float epsilon;
float log_epsilon; // 1.0 
float falloff_param; // units km.
float model_weight_cutoff;
float map_scale;

int nboundary;

float deep_lon;
float deep_lat;
float sand_lon;
float sand_lat;

int n_local_models;
int depth_radius;
int n_depth_bins;
float lsb_spread;

float sounding_error;
float constraint_error;
float local_constraint_error;

float noise_threshold;
float min_radiance_threshold;
float max_radiance_threshold;

int bands_indexes[MAX_GRIDS];
int spectral_indexes[MAX_GRIDS];
int land_index;
int shallow_index;

int n_model_bands, **model_origins, NAlpha[MAX_GRIDS]; 
float thetas[MAX_GRIDS], model_wavelengths[MAX_GRIDS], Lsm[MAX_GRIDS], Lsb[MAX_GRIDS], Lsb_sigma[MAX_GRIDS], 
  Alpha[MAX_GRIDS], Alpha_sigma[MAX_GRIDS], Rrs_sigma[MAX_GRIDS], Hj[MAX_GRIDS], Lsm_sigma[MAX_GRIDS], *mparams;

float jerlov_water_type; 
int n_water_type_estimates;

float K_ratios[MAX_GRIDS][4];

float slope_line[LINE_MAX_POINTS];

float depthmin;
float depthmax;

float line_xa, line_ya, line_xb, line_yb;

float bound_error, coeff_size_preference, regression_error;

float weight_soundings, weight_constraints, weight_local_constraints, 
  weight_coeff_size_preference, weight_bound_error, weight_regression;

float rms_error;

int model_optimisation;

int samodel_n_smoothing_radius, samodel_n_spatial, samodel_n_bottoms;



#endif

