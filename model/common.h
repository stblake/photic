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

// bool is an int (not _Bool): some code stores non-0/1 values in it. C23 makes
// bool a keyword, so compile with -std=gnu17 or earlier.
#if defined(__STDC_VERSION__) && __STDC_VERSION__ > 201710L
#error "photic must be compiled as C17 or earlier (e.g. -std=gnu17)."
#endif
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
	double *Y_scene; // Fitted backscattering slope per scene (Y_MODEL 3).
	double **R_inf;  // Deep-water Rrs per scene and band (SCENE RINF), 0 if not given.
	double *prev;
	double ***Rrs_modelled; 
	double ***rrs_modelled;
	double ***rrs_bottom;
	double ***rrs_dp;
	double *iop_cache;    // Water optics per scene and band for the last P, G, X (samodel_Rrs), or NULL.
	double **K;
	double ***rho;
	double depth;
	double h0;            // Physics-based first estimate of depth (samodel_h0), -1 if none.
	bool optically_deep;  // The bottom was not detected (EXPERIMENTS.md C2).
	double deep_chi2;     // Fit improvement of the bottom over an optically deep model.
	double ident_chi2;    // Misfit gap to the best clearly different depth (EXPERIMENTS.md B20).
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
bool grid_matches(int n, const char *name);
bool grid_name_taken(int n, const char *name);
void release_orphan_grids();

#if PLOTTING
extern int n_fn_evals;
extern double Error_path[6000];
extern double P_path[6000]; 
extern double G_path[6000]; 
extern double X_path[6000];
extern double B_path[6000];
extern double H_path[6000];
extern double Delta_path[6000];
#endif


// Globals (defined in common.c). 
extern int n_grids_in_use;
extern int n_images_in_use;
extern float meminuse; 

extern char command_file[MAX_STRING_LEN];
extern bool allocated_grids[MAX_GRIDS];
extern bool allocated_images[MAX_IMAGES];
extern char grid_names[MAX_GRIDS][MAX_STRING_LEN];
extern char image_names[MAX_IMAGES][MAX_STRING_LEN];
extern geogrid gridded_data[MAX_GRIDS];
extern image image_data[MAX_IMAGES];
extern line_segment line_data[MAX_GRIDS];
extern char parsed[MAX_ARGS][MAX_ARG_SIZE];


extern scene scene_data[MAX_SCENES];
extern bool allocated_scenes[MAX_SCENES];
extern char scene_names[MAX_SCENES][MAX_STRING_LEN];
extern int n_scenes_in_use;

extern char prev1[MAX_STRING_LEN];
extern char prev2[MAX_STRING_LEN];
extern char prev3[MAX_STRING_LEN];
extern char prev4[MAX_STRING_LEN];
extern char prev5[MAX_STRING_LEN];
extern char prev6[MAX_STRING_LEN];
extern char prev7[MAX_STRING_LEN];
extern char prev8[MAX_STRING_LEN];
extern char prev9[MAX_STRING_LEN];

// Depth soundings. 

extern bool read_soundings;
extern bool first_model_error_soundings;
extern int nsoundings;
extern float **soundings;
extern float *model_depths, *sounding_depths, *sounding_errors, *sounding_weights, **radiance_at_soundings, 
  *dz_at_soundings, *lsb_at_soundings;

// Bottom spectra. 

extern float **bottom_spectra;
extern int n_bottom_spectra;

// Graphics defaults. 

extern int background, palette;
extern float pagesize, textsize, pointsize;
extern int line_width;

// Satellite metadata. 

extern float sun_azimuth;
extern float sun_elevation;
extern float earth_sun_distance;
extern float radiance_mult_band[MAX_GRIDS];
extern float radiance_add_band[MAX_GRIDS];
extern float quantize_cal_min_band[MAX_GRIDS];
extern float quantize_cal_max_band[MAX_GRIDS];
extern float radiance_minimum_band[MAX_GRIDS];
extern float radiance_maximum_band[MAX_GRIDS];

// Model defaults.

#define LINEAR 0
#define LSB_SIGMA 0

extern float H0;

extern bool relative_error;

extern float (*zdepth)(int, int);

extern char model_name[MAX_STRING_LEN];

extern bool pan_normalise;

extern float extent_eps;

extern float epsilon;
extern float log_epsilon; // 1.0 
extern float falloff_param; // units km.
extern float model_weight_cutoff;
extern float map_scale;

extern int nboundary;

extern float deep_lon;
extern float deep_lat;
extern float sand_lon;
extern float sand_lat;

extern int n_local_models;
extern int depth_radius;
extern int n_depth_bins;
extern float lsb_spread;

extern float sounding_error;
extern float constraint_error;
extern float local_constraint_error;

extern float noise_threshold;
extern float min_radiance_threshold;
extern float max_radiance_threshold;

extern int bands_indexes[MAX_GRIDS];
extern int spectral_indexes[MAX_GRIDS];
extern int land_index;
extern int shallow_index;

extern int n_model_bands, **model_origins, NAlpha[MAX_GRIDS]; 
extern float thetas[MAX_GRIDS], model_wavelengths[MAX_GRIDS], Lsm[MAX_GRIDS], Lsb[MAX_GRIDS], Lsb_sigma[MAX_GRIDS], 
  Alpha[MAX_GRIDS], Alpha_sigma[MAX_GRIDS], Rrs_sigma[MAX_GRIDS], Hj[MAX_GRIDS], Lsm_sigma[MAX_GRIDS], *mparams;

extern float jerlov_water_type; 
extern int n_water_type_estimates;

extern float K_ratios[MAX_GRIDS][4];

extern float slope_line[LINE_MAX_POINTS];

extern float depthmin;
extern float depthmax;

extern float line_xa, line_ya, line_xb, line_yb;

extern float bound_error, coeff_size_preference, regression_error;

extern float weight_soundings, weight_constraints, weight_local_constraints, 
  weight_coeff_size_preference, weight_bound_error, weight_regression;

extern float rms_error;

extern int model_optimisation;

extern int samodel_n_smoothing_radius, samodel_n_spatial, samodel_n_bottoms;



#endif

