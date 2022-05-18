
#include "common.h"
#include "nc.h"
#include "asa047.h"
#include "samodelgraphics.h"


void samodel(
		// inputs
		scene scene_data[], geogrid gridded_data[], int *scene_indexes, int nscenes,
		bool empirical_depth_present, geogrid empirical_depths,  
		// parameterisation
		int n_smoothing_radius, int n_spatial, int n_bottoms,
		// outputs
		float **depth, float **depth_sigma, float **model_error, float **bottom_albedo, 
		float **bottom_sand, float **bottom_seagrass, float **bottom_coral, 
		float **K_min, float **bottom_type, float **index_optical_depth,
		// plotting parameters
		float pagesize, int background, int linewidth);

void samodel_optimise(model_data *md);

void samodel_optimise_one_bottom_combination(model_data *md, double *params);

double samodel_error(double *params, model_data *md);

void samodel_Rrs(
                 double P, double G, double X, double B[], double q[], double H, double Delta, 
		 model_data *md, int k_scene, int k_band, int k_region);

void samodel_graphics(model_data *md, float pagesize, int linewidth);

void extract_Rrs_data(int i, int j, scene scene_data[], geogrid gridded_data[], int *scene_indexes, 
	int nscenes, int n_spatial, int n_smoothing_radius, int nrows, int ncols, model_data *md, 
	float n_sigma);

void polcyn_depth(model_data *md, double *log_ratio_depth, double *log_ratio_depth_error, bool printing);

void samodel_display(model_data *md, double *min);

float model_spectral_diff(double ***Rrs_measured, double ***lut_modelled_Rrs, model_data *md);

float model_spectral_diff_mean(double ***Rrs_measured, double ***lut_modelled_Rrs, model_data *md);

void partial_copy_md(model_data *md, model_data *md_out, bool allocate);
