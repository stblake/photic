

#include "common.h"
#include "io.h"
#include "interp.h"
#include "asa047.h"

void run_align();
double error_aligned(double *params, double *z_model, double *z_soundings, double *sounding_weights, int n_soundings, double spval);
double error_weighted_sum(double *params, double **z_model, double *z_soundings, double *sounding_weights, int n_grids, int n_soundings, double spval);
void align_compute_weights(double **weights, double *soundings, int nsoundings);