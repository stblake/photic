
#include "common.h"

void nelmin ( double fn ( double x[], model_data *md ), model_data *md, 
	      int n, double start[], double xmin[], 
	      double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
	      int *icount, int *numres, int *ifault );

void nelmin2 ( double fn ( double x[], double *z_model, double *z_soundings, double *sounding_weights, int n_soundings, double spval),
	       double *z_model, double *z_soundings, double *sounding_weights, int n_soundings, double spval, 
	      int n, double start[], double xmin[],
	      double *ynewlo, double reqmin, double step[], int konvge, int kcount,
	      int *icount, int *numres, int *ifault );

void nelmin3 ( double fn ( double x[], double **z_model, double *z_soundings, double *sounding_weights, int n_grids, int n_soundings, double spval), 
  double **z_model, double *z_soundings, double *sounding_weights, int n_grids, int n_soundings, double spval, 
  int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault );

void timestamp ( );
