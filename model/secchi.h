

#include "common.h"

void Kd_LS8(
	float **coastal, float **blue, float **green, float **red, float **kd,
	int nrows, int ncols, 
	float coastal_spv, float blue_spv, float green_spv, float red_spv, float theta_s);

float Lee_kd_LS8(float coastal, float blue, float green, float red, float theta_s);

void secchi_disk_depth(
	float **coastal, float **blue, float **green, float **red, float **zsd,
	int nrows, int ncols, 
	float coastal_spv, float blue_spv, float green_spv, float red_spv, float theta_s);

void Lee_Kd_LS8(float Rrs443, float Rrs481, float Rrs554, float Rrs656, float theta_s, 
	float *Kd443, float *Kd481, float *Kd554, float *Kd656);

float Lee_secchi_LS8(float coastal, float blue, float green, float red, float theta_s);