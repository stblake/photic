

#include "common.h"

#define NWAVELENS 13
#define NJERLOVTYPES 10


bool jerlov(float wlen_i, float wlen_j, float Lsmi, float Lsmj, 
	float *Li, float *Lj, int npoints, float *ki, float *kj, 
	float *m, float *c, float *r, float *water_type, float manual_ratio);
bool interp_jerlov_wavelength(float jdata[NWAVELENS], float wl);
bool interp_jerlov_ratio(float ratio, float ratios[NWAVELENS], float jerlov_i[NWAVELENS], float *ki);
void jerlov_water_type_str(float wtype, char *str_type, int str_len);
bool estimate_water_type(float ratio, float ratios[NJERLOVTYPES], float *water_type);
void compute_k_from_jerlov(float water_type, float *alphas, int *spectral_indexes, float *wavelengths, int nspec);
float compute_k(float water_type, float wlen);
bool compute_k_from_ratio(float ratio, float wlen_i, float wlen_j, float *water_type, float *k, float *wavelengths, float n_wlens);