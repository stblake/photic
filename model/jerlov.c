//
// Two-way interpolation of Jerlov's table of spectral attenuation coefficients
//

// Ref: Jerlov, N, Marine Optics, 2nd Edition, Table XXVII, pp. 135, 1976.


#include "jerlov.h"

float jerlov_data[NJERLOVTYPES][NWAVELENS];

int wavelengths[] = {400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700};

// Original data from Jerlov. 

// float Jerlov_I[] = {0.027, 0.022, 0.019, 0.018, 0.027, 0.043, 0.063, 0.089, 0.235, 0.305, 0.36, 0.42, 0.56};

// float Jerlov_IA[] = {0.038, 0.031, 0.026, 0.025, 0.032, 0.048, 0.067, 0.094, 0.24, 0.31, 0.37, 0.43, 0.57};

// float Jerlov_IB[] = {0.051, 0.042, 0.036, 0.033, 0.042, 0.054, 0.072, 0.099, 0.245, 0.315, 0.375, 0.435, 0.58};

// float Jerlov_II[] = {0.096, 0.081, 0.068, 0.062, 0.07, 0.076, 0.089, 0.115, 0.26, 0.335, 0.4, 0.465, 0.61};

// float Jerlov_III[] = {0.185, 0.16, 0.135, 0.116, 0.115, 0.116, 0.12, 0.148, 0.295, 0.375, 0.445, 0.52, 0.66};

// float Jerlov_C1[] = {0.51, 0.36, 0.25, 0.17, 0.14, 0.13, 0.12, 0.15, 0.3, 0.37, 0.45, 0.51, 0.65};

float Jerlov_C3[] = {0.78, 0.54, 0.39, 0.29, 0.22, 0.2, 0.19, 0.21, 0.33, 0.4, 0.46, 0.56, 0.71};

float Jerlov_C5[] = {1.1, 0.78, 0.56, 0.43, 0.36, 0.31, 0.3, 0.33, 0.4, 0.48, 0.54, 0.65, 0.8};

float Jerlov_C7[] = {1.6, 1.2, 0.89, 0.71, 0.58, 0.49, 0.46, 0.46, 0.48, 0.54, 0.63, 0.78, 0.92};

float Jerlov_C9[] = {2.4, 1.9, 1.6, 1.23, 0.99, 0.78, 0.63, 0.58, 0.6, 0.65, 0.76, 0.92, 1.1};


// Updates from Austin & Petzold, 1984. 

float Jerlov_I[] = {0.0217, 0.0185, 0.0176, 0.0184, 0.0280, 0.0504, 0.0640, 0.0931, 0.2408, 0.3174, 0.3559, 0.4372, 0.6513};

float Jerlov_IA[] = {0.0316, 0.0280, 0.0257, 0.0250, 0.0332, 0.0545, 0.0674, 0.0960, 0.2437, 0.3206, 0.3601, 0.4410, 0.6530};

float Jerlov_IB[] = {0.0438, 0.0395, 0.0355, 0.0330, 0.0396, 0.0596, 0.0715, 0.0995, 0.2471, 0.3245, 0.3652, 0.4457, 0.6550};

float Jerlov_II[] = {0.0878, 0.0814, 0.0714, 0.0620, 0.0627, 0.0779, 0.0863, 0.1122, 0.2595, 0.3389, 0.3837, 0.4626, 0.6623};

float Jerlov_III[] = {0.1697, 0.1594, 0.1381, 0.1160, 0.1056, 0.1120, 0.1139, 0.1359, 0.2826, 0.3655, 0.4181, 0.4942, 0.6760};

float Jerlov_C1[] = {0.2516, 0.2374, 0.2048, 0.1700, 0.1486, 0.1461, 0.1415, 0.1596, 0.3057, 0.3922, 0.4525, 0.5257, 0.6896};

/* This routine takes as input the following:

	wlen_i   - wavelength of the first spectral band. 
	wlen_j   - wavelength of the second spectral band.
	Lsmi    - Deep water signal for band i.
	Lsmj    - Deep water signal for band j.
	Li      - vector of radiances for the first spectral band. This
			  represents a line from optically shallow water to 
			  optically deep water. Ideally over the same bottom type. 
	Lj      - vector of radiances for the second spectral band. 
	npoints - number of points in each line. 

	And returns an approximation to the attenuation coefficients for each spectral band. 

	Additional parameters are: 

	m 			- slope of line of regression line. 
	c 			- y-intercept of line of best fit. 
	water_type 	- Jerlov water type. 
*/


bool jerlov(float wlen_i, float wlen_j, float Lsmi, float Lsmj, 
	float *Li, float *Lj, int npoints,
	float *ki, float *kj, float *m, float *c, float *r, float *water_type, 
	float manual_ratio) {

	int n, i, k, n_too_deep = 0, n_shallow;
	float jerlov_i[NJERLOVTYPES], jerlov_j[NJERLOVTYPES], ratios[NJERLOVTYPES], 
	*Xi, *Xj;
	bool success;
	char str_water_type[32];

// Create table of attenuation coefficients. 

	for (i = 0; i < NWAVELENS; i++) {
		jerlov_data[0][i] = Jerlov_I[i];
		jerlov_data[1][i] = Jerlov_IA[i];
		jerlov_data[2][i] = Jerlov_IB[i];
		jerlov_data[3][i] = Jerlov_II[i];
		jerlov_data[4][i] = Jerlov_III[i];	
		jerlov_data[5][i] = Jerlov_C1[i];
		jerlov_data[6][i] = Jerlov_C3[i];
		jerlov_data[7][i] = Jerlov_C5[i];
		jerlov_data[8][i] = Jerlov_C7[i];
		jerlov_data[9][i] = Jerlov_C9[i];
	}

// Compute Xi = log(Li - Lsmi) and Xj = log(Lj - Lsmj).

	Xi = malloc(npoints*sizeof(float));
	Xj = malloc(npoints*sizeof(float));
	
	k = 0;

	for (n = 0; n < npoints; n++) {
		if (Li[n] > Lsmi + 1.0 && Lj[n] > Lsmj + 1.0) {
			Xi[k] = log(Li[n] - Lsmi);
			Xj[k] = log(Lj[n] - Lsmj);
			k++;
		} else {
			n_too_deep++;
		}
	}

	if (n_too_deep > 0) {
		float n_percent_too_deep = 100.0*((float) n_too_deep)/((float) npoints);
		printf("\nWARNING: %.1f(%%) of points below deep water signal.\n", n_percent_too_deep);
	}

	n_shallow = k;

	printf("\nNumber of optically shallow points = %d\n", n_shallow);

// Compute linear fit to the data. 

	success = linear_fit(Xi, Xj, n_shallow, m, c, r);
	if (! success) {
		free(Xi);
		free(Xj);
		return false;
	}

	printf("\nm, c, r**2 = %.3f, %.3f, %.3f\n", *m, *c, pow(*r, 2));

// Manual set the ratio? 

	if (! approx_equal(manual_ratio, 0.0, 1.0e-4)) {
		*m = manual_ratio;
	}

// Interpolate Jerlov data at the two wavelengths.  

	success = interp_jerlov_wavelength(jerlov_i, wlen_i);
	if (! success) {
		free(Xi);
		free(Xj);
		return false;
	}

	success = interp_jerlov_wavelength(jerlov_j, wlen_j);
	if (! success) {
		free(Xi);
		free(Xj);
		return false;
	}

	printf("\n type         =    OI   OIA   OIB   OII  OIII    C1    C3    C5    C7    C9");
	printf("\n jerlov @ %d = ", ((int) wlen_i));	
	for (i = 0; i < NJERLOVTYPES; i++) 
		printf("%3.3f ", jerlov_i[i]);

	printf("\n jerlov @ %d = ", ((int) wlen_j));
	for (i = 0; i < NJERLOVTYPES; i++) 
		printf("%3.3f ", jerlov_j[i]);

// Compute ratios. 

	for (i = 0; i < NJERLOVTYPES; i++) 
		ratios[i] = jerlov_j[i]/jerlov_i[i];

	printf("\n ratios       = ");
	for (i = 0; i < NJERLOVTYPES; i++) 
		printf("%3.3f ", ratios[i]);
	printf("\n");

// Interpolate from ratios to estimate attenuation coefficients. 

	success = interp_jerlov_ratio(*m, ratios, jerlov_i, ki);
	if (! success) {
		free(Xi);
		free(Xj);
		return false;
	}

	success = interp_jerlov_ratio(*m, ratios, jerlov_j, kj);
	if (! success) {
		free(Xi);
		free(Xj);
		return false;
	}

// Compute water type. 

	success = estimate_water_type(*m, ratios, water_type);
	if (! success) {
		free(Xi);
		free(Xj);
		return false;
	}

	printf("\nWater index = %.3f", *water_type);
	printf("\nKi,Kj = %.3f, %.3f\n", *ki, *kj);


// Free memory. 

	free(Xi);
	free(Xj);

	return true;
}



bool compute_k_from_ratio(float ratio, float wlen_i, float wlen_j, float *water_type, float *k, float *wavelengths, float n_wlens) {

	int i;
	float jerlov_i[NJERLOVTYPES], jerlov_j[NJERLOVTYPES], ratios[NJERLOVTYPES], *Xi, *Xj;
	bool success;

// Interpolate Jerlov data at the two wavelengths.  

	success = interp_jerlov_wavelength(jerlov_i, wlen_i);
	if (! success) {
		return false;
	}

	success = interp_jerlov_wavelength(jerlov_j, wlen_j);
	if (! success) {
		return false;
	}

	printf("\n type         =    OI   OIA   OIB   OII  OIII    C1    C3    C5    C7    C9");
	printf("\n jerlov @ %d = ", ((int) wlen_i));	
	for (i = 0; i < NJERLOVTYPES; i++) {
		printf("%3.3f ", jerlov_i[i]);
	}

	printf("\n jerlov @ %d = ", ((int) wlen_j));
	for (i = 0; i < NJERLOVTYPES; i++) {
		printf("%3.3f ", jerlov_j[i]);
	}

// Compute ratios. 

	for (i = 0; i < NJERLOVTYPES; i++) 
		ratios[i] = jerlov_j[i]/jerlov_i[i];

	printf("\n ratios       = ");
	for (i = 0; i < NJERLOVTYPES; i++) 
		printf("%3.3f ", ratios[i]);
	printf("\n");

// Compute water type. 

	success = estimate_water_type(ratio, ratios, water_type);
	if (! success) {
		return false;
	}

	printf("\nWater index = %.3f", *water_type);

// Compute attenuation coefficients at all wavelengths. 

	for (i = 0; i < n_wlens; i++) {
		k[i] = compute_k(*water_type, wavelengths[i]);
	}

	return true;
}



void compute_k_from_jerlov(float water_type, float *alphas, int *spectral_indexes, float *wavelengths, int nspec) {

	int k;

	for (k = 0; k < nspec; k++) {
		alphas[spectral_indexes[k]] = compute_k(water_type, wavelengths[k]);
	}
}

float compute_k(float water_type, float wlen) {

	int i;
	float w, alpha = 0.0, jerlov_all_types[NJERLOVTYPES];
	bool success;

// Create table of attenuation coefficients. 

	for (i = 0; i < NWAVELENS; i++) {
		jerlov_data[0][i] = Jerlov_I[i];
		jerlov_data[1][i] = Jerlov_IA[i];
		jerlov_data[2][i] = Jerlov_IB[i];
		jerlov_data[3][i] = Jerlov_II[i];
		jerlov_data[4][i] = Jerlov_III[i];	
		jerlov_data[5][i] = Jerlov_C1[i];
		jerlov_data[6][i] = Jerlov_C3[i];
		jerlov_data[7][i] = Jerlov_C5[i];
		jerlov_data[8][i] = Jerlov_C7[i];
		jerlov_data[9][i] = Jerlov_C9[i];
	}


	success = interp_jerlov_wavelength(jerlov_all_types, wlen);
	if (! success) 
		return 0.0;

	i = (int) floor(water_type);
	w = water_type - floor(water_type);

	alpha = jerlov_all_types[i]*(1.0 - w) + w*jerlov_all_types[i + 1];

	return alpha;
}


// Estimate the Jerlov water type from the ratios of attenuation coefficients 
// and the slope of the line of best fit. 

bool estimate_water_type(float ratio, float ratios[NJERLOVTYPES], float *water_type) {

	int i;

	for (i = 1; i < NJERLOVTYPES; i++) {
		if (ratio > ratios[i - 1] && ratio < ratios[i]) {
			*water_type = ((float) i - 1) + (ratio - ratios[i - 1])/(ratios[i] - ratios[i - 1]);		
			return true;
		}

		if (ratio < ratios[i - 1] && ratio > ratios[i]) {
			*water_type = ((float) i - 1) + (ratio - ratios[i])/(ratios[i - 1] - ratios[i]);
			return true;
		}
	}

	return false;
}


bool interp_jerlov_ratio(float ratio, float ratios[NJERLOVTYPES], float jerlov_i[NJERLOVTYPES], float *ki) {

	int i, indx;
	float alpha;
	bool present = false; 

// Do we need to extrapolate? (Risky, yet the coastal band of LANDSAT often has ratios below Jerlov.)

	if ((ratios[0] > ratios[NJERLOVTYPES - 1] && ratio > ratios[0]) || 
			(ratios[0] < ratios[NJERLOVTYPES - 1] && ratio < ratios[0])) {
		// printf("\nWARNING: using linear extrapolation to estimate attenuation coefficient.\n\n");
		// *ki = jerlov_i[0] + (ratio - ratios[0])/(ratios[1] - ratios[0])*(jerlov_i[1] - jerlov_i[0]);
		// *ki = jerlov_i[0];
		*ki = 0.0;
		return false;
	}

	if ((ratios[0] < ratios[NJERLOVTYPES - 1] && ratio > ratios[NJERLOVTYPES - 1]) || 
			(ratios[0] > ratios[NJERLOVTYPES - 1] && ratio < ratios[NJERLOVTYPES - 1])) {
		// printf("\nWARNING: using linear extrapolation to estimate attenuation coefficient.\n\n");
		// *ki = jerlov_i[NJERLOVTYPES - 1] + (ratio - ratios[NJERLOVTYPES - 1])/(ratios[NJERLOVTYPES - 2] - ratios[NJERLOVTYPES - 1])*(jerlov_i[NJERLOVTYPES - 2] - jerlov_i[NJERLOVTYPES - 1]);
		// *ki = jerlov_i[NJERLOVTYPES - 1];
		*ki = 0.0;
		return false;
	}

// Determine ratio index, indx, such that ratios[indx] < ratio < ratios[indx + 1].

	for (i = 1; i < NJERLOVTYPES; i++) {
		if ((ratio <= ratios[i - 1] && ratio >= ratios[i]) || (ratio >= ratios[i - 1] && ratio <= ratios[i])) {
			indx = i - 1;
			present = true;
			break ;
		}
	}

	if (! present) {
		printf("\nERROR: ratio not in Jerlov's table of spectral attenuation coefficients.\n\n");
		return false;
	}

// Linear interpolation. 

	alpha = (ratio - ratios[indx])/(ratios[indx + 1] - ratios[indx]);
	*ki = alpha*jerlov_i[indx] + (1 - alpha)*jerlov_i[indx + 1];

	return true; 
}


bool interp_jerlov_wavelength(float jdata[NJERLOVTYPES], float wl) {

	int i, indx;
	float alpha;

// Check boundary cases. 

	if (wl < wavelengths[0] || wl > wavelengths[NWAVELENS - 1]) {
		printf("\nERROR: wavelength (%.0f) not in Jerlov's table of spectral attenuation coefficients.\n\n", wl);
		return false;
	}

// Determine wavelength index, indx, such that wavelengths[indx] < wl < wavelengths[indx + 1]

	for (i = 1; i < NWAVELENS; i++) {
		if (wl <= wavelengths[i]) {
			indx = i - 1;
			break ;
		}
	}

// Linear interpolation. 

	alpha = (wl - wavelengths[indx])/(wavelengths[indx + 1] - wavelengths[indx]);

	for (i = 0; i < NJERLOVTYPES; i++) {
		jdata[i] = (1 - alpha)*jerlov_data[i][indx] + alpha*jerlov_data[i][indx + 1];
	}

	return true;
}



void jerlov_water_type_str(float wtype, char *str_type, int str_len) {

	if (wtype < 1.0) {
		sprintf(str_type, "OI + %.2f", wtype);
	} else if (wtype < 2.0) {
		sprintf(str_type, "OIA + %.2f", wtype - 1.0);
	} else if (wtype < 3.0) {
		sprintf(str_type, "OIB + %.2f", wtype - 2.0);
	} else if (wtype < 4.0) {
		sprintf(str_type, "OII + %.2f", wtype - 3.0);
	} else if (wtype < 5.0) {
		sprintf(str_type, "OIII + %.2f", wtype - 4.0);
	} else if (wtype < 6.0) {
		sprintf(str_type, "C1 + %.2f", wtype - 5.0);
	} else if (wtype < 7.0) {
		sprintf(str_type, "C3 + %.2f", wtype - 6.0);
	} else if (wtype < 8.0) {
		sprintf(str_type, "C5 + %.2f", wtype - 7.0);
	} else if (wtype < 9.0) {
		sprintf(str_type, "C7 + %.2f", wtype - 8.0);
	} else if (wtype < 10.0) {
		sprintf(str_type, "C9 + %.2f", wtype - 9.0);
	} 

}
