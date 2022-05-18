//
//	Secchi Disk and Diffuse Attenuation Coefficients
//	for LANDSAT 8 Imagery. 
//

// Ref: Zhongping Lee, Shaoling Shang, Lin Qi, Jing Yan, Gong Lin, 
//		A semi-analytical scheme to estimate Secchi-disk depth from Landsat-8 measurements, 
//		Remote Sensing of Environment 177 (2016) 101–106


#include "secchi.h"

void Kd_LS8(
	float **coastal, float **blue, float **green, float **red, float **kd,
	int nrows, int ncols, 
	float coastal_spv, float blue_spv, float green_spv, float red_spv, float theta_s) {

	int i, j;

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			if (coastal[i][j] != coastal_spv && blue[i][j] != blue_spv && 
				green[i][j] != green_spv && red[i][j] != red_spv) {
				kd[i][j] = Lee_kd_LS8(
							coastal[i][j], 
							blue[i][j], 
							green[i][j], 
							red[i][j], 
							theta_s);
				// printf("\nZsd = %.2f", zsd[i][j]);
			} else {
				kd[i][j] = coastal_spv;
			}
		}
	}
}


float Lee_kd_LS8(float coastal, float blue, float green, float red, float theta_s) {

	float Kd[5], Kd_coastal, Kd_blue, Kd_green, Kd_red, Kd_min;

	Lee_Kd_LS8(coastal, blue, green, red, theta_s, 
		&Kd_coastal, &Kd_blue, &Kd_green, &Kd_red);

	Kd[0] = Kd_coastal;
	Kd[1] = Kd_blue;
	Kd[2] = 0.20*Kd_blue + 0.75*Kd_green; // Kd(530) per Lee 
	Kd[3] = Kd_green;
	Kd[4] = Kd_red;

	Kd_min = vec_min(Kd, 5);

	return Kd_min;
}


void secchi_disk_depth(
	float **coastal, float **blue, float **green, float **red, float **zsd,
	int nrows, int ncols, 
	float coastal_spv, float blue_spv, float green_spv, float red_spv, float theta_s) {

	int i, j;

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			if (coastal[i][j] != coastal_spv && blue[i][j] != blue_spv && 
				green[i][j] != green_spv && red[i][j] != red_spv) {
				zsd[i][j] = Lee_secchi_LS8(
							coastal[i][j], 
							blue[i][j], 
							green[i][j], 
							red[i][j], 
							theta_s);
				// printf("\nZsd = %.2f", zsd[i][j]);
			} else {
				zsd[i][j] = coastal_spv;
			}
		}
	}
}


// Secchi disk depth for LS8 imagery

float Lee_secchi_LS8(float coastal, float blue, float green, float red, float theta_s) {

	float Rrs[4], Kd[5], Kd_coastal, Kd_blue, Kd_green, Kd_red, Rrs_max, Kd_min;

	Rrs[0] = coastal;
	Rrs[1] = blue; 
	Rrs[2] = green; 
	Rrs[3] = red; 

	Rrs_max = vec_max(Rrs, 4);

	Lee_Kd_LS8(coastal, blue, green, red, theta_s, 
		&Kd_coastal, &Kd_blue, &Kd_green, &Kd_red);

	Kd[0] = Kd_coastal;
	Kd[1] = Kd_blue;
	Kd[2] = 0.20*Kd_blue + 0.75*Kd_green; // Kd(530) per Lee 
	Kd[3] = Kd_green;
	Kd[4] = Kd_red;

	// printf("\nRrs443,481,554,656 = %.7f, %.7f, %.7f, %.7f", coastal, blue, green, red);
	// printf("\nKd443,481,554,656 = %.5f, %.5f, %.5f, %.5f", Kd[0], Kd[1], Kd[2], Kd[3]);

	Kd_min = vec_min(Kd, 5);

	return log(fabs(0.14 - Rrs_max)/0.013)/(2.5*Kd_min);
}


// Diffuse Attenuation Coefficients of Downwelling Radiance. 

void Lee_Kd_LS8(float Rrs443, float Rrs481, float Rrs554, float Rrs656, float theta_s, 
	float *Kd443, float *Kd481, float *Kd554, float *Kd656) {

	double rrs443, rrs481, rrs554, rrs656, 
		g0, g1, aw,
		u443, u481, u554, u656, 
		bbw443, bbw481, bbw554, bbw656,
		h0, h1, h2, chi, eta, 
		a443, a481, a554, a656, 
		bb443, bb481, bb554, bb656,
		bbp443, bbp481, bbp554, bbp656,
		m0, m1, m2, m3, gamma,
		kk1, kk2, kk3;

	// Compute rrs from Rrs. 

	rrs443 = Rrs443/(0.52 + 1.7*Rrs443);
	rrs481 = Rrs481/(0.52 + 1.7*Rrs481);
	rrs554 = Rrs554/(0.52 + 1.7*Rrs554);
	rrs656 = Rrs656/(0.52 + 1.7*Rrs656);

	// Model constants per Lee, 2002.

	g0 = 0.0895; 
	g1 = 0.1247; 

	// Compute u. Per model of Gordon et al, 1988. 

	u443 = (-g0 + sqrt(g0*g0 + 4.0*g1*rrs443))/(2.0*g1); 
	u481 = (-g0 + sqrt(g0*g0 + 4.0*g1*rrs481))/(2.0*g1); 
	u554 = (-g0 + sqrt(g0*g0 + 4.0*g1*rrs554))/(2.0*g1); 
	u656 = (-g0 + sqrt(g0*g0 + 4.0*g1*rrs656))/(2.0*g1);

	// Compute bbw, via interpolation from Morel, 1974. 

	bbw443 = 0.00244761;
	bbw481 = 0.00171397;
	bbw554 = 0.000931339;
	bbw656 = 0.000448682;

	// Compute aw (at reference wavelength 554 nm) via interpolation from Pope & Fry, 1997. 

	aw = 0.05866;

	// Compute a (at reference wavelength 554 nm).

	h0 = -1.146;
	h1 = -1.366;
	h2 = 0.469;

	chi = log10((rrs443 + rrs481)/(rrs554 + 5.0*rrs656*rrs656/rrs481));

	a554 = aw + pow(10.0, h0 + h1*chi + h2*chi*chi);

	// Compute bb (at reference wavelength 554 nm).

	bb554 = (-a554*g0 + 2.0*a554*rrs554 + a554*sqrt(g0*g0 + 4.0*g1*rrs554))/(2.0*(g0 + g1 - rrs554));

	// Compute bbp (at reference wavelength 554 nm).	

	bbp554 = (u554*a554)/(1 - u554) - bbw554;

	// Compute bbp. 

	eta = 2.0*(1.0 - 1.2*exp(-0.9*rrs443/rrs554));

	bbp443 = bbp554*pow(554.0/443.0, eta);
	bbp481 = bbp554*pow(554.0/481.0, eta);
	bbp656 = bbp554*pow(554.0/656.0, eta);

	// Compute a. 

	a443 = (1.0 - u443)*(bbw443 + bbp443)/u443;
	a481 = (1.0 - u481)*(bbw481 + bbp481)/u481;
	a656 = (1.0 - u656)*(bbw656 + bbp656)/u656;

	// Compute bb. 

	bb443 = (-a443*g0 + 2.0*a443*rrs443 + a443*sqrt(g0*g0 + 4.0*g1*rrs443))/(2.0*(g0 + g1 - rrs443));
	bb481 = (-a481*g0 + 2.0*a481*rrs481 + a481*sqrt(g0*g0 + 4.0*g1*rrs481))/(2.0*(g0 + g1 - rrs481));
	bb656 = (-a656*g0 + 2.0*a656*rrs656 + a656*sqrt(g0*g0 + 4.0*g1*rrs656))/(2.0*(g0 + g1 - rrs656));

	// Compute Kd. 

	m0 = 0.005;
	m1 = 4.26;
	m2 = 0.52;
	m3 = 10.8;
	gamma = 0.265;

	// Kd_443

	kk1 = (1.0 + m0*theta_s)*a443; 
	kk2 = m1*(1.0 - gamma*bbw443/bb443);
	kk3 = (1.0 - m2*exp(-m3*a443))*bb443; 
	*Kd443 = kk1 + kk2*kk3;

	// Kd_481

	kk1 = (1.0 + m0*theta_s)*a481; 
	kk2 = m1*(1.0 - gamma*bbw481/bb481);
	kk3 = (1.0 - m2*exp(-m3*a481))*bb481; 
	*Kd481 = kk1 + kk2*kk3;

	// Kd_554

	kk1 = (1.0 + m0*theta_s)*a554; 
	kk2 = m1*(1.0 - gamma*bbw554/bb554);
	kk3 = (1.0 - m2*exp(-m3*a554))*bb554; 
	*Kd554 = kk1 + kk2*kk3;

	// Kd_656

	kk1 = (1.0 + m0*theta_s)*a656; 
	kk2 = m1*(1.0 - gamma*bbw656/bb656);
	kk3 = (1.0 - m2*exp(-m3*a656))*bb656; 
	*Kd656 = kk1 + kk2*kk3;
}


















