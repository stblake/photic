//
//	photic -- A Multi-spectral, Multi-temporal, Multi-spatial, Remote Sensing Model for Shallow Waters
//	


// Written by Sam Blake.

// Started 4 July 2017.

// An implementation of the HOPE (Hyperspectral Optimization Process Exemplar) model 
// of Lee et al based on publications of Lee et al. This implementation is generalised
// for the use of multi-spectral imagery. In particular, we compensate for the lack of 
// spectral resolution (in comparison to hyperspectral imagery) as follows: 

// Simultaneously processing multiple scenes -- across multiple scenes we assume the 
// depths are static (up to tidal variations) and the bottom reflectance is unchanged. 
// However, the attenuation coefficients (k) do vary across scenes. 

// Simultaneously processing multiple pixels -- we further increase the data available to 
// the model by introducing a _region_ around the given pixel. Within 
// this region (of each scene) we assume that k is static. 

// To accomodate multiple bottom types we use the approach of Gege (2006) 
// where R_b = sum B_i*q_i*rho_i where sum q_i = 1. We also employ the approach of Wettel (2006) 
// to process all partial combinations of bottom types (however this is computationally expensive).

// Note that this model assumes the input spectra has been atmospherically, adjacency and 
// sun glitter corrected. It also assumes that land and clouds have been masked from the dataset. 

// In regards to optimisation, the simplex algorithm takes ~20m^2 iterations to converge. The 
// original Lee algorithm for hyperspectral data only had 5 parameters (P, G, X, B, H). 
// Incorporating multiple bottom types requires rho_i (i = 1 .. n_bottoms). The _photic_
// model models a region around a given pixel (of n_regions pixels) and multiple 
// scenes (n_scenes). This increases the number of unknown parameters in the model to 
// n_regions + n_regions*n_bottoms + 4*n_scenes. Consequently, this model requires significantly 
// more iterations to converge. In order to efficiently model large scenes we make the following 
// optimisations: 
//
//     * We use the statistical model of Lyzenga et al (2006) to provide an initial depth approximation. 
//
//     * We model the pixels from deepest to shallow. (Using the statistical depth approximation.)
//
//     * We start the K parameters (P,G,X) in the optimisation process at the previously derived 
//           optima. 
//
//     * Prior to running the optimisation algorithm, we check a lookup table of previously computed 
//           spectra, if these are "close" to the current pixel we use the lookup table results. 


// References:
//
//		Lee et al, Hyperspectral remote sensing for shallow waters. I. A semi-analytical model, 
//			Vol. 37, No. 27, APPLIED OPTICS, 20 September 1998
//
//		Lee et al, Hyperspectral remote sensing for shallow waters: 2. Deriving bottom depths 
//			and water properties by optimization, Vol. 38, No. 18, APPLIED OPTICS, 20 June 1999
//
//		Nelder, J.A., and Mead, A simplex method for function minimisation, R. 1965, 
//			Computer Journal, vol. 7, pp. 308–313
//
//		Gege P, Albert A, Inversion of irradiance and remote sensing reflectance in shallow water 
//      	between 400 and 800 nm for calculations of water and bottom properties, APPLIED OPTICS, 
//      	vol. 45, no. 10, April 2006
//
//		Wettle M, Brando V, SAMBUCA - Semi-Analytical Model for Bathymetry, Un-mixing, and 
//			Concentration Assessment, Environmental Remote Sensing Group CSIRO Land and Water, 
//			July 2006  
//
//		Lyzenga D, Malinas N, Tanis F, Multispectral Bathymetry Using a Simple Physically 
//			Based Algorithm, IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, 
//			vol. 44, no. 8, August 2006
//

// Note on performance: Lee's HOPE implementation could process ~15 pixels per second in 
// 2001 on a 400Mhz processor using a "predictor-corrector" optimisation scheme. In 2009 
// It could process ~157 pixels per second. 

// In 2019, _photic_ can process ~600 pixels per second per core, where on average 96% of pixels 
// are modelled using the LUT (using the configuration n_regions = 9, n_scenes = 2, n_bottoms = 2.) 
// 


#include "samodel.h"
#if _OPENMP
#include <omp.h>
#endif

#define DEBUG 0
#define WRITE_EPGXBHD 0

#define LOCALLY_AVERAGE 1

#define ALL_BOTTOMS 1
#define MAX_BOTTOMS 12

#define SAM 0
#define LUT 1
#define RESTART 0
#define DELTA 0
#define K_PRIOR 0

// Backscattering slope Y of b_bp = X (440/lambda)^Y:
//   0  from chi = Rrs440/Rrs490 of the measured spectrum (Lee 1999), as originally;
//   1  self-consistent: QAA's Y = 2[1 - 1.2 exp(-0.9 rrs_dp(443)/rrs_dp(555))] from the
//      model's own deep-water reflectance (no bottom contamination);
//   2  self-consistent: Lee's chi relation applied to the model's rrs_dp(440)/rrs_dp(490);
//   3  fitted per scene, in [Y_MIN, Y_MAX].

// Optically deep test (EXPERIMENTS.md C2): the bottom counts as detected only if the best fit
// beats an optically deep fit (no bottom; the water re-fitted) by
//   chi2 = n (E_deep^2 - E_best^2)/(100 sigma_rel)^2 > DEEP_PENALTY k ln(n),
// a multiple of the Bayesian information criterion's penalty for the k depth and bottom
// parameters the shallow model adds. E is the Rrs error (RMS %), n the number of residuals and
// sigma_rel^2 = (sigma_noise/mean Rrs)^2 + DEEP_MODEL_ERROR^2 (the forward model's own error,
// about 2-3% against exact radiative transfer, tier 2). DEEP_PENALTY was chosen on clear and
// moderate tier 3 pixels and checked on turbid ones. Pixels failing the test get no depth.

#ifndef OPTICALLY_DEEP
#define OPTICALLY_DEEP 1
#endif
#ifndef DEEP_MODEL_ERROR
#define DEEP_MODEL_ERROR 0.03
#endif
#ifndef DEEP_PENALTY
#define DEEP_PENALTY 2.0
#endif

// Coefficients of Lee's shallow-water model (Lee et al. 1998, 1999), fitted to Hydrolight:
// rrs_dp = (G0 + G1 u) u; D_uC = C0 sqrt(1 + C1 u); D_uB = B0 sqrt(1 + B1 u). EXPERIMENTS.md B10
// refits them to exact radiative transfer (OSOAA).

#ifndef LEE_G0
#define LEE_G0 0.084
#define LEE_G1 0.170
#define LEE_C0 1.03
#define LEE_C1 2.4
#define LEE_B0 1.04
#define LEE_B1 5.4
#endif

// Bottom-surface interreflection (EXPERIMENTS.md B11). Light reflected by the bottom is partly
// reflected back down by the water-air interface and the water column, and re-reflected by the
// bottom, so the bottom term is rho/pi exp(-M K H)/(1 - c rho), where the reflectance seen
// from the bottom looking up is
//   c = R_INTERNAL exp(-2 D_DIFFUSE K H) + pi rrs_dp (1 - exp(-2 D_DIFFUSE K H)):
// the interface's internal reflectance for diffuse light over the round trip, plus the water
// column's own reflectance. R_INTERNAL and D_DIFFUSE were fitted to OSOAA (0.464, theory about
// 0.48; 2.07, about 1/mu for diffuse light). Lee's model (INTERREFLECTION 0) is linear in rho.

#ifndef INTERREFLECTION
#define INTERREFLECTION 0
#endif
#define R_INTERNAL 0.464
#define D_DIFFUSE 2.07

// Seabed continuity (EXPERIMENTS.md B4): weak prior residuals on the second differences of
// depth along the rows and columns of the spatial window (a plane is free), each divided by
// CONT_SIGMA_0 + CONT_SIGMA_H*mean depth (m), in the same chi units as the Rrs residuals.

#ifndef CONTINUITY
#define CONTINUITY 0
#endif
#ifndef CONT_SIGMA_0
#define CONT_SIGMA_0 0.1
#endif
#ifndef CONT_SIGMA_H
#define CONT_SIGMA_H 0.05
#endif

// Water tied between scenes (EXPERIMENTS.md B12): prior residuals (ln W_s - mean ln W)/TIE_SIGMA
// for W = P, G, X, in the same chi units as the Rrs residuals. With the water tied, the scenes'
// different light paths (sun angle, tide) constrain depth independently of the bottom albedo.

#ifndef WATER_TIE
#define WATER_TIE 0
#endif
#ifndef TIE_SIGMA
#define TIE_SIGMA 0.25
#endif

#ifndef Y_MODEL
#define Y_MODEL 0
#endif
#define Y_MIN 0.0
#define Y_MAX 2.5

#if Y_MODEL == 3
#define N_WATER 4 // ln P, ln G, ln X, Y per scene.
#else
#define N_WATER 3 // ln P, ln G, ln X per scene.
#endif // Penalise low K in shallow water and very high K (off: it excludes the truth, EXPERIMENTS.md A-2).

#if DELTA
#error "DELTA is not supported by the bounded parameterisation"
#endif

// Bounds of the water parameters (m^-1): phytoplankton absorption P = a_phi(440), gelbstoff
// absorption G = a_g(440) and particle backscattering X = b_bp(440).

#define P_MIN 0.0005
#define P_MAX 1.0
#define G_MIN 0.0005
#define G_MAX 2.0
#define X_MIN 0.0001
#define X_MAX 0.5

// Bounds of depth (m), log-uniform so the search keeps its resolution in shallow water.

#define H_MIN 0.05
#define H_MAX 60.0

// Starting bottom albedo.

#define ALBEDO_START 0.3

// Reference wavelengths. 

int n_ref_wlens = 41;

double ref_wlens[] = {400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 
	530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 
	670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800};

//
//	Spectral Weights
//

// These are ~ 1/K where K is the average of the ocean-type waters (I, IA, IB, II, III) from 
// Jerlov's table of spectral attenuation coefficients. 

// Ref: Jerlov, N, Marine Optics, 2nd Edition, Table XXVII, pp. 135, 1976.

double spectral_weights[] = {12.56,13.40,14.35,15.36,16.40,17.61,18.38,19.23,19.20,18.30,
	17.48,16.32,15.30,14.21,13.11,12.17,10.76,9.649,7.236,5.086,3.922,3.519,3.191,2.938,
	2.738,2.564,2.406,2.267,2.073,1.855,1.678,1.500,1.500,1.500,1.500,1.500,1.500,1.500,
	1.500,1.500,1.500};

//
//	Absorption Coefficients of Pure Water
//

// Ref: R. Pope and E. Fry, “Absorption spectrum 380-700nm of pure waters: II. Integrating 
//		cavity measurements,” Appl. Opt. 36, pp. 8710–8723, 1997.

double aw_Pope_Fry1997[] = {0.00663,0.00473,0.00454,0.00495,0.00635,0.00922,0.00979,0.0106,0.0127,
	0.015,0.0204,0.0325,0.0409,0.0434,0.0474,0.0565,0.0619,0.0695,0.0896,0.1351,0.2224,0.2644,
	0.2755,0.2916,0.3108,0.34,0.41,0.439,0.465,0.516,0.624,0.827,1.231,1.65768,2.38147,2.47014,
	2.5045,2.47509,2.35431,2.17654,2.01226};

//
//	Backscattering Coefficients of Pure Water
//

// Ref: A. Morel, “Optical properties of pure water and pure sea waters,” in Optical Aspects of Oceanography, 
//		N. G. Jerlov and E. S. Nielsen, eds. Academic, pp. 1–24, New York, 1974

double bbw_Morel_1974[] = {0.0038,0.00341552,0.00307784,0.00278035,0.00251749,0.00228457,0.00207763,
	0.0018933,0.0017287,0.00158138,0.00144921,0.00133039,0.00122334,0.00112671,0.0010393,0.000960099,
	0.000888199,0.000822817,0.000763262,0.000708928,0.000659279,0.000613844,0.000572204,0.000533988,
	0.000498868,0.000466549,0.00043677,0.000409298,0.000383923,0.000360458,0.000338734,0.000318601,
	0.000299921,0.000282571,0.000266441,0.000251431,0.000237448,0.00022441,0.000212243,0.000200879,
	0.000190254};

// 
//	a_0(lambda), a_1(lambda) as given in Lee, 1994
//

// Ref: Z. P. Lee, “Visible-infrared remote-sensing model and applications for ocean waters,” 
//		Ph.D. dissertation  Department of Marine Science, University of South Florida, St. Petersburg, Fla., 1994 .

double a0_Lee[] = {0.684266,0.778234,0.86369,0.960272,1,0.96339,0.931107,0.869651,0.788982,
	0.755811,0.733335,0.691062,0.632722,0.568079,0.504551,0.426236,0.343342,0.294978,0.278405,
	0.259509,0.238935,0.274531,0.319712,0.342135,0.333133,0.350196,0.560974,0.843456,0.74852,
	0.4,0.15,0.06,0.026,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013};

double a1_Lee[] = {0.0205345,0.0129366,0.0064013,0.00168429,0,0.00602146,0.0109083,0.0157169,0.0152315,
	0.0255969,0.0559074,0.0865174,0.098141,0.0968975,0.0899322,0.0780758,0.0658936,0.0596445,0.0581118,
	0.0539721,0.0494872,0.057845,0.0673832,0.0717643,0.0684625,0.0713488,0.112818,0.159471,0.138835,
	0.0811786,0.03,0.012,0.0054,0.0028,0.0028,0.0028,0.0028,0.0028,0.0028,0.0028,0.0028};

//
//	Bottom Albedo Spectrums
//

// Ref: Lee et al, Properties of the water column and bottom derived from Airborne Visible Infrared 
//		Imaging Spectrometer (AVIRIS) data, JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 106, 
//		NO. C6, PAGES 11639-11651, JUNE 15,2001

double bottom_sand_Lee[] = {0.583466,0.596165,0.623453,0.659155,0.696085,0.721449,0.752614,0.784261,
	0.810131,0.838364,0.870429,0.898984,0.930959,0.955174,0.976314,1,1.02848,1.05523,1.08057,
	1.11548,1.13863,1.15436,1.16485,1.18209,1.19692,1.20153,1.17874,1.15904,1.18503,1.25993,
	1.34427,1.39835,1.43151,1.46162,1.49477,1.51586,1.52225,1.54963,1.56201,1.57216,1.57973};

double bottom_seagrass_Lee[] = {0.271406,0.277846,0.274747,0.274747,0.271616,0.267027,0.275068,
	0.277633,0.288372,0.31707,0.402691,0.540072,0.697163,0.817698,0.890198,0.973223,1.06845,
	1.109,1.03897,0.948927,0.895306,0.857696,0.802312,0.722199,0.606007,0.489014,0.445886,
	0.451764,0.570228,1.08633,2.52329,4.4967,6.0473,6.69215,6.91428,7.02582,7.09853,7.1631,
	7.23864,7.29724,7.35732};

double a0_Blake[] = {0.693224,0.898908,0.992684,1.02963,1.,0.792275,0.748701,0.776112,0.726627,0.632112,
		     0.515369,0.407641,0.328754,0.27633,0.236208,0.197239,0.168193,0.154953,0.157945,
		     0.164272,0.160942,0.164892,0.178776,0.199702,0.241388,0.264834,0.314421,0.49153,
		     0.464896,0.19484,0.0780387,0.0442571,0.0320876,0.0209004,0.0167131,0.0126202};

double a1_Blake[] = {0.0103484,0.0227841,0.0237051,0.0138089,0.0,-0.0229181,-0.0227291,-0.00745802,
		     -0.0000655587,-0.00108625,0.000734957,0.00391282,0.00617119,0.008739,0.0100529,0.00974627,
		     0.0104398,0.0107561,0.0110176,0.0105372,0.0106184,0.0107536,0.0107262,0.0132929,0.0212209,
		     0.0258873,0.023263,0.0317735,0.0294253,0.0100652,0.00588207,0.00401144,0.0036124,0.00252763,
		     0.00382846,0.00287507};


// Ref: Mapping seaweed forests with IKONOS image based on bottom surface reflectance, Remote 
//		Sensing of the Marine Environment II, edited by Robert J. Frouin, Naoto Ebuchi, Delu Pan, 
//		Toshiro Saino, Proc. of SPIE Vol. 8525

double bottom_rock_Sagawa[] = {0.446844,0.463106,0.489994,0.520033,0.547332,0.577717,0.613093,0.640429,
	0.664388,0.686469,0.716171,0.746174,0.789359,0.854467,0.937349,1.,1.06547,1.10031,1.11612,
	1.12155,1.11591,1.09913,1.07528,1.06594,1.06857,1.07232,1.05406,1.03511,1.05512,1.15098,
	1.28489,1.38377,1.45716,1.5079,1.54188,1.55929,1.57666,1.59186,1.60415,1.61531,1.6261};

double bottom_brown_alga_Sagawa[] = {1.02407,1.01035,0.95028,0.94553,0.942414,0.915504,0.874995,
	0.847023,0.85149,0.855956,0.864693,0.875119,0.885544,0.913803,0.951189,1.,1.12001,
	1.24112,1.19843,1.19843,1.22309,1.10903,0.981867,0.926624,1.00622,0.947069,0.817357,
	0.803944,0.979778,1.35551,3.15296,7.23477,10.3982,11.0253,10.7855,10.9471,11.1577,
	11.1373,11.3919,11.8403,12.2725};

double bottom_red_seaweed_Sagawa[] = {0.659643,0.74779,0.753905,0.759553,0.734172,0.735914,0.743928,
	0.752305,0.762707,0.784911,0.828281,0.86661,0.924044,0.982859,1.00066,1.,0.980934,0.963559,
	0.966101,0.967321,0.958482,0.967321,0.992428,1.00604,1.01299,1.05706,1.07768,1.13556,1.17604,
	1.25306,1.4027,1.78951,2.26202,2.54165,2.64516,2.72694,2.79462,2.86691,2.92661,2.98747,3.04329};

// Ref: Nicole Pinnel, Bolivar, South Australia 

double bottom_sand_Bolivar[] = {0.831952,0.830508,0.834468,0.831049,0.845008,0.862618,0.878175,0.893138,
	0.911949,0.917074,0.927628,0.941394,0.951449,0.96914,0.982218,1.,1.01706,1.0287,1.03817,
	1.03768,1.0442,1.03516,1.02352,1.01749,1.02367,1.01656,0.971376,0.924693,0.940348,1.02032,
	1.13057,1.19403,1.20302,1.21271,1.22576,1.23423,1.23815,1.24207,1.24403,1.24558,1.25408};

// Ref: Dekker, Imaging Spectronomy of Waters, 2002

double bottom_coral_dekker[] = {0.383243, 0.410755, 0.393293, 0.401978, 0.44859, 0.503933, 0.533888,
	0.575502, 0.639226, 0.695987, 0.734558, 0.797606, 0.86595, 0.905876, 
	0.939073, 1., 1.0991, 1.24556, 1.33291, 1.35833, 1.44966, 1.59142, 
	1.53355, 1.47955, 1.44494, 1.43629, 1.28982, 1.03614, 1.08886, 1.72366, 2.99676, 3.1, 3.2, 3.3, 
	3.4, 3.5, 3.6, 3.7, 3.8, 3.9};

// Averaged bottom types from HYDROLIGHT. 

int n_bottom_table = 8;

double bottom_type_sand[] = {0.634432,0.654965,0.682443,0.712079,0.74123,0.769656,0.802245,0.827992,
	0.854624,0.878036,0.89898,0.924542,0.94899,0.965671,0.982481,1.,1.01237,1.02994,1.04658,1.07139,
	1.08139,1.08884,1.09447,1.09916,1.10557,1.10758,1.101,1.08464,1.08658,1.12477,1.17063,1.20125,
	1.21407,1.22691,1.23186,1.24098,1.24708,1.25257,1.25607,1.26361,1.26789};

double bottom_type_seagrass[] = {0.542201,0.538755,0.543345,0.550218,0.550218,0.562483,0.606813,0.607225,
	0.63619,0.650873,0.670997,0.735138,0.808387,0.887508,0.949743,1.,1.03488,1.04312,1.03199,1.00882,
	0.95964,0.903696,0.85915,0.848253,0.830917,0.807567,0.785207,0.839027,1.09213,1.68562,2.5015,
	4.15478,5.11731,5.77383,6.02484,6.11441,6.16467,6.22648,6.27581,6.29811,6.34481};

double bottom_type_coral[] = {0.437155,0.381353,0.353904,0.339709,0.358934,0.382172,0.390592,0.407088,
	0.462696,0.522577,0.60409,0.707443,0.764575,0.826841,0.892756,1.,1.15687,1.3644,1.44365,1.40967,
	1.48473,1.53551,1.41572,1.32706,1.25349,1.20397,0.955824,0.643935,0.669901,1.38703,2.42848,
	3.3846,3.6429,3.89712,4.10336,4.23733,4.32766,4.37486,4.43428,4.47756,4.52468};

double bottom_type_macrophytes[] = {0.443737,0.401082,0.382427,0.382427,0.382427,0.382427,0.402497,
	0.409674,0.443211,0.44772,0.494445,0.620241,0.755823,0.882367,0.973539,1.,1.09947,1.09321,
	1.11723,1.09127,1.06283,1.03064,0.919521,0.896435,0.839186,0.769511,0.645607,0.60577,
	0.743085,1.1217,2.17169,3.27124,3.47565,3.3989,3.48151,3.59713,3.65614,3.65521,3.67159,3.6131,3.58176};

double bottom_type_dark_sediment[] = {0.57054,0.564806,0.583536,0.596254,0.616995,0.640224,0.669295,
	0.691865,0.722327,0.740098,0.765618,0.801745,0.850374,0.901024,0.943913,1.,1.05834,1.10294,
	1.13346,1.17132,1.20427,1.20038,1.16874,1.16074,1.14623,1.14875,1.03881,0.903713,0.894767,
	1.05454,1.30921,1.49551,1.58863,1.64064,1.65234,1.66933,1.68055,1.69243,1.70471,1.72121,1.72876};

double bottom_type_coral_sand[] = {0.616201,0.636915,0.659543,0.685341,0.71562,0.742639,0.769601,
	0.792999,0.818446,0.840791,0.86273,0.890834,0.912898,0.940877,0.970795,1.,1.03073,1.05991,
	1.08887,1.1141,1.13666,1.15865,1.18297,1.20189,1.22277,1.23953,1.25598,1.27532,1.29161,
	1.31127,1.3302,1.348,1.36099,1.37632,1.38579,1.39478,1.40298,1.41371,1.42261,1.43402,1.44206};

double bottom_type_green_algae[] = {0.16848,0.188393,0.188646,0.17936,0.170655,0.190607,0.220843,0.241126,
	0.262068,0.294282,0.357834,0.516557,0.717781,0.898626,0.975598,1.,1.00009,0.946152,0.874696,0.829798,
	0.807546,0.763138,0.701775,0.673854,0.606387,0.524554,0.426766,0.3044,0.284069,0.48541,0.912011,
	1.33421,1.52307,1.61754,1.65643,1.66741,1.66741,1.66741,1.66741,1.66741,1.66298};

double bottom_type_red_algae[] = {1.25032,1.13331,1.06326,1.00005,0.925705,0.983453,1.,1.,0.975044,
	0.877564,0.877766,0.987416,1.04918,1.03747,1.,1.,1.01387,1.16916,1.56137,2.08995,2.44941,
	2.40749,2.30252,2.33337,2.48719,2.53279,2.21543,1.72598,1.80604,2.72284,4.35056,5.67543,
	6.32734,6.72766,6.90232,6.98939,7.04997,7.09008,7.11525,7.09872,7.04973};

double bottom_type_sand66_green_algea33[] = {0.478636,0.498941,0.517326,0.533971,0.550487,0.576063,
	0.607836,0.631737,0.656448,0.682768,0.717879,0.787759,0.871049,0.942379,0.979206,0.999,
	1.00727,1.00101,0.988295,0.98987,0.989118,0.979293,0.962609,0.956435,0.938234,0.912327,
	0.875376,0.823737,0.818259,0.910737,1.08334,1.24433,1.31575,1.35576,1.37201,1.38174,
	1.3858,1.38946,1.39179,1.39681,1.39819};

double bottom_type_sand50_green_algea50[] = {0.401456,0.421679,0.435545,0.445719,0.455943,0.480132,
	0.511544,0.534559,0.558346,0.586159,0.628407,0.720549,0.833386,0.932149,0.97904,1.,1.00623,
	0.988044,0.960637,0.950595,0.944467,0.92599,0.898123,0.886508,0.855977,0.816069,0.763881,
	0.694521,0.685326,0.805089,1.04132,1.26773,1.36857,1.42222,1.44415,1.45419,1.45724,1.45999,
	1.46174,1.46551,1.46544};


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
		float pagesize, int background, int linewidth) {
  
	int thread_id, n_done;

#if _OPENMP
  #pragma omp parallel private(thread_id, n_done)
#endif
  {

#if _OPENMP
    thread_id = omp_get_thread_num();
    
    if (thread_id == 0) {
      printf("\n\nOpenMP number of threads = %d\n", omp_get_num_threads());
    }
#else
    thread_id = 0;
#endif

	int i, j, k, ii, jj, iii, jjj, ir, jr, k_region, k_scene, k_band, k_bottom, k_lut, *n_bands, *n_raw_bands, max_n_bands,
		n_regions, max_n_regions, nrows, ncols, n_shallow, n_complete, nprev, nperc, 
		n_converged, n_not_converged, k_sample, k_depth, k_trial, n_samples, n_trials, n_depth_intervals, 
		n_pertubations, n_total_per_interval, k_coastal, k_blue, k_green, k_red, **indexes, 
		**shallow_indexes, n_shallow_pre, n_lut, lut_position, *lut_index, oldest_position, oldest_index, 
		total_lut, max_bottom_types, n_modelled_optimisation;

	float ***K_coastal, ***K_blue, ***K_green, ***K_red, ***P, ***G, ***X, ***D, **spectral_index, 
		blue_sum, green_sum, red_sum, blue, green, red;

	double ***Rrs_measured, ***Rrs_modelled, ***rrs_modelled, ***rrs_bottom, ***rrs_dp, ***rho,
		**wavelengths, **raw_wavelengths, *start, *min, *step, **bottom_types, ***bottom_reflectance,
		**a_0, **a_1, **a_w, **b_bw, **K, **rrs_noise, *theta_view, *theta_sun, *sec_theta_view, *sec_theta_sun,
		**Rrs440, **Rrs490, **Rrs550, **Rrs640, **Rrs750, *h_tide, 
		a_w640, pixels_per_sec, elapsed, refl, total_start_depth,
		total_depth, n_modelled, total_bottom_albedo, total_Rrs_error, 
		total_depth_error, total_bottom_error, total_model_error, 
		total_kmin, total_n_iterations, total_B_sand, total_B_seagrass, 
		total_B_coral, total_index_optical_depth, d, depth_interval, max_depth_reached, 
		*trial_depths, *depth_sigma_table, mean_trial_depths, n_sigma, time_remaining, 
		*PP, *GG, *XX, *DD, *prev, total_sam_error, total_K_error, total_spectral_diff, 
		spectral_diff, n_skipped, lut_model_error, best_lut_spectral_match, lut_match_threshold, 
		n_restarted, n_obs, lut_min_H, lut_mean_H, lut_max_H, lut_mean_K_min, lut_min_K_min, 
		lut_max_K_min, lut_mean_age, continuity_threshold, skipped_ratio, max_error_before_restart, 
		min_preferable_skip, max_preferable_skip, min_lut_match_threshold, max_lut_match_threshold, 
		blue_noise, green_noise, red_noise, spectral_angle, *error_list, lut_entry_threshold, 
		error_mean, error_stddev, total_P, total_G, total_X, mean_P, mean_G, mean_X; 
	model_data *md, *md_lut, md_copy;
	bool nodata, neg_rrs, within_interval, *lut_in_use, rrs_noise_present;
	char string[4];
	clock_t t0, t1, t2;

	// n_smoothing_radius = 1; // TODO -- these should be set at runtime! 
	// n_spatial = 2; // 2, 1 
	// n_bottoms = 3;

	printf("\nn_smoothing_radius = %d\nn_spatial = %d, \nn_bottoms = %d\n", 
		n_smoothing_radius, n_spatial, n_bottoms);

	// Indexes for named spectral bands. 

	k_coastal = 0;
	k_blue    = 1;
	k_green   = 2;
	k_red     = 3;

#if DEBUG
	printf("\nENTER samodel\n");
#endif

	// Set random seed. 

	srand(time(NULL));

	// Allocate memory for the model data struct. 

	md = (struct model_data_struct*) malloc(sizeof(struct model_data_struct));

	// Maximum number of regions. 

	max_n_regions = pow(2*n_spatial + 1, 2);
	md->max_n_regions = max_n_regions;

	// Allocate memory for the number of spectral bands in each scene. 

	n_bands = (int*) malloc(nscenes*sizeof(int));
	md->n_bands = n_bands;

	// Array dimensions. 

	nrows = scene_data[scene_indexes[0]].nrows;
	ncols = scene_data[scene_indexes[0]].ncols;

	printf("\nnrows,ncols = %d,%d\n", nrows, ncols);

	// Number of bands in each scene. 

	md->n_scenes = nscenes;

	n_raw_bands = (int*) malloc(nscenes*sizeof(int));
	md->n_raw_bands = n_raw_bands;

	max_n_bands = 0;

	for (i = 0; i < nscenes; i++) {
		md->n_bands[i] = scene_data[scene_indexes[i]].n_bands;
		md->n_raw_bands[i] = md->n_bands[i];
		if (md->n_bands[i] > max_n_bands) {
			max_n_bands = md->n_bands[i];
		}
	}

	md->max_n_bands = max_n_bands;

	printf("\nn_scenes = %d, max_n_bands = %d", nscenes, max_n_bands);
	fflush(stdout);

	// Allocate memory for fixed model parameters. 

	allocate_double_array_3d(&Rrs_measured, max_n_regions, nscenes, max_n_bands);
	allocate_double_array_3d(&Rrs_modelled, max_n_regions, nscenes, max_n_bands);
	allocate_double_array_3d(&rrs_modelled, max_n_regions, nscenes, max_n_bands);
	allocate_double_array_3d(&rrs_bottom,   max_n_regions, nscenes, max_n_bands);
	allocate_double_array_3d(&rho,          max_n_regions, nscenes, max_n_bands);
	allocate_double_array_3d(&rrs_dp,       max_n_regions, nscenes, max_n_bands);

	allocate_double_array_2d(&wavelengths, nscenes, max_n_bands);
	allocate_double_array_2d(&raw_wavelengths, nscenes, max_n_bands);
	allocate_double_array_2d(&a_0, nscenes, max_n_bands);
	allocate_double_array_2d(&a_1, nscenes, max_n_bands);
	allocate_double_array_2d(&a_w, nscenes, max_n_bands);
	allocate_double_array_2d(&b_bw, nscenes, max_n_bands);
	allocate_double_array_2d(&K, nscenes, max_n_bands);
	allocate_double_array_2d(&rrs_noise, nscenes, max_n_bands);

	theta_view     = (double*) malloc(nscenes*sizeof(double));
	theta_sun      = (double*) malloc(nscenes*sizeof(double));
	sec_theta_view = (double*) malloc(nscenes*sizeof(double));
	sec_theta_sun  = (double*) malloc(nscenes*sizeof(double));
	h_tide         = (double*) malloc(nscenes*sizeof(double));

	md->n_bottoms       = n_bottoms;
	md->Rrs_measured    = Rrs_measured;
	md->Rrs_modelled    = Rrs_modelled;
	md->rrs_modelled    = rrs_modelled;
	md->rrs_bottom      = rrs_bottom;
	md->rrs_dp          = rrs_dp;
	md->rho             = rho;
	md->rrs_noise       = rrs_noise;
	md->wavelengths     = wavelengths;
	md->raw_wavelengths = raw_wavelengths;
	md->a_0             = a_0;
	md->a_1             = a_1;
	md->a_w             = a_w;
	md->b_bw            = b_bw;
	md->K               = K;
	md->theta_view      = theta_view;
	md->theta_sun       = theta_sun;
	md->sec_theta_view  = sec_theta_view;
	md->sec_theta_sun   = sec_theta_sun;
	md->H_tide          = h_tide; 
	md->empirical_depth_present = empirical_depth_present;

	// Set the table of bottom reflectances. 

	allocate_double_array_2d(&bottom_types, n_bottom_table, n_ref_wlens);

	for (i = 0; i < n_ref_wlens; i++) {
		bottom_types[0][i] = bottom_type_sand[i];
		// bottom_types[0][i] = bottom_sand_Lee[i];
		// bottom_types[0][i] = bottom_type_sand66_green_algea33[i];
		// bottom_types[0][i] = bottom_type_sand50_green_algea50[i];
		bottom_types[1][i] = bottom_type_seagrass[i];
		bottom_types[2][i] = bottom_type_coral[i];
		bottom_types[3][i] = bottom_type_macrophytes[i];
		bottom_types[4][i] = bottom_type_dark_sediment[i];
		bottom_types[5][i] = bottom_type_coral_sand[i];
		bottom_types[6][i] = bottom_type_green_algae[i];
		bottom_types[7][i] = bottom_type_red_algae[i];
	}

	strcpy(md->B_type_name[0], "sand");
	strcpy(md->B_type_name[1], "seagrass");
	strcpy(md->B_type_name[2], "coral");
	strcpy(md->B_type_name[3], "macrophytes");
	strcpy(md->B_type_name[4], "dark sediment");
	strcpy(md->B_type_name[5], "coral_sand");
	strcpy(md->B_type_name[6], "green_algae");
	strcpy(md->B_type_name[7], "red algae");

	// Allocate memory for variable model parameters. 

	allocate_double_array_2d(&Rrs440, max_n_regions, md->n_scenes);
	allocate_double_array_2d(&Rrs490, max_n_regions, md->n_scenes);
	allocate_double_array_2d(&Rrs550, max_n_regions, md->n_scenes);
	allocate_double_array_2d(&Rrs640, max_n_regions, md->n_scenes);
	allocate_double_array_2d(&Rrs750, max_n_regions, md->n_scenes);

	md->Rrs440 = Rrs440; 
	md->Rrs490 = Rrs490;
	md->Rrs550 = Rrs550;
	md->Rrs640 = Rrs640;
	md->Rrs750 = Rrs750;

	// Compute a_w(640).

	md->a_w640 = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, 640.0);

	printf("\na_w(640) = %.5f\n", md->a_w640);

	// Set the tidal heights. 

	printf("\nh_tide = \n\t");

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		h_tide[k_scene] = scene_data[scene_indexes[k_scene]].H_tide;
		printf("%.2f\t", h_tide[k_scene]);
	}

	// Set wavelengths from the input model bands. 

	printf("\nlambda = \n");

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		printf("\t");
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
			wavelengths[k_scene][k_band] = scene_data[scene_indexes[k_scene]].wavelengths[k_band];
			raw_wavelengths[k_scene][k_band] = wavelengths[k_scene][k_band];
			printf("%.2f\t", wavelengths[k_scene][k_band]);
		}
		printf("\n");
	}

	// Set the subsurface solar zenith angle (theta_sun), and the subsurface viewing angle from nadir (theta_view).

	printf("\ntheta_view = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		md->theta_view[k_scene] = scene_data[scene_indexes[k_scene]].theta_v;
		printf("%.3f\t", md->theta_view[k_scene]);
		md->theta_view[k_scene] *= PI/180.0; // degrees 2 radians.
		md->sec_theta_view[k_scene] = 1.0/cos(md->theta_view[k_scene]);
	}

	printf("\ntheta_sun = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		md->theta_sun[k_scene] = scene_data[scene_indexes[k_scene]].theta_w;
		printf("%.3f\t", md->theta_sun[k_scene]);
		md->theta_sun[k_scene] *= PI/180.0; // degrees 2 radians.
		md->sec_theta_sun[k_scene] = 1.0/cos(md->theta_sun[k_scene]);
	}

	// Compute a_0, a_1, a_w, b_bw via interpolation at the spectral wavelengths.

	printf("\na_0 = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
		  a_0[k_scene][k_band] = interp_1d(ref_wlens, a0_Lee, n_ref_wlens, wavelengths[k_scene][k_band]);
		  // a_0[k_scene][k_band] = interp_1d(ref_wlens, a0_Blake, n_ref_wlens, wavelengths[k_scene][k_band]);
			printf("%.5f\t", a_0[k_scene][k_band]);
		}
		printf("\n\t");
	}

	printf("\na_1 = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
		  a_1[k_scene][k_band] = interp_1d(ref_wlens, a1_Lee, n_ref_wlens, wavelengths[k_scene][k_band]);
		  // a_1[k_scene][k_band] = interp_1d(ref_wlens, a1_Blake, n_ref_wlens, wavelengths[k_scene][k_band]);
			printf("%.5f\t", a_1[k_scene][k_band]);
		}
		printf("\n\t");
	}

	printf("\nb_bw = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
			b_bw[k_scene][k_band] = interp_1d(ref_wlens, bbw_Morel_1974, n_ref_wlens, wavelengths[k_scene][k_band]);
			printf("%.5f\t", b_bw[k_scene][k_band]);
		}
		printf("\n\t");
	}

	printf("\na_w = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
			a_w[k_scene][k_band] = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, wavelengths[k_scene][k_band]);
			printf("%.5f\t", a_w[k_scene][k_band]);
		}
		printf("\n\t");
	}

	rrs_noise_present = true;
	printf("\nrrs_noise = \n\t");
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
			rrs_noise[k_scene][k_band] = scene_data[scene_indexes[k_scene]].R_sigma[k_band];
			printf("%.7f\t", rrs_noise[k_scene][k_band]);
			if (rrs_noise[k_scene][k_band] < 1.0e-6) {
				rrs_noise_present = false;
			}
		}
		printf("\n\t");
	}

	// Allocate memory for the different bottom types. 

	max_bottom_types = md->n_bottoms; 
	allocate_double_array_3d(&bottom_reflectance, max_bottom_types, md->n_scenes, md->max_n_bands);

	// Interpolate bottom types across input spectrum for each scene. 

	for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
		for (i = 0; i < md->n_scenes; i++) {
			for (j = 0; j < md->n_bands[i]; j++) {
				bottom_reflectance[k_bottom][i][j] = interp_1d(ref_wlens, bottom_types[k_bottom], n_ref_wlens, md->wavelengths[i][j]);
			}
		}
	}

	md->bottom_reflectance = bottom_reflectance;

	// Allocate memory for model derived gridded datasets. 

	allocate_float_array_3d(&K_coastal, md->n_scenes, nrows, ncols);
	allocate_float_array_3d(&K_blue,    md->n_scenes, nrows, ncols);
	allocate_float_array_3d(&K_green,   md->n_scenes, nrows, ncols);
	allocate_float_array_3d(&K_red,     md->n_scenes, nrows, ncols);

	allocate_float_array_3d(&P, md->n_scenes, nrows, ncols);
	allocate_float_array_3d(&G, md->n_scenes, nrows, ncols);
	allocate_float_array_3d(&X, md->n_scenes, nrows, ncols);
	allocate_float_array_3d(&D, md->n_scenes, nrows, ncols);

	PP = malloc(md->n_scenes*sizeof(double));
	GG = malloc(md->n_scenes*sizeof(double));
	XX = malloc(md->n_scenes*sizeof(double));
	DD = malloc(md->n_scenes*sizeof(double));

	md->P = PP;
	md->G = GG;
	md->X = XX;
	md->D = DD;

	prev = malloc(4*md->n_scenes*sizeof(double)); // P,G,X,D for each scene.
	md->prev = prev;
	md->Y_scene = calloc(md->n_scenes, sizeof(double));

	printf("\n\n");

	// Initialise lookup table. 

#if LUT
	n_lut = 16000; // 4096, 2048, 256
	printf("\nLookup table size = %d\n", n_lut);
	md_lut = (struct model_data_struct*) malloc(n_lut*sizeof(struct model_data_struct));
	lut_in_use = (bool*) malloc(n_lut*sizeof(bool));
	lut_index = (int*) malloc(n_lut*sizeof(int));

	for (i = 0; i < n_lut; i++) {
		lut_in_use[i] = false;
		lut_index[i] = 0;
	}
#endif

	// Count number of optically shallow pixels with positive rrs. 

	n_shallow_pre = 0;
	n_complete = 0;
	nprev = 0;
	n_converged = 0;
	n_not_converged = 0;
	n_converged = 0;
	n_not_converged = 0; 

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {

			// Check for no data in each spectral band. 

			nodata = false;

			for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
				for (k_band = 0; k_band < md->n_raw_bands[k_scene]; k_band++) {
					k = scene_data[scene_indexes[k_scene]].band_indexes[k_band];
					refl = gridded_data[k].array[i][j];
					if (approx_equal(refl, gridded_data[k].nodata_value, 1.0e-6) || refl < 0.0) {
						nodata = true;
						break ;
					}
				}
				if (nodata) break ;
			}

			if (! nodata) {
				n_shallow_pre++;
			}
		}
	}

	printf("\nNumber of optically shallow pixels = %d\n", n_shallow_pre);

	// Determine the order to model pixels. 

	printf("\nSorting depths...");

	allocate_int_array_2d(&indexes, nrows, ncols);
	allocate_float_array_2d(&spectral_index, nrows, ncols);

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {

			if (empirical_depth_present) {
				spectral_index[i][j] = empirical_depths.array[i][j];
			} else {
				green_sum = 0.0;
				blue_sum = 0.0;
				red_sum = 0.0;

				for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {

					k = scene_data[scene_indexes[k_scene]].band_indexes[k_blue];				
					// blue = gridded_data[k].array[i][j];
					blue = smoothed_index_2d(gridded_data[k].array, i, j, nrows, ncols, 
						n_smoothing_radius, gridded_data[k].nodata_value);
					if (! approx_equal(blue, gridded_data[k].nodata_value, 1.0e-6)) {
						blue_sum += blue;
					}

					k = scene_data[scene_indexes[k_scene]].band_indexes[k_green];
					// green = gridded_data[k].array[i][j];
					green = smoothed_index_2d(gridded_data[k].array, i, j, nrows, ncols, 
						n_smoothing_radius, gridded_data[k].nodata_value);				
					if (! approx_equal(green, gridded_data[k].nodata_value, 1.0e-6)) {
						green_sum += green;
					}

					k = scene_data[scene_indexes[k_scene]].band_indexes[k_red];
					// red = gridded_data[k].array[i][j];
					red = smoothed_index_2d(gridded_data[k].array, i, j, nrows, ncols, 
						n_smoothing_radius, gridded_data[k].nodata_value);
					if (! approx_equal(red, gridded_data[k].nodata_value, 1.0e-6)) {
						red_sum += red;
					}
				}

				// Norm. 

				spectral_index[i][j] = sqrt(pow(blue_sum, 2) + pow(green_sum, 2) + pow(red_sum, 2));

				// Spectral angle. 

				if (rrs_noise_present) {

					blue_noise  = 0.0;
					green_noise = 0.0;
					red_noise   = 0.0;

					for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
						blue_noise  += rrs_noise[k_scene][k_blue];
						green_noise += rrs_noise[k_scene][k_green];
						red_noise   += rrs_noise[k_scene][k_red];
					}

					spectral_angle = blue_sum*blue_noise + green_sum*green_noise + red_sum*red_noise;
					spectral_angle /= spectral_index[i][j]*sqrt(pow(blue_noise, 2) + pow(green_noise, 2) + pow(red_noise, 2));
					spectral_angle = acos(spectral_angle);
					spectral_index[i][j] *= spectral_angle;
				}
			}
		}
	}

	sort_indexes_2d(spectral_index, nrows, ncols, indexes);

	free_float_array_2d(spectral_index, nrows);
	printf("...finished sorting depths.\n");

	// Exclude optically deep indexes. 

	allocate_int_array_2d(&shallow_indexes, n_shallow_pre, 2);

	int k_shallow = 0;
	for (ir = 0; ir < nrows; ir++) {
		for (jr = 0; jr < ncols; jr++) {

			i = floor(indexes[ir][jr]/ncols);
			j = indexes[ir][jr]%ncols;

			// Check for no data in each spectral band. 

			nodata = false;

			for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
				for (k_band = 0; k_band < md->n_raw_bands[k_scene]; k_band++) {
					k = scene_data[scene_indexes[k_scene]].band_indexes[k_band];
					refl = gridded_data[k].array[i][j];
					if (approx_equal(refl, gridded_data[k].nodata_value, 1.0e-6) || refl < 0.0) {
						nodata = true;
						break ;
					}
				}
				if (nodata) break ;
			}

			if (! nodata) {
				shallow_indexes[k_shallow][0] = i;
				shallow_indexes[k_shallow][1] = j;
				k_shallow++;
			}

			if (k_shallow > n_shallow_pre) {
				break ;
			}
		}
	}

	n_shallow = k_shallow;
	// printf("\nn_shallow = %d", n_shallow);

	// Set defaults. 

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			depth[i][j]               = 0.0;
			model_error[i][j]         = 0.0;
			bottom_albedo[i][j]       = 0.0;
			bottom_sand[i][j]         = -9999.0;
			bottom_seagrass[i][j]     = -9999.0;
			bottom_coral[i][j]        = -9999.0;
			K_min[i][j]               = 0.0;
			bottom_type[i][j]         = -9999.0;
			index_optical_depth[i][j] = 0.0;
	  		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
  				K_coastal[k_scene][i][j] = 0.0;
  				K_blue[k_scene][i][j]    = 0.0;
  				K_green[k_scene][i][j]   = 0.0;
  				K_red[k_scene][i][j]     = 0.0;
  				P[k_scene][i][j]         = 0.0;
  				G[k_scene][i][j]         = 0.0;
  				X[k_scene][i][j]         = 0.0;
  				D[k_scene][i][j]         = 0.0;
  			}
		}
	}

	// Run model at each optically shallow pixel. 

	t0 = clock();
	t1 = clock();

	total_model_error    = 0.0;
	total_depth_error    = 0.0;
	total_bottom_error   = 0.0;
	total_Rrs_error      = 0.0;
	total_sam_error      = 0.0;
	total_K_error        = 0.0;
	total_depth          = 0.0;
	total_kmin           = 0.0;
	total_start_depth    = 0.0;
	total_bottom_albedo  = 0.0;
	n_modelled           = 0.0;
	n_skipped            = 0.0;
	n_restarted          = 0.0;
	total_n_iterations   = 0.0;
	total_B_sand         = 0.0;
	total_B_seagrass     = 0.0;
	total_B_coral        = 0.0;
	total_P              = 0.0;
	total_G              = 0.0;
	total_X              = 0.0;
	total_index_optical_depth = 0.0;
	total_spectral_diff       = 0.0;

	md->start_at_previous = false;
	md->depth_prev  = 0.0;
	md->K_min_prev  = 0.0;
	md->model_error = 100.0;
	total_lut       = 0;
	n_done          = 0;
	n_modelled_optimisation = 0;

	min_preferable_skip      = 0.95; // Skip ratio threshold before increasing LUT match. 
	max_preferable_skip      = 0.9925; // Skip ratio threshold before decreasing LUT match. 
	skipped_ratio            = (min_preferable_skip + max_preferable_skip)/2.0; // Initialisation. 
	// max_lut_match_threshold  = 2.5 + MIN(5.0, 2.5*((double) md->n_scenes));
	// max_lut_match_threshold  = 7.5*((double) md->n_scenes);

	max_lut_match_threshold = 5.0*((double) md->n_scenes);

#if SAM
	// min_lut_match_threshold  = 0.25*sqrt((double) md->n_scenes);
	min_lut_match_threshold  = 0.25;
#else 
	// min_lut_match_threshold  = 2.0*sqrt((double) md->n_scenes);
	min_lut_match_threshold  = 1.5;
#endif

	max_error_before_restart = 5.0;
	lut_match_threshold      = MAX(1.5, 1.125*((double) md->n_scenes)); // RMS relative % difference.
	lut_entry_threshold      = lut_match_threshold;

	error_list = malloc(nrows*ncols*sizeof(double));

#if _OPENMP
	#pragma omp for schedule(dynamic)
#endif
	for (k_shallow = 0; k_shallow < n_shallow; k_shallow++) {

			// printf("\nthread_id = %d, complete = %d", thread_id, n_done++);

			i = shallow_indexes[k_shallow][0];
			j = shallow_indexes[k_shallow][1];

			md->i = i;
			md->j = j;

			// printf("\ni,j = %d, %d\n", i, j);

			// Set defaults. 

			for (k_region = 0; k_region < max_n_regions; k_region++) {
				for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
					for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
						Rrs_measured[k_region][k_scene][k_band] = 0.0;
						Rrs_modelled[k_region][k_scene][k_band] = 0.0;
						rrs_modelled[k_region][k_scene][k_band] = 0.0;
						rrs_bottom[k_region][k_scene][k_band]   = 0.0;
						rrs_dp[k_region][k_scene][k_band]       = 0.0;
						rho[k_region][k_scene][k_band]          = 0.0;
					}
				}
			}

			// Check for no data in each spectral band. 

			nodata = false;

			for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
				for (k_band = 0; k_band < md->n_raw_bands[k_scene]; k_band++) {
					k = scene_data[scene_indexes[k_scene]].band_indexes[k_band];
					refl = gridded_data[k].array[i][j];
					if (approx_equal(refl, gridded_data[k].nodata_value, 1.0e-6) || refl < 0.0) {
						nodata = true;
						break ;
					}
				}
				if (nodata) break ;
			}

			if (nodata) continue ;

			// Extract Rrs within the region of each band of each scene. 

			extract_Rrs_data(i, j, scene_data, gridded_data, scene_indexes, 
				nscenes, n_spatial, n_smoothing_radius, nrows, ncols, md, 0.0);

			if (md->n_regions == 0) {
				continue ;
			}

			// Empirical depth. 

			if (empirical_depth_present) {
				if (! approx_equal(empirical_depths.array[i][j], empirical_depths.nodata_value, 1.0e-6)) {
					md->empirical_depth_present = true;
					if (empirical_depths.array[i][j] > -1.0) { 
						md->h_empirical = 1.0; // Very shallow. 
					} else {
						md->h_empirical = fabs(empirical_depths.array[i][j]);
					}
				} else {
					// continue ; // Temporary.
					md->empirical_depth_present = false;
					md->start_at_previous = false;
					md->h_empirical = 0.0; // Unknown: fit all bottom types.
				}
			} else {
			  // continue ; // Temporary. 
				md->empirical_depth_present = false;
				md->h_empirical = 0.0; // Unknown: fit all bottom types.
			}

			// Check LUT. 

#if LUT
			continuity_threshold = 0.1; // 10% relative difference. 
			best_lut_spectral_match = 1.0e4;
			lut_model_error = 1.0e4;
			lut_position = 0;

			for (k_lut = 0; k_lut < n_lut; k_lut++) {
				if (lut_in_use[k_lut]) {
					if (md_lut[k_lut].n_regions == md->n_regions && 
						fabs(md_lut[k_lut].depth - md->depth)/(md->depth + epsilon) < continuity_threshold && 
						fabs(md_lut[k_lut].K_min - md->K_min)/(md->K_min + epsilon) < continuity_threshold) {
						spectral_diff = model_spectral_diff_mean(Rrs_measured, md_lut[k_lut].Rrs_modelled, md);
						if (spectral_diff < best_lut_spectral_match) {
							best_lut_spectral_match = spectral_diff;
							lut_model_error = md_lut[k_lut].model_error;
							lut_position = k_lut;
						}
					}
				}
			}

			// printf("\nbest_lut_spectral_match = %.2f, lut_position = %d", 
			//		best_lut_spectral_match, lut_position);
			// fflush(stdout);
#endif
			
			// Decide to model or skip. 

			if (best_lut_spectral_match > lut_match_threshold) {

				// Call the core modelling routine. 

				samodel_optimise(md);

#if RESTART
				if (md->start_at_previous && md->model_error > max_error_before_restart) { 
					md->start_at_previous = false;

					samodel_optimise(md);

					if (md->model_error > max_error_before_restart) {
						md->start_at_previous = false;
					} else {
						md->start_at_previous = true;
					}
				}

				error_list[n_modelled_optimisation++] = md->model_error;

				if (! md->start_at_previous) {
					n_restarted += 1.0;
				}
#else
				error_list[n_modelled_optimisation++] = md->model_error;
#endif
				// Enter model output into oldest entry in LUT. 
#if LUT
				if (md->model_error < lut_entry_threshold && md->depth < 32.0) {
					total_lut += 1;
					oldest_position = INT_MAX;
					oldest_index = 0;
					for (k_lut = 0; k_lut < n_lut; k_lut++) {
						if (! lut_in_use[k_lut]) {
							oldest_index = k_lut;
							break ;
						}
						if (lut_index[k_lut] < oldest_position) {
							oldest_position = lut_index[k_lut];
							oldest_index = k_lut;
						}
					}
					// printf("\noldest_index = %d", oldest_index);
					// fflush(stdout);

					if (lut_in_use[oldest_index]) {
						partial_copy_md(md, &md_lut[oldest_index], false);
					} else {
						// printf("\nallocating index %d\n", oldest_index);
						partial_copy_md(md, &md_lut[oldest_index], true);
					}
					
					lut_index[oldest_index] = k_shallow;
					lut_in_use[oldest_index] = true;
				}
#endif
			} else { 
#if LUT
				// Use LUT (skip re-modelling).

				md->depth_prev = md->depth;
				md->K_min_prev = md->K_min;
				partial_copy_md(&md_lut[lut_position], md, false);
				md->model_error = best_lut_spectral_match;
				n_skipped += 1.0;
#else

				// No LUT. 

				samodel_optimise(md);

#if RESTART
				if (md->start_at_previous && md->model_error > max_error_before_restart) {
					md->start_at_previous = false;

					samodel_optimise(md);

					if (md->model_error > max_error_before_restart) {
						md->start_at_previous = false;
					} else {
						md->start_at_previous = true;
					}

					error_list[n_modelled_optimisation++] = md->model_error;
				}
				
				if (! md->start_at_previous) {
					n_restarted += 1.0;
				}
#endif
#endif
			}

			// Decide whether to restart model initialisation parameters. 

			if (fabs(md->depth - md->depth_prev)/(epsilon + md->depth_prev) < continuity_threshold && 
				fabs(md->K_min - md->K_min_prev)/(epsilon + md->K_min_prev) < continuity_threshold) {
				md->start_at_previous = true;
			} else {
				md->start_at_previous = false;
			}

			// Remember last modelled. 

			md->prev_i = i;
			md->prev_j = j;
			md->depth_prev = md->depth;
			md->K_min_prev = md->K_min;

			// Store model data. 

			depth[i][j]               = (float) md->depth;  
			model_error[i][j]         = (float) md->Rrs_error;
			bottom_albedo[i][j]       = (float) md->bottom_albedo;
			bottom_sand[i][j]         = (float) md->B_type_percent[0];
			bottom_seagrass[i][j]     = (float) md->B_type_percent[1];
			bottom_coral[i][j]        = (float) md->B_type_percent[2];
			K_min[i][j]               = (float) md->K_min;
			index_optical_depth[i][j] = (float) md->index_optical_depth;
			bottom_type[i][j]         = (float) md->bottom_type; 

			total_model_error         += md->model_error;
			total_Rrs_error           += md->Rrs_error;
			total_sam_error           += md->SAM_error;
			total_K_error             += md->K_error;
			total_depth_error         += md->depth_error;
			total_bottom_error        += md->bottom_error;
			total_depth               += md->depth;
			total_kmin                += md->K_min;
			total_start_depth         += md->h_empirical;
			total_bottom_albedo       += md->bottom_albedo;
			total_n_iterations        += (double) md->n_iterations;
			total_B_sand              += (double) md->B_type_percent[0];
			total_B_seagrass          += (double) md->B_type_percent[1];
			total_B_coral             += (double) md->B_type_percent[2];
			total_index_optical_depth += md->index_optical_depth;
			if (best_lut_spectral_match < 100.0) {
				total_spectral_diff       += best_lut_spectral_match;
			}

			mean_P = 0.0;
			mean_G = 0.0;
			mean_X = 0.0;
  			for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
  				K_coastal[k_scene][i][j] = (float) md->K[k_scene][k_coastal];
  				K_blue[k_scene][i][j]    = (float) md->K[k_scene][k_blue];
  				K_green[k_scene][i][j]   = (float) md->K[k_scene][k_green];
  				K_red[k_scene][i][j]     = (float) md->K[k_scene][k_red];
  				P[k_scene][i][j]         = (float) md->P[k_scene];
  				G[k_scene][i][j]         = (float) md->G[k_scene];
  				X[k_scene][i][j]         = (float) md->X[k_scene];
  				D[k_scene][i][j]         = (float) md->D[k_scene];
  				mean_P += md->P[k_scene];
  				mean_G += md->G[k_scene];
  				mean_X += md->X[k_scene];
  			}

  			mean_P /= (double) md->n_scenes;
  			mean_G /= (double) md->n_scenes;
  			mean_X /= (double) md->n_scenes;

			total_P += mean_P;
			total_G += mean_G;
			total_X += mean_X;

			n_modelled += 1.0;

			if (md->converged) {
				n_converged++;
			} else {
				n_not_converged++;
			}

			// Update LUT adaptive threshold.  

			nperc = (int) 1000.0*((double) ++n_complete)/((double) n_shallow);

			if (nperc > nprev) {

				if (n_modelled_optimisation > 0) {
					error_mean   = vec_mean2_double(error_list, n_modelled_optimisation, 0.0);
					error_stddev = vec_stddev_double(error_list, n_modelled_optimisation, 0.0, error_mean); 
					// printf("\n\nn_modelled_optimisation, error_mean, error_stddev = %d\t%f\t%f\n\n", n_modelled_optimisation, error_mean, error_stddev);

					lut_match_threshold += 0.75*error_mean + 0.25*total_model_error/n_modelled;
					lut_match_threshold /= 2.0;
					
					// lut_entry_threshold = error_mean - 1.645*error_stddev; // 1.645 - 90% CI

					lut_entry_threshold = lut_match_threshold;
				}	

				skipped_ratio = n_skipped/n_modelled;
				if (skipped_ratio > max_preferable_skip) {
					lut_match_threshold *= 0.75;
					lut_entry_threshold *= 0.75;
				}

				lut_match_threshold = MAX(min_lut_match_threshold, lut_match_threshold);
				lut_entry_threshold = MAX(min_lut_match_threshold, lut_entry_threshold);

				// lut_match_threshold = MIN(max_lut_match_threshold, lut_match_threshold);
				// lut_entry_threshold = MIN(max_lut_match_threshold, lut_entry_threshold);

				// Display summary of model performance. 

				nprev = nperc;
				nperc /= 10.0;
				t2 = clock();
				elapsed = ((double) t2 - t1)/CLOCKS_PER_SEC; // elapsed time since last display.
				pixels_per_sec = ((double) n_modelled)/elapsed;
				elapsed = ((double) t2 - t0)/CLOCKS_PER_SEC;  // total elapsed time.
				time_remaining = elapsed*(((double) n_shallow)/((double) n_complete) - 1.0)/3600.0; 
				printf("\nPERCENTAGE COMPLETE = %3d%% (%6.2f px/sec) [%.2f hrs remaining]\n\t\
Mean Error          = %7.2f (%%)\n\t\
Mean Error Rrs      = %7.2f (%%)\n\t\
Mean Error SAM      = %7.2f     \n\t\
Mean Error K        = %7.2f     \n\t\
Mean Error H        = %7.2f (%%)\n\t\
Mean Error B        = %7.2f (%%)\n\t\
Mean H              =  \033[1m\033[31m %5.2f\033[0m (m)\n\t\
Mean H start        =   %5.2f\n\t\
Mean P              =   %6.4f \n\t\
Mean G              =   %6.4f \n\t\
Mean X              =   %6.4f \n\t\
Mean B              =   %4.3f \n\t\
sand                =   %5.2f (%%)\n\t\
seagrass            =   %5.2f (%%)\n\t\
coral               =   %5.2f (%%)\n\t\
Mean K min          =   %5.4f \n\t\
Mean IOD            =   %5.2f (%%)\n\t\
Mean iterations     =   %5.0f \n\t\
diverged            =      %2.0f (%%)\n\t\
Mean spectral diff  =   %5.2f (%%)\n\t\
skipped             =   %5.2f (%%)\n\t\
restarted           =   %5.2f (%%)\n\t", 
					nperc, 
					pixels_per_sec,
					time_remaining,
					total_model_error/n_modelled,
					total_Rrs_error/n_modelled,
					total_sam_error/n_modelled,
					total_K_error/n_modelled,
					total_depth_error/n_modelled,
					total_bottom_error/n_modelled,
					total_depth/n_modelled,
					total_start_depth/n_modelled,
					total_P/n_modelled,
					total_G/n_modelled,
					total_X/n_modelled,
					total_bottom_albedo/n_modelled, 
					total_B_sand/n_modelled, 
					total_B_seagrass/n_modelled, 
					total_B_coral/n_modelled, 
					total_kmin/n_modelled,
					total_index_optical_depth/n_modelled,
					total_n_iterations/n_modelled,
					100.0*((double) n_not_converged)/((double) (n_converged + n_not_converged)), 
					total_spectral_diff/n_modelled,
					100.0*n_skipped/n_modelled,
					(((int) n_restarted == 0) || ((int) n_modelled) == ((int) n_skipped)) ? 0.0 : 100.0*n_restarted/((double) n_modelled - n_skipped));
#if LUT
				// Display LUT summary. 

				n_obs = 0.0;
				lut_mean_age = 0.0;
				lut_mean_H = 0.0;
				lut_min_H = 100.0;
				lut_max_H = 0.0;
				lut_mean_K_min = 0.0;
				lut_min_K_min = 100.0;
				lut_max_K_min = 0.0;

				for (k_lut = 0; k_lut < n_lut; k_lut++) {
					if (lut_in_use[k_lut]) {

						lut_mean_H += md_lut[k_lut].depth;
						lut_mean_K_min += md_lut[k_lut].K_min;
						lut_mean_age += (double) k_shallow - lut_index[k_lut];
						n_obs += 1.0;

						if (md_lut[k_lut].depth < lut_min_H) {
							lut_min_H = md_lut[k_lut].depth;
						}
						if (md_lut[k_lut].depth > lut_max_H) {
							lut_max_H = md_lut[k_lut].depth;
						}
						if (md_lut[k_lut].K_min < lut_min_K_min) {
							lut_min_K_min = md_lut[k_lut].K_min;
						}
						if (md_lut[k_lut].K_min > lut_max_K_min) {
							lut_max_K_min = md_lut[k_lut].K_min;
						}
					}
				}

				if (n_obs > 0.5) {
					printf(
"\n\tLUT entries         = %d\n\t\
LUT match threshold = %.2f\n\t\
LUT entry threshold = %.2f\n\t\
mean age            = %.0f\n\t\
min/mean/max H      = %.1f %.1f %.1f\n\t\
min/mean/max min(K) = %.4f %.4f %.4f\n", 
						total_lut, 
						lut_match_threshold,
						lut_entry_threshold,
						lut_mean_age/n_obs,
						lut_min_H,
						lut_mean_H/n_obs, 
						lut_max_H,
						lut_min_K_min,
						lut_mean_K_min/n_obs,
						lut_max_K_min);
				}
#endif
				fflush(stdout);

				total_model_error = 0.0;
				total_Rrs_error = 0.0;
				total_sam_error = 0.0;
				total_K_error = 0.0;
				total_depth_error = 0.0;
				total_bottom_error = 0.0;
				total_depth = 0.0;
				total_kmin = 0.0;
				total_start_depth = 0.0;
				total_bottom_albedo = 0.0;
				total_n_iterations = 0.0;
				total_B_sand = 0.0;
				total_B_seagrass = 0.0;
				total_B_coral = 0.0;
				total_P = 0.0;
				total_G = 0.0;
				total_X = 0.0;
				total_index_optical_depth = 0.0;
				n_modelled = 0.0;
				n_converged = 0;
				n_not_converged = 0;
				total_spectral_diff = 0.0;
				n_skipped = 0.0;
				n_restarted = 0.0;
				n_modelled_optimisation = 0;
				t1 = clock();
			}

			// Plotting (debugging only)

#if PLOTTING
			samodel_graphics(md, pagesize, linewidth);

			printf("\nPress any key to continue...");
			fgets(string, 4, stdin); 
			printf("\n");
			exit(1);
			}
#endif
	}

	//
	//	Compute depth error estimate. 
	// 

	// This is an hoc estimate of the depth error. A thorough estimate would require
	// modelling each pixel many times with an appropriate pertubation to the reflectance
	// as dictated by the variance of the reflectance over optically deep water. 

	printf("\nComputing depth error estimates...");

	n_samples = 128;
	n_pertubations = 16;
	n_total_per_interval = n_samples*n_pertubations;
	n_trials = (int) sqrt(nrows*ncols);
	depth_interval = 0.25;

	max_depth_reached = array_max2(depth, nrows, ncols, 0.0);
	max_depth_reached = depth_interval*((int) max_depth_reached/depth_interval);
	max_depth_reached = MIN(max_depth_reached, 30.0);
	n_depth_intervals = (int) max_depth_reached/depth_interval;

	trial_depths = malloc(n_samples*sizeof(double));
	depth_sigma_table = malloc(n_depth_intervals*sizeof(double));

	// Iterate over all depth intervals. 

	md->start_at_previous = false;
	k_depth = 0;
	for (d = 0.0; d < max_depth_reached; d += depth_interval) {

		printf("\nDepth interval %.2f - %.2f\n", d, d + depth_interval);

		// Create samples within the interval. 

		for (k_sample = 0; k_sample < n_samples; k_sample++) {

			// Randomly find a depth within the interval. 

			within_interval = false;

			for (k_trial = 0; k_trial < n_trials; k_trial++) {
				i = random_in_range(0, nrows);
				j = random_in_range(0, ncols);

				if (depth[i][j] > d && depth[i][j] < d + depth_interval) {
					within_interval = true; 
					break ;
				}
			}

			if (! within_interval) {
				trial_depths[k_sample] = 0.0;
				continue ;
			}

			// Pertubate and remodel. 

			n_sigma = frand2(1.0);

			extract_Rrs_data(i, j, scene_data, gridded_data, scene_indexes, 
				nscenes, n_spatial, n_smoothing_radius, nrows, ncols, md, n_sigma);

			if (md->n_regions == 0) {
				continue ;
			}

			// Empirical depth. 

			if (empirical_depth_present) {
				if (! approx_equal(empirical_depths.array[i][j], empirical_depths.nodata_value, 1.0e-6)) {
					if (empirical_depths.array[i][j] > -1.0) { 
						md->h_empirical = 1.0; // Very shallow. 
					} else {
						md->h_empirical = fabs(empirical_depths.array[i][j]);
					}
				} else {
					trial_depths[k_sample] = 0.0;
					continue ;
				}
			} else {
				md->h_empirical = 5.0;
			}

			// Call the core modelling routine. 

			samodel_optimise(md);
			md->start_at_previous = true;

			trial_depths[k_sample] = md->depth;
		}

		mean_trial_depths = vec_mean2_double(trial_depths, n_samples, 0.0);
		depth_sigma_table[k_depth++] = vec_stddev_double(trial_depths, n_samples, 0.0, mean_trial_depths);
	}

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			depth_sigma[i][j] = 0.0;
			if (depth[i][j] > 0.0) {
				k_depth = 0;
				for (d = 0.0; d < max_depth_reached; d += depth_interval) {
					if (depth[i][j] > d && depth[i][j] <= d + depth_interval) {
						depth_sigma[i][j] = depth_sigma_table[k_depth];  
						break ;
					}
					k_depth++;
				}
			}
		}
	}

	printf("...finished.\n");

	free(trial_depths);
	free(depth_sigma_table);

	// NRS uses negatives depths below MSL.

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			depth[i][j] *= -1.0; 
		}
	}

	// Write model data to file. 

	printf("\nWriting model data to file...");

	int n, k_index;
	float *lons, *lats;
	char *scene_name, file_name[64];

  	k_index = scene_data[scene_indexes[0]].band_indexes[0];

  	lons = (float*) malloc(gridded_data[k_index].ncols*sizeof(float));
  	lats = (float*) malloc(gridded_data[k_index].nrows*sizeof(float));

  	for (n = 0; n < gridded_data[k_index].ncols; n++) {
  		lons[n] = gridded_data[k_index].wlon + ((float) n)*gridded_data[k_index].cellsize;
  	}

  	for (n = 0; n < gridded_data[k_index].nrows; n++) {
  		lats[n] = gridded_data[k_index].slat + ((float) n)*gridded_data[k_index].cellsize;
  	}

  	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
  		scene_name = scene_data[scene_indexes[k_scene]].scene_name;

  		sprintf(file_name, "%s_K_coastal.nc", trim(scene_name));

  		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			K_coastal[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

  		sprintf(file_name, "%s_K_blue.nc", trim(scene_name));

  		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			K_blue[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

  		sprintf(file_name, "%s_K_green.nc", trim(scene_name));

  		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			K_green[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

  		sprintf(file_name, "%s_K_red.nc", trim(scene_name));

  		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			K_red[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

  		sprintf(file_name, "%s_P.nc", trim(scene_name));

  		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			P[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

  		sprintf(file_name, "%s_G.nc", trim(scene_name));

  		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			G[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

		sprintf(file_name, "%s_X.nc", trim(scene_name));

		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			X[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);

		sprintf(file_name, "%s_D.nc", trim(scene_name));

#if DELTA
		printf("\nWriting file %s...", file_name);
		write_nc(file_name, 
			D[k_scene], 
			gridded_data[k_index].ncols,
			gridded_data[k_index].nrows,
			lons, lats, 
			gridded_data[k_index].nodata_value);
#endif
	}

	sprintf(file_name, "modelled_H.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		depth, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_H_sigma.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		depth_sigma, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_error.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		model_error, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_albedo.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		bottom_albedo, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_bottom_sand.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		bottom_sand, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_bottom_seagrass.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		bottom_seagrass, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_bottom_coral.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		bottom_coral, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_min_K.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		K_min, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_bottom_type.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		bottom_type, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	sprintf(file_name, "modelled_index_optical_depth.nc");
	printf("\nWriting file %s...", file_name);
	write_nc(file_name, 
		index_optical_depth, 
		gridded_data[k_index].ncols,
		gridded_data[k_index].nrows,
		lons, lats, 
		gridded_data[k_index].nodata_value);

	printf("\n... finished.\n");

	free(lons);
	free(lats);
	free_float_array_3d(K_coastal, nscenes, nrows);
	free_float_array_3d(K_blue,    nscenes, nrows);
	free_float_array_3d(K_green,   nscenes, nrows);
	free_float_array_3d(K_red,     nscenes, nrows);
	
	free_int_array_2d(indexes, nrows);
	free_int_array_2d(shallow_indexes, n_shallow_pre);

	free(PP);
	free(GG);
	free(XX);
	free(DD);

	free(prev);
	free(error_list);

	free_float_array_3d(P, nscenes, nrows);
	free_float_array_3d(G, nscenes, nrows);
	free_float_array_3d(X, nscenes, nrows);
	free_float_array_3d(D, nscenes, nrows);

	printf("\n\n");

	free(md);
#if LUT
	free(md_lut);
	free(lut_index);
	free(lut_in_use);
#endif

	free(n_bands);
	free(n_raw_bands);
	free(theta_view);
	free(theta_sun);
	free(sec_theta_view);
	free(sec_theta_sun);
	free(h_tide);

	free_double_array_2d(Rrs440, max_n_regions);
	free_double_array_2d(Rrs490, max_n_regions);
	free_double_array_2d(Rrs550, max_n_regions);
	free_double_array_2d(Rrs640, max_n_regions);
	free_double_array_2d(Rrs750, max_n_regions);

	free_double_array_3d(Rrs_measured, max_n_regions, nscenes);
	free_double_array_3d(Rrs_modelled, max_n_regions, nscenes);
	free_double_array_3d(rrs_modelled, max_n_regions, nscenes);
	free_double_array_3d(rrs_bottom,   max_n_regions, nscenes);
	free_double_array_3d(rho,          max_n_regions, nscenes);
	free_double_array_3d(rrs_dp,       max_n_regions, nscenes);

	free_double_array_2d(wavelengths,      nscenes);
	free_double_array_2d(raw_wavelengths,  nscenes);

	free_double_array_2d(a_0, nscenes);
	free_double_array_2d(a_1, nscenes); 
	free_double_array_2d(a_w, nscenes);
	free_double_array_2d(b_bw, nscenes);
	free_double_array_2d(K, nscenes);
	free_double_array_2d(rrs_noise, nscenes);

	free_double_array_2d(bottom_types, n_bottom_table);
	free_double_array_3d(bottom_reflectance, max_bottom_types, nscenes); 

  // End of parallel region.
  }

#if DEBUG
	printf("\nEXIT samodel\n");
#endif
}




// The parameter vector, each element within a box (samodel_bounds):
//
//   [H (n_regions) | B (n_regions*n_bottoms) | ln P, ln G, ln X per scene]
//
// H is depth (m) in [H_MIN, H_MAX]. B_k in [0, 1] is the albedo of bottom type k in the
// region, so the bottom reflectance is the linear mixture sum_k B_k rho_k(lambda); a pure
// bottom lies on the bounds. The region's albedo is sum_k B_k and its fractions are
// B_k/sum_k B_k. P, G and X lie within their bounds. Values are clamped to the box when
// unpacked.

static double clamp(double v, double lo, double hi) {
	return v < lo ? lo : (v > hi ? hi : v);
}

int samodel_n_params(model_data *md) {
	return md->n_regions*(1 + md->n_bottoms) + N_WATER*md->n_scenes;
}

void samodel_bounds(model_data *md, double *lo, double *hi) {

	int k, k_scene, k_band, k_bottom, n_regions = md->n_regions, offset = md->n_regions*(1 + md->n_bottoms);
	double shape_max;

	for (k = 0; k < n_regions; k++) {
		lo[k] = H_MIN;
		hi[k] = H_MAX;
	}
	// A bottom type's albedo is bounded so that its reflectance B_k rho_k(lambda) is at most 1 at
	// every fitted band: the tabulated shapes are normalised near 550 nm, and the red edge of
	// seagrass and coral exceeds 1 (3.3 and 2.9 at 705 nm).

	for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
		shape_max = 1.0;
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
				shape_max = MAX(shape_max, md->bottom_reflectance[md->bottom_type_indexes[k_bottom]][k_scene][k_band]);
			}
		}
		for (k = 0; k < n_regions; k++) {
			lo[n_regions + k*md->n_bottoms + k_bottom] = 0.0;
			hi[n_regions + k*md->n_bottoms + k_bottom] = 1.0/shape_max;
		}
	}
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		lo[offset + N_WATER*k_scene]     = log(P_MIN);
		hi[offset + N_WATER*k_scene]     = log(P_MAX);
		lo[offset + N_WATER*k_scene + 1] = log(G_MIN);
		hi[offset + N_WATER*k_scene + 1] = log(G_MAX);
		lo[offset + N_WATER*k_scene + 2] = log(X_MIN);
		hi[offset + N_WATER*k_scene + 2] = log(X_MAX);
#if Y_MODEL == 3
		lo[offset + N_WATER*k_scene + 3] = Y_MIN;
		hi[offset + N_WATER*k_scene + 3] = Y_MAX;
#endif
	}
}

double samodel_depth(double *params, int k_region) {
	return clamp(params[k_region], H_MIN, H_MAX);
}

void samodel_unpack_region(double *params, model_data *md, int k_region, double *H, double *albedo, double q[]) {

	int k, n_bottoms = md->n_bottoms;
	double sum;

	*H = samodel_depth(params, k_region);

	sum = 0.0;
	for (k = 0; k < n_bottoms; k++) {
		q[k] = clamp(params[md->n_regions + k_region*n_bottoms + k], 0.0, 1.0);
		sum += q[k];
	}
	*albedo = sum;

	for (k = 0; k < n_bottoms; k++) {
		q[k] = sum > 0.0 ? q[k]/sum : 1.0/n_bottoms;
	}
}

void samodel_unpack_scene(double *params, model_data *md, int k_scene, double *P, double *G, double *X) {

	int offset = md->n_regions*(1 + md->n_bottoms);

	*P = exp(clamp(params[offset + N_WATER*k_scene],     log(P_MIN), log(P_MAX)));
	*G = exp(clamp(params[offset + N_WATER*k_scene + 1], log(G_MIN), log(G_MAX)));
	*X = exp(clamp(params[offset + N_WATER*k_scene + 2], log(X_MIN), log(X_MAX)));
#if Y_MODEL == 3
	md->Y_scene[k_scene] = clamp(params[offset + N_WATER*k_scene + 3], Y_MIN, Y_MAX);
#endif
}

// Pack physical values into a parameter vector. q holds n_bottoms fractions per region.

void samodel_pack(model_data *md, double *H, double *albedo, double *q, double *P, double *G, double *X,
	double *params) {

	int k_region, k_scene, k, n_regions = md->n_regions, n_bottoms = md->n_bottoms,
		offset = md->n_regions*(1 + md->n_bottoms);
	for (k_region = 0; k_region < n_regions; k_region++) {
		params[k_region] = clamp(H[k_region], H_MIN, H_MAX);
		for (k = 0; k < n_bottoms; k++) {
			params[n_regions + k_region*n_bottoms + k] = clamp(albedo[k_region]*q[k_region*n_bottoms + k], 0.0, 1.0);
		}
	}

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		params[offset + N_WATER*k_scene]     = log(clamp(P[k_scene], P_MIN, P_MAX));
		params[offset + N_WATER*k_scene + 1] = log(clamp(G[k_scene], G_MIN, G_MAX));
		params[offset + N_WATER*k_scene + 2] = log(clamp(X[k_scene], X_MIN, X_MAX));
#if Y_MODEL == 3
		params[offset + N_WATER*k_scene + 3] = 1.0; // Starting Y.
#endif
	}
}

// Starting point for depth h and bottom albedo: equal bottom fractions, and P, G, X estimated
// from each scene's mean Rrs(490), Rrs(550) and Rrs(640). is_free[i] is set non-zero for every
// parameter (callers clear the ones they hold fixed). Requires md->Rrs490, Rrs550 and
// Rrs640 (set in samodel_optimise).

void samodel_initial(model_data *md, double h, double albedo, double *start, double *is_free) {

	int k_region, k_scene, k, n_regions = md->n_regions, n_bottoms = md->n_bottoms;
	double *H, *A, *q, P[MAX_SCENES], G[MAX_SCENES], X[MAX_SCENES], mean_Rrs490, mean_Rrs550, mean_Rrs640;

	H = malloc(n_regions*sizeof(double));
	A = malloc(n_regions*sizeof(double));
	q = malloc(n_regions*n_bottoms*sizeof(double));

	for (k_region = 0; k_region < n_regions; k_region++) {
		H[k_region] = h;
		A[k_region] = albedo;
		for (k = 0; k < n_bottoms; k++) {
			q[k_region*n_bottoms + k] = 1.0/n_bottoms;
		}
	}

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		mean_Rrs490 = 0.0;
		mean_Rrs550 = 0.0;
		mean_Rrs640 = 0.0;
		for (k_region = 0; k_region < n_regions; k_region++) {
			mean_Rrs490 += md->Rrs490[k_region][k_scene];
			mean_Rrs550 += md->Rrs550[k_region][k_scene];
			mean_Rrs640 += md->Rrs640[k_region][k_scene];
		}
		mean_Rrs490 /= (double) n_regions;
		mean_Rrs550 /= (double) n_regions;
		mean_Rrs640 /= (double) n_regions;

		P[k_scene] = 0.072*pow(mean_Rrs490/mean_Rrs550, -1.7);
		G[k_scene] = 1.5*P[k_scene];
		X[k_scene] = 30.0*md->a_w640*mean_Rrs640;
	}

	samodel_pack(md, H, A, q, P, G, X, start);

	free(H);
	free(A);
	free(q);

	for (k = 0; k < samodel_n_params(md); k++) {
		is_free[k] = 1.0;
	}
}

// Residuals of the modelled Rrs, scaled so that their sum of squares is the square of the
// Rrs error in samodel_error (the RMS difference as a percentage of the mean measured Rrs).
// Residuals are ordered by region, then scene, then band. Returns the number of residuals
// (n_regions*sum(n_bands)).

static int residuals_per_region(model_data *md) {
	int k_scene, n = 0;
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) n += md->n_bands[k_scene];
	return n;
}

static double residual_scale(model_data *md) {
	int k_region, k_scene, k_band, n = 0;
	double mean = 0.0;
	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
				mean += md->Rrs_measured[k_region][k_scene][k_band];
				n++;
			}
		}
	}
	mean /= (double) n;
	return 100.0/(mean*sqrt((double) n));
}

// Side of the square spatial window, or 0 if the regions do not form a full square (nodata).

static int window_side(model_data *md) {
	int side = (int) floor(sqrt((double) md->n_regions) + 0.5);
	return (side*side == md->n_regions && side >= 3) ? side : 0;
}

static int continuity_count(model_data *md) {
#if CONTINUITY
	int side = window_side(md);
	return side ? 2*side*(side - 2) : 0;
#else
	return 0;
#endif
}

// Rrs residual units per chi unit: 100 sigma_rel/sqrt(n), with sigma_rel as in the optically
// deep test.

static double chi_unit(model_data *md) {
	int k_scene, k_band, k_region, n = 0;
	double noise = 0.0, mean = 0.0, sigma_rel;
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
			noise += md->rrs_noise[k_scene][k_band]*md->rrs_noise[k_scene][k_band];
			for (k_region = 0; k_region < md->n_regions; k_region++) {
				mean += md->Rrs_measured[k_region][k_scene][k_band];
			}
			n++;
		}
	}
	mean /= (double) (n*md->n_regions);
	noise = sqrt(noise/n);
	sigma_rel = sqrt(pow(noise/mean, 2.0) + DEEP_MODEL_ERROR*DEEP_MODEL_ERROR);
	return 100.0*sigma_rel/sqrt((double) (n*md->n_regions));
}

static int tie_count(model_data *md) {
#if WATER_TIE
	return md->n_scenes > 1 ? 3*md->n_scenes : 0;
#else
	(void) md;
	return 0;
#endif
}

// Tie residuals (tie_count values) for the water parameters in params.

static void tie_residuals(double *params, model_data *md, double unit, double *r) {
#if WATER_TIE
	int k_scene, w, offset = md->n_regions*(1 + md->n_bottoms);
	double mean;
	if (md->n_scenes < 2) return;
	for (w = 0; w < 3; w++) {
		mean = 0.0;
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) mean += params[offset + N_WATER*k_scene + w];
		mean /= (double) md->n_scenes;
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			r[w*md->n_scenes + k_scene] = unit*(params[offset + N_WATER*k_scene + w] - mean)/TIE_SIGMA;
		}
	}
#else
	(void) params; (void) md; (void) unit; (void) r;
#endif
}

// Continuity residuals (continuity_count values) for the depths in params.

static void continuity_residuals(double *params, model_data *md, double unit, double *r) {
#if CONTINUITY
	int side = window_side(md), row, col, n = 0;
	double H_mean = 0.0, sigma;
	if (! side) return;
	for (row = 0; row < md->n_regions; row++) H_mean += samodel_depth(params, row);
	H_mean /= (double) md->n_regions;
	sigma = CONT_SIGMA_0 + CONT_SIGMA_H*H_mean;
	for (row = 0; row < side; row++) {
		for (col = 1; col < side - 1; col++) {
			r[n++] = unit*(samodel_depth(params, row*side + col - 1) - 2.0*samodel_depth(params, row*side + col)
				+ samodel_depth(params, row*side + col + 1))/sigma;
			r[n++] = unit*(samodel_depth(params, (col - 1)*side + row) - 2.0*samodel_depth(params, col*side + row)
				+ samodel_depth(params, (col + 1)*side + row))/sigma;
		}
	}
#else
	(void) params; (void) md; (void) unit; (void) r;
#endif
}

// The residuals of one region, written to r (residuals_per_region values).

static void region_residuals(double *params, model_data *md, int k_region, double scale, double *r) {

	int k_scene, k_band, k_bottom, n = 0;
	double P, G, X, H, albedo, B[MAX_BOTTOMS], q[MAX_BOTTOMS];

	samodel_unpack_region(params, md, k_region, &H, &albedo, q);
	for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
		B[k_bottom] = albedo;
	}

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		samodel_unpack_scene(params, md, k_scene, &P, &G, &X);
		for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
			samodel_Rrs(P, G, X, B, q, H, 0.0, md, k_scene, k_band, k_region);
			r[n++] = scale*(md->Rrs_modelled[k_region][k_scene][k_band] - md->Rrs_measured[k_region][k_scene][k_band]);
		}
	}
}

int samodel_residuals(double *params, model_data *md, double *r) {

	int k_region, n_region = residuals_per_region(md);
	double scale = residual_scale(md);

	for (k_region = 0; k_region < md->n_regions; k_region++) {
		region_residuals(params, md, k_region, scale, &r[k_region*n_region]);
	}
	continuity_residuals(params, md, chi_unit(md), &r[md->n_regions*n_region]);
	tie_residuals(params, md, chi_unit(md), &r[md->n_regions*n_region + continuity_count(md)]);

	return md->n_regions*n_region + continuity_count(md) + tie_count(md);
}

// Bounded Levenberg-Marquardt on samodel_residuals, over the parameters with is_free[i] != 0
// (the others are held at their values in x). The Jacobian is by forward differences, taken
// inwards at the bounds. Steps are projected onto the box, and a parameter at a bound whose
// gradient points out of the box is held there for that iteration. Damping is Marquardt's
// (scaled by the diagonal of J'J). On return x holds the solution; returns the Rrs error
// (the root sum of squares of the residuals) and the number of residual evaluations. Stops
// after max_iterations, or when the relative fall in the cost is below ftol.

static int cholesky_solve(double *A, double *b, int n) {

	// Solves A y = b in place (y in b) for symmetric positive definite A; A is overwritten.

	int i, j, k;
	double s;

	for (j = 0; j < n; j++) {
		s = A[j*n + j];
		for (k = 0; k < j; k++) s -= A[j*n + k]*A[j*n + k];
		if (s <= 0.0) return 0;
		A[j*n + j] = sqrt(s);
		for (i = j + 1; i < n; i++) {
			s = A[i*n + j];
			for (k = 0; k < j; k++) s -= A[i*n + k]*A[j*n + k];
			A[i*n + j] = s/A[j*n + j];
		}
	}
	for (i = 0; i < n; i++) {
		s = b[i];
		for (k = 0; k < i; k++) s -= A[i*n + k]*b[k];
		b[i] = s/A[i*n + i];
	}
	for (i = n - 1; i >= 0; i--) {
		s = b[i];
		for (k = i + 1; k < n; k++) s -= A[k*n + i]*b[k];
		b[i] = s/A[i*n + i];
	}
	return 1;
}

double samodel_lm(model_data *md, double *x, double *is_free, int max_iterations, double ftol,
	int *n_evaluations) {

	int i, j, k, m, n, n_free, iteration, accepted, pass, solved, new_bounds, *index, *active, *bound,
		n_region, offset, k_region;
	double *lo, *hi, *r, *r_new, *J, *JTJ, *A, *g, *delta, *x_new, *diag,
		cost, cost_new, lambda, nu, h, step_norm, x_norm, predicted, gain, *fixed_step, scale, unit;

	n = samodel_n_params(md);
	n_region = residuals_per_region(md);
	m = md->n_regions*n_region + continuity_count(md) + tie_count(md);
	scale = residual_scale(md);
	unit = chi_unit(md);
	offset = md->n_regions*(1 + md->n_bottoms);

	lo     = malloc(n*sizeof(double));
	hi     = malloc(n*sizeof(double));
	index  = malloc(n*sizeof(int));
	active = malloc(n*sizeof(int));
	r      = malloc(m*sizeof(double));
	r_new  = malloc(m*sizeof(double));
	J      = malloc(m*n*sizeof(double));
	JTJ    = malloc(n*n*sizeof(double));
	A      = malloc(n*n*sizeof(double));
	g      = malloc(n*sizeof(double));
	delta  = malloc(n*sizeof(double));
	diag   = malloc(n*sizeof(double));
	x_new  = malloc(n*sizeof(double));
	bound  = malloc(n*sizeof(int));
	fixed_step = malloc(n*sizeof(double));

	samodel_bounds(md, lo, hi);

	n_free = 0;
	for (i = 0; i < n; i++) {
		x[i] = clamp(x[i], lo[i], hi[i]);
		if (is_free[i] != 0.0) index[n_free++] = i;
	}

	samodel_residuals(x, md, r);
	*n_evaluations = 1;
	cost = 0.0;
	for (k = 0; k < m; k++) cost += r[k]*r[k];

	lambda = 1.0e-3;
	nu = 2.0;

	for (iteration = 0; iteration < max_iterations && n_free > 0; iteration++) {

		// Jacobian of the free parameters (column j of J is parameter index[j]). A region's
		// depth and bottom albedos affect only that region's residuals, so only that region is
		// re-evaluated; the water parameters affect every region.

		copy_double_vec(x, x_new, n);
		for (j = 0; j < n_free; j++) {
			i = index[j];
			h = 1.0e-6*MAX(fabs(x[i]), 1.0);
			if (x[i] + h > hi[i]) h = -h;
			x_new[i] = x[i] + h;
			if (i < offset) {
				k_region = i < md->n_regions ? i : (i - md->n_regions)/md->n_bottoms;
				for (k = 0; k < m; k++) J[k*n_free + j] = 0.0;
				region_residuals(x_new, md, k_region, scale, &r_new[k_region*n_region]);
				for (k = k_region*n_region; k < (k_region + 1)*n_region; k++) {
					J[k*n_free + j] = (r_new[k] - r[k])/h;
				}
				if (i < md->n_regions && continuity_count(md) > 0) {
					continuity_residuals(x_new, md, unit, &r_new[md->n_regions*n_region]);
					for (k = md->n_regions*n_region; k < md->n_regions*n_region + continuity_count(md); k++) {
						J[k*n_free + j] = (r_new[k] - r[k])/h;
					}
				}
			} else {
				samodel_residuals(x_new, md, r_new);
				for (k = 0; k < m; k++) J[k*n_free + j] = (r_new[k] - r[k])/h;
			}
			(*n_evaluations)++;
			x_new[i] = x[i];
		}

		// Gradient J'r and J'J.

		for (j = 0; j < n_free; j++) {
			g[j] = 0.0;
			for (k = 0; k < m; k++) g[j] += J[k*n_free + j]*r[k];
			for (i = 0; i <= j; i++) {
				JTJ[j*n_free + i] = 0.0;
				for (k = 0; k < m; k++) JTJ[j*n_free + i] += J[k*n_free + j]*J[k*n_free + i];
				JTJ[i*n_free + j] = JTJ[j*n_free + i];
			}
			diag[j] = MAX(JTJ[j*n_free + j], 1.0e-12);
		}

		// Hold parameters at a bound whose descent direction (-g) points out of the box.

		for (j = 0; j < n_free; j++) {
			i = index[j];
			active[j] = (x[i] <= lo[i] && g[j] > 0.0) || (x[i] >= hi[i] && g[j] < 0.0);
		}

		// Try damped steps until the cost falls.

		accepted = 0;
		while (lambda < 1.0e12) {

			// Solve the damped system for the free parameters. Any parameter the step would take
			// out of the box is put on its bound and held there, and the rest are solved again
			// (bound[j]: 0 free, -1 at its lower bound, +1 at its upper bound, 2 held by the
			// gradient test above).

			for (j = 0; j < n_free; j++) {
				bound[j] = active[j] ? 2 : 0;
			}

			solved = 1;
			for (pass = 0; pass <= n_free; pass++) {

				for (j = 0; j < n_free; j++) {
					i = index[j];
					fixed_step[j] = bound[j] == -1 ? lo[i] - x[i] : (bound[j] == 1 ? hi[i] - x[i] : 0.0);
				}

				for (j = 0; j < n_free; j++) {
					for (i = 0; i < n_free; i++) {
						A[j*n_free + i] = (bound[j] || bound[i]) ? 0.0 : JTJ[j*n_free + i];
					}
					if (bound[j]) {
						A[j*n_free + j] = 1.0;
						delta[j] = 0.0;
					} else {
						A[j*n_free + j] = JTJ[j*n_free + j] + lambda*diag[j];
						delta[j] = -g[j];
						for (i = 0; i < n_free; i++) {
							if (bound[i]) delta[j] -= JTJ[j*n_free + i]*fixed_step[i];
						}
					}
				}

				if (! cholesky_solve(A, delta, n_free)) {
					solved = 0;
					break ;
				}

				for (j = 0; j < n_free; j++) {
					if (bound[j]) delta[j] = fixed_step[j];
				}

				new_bounds = 0;
				for (j = 0; j < n_free; j++) {
					i = index[j];
					if (! bound[j] && x[i] + delta[j] < lo[i]) {
						bound[j] = -1;
						new_bounds++;
					} else if (! bound[j] && x[i] + delta[j] > hi[i]) {
						bound[j] = 1;
						new_bounds++;
					}
				}
				if (new_bounds == 0) break ;
			}

			if (! solved) {
				lambda *= nu;
				nu *= 2.0;
				continue ;
			}

			copy_double_vec(x, x_new, n);
			step_norm = 0.0;
			x_norm = 0.0;
			for (j = 0; j < n_free; j++) {
				i = index[j];
				x_new[i] = clamp(x[i] + delta[j], lo[i], hi[i]);
				delta[j] = x_new[i] - x[i]; // The projected step.
				step_norm += delta[j]*delta[j];
				x_norm += x[i]*x[i];
			}

			// Reduction predicted by the linear model for the projected step:
			// -(2 g'd + d'J'Jd).

			predicted = 0.0;
			for (j = 0; j < n_free; j++) {
				h = 0.0;
				for (i = 0; i < n_free; i++) h += JTJ[j*n_free + i]*delta[i];
				predicted -= 2.0*g[j]*delta[j] + delta[j]*h;
			}

			samodel_residuals(x_new, md, r_new);
			(*n_evaluations)++;
			cost_new = 0.0;
			for (k = 0; k < m; k++) cost_new += r_new[k]*r_new[k];

			if (cost_new < cost) {
				accepted = 1;
				break ;
			}
			lambda *= nu;
			nu *= 2.0;
			if (sqrt(step_norm) < 1.0e-10*(sqrt(x_norm) + 1.0e-10)) break ;
		}

		if (! accepted) break ;

		// Nielsen's update of the damping from the gain ratio.

		gain = predicted > 0.0 ? (cost - cost_new)/predicted : 0.0;
		lambda *= MAX(1.0/3.0, 1.0 - pow(2.0*gain - 1.0, 3.0));
		lambda = MAX(lambda, 1.0e-12);
		nu = 2.0;

		copy_double_vec(x_new, x, n);
		copy_double_vec(r_new, r, m);

		// Converged: the cost is negligible, has stopped falling (relative change below ftol)
		// or the step is negligible.

		if (cost_new < 1.0e-12 || cost - cost_new < ftol*cost ||
				sqrt(step_norm) < 1.0e-8*(sqrt(x_norm) + 1.0e-8)) {
			cost = cost_new;
			break ;
		}
		cost = cost_new;
	}

	free(lo);
	free(hi);
	free(index);
	free(active);
	free(r);
	free(r_new);
	free(J);
	free(JTJ);
	free(A);
	free(g);
	free(delta);
	free(diag);
	free(x_new);
	free(bound);
	free(fixed_step);

	return sqrt(cost);
}



void samodel_optimise(model_data *md) {

	int i, j, k_region, k_scene, k_band, k_bottom, max_n_params, max_bottom_types, best_single_bottom, offset,
	best_double_bottom_A, best_double_bottom_B, best_single_iterations, best_double_iterations; 
	double max_bottom_retrieval_depth, error, lowest_single_error, lowest_double_error, 
	origin_weight, H, B_type, extinction_depth, nobs, env_system_noise,
	*min, *best, *best_single, *best_double, largest_percentage, albedo, q[MAX_BOTTOMS];
	bool best_single_convergence, best_double_convergence;

#if DEBUG
	printf("\nENTER samodel_optimise\n");
#endif

	max_bottom_retrieval_depth = 8.0; // TODO -- should be set at runtime! 

	// Estimate Rrs(440), Rrs(550), Rrs(640) for each scene. 

	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			md->Rrs440[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], md->n_bands[k_scene], 440.0);
			md->Rrs490[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], md->n_bands[k_scene], 490.0);
			md->Rrs550[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], md->n_bands[k_scene], 550.0);
			md->Rrs640[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], md->n_bands[k_scene], 640.0);
			md->Rrs750[k_region][k_scene] = interp_1d(md->wavelengths[k_scene], md->Rrs_measured[k_region][k_scene], md->n_bands[k_scene], 750.0);

			// Prevent negative Rrs from extrapolation with Sentinel-2 scenes. 
		
			if (md->Rrs440[k_region][k_scene] < 0.0) {
				md->Rrs440[k_region][k_scene] = 0.0001; 
			}
		}
	}

//	printf("\nRrs(440,490,550,640,750) = %.5f, %.5f, %.5f, %.5f, %.5f\n", 
//		md->Rrs440, md->Rrs490, md->Rrs550, md->Rrs640, md->Rrs750);

	// Max number of model parameters. 

	max_bottom_types = md->n_bottoms;

	max_n_params  = md->n_regions; // H
	max_n_params += md->n_regions*md->n_bottoms; // B0, B1, ...
	max_n_params += md->n_regions*md->n_bottoms; // q0, q1, ...
#if DELTA
	max_n_params += 4*md->n_scenes; // P,G,X,D
#else
	max_n_params += N_WATER*md->n_scenes; // P,G,X (and Y)
#endif

	min  = malloc(max_n_params*sizeof(double));
	best = malloc(max_n_params*sizeof(double));
	best_single = malloc(max_n_params*sizeof(double));
	best_double = malloc(max_n_params*sizeof(double));

	// All bottom types simultaneously. 

#if ALL_BOTTOMS

	if (md->h_empirical > max_bottom_retrieval_depth) {
		
		// Too deep to retrieve bottom type - sand only. 

		md->n_bottoms = 1;
		md->bottom_type_indexes[0] = 0;

		md->n_params = samodel_n_params(md);

	} else {

		for (i = 0; i < md->n_bottoms; i++) {
			md->bottom_type_indexes[i] = i;
		}

		md->n_params = samodel_n_params(md);
	}

	samodel_optimise_one_bottom_combination(md, best);
	// printf("\n... error is %.4f", md->model_error);

#else

	// Iterate over all individual bottom types and all pairs.  

	if (md->h_empirical > max_bottom_retrieval_depth || max_bottom_types == 1) {
		
		// Too deep to retrieve bottom type - sand only. 

		md->n_bottoms = 1;
		md->bottom_type_indexes[0] = 0;

		md->n_params = samodel_n_params(md);

		samodel_optimise_one_bottom_combination(md, best);
	} else {

		//	Individual bottom types.

		lowest_single_error = 1.0e4;
		for (i = 0; i < max_bottom_types; i++) {
			// printf("\nTrying %s...", md->B_type_name[i]);

			md->n_bottoms = 1;
			md->bottom_type_indexes[0] = i;

			md->n_params = samodel_n_params(md);

			samodel_optimise_one_bottom_combination(md, min);
			// printf("\n... error is %.4f", md->model_error);

			if (md->model_error < lowest_single_error) {
				lowest_single_error = md->model_error;
				best_single_bottom = i;
				copy_double_vec(min, best_single, md->n_params);

				best_single_iterations = md->n_iterations;
				best_single_convergence = md->converged;
			}
		}

		// All pairs of bottom types. 

		lowest_double_error = 1.0e4;
		for (i = 0; i < max_bottom_types - 1; i++) {
			for (j = i + 1; j < max_bottom_types; j++) {
				// printf("\nTrying %s and %s...", md->B_type_name[i], md->B_type_name[j]);

				md->n_bottoms = 2;
				md->bottom_type_indexes[0] = i;
				md->bottom_type_indexes[1] = j;

				md->n_params = samodel_n_params(md);

				samodel_optimise_one_bottom_combination(md, min);
				// printf("\n... error is %.4f", md->model_error);

				if (md->model_error < lowest_double_error) {
					lowest_double_error = md->model_error;
					best_double_bottom_A = i;
					best_double_bottom_B = j;
					copy_double_vec(min, best_double, md->n_params);

					best_double_iterations = md->n_iterations;
					best_double_convergence = md->converged;
				}
			}
		}

		// Recompute spectrum at optimum. 

		if (lowest_single_error < lowest_double_error) {
			md->n_bottoms = 1;
			md->bottom_type_indexes[0] = best_single_bottom;

			md->n_params = samodel_n_params(md);
			copy_double_vec(best_single, best, md->n_params);

			md->n_iterations = best_single_iterations;
			md->converged = best_single_convergence;
		} else {
			md->n_bottoms = 2;
			md->bottom_type_indexes[0] = best_double_bottom_A;
			md->bottom_type_indexes[1] = best_double_bottom_B;

			md->n_params = samodel_n_params(md);
			copy_double_vec(best_double, best, md->n_params);	

			md->n_iterations = best_double_iterations;
			md->converged = best_double_convergence;					
		}

		error = samodel_error(best, md);
		// printf("\n\nerror = %f\n", error);
	}
#endif

	// Compute min(K) across wavelengths and scenes. 

	md->K_min = array_min_double2(md->K, md->n_scenes, md->max_n_bands, 0.0);
	extinction_depth = 1.7/md->K_min;

	// Estimate depth from model.  

#if 1
	md->depth = 0.0;
	origin_weight = sqrt((double) md->n_regions);

	for (i = 0; i < md->n_regions; i++) {
		H = samodel_depth(best, i);
#if 0
		if (H > extinction_depth) {  
			H = extinction_depth;
		}
#endif
		if (i == md->origin) {
			md->depth += origin_weight*H;
		} else {
			md->depth += H;
		}
	}

	md->depth /= origin_weight + ((double) md->n_regions) - 1.0;

	// No depth where the bottom is not detected (optically deep).

	if (md->optically_deep) {
		md->depth = 0.0;
	}
#else
	md->depth = samodel_depth(best, md->origin);
#endif

	// Store percentages of each bottom type. 

	for (k_bottom = 0; k_bottom < 12; k_bottom++) {
		md->B_type_percent[k_bottom] = 0.0;
	}

	samodel_unpack_region(best, md, md->origin, &H, &albedo, q);

	for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
		md->B_type_percent[md->bottom_type_indexes[k_bottom]] = 100.0*q[k_bottom];
	}

	// Compute bottom type. 

	largest_percentage = 0.0;

	for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
		if (md->B_type_percent[md->bottom_type_indexes[k_bottom]] > largest_percentage) {
			largest_percentage = md->B_type_percent[md->bottom_type_indexes[k_bottom]];
			md->bottom_type = 1.0 + md->bottom_type_indexes[k_bottom];
		} 
	}

	// Compute index of optical depth.  

	md->index_optical_depth = 0.0;
	nobs = 0.0;

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
				env_system_noise = md->rrs_noise[k_scene][k_band];
				md->index_optical_depth += md->rrs_bottom[k_region][k_scene][k_band]/md->rrs_modelled[k_region][k_scene][k_band];
				nobs += 1.0;
			}
		}
	}

	md->index_optical_depth = 100.0*md->index_optical_depth/nobs; 

	// Store P, G, X, and the transformed values in prev (start_at_previous).

	offset = md->n_regions*(1 + md->n_bottoms);
	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		samodel_unpack_scene(best, md, k_scene, &md->P[k_scene], &md->G[k_scene], &md->X[k_scene]);
		for (i = 0; i < N_WATER; i++) {
			md->prev[N_WATER*k_scene + i] = best[offset + N_WATER*k_scene + i];
		}
	}

	// Display all model outputs. 

	// samodel_display(md, best);

	// Restore the total number of bottom types. 

	md->n_bottoms = max_bottom_types;

	free(min);
	free(best);
	free(best_single);
	free(best_double);

#if DEBUG
	printf("\nEXIT samodel_optimise\n");
#endif
}



// Search for the best parameters for one combination of bottom types, by bounded
// Levenberg-Marquardt (samodel_lm) on the Rrs misfit.
//
// Without an empirical depth, the search is in two stages:
//   1. Profile: at each depth on a coarse grid, hold every regional depth there and fit the
//      bottom and water parameters.
//   2. Refine: from the N_REFINE best depths, fit all parameters.
// With an empirical depth, only stage 2 is run, from that depth. After a good fit at the
// previous pixel, the water parameters start from its optimum.

#define N_REFINE 2

// Levenberg-Marquardt limits: the profile only ranks depths, so it stops earlier.

#define LM_PROFILE_ITERATIONS 30
#define LM_PROFILE_FTOL 1.0e-4
#define LM_ITERATIONS 100
#define LM_FTOL 1.0e-8

void samodel_optimise_one_bottom_combination(model_data *md, double *params) { 

	int i, k, k_h, n_h, offset, n_eval, total_count;
	double error, lowest_error, *start, *is_free, *best, *profile, profile_error[N_REFINE];

	double H_grid[] = {0.5, 1.0, 2.0, 3.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0};
	int n_grid = sizeof(H_grid)/sizeof(H_grid[0]);
	double H_refine[N_REFINE];

#if DEBUG
	printf("\nENTER samodel_optimise_one_bottom_combination\n");
#endif

	start   = malloc(md->n_params*sizeof(double));
	is_free = malloc(md->n_params*sizeof(double));
	best    = malloc(md->n_params*sizeof(double));
	profile = malloc(N_REFINE*md->n_params*sizeof(double));
	offset  = md->n_regions*(1 + md->n_bottoms);
	total_count = 0;

	// Stage 1: profile over depth (keeps the N_REFINE lowest errors, in order).

	for (k = 0; k < N_REFINE; k++) {
		profile_error[k] = 1.0e30;
	}

	if (md->empirical_depth_present) {
		n_h = 1;
		H_refine[0] = md->h_empirical;
		samodel_initial(md, md->h_empirical, ALBEDO_START, profile, is_free);
		if (md->start_at_previous) {
			for (i = 0; i < N_WATER*md->n_scenes; i++) {
				profile[offset + i] = md->prev[i];
			}
		}
		profile_error[0] = 0.0;
	} else {
		for (k_h = 0; k_h < n_grid; k_h++) {
			samodel_initial(md, H_grid[k_h], ALBEDO_START, start, is_free);
			if (md->start_at_previous) {
				for (i = 0; i < N_WATER*md->n_scenes; i++) {
					start[offset + i] = md->prev[i];
				}
			}
			for (i = 0; i < md->n_regions; i++) {
				is_free[i] = 0.0;
			}
			error = samodel_lm(md, start, is_free, LM_PROFILE_ITERATIONS, LM_PROFILE_FTOL, &n_eval);
			total_count += n_eval;

			// Insert into the sorted list of the best depths.

			for (k = 0; k < N_REFINE; k++) {
				if (error < profile_error[k]) {
					for (i = N_REFINE - 1; i > k; i--) {
						profile_error[i] = profile_error[i - 1];
						H_refine[i] = H_refine[i - 1];
						copy_double_vec(&profile[(i - 1)*md->n_params], &profile[i*md->n_params], md->n_params);
					}
					profile_error[k] = error;
					H_refine[k] = H_grid[k_h];
					copy_double_vec(start, &profile[k*md->n_params], md->n_params);
					break ;
				}
			}
		}
		n_h = MIN(N_REFINE, n_grid);
	}

	// Stage 2: refine all parameters from the best depths.

	lowest_error = 1.0e30;
	copy_double_vec(profile, best, md->n_params);

	for (k = 0; k < n_h; k++) {
		if (profile_error[k] > 1.0e29) continue ;
		copy_double_vec(&profile[k*md->n_params], start, md->n_params);
		for (i = 0; i < md->n_params; i++) {
			is_free[i] = 1.0;
		}
		error = samodel_lm(md, start, is_free, LM_ITERATIONS, LM_FTOL, &n_eval);
		total_count += n_eval;
		if (error < lowest_error) {
			lowest_error = error;
			copy_double_vec(start, best, md->n_params);
		}
	}

	// Optically deep test: fit the water with every region at H_MAX and no bottom.

	md->optically_deep = false;
	md->deep_chi2 = -1.0;

#if OPTICALLY_DEEP
	{
		int k_scene, k_band, n_resid = 0;
		double deep_error, noise, mean = 0.0, sigma_rel;

		copy_double_vec(best, start, md->n_params);
		for (i = 0; i < offset; i++) {
			start[i] = i < md->n_regions ? H_MAX : 0.0;
			is_free[i] = 0.0;
		}
		for (i = offset; i < md->n_params; i++) {
			is_free[i] = 1.0;
		}
		deep_error = samodel_lm(md, start, is_free, LM_ITERATIONS, LM_FTOL, &n_eval);
		total_count += n_eval;

		noise = 0.0;
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
				noise += md->rrs_noise[k_scene][k_band]*md->rrs_noise[k_scene][k_band];
				for (k = 0; k < md->n_regions; k++) {
					mean += md->Rrs_measured[k][k_scene][k_band];
				}
				n_resid++;
			}
		}
		mean /= (double) (n_resid*md->n_regions);
		noise = sqrt(noise/n_resid);
		n_resid *= md->n_regions;

		sigma_rel = sqrt(pow(noise/mean, 2.0) + DEEP_MODEL_ERROR*DEEP_MODEL_ERROR);
		md->deep_chi2 = n_resid*(deep_error*deep_error - lowest_error*lowest_error)/pow(100.0*sigma_rel, 2.0);
		md->optically_deep = md->deep_chi2 < DEEP_PENALTY*offset*log((double) n_resid);
	}
#endif

	md->n_iterations = total_count;
	md->converged = true;

	// Recompute all model-derived data at the optimum. 

	error = samodel_error(best, md);

	for (i = 0; i < md->n_params; i++) {
		params[i] = best[i];
	}

	free(start);
	free(is_free);
	free(best);
	free(profile);

#if DEBUG
	printf("\nEXIT samodel_optimise_one_bottom_combination\n");
#endif
}




double samodel_error(double *params, model_data *md) {

	int i, k_region, k_scene, k_band, k_bottom, n_total_bands;
	double P, G, X, H, B[12], q[12], Delta, error_Rrs, total_Rrs, error_depth, 
		Rrs_modelled, Rrs_measured, Rrs_measured_mean,
		depth_threshold, depth_mean, n_outliers, error_model, error_bottom, 
		bottom, bottom_mean, bottom_total, bottom_threshold, 
		Rrs_modelled_sum_sq, Rrs_measured_sum_sq, sam_error, 
		error_spectrum, K_min, K_error, albedo, H_region, q_region[12];

#if WRITE_EPGXBHD
	char filename[64];
	FILE *fp;
#endif

#if DEBUG
	printf("\nENTER samodel_error\n");
#endif

	// Compute the modelled Rrs at each wavelength for each scene.

	error_Rrs = 0.0;
	total_Rrs = 0.0;
	n_total_bands = 0;

#if SAM
	Rrs_modelled_sum_sq = 0.0;
	Rrs_measured_sum_sq = 0.0;
	sam_error = 0.0;
#endif

	for (k_region = 0; k_region < md->n_regions; k_region++) {

		samodel_unpack_region(params, md, k_region, &H, &albedo, q);
		for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
			B[k_bottom] = albedo;
		}

		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {

			samodel_unpack_scene(params, md, k_scene, &P, &G, &X);
			Delta = 0.0;

#if WRITE_EPGXBHD
			// P
			sprintf(filename, "P_i%d_j%d_r%d_s%d.csv", md->i, md->j, k_region, k_scene);
			fp = fopen(filename, "a");
			fprintf(fp, "%.6f,", P);
			fclose(fp);

			// G
			sprintf(filename, "G_i%d_j%d_r%d_s%d.csv", md->i, md->j, k_region, k_scene);
			fp = fopen(filename, "a");
			fprintf(fp, "%.6f,", G);
			fclose(fp);

			// X
			sprintf(filename, "X_i%d_j%d_r%d_s%d.csv", md->i, md->j, k_region, k_scene);
			fp = fopen(filename, "a");
			fprintf(fp, "%.6f,", X);
			fclose(fp);

			for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
				// B
				sprintf(filename, "B_%d_i%d_j%d_r%d_s%d.csv", k_bottom, md->i, md->j, k_region, k_scene);
				fp = fopen(filename, "a");
				fprintf(fp, "%.6f,", B[k_bottom]);
				fclose(fp);

				// q
				sprintf(filename, "q_%d_i%d_j%d_r%d_s%d.csv", k_bottom, md->i, md->j, k_region, k_scene);
				fp = fopen(filename, "a");
				fprintf(fp, "%.6f,", q[k_bottom]);
				fclose(fp);
			}

			// H
			sprintf(filename, "H_i%d_j%d_r%d_s%d.csv", md->i, md->j, k_region, k_scene);
			fp = fopen(filename, "a");
			fprintf(fp, "%.6f,", H);
			fclose(fp);

#if DELTA
			// Delta
			sprintf(filename, "D_i%d_j%d_r%d_s%d.csv", md->i, md->j, k_region, k_scene);
			fp = fopen(filename, "a");
			fprintf(fp, "%.8f,", Delta);
			fclose(fp);
#endif
#endif
			for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {

				if (true || !(md->wavelengths[k_scene][k_band] > 675.0 && md->wavelengths[k_scene][k_band] < 750.0)) {

				  samodel_Rrs(P, G, X, B, q, H, Delta, md, k_scene, k_band, k_region);

					Rrs_modelled = md->Rrs_modelled[k_region][k_scene][k_band];
					Rrs_measured = md->Rrs_measured[k_region][k_scene][k_band];

					// RMS relative error.

					error_Rrs += pow(Rrs_modelled - Rrs_measured, 2.0);
					total_Rrs += Rrs_measured;

					// error_Rrs += pow((Rrs_modelled - Rrs_measured)/Rrs_measured, 2.0);
					n_total_bands++;

					// SAM error.
#if SAM
					Rrs_modelled_sum_sq += pow(Rrs_modelled, 2.0);
					Rrs_measured_sum_sq += pow(Rrs_measured, 2.0);
					sam_error += Rrs_modelled*Rrs_measured;
#endif
				} 
			}
		}
	}

	// Compute RMS relative percentage error of Rrs modelled vs Rrs measured.

	Rrs_measured_mean = total_Rrs/((double) n_total_bands);
	error_Rrs = 100.0*sqrt(error_Rrs/((double) n_total_bands))/Rrs_measured_mean;

	// error_Rrs = 100.0*sqrt(error_Rrs/((double) n_total_bands));
	md->Rrs_error = error_Rrs;

	// Compute the SAM (Spectral angle mapper) error. 

#if SAM
	sam_error = acos(sam_error/(sqrt(Rrs_modelled_sum_sq)*sqrt(Rrs_measured_sum_sq)));
	sam_error /= PI;
	md->SAM_error = sam_error;
#else
	sam_error = 1.0;
	md->SAM_error = sam_error;
#endif

	// Model spectrum error. 

	error_spectrum = error_Rrs*sam_error;

	// Compute the RMS relative percentage depth error (across the region). 

	// depth_threshold = 0.15 // 15%
	error_depth = 0.0;
	depth_mean = 0.0;
	n_outliers = 0.0;

	for (k_region = 0; k_region < md->n_regions; k_region++) {
		depth_mean += samodel_depth(params, k_region);
	}
	depth_mean /= (double) md->n_regions;

	if (depth_mean < 4.0) {
		depth_threshold = 0.4;
	} else if (depth_mean < 8.0) {
		depth_threshold = 0.2;
	} else if (depth_mean < 12.0) {
		depth_threshold = 0.1;
	} else {
		depth_threshold = 0.05;
	}

	for (k_region = 0; k_region < md->n_regions; k_region++) {
		H_region = samodel_depth(params, k_region);
		if (H_region < (1.0 - depth_threshold)*depth_mean || H_region > (1.0 + depth_threshold)*depth_mean) {
			error_depth += pow(H_region - depth_mean, 2.0);
			n_outliers += 1.0;
		}
	}

	if (n_outliers > 0.5) {
		error_depth = 100.0*sqrt(error_depth/n_outliers)/depth_mean;
	}

	md->depth_error = error_depth; 

	// Compute the RMS relative percentage bottom difference (across the region). 

	if (depth_mean < 5.0) {
		bottom_threshold = 0.25;
	} else if (depth_mean < 10.0) {
		bottom_threshold = 0.1;
	} else if (depth_mean < 15.0) {
		bottom_threshold = 0.05;
	} else {
		bottom_threshold = 0.01;
	}

	error_bottom = 0.0;
	n_outliers = 0.0;
	bottom_total = 0.0;

	for (k_bottom = 0; k_bottom < md->n_bottoms; k_bottom++) {
		
		// Compute mean q*B. 

		bottom_mean = 0.0;
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			samodel_unpack_region(params, md, k_region, &H_region, &albedo, q_region);
			bottom_mean += albedo*q_region[k_bottom];
		}

		bottom_mean /= (double) md->n_regions;
		bottom_total += bottom_mean;

		// Compute deviation from mean. 

		for (k_region = 0; k_region < md->n_regions; k_region++) {

			samodel_unpack_region(params, md, k_region, &H_region, &albedo, q_region);
			bottom = albedo*q_region[k_bottom];

			if (bottom < (1.0 - bottom_threshold)*bottom_mean || bottom > (1.0 + bottom_threshold)*bottom_mean) {
				error_bottom += pow(bottom - bottom_mean, 2.0);
				n_outliers += 1.0;
			}
		}
	}

	if (n_outliers > 0.5) {
		bottom_mean = bottom_total/((double) md->n_bottoms);
		error_bottom = 100.0*sqrt(error_bottom/n_outliers)/bottom_mean;
	}

	md->bottom_error = error_bottom;

	K_error = 0.0;

#if K_PRIOR
	// Move away from noisey K in shallow water. That is, very low K and shallow H. 

	double min_mean_coastal_K = 0.275;
	double min_min_coastal_K = 0.185;
	H = samodel_depth(params, md->origin);
	K_error = 0.0;

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {

		// Compute min(K) in the scene. 
		K_min = 1.0e4;
		for (k_band = 0; k_band < md->max_n_bands; k_band++) {
			if (! approx_equal(md->K[k_scene][k_band], 0.0, 1.0e-6) && md->K[k_scene][k_band] < K_min) {
				K_min = md->K[k_scene][k_band];
			}
		}

		if (H < 1.0 && K_min < min_mean_coastal_K) {
		  K_error += 100.0*pow(1.0/(0.01 + K_min) - 1.0/(0.01 + min_min_coastal_K), 2); 
		} else if (H < 2.0 && K_min < 0.5*(1.5*min_min_coastal_K + 0.5*min_mean_coastal_K)) {
		  K_error += 100.0*pow(1.0/(0.01 + K_min) - 1.0/(0.01 + 0.5*(1.5*min_min_coastal_K + 0.5*min_mean_coastal_K)), 2); 
		} else if (H < 3.0 && K_min < 0.5*(1.25*min_min_coastal_K + 0.75*min_mean_coastal_K)) {
		  K_error += 100.0*pow(1.0/(0.01 + K_min) - 1.0/(0.01 + 0.5*(1.25*min_min_coastal_K + 0.75*min_mean_coastal_K)), 2);
		} else if (H < 4.0 && K_min < 0.5*(1.5*min_min_coastal_K + 0.5*min_mean_coastal_K)) {
		  K_error += 100.0*pow(1.0/(0.01 + K_min) - 1.0/(0.01 + 0.5*(1.5*min_min_coastal_K + 0.5*min_mean_coastal_K)), 2); 
		} else if (H < 5.0 && K_min < 0.5*(1.75*min_min_coastal_K + 0.25*min_mean_coastal_K)) {
		  K_error += 100.0*pow(1.0/(0.01 + K_min) - 1.0/(0.01 + 0.5*(1.75*min_min_coastal_K + 0.25*min_mean_coastal_K)), 2);
		}
	}

	md->K_error = K_error;

	// Move away from extremely high K. 

	if (K_min > 0.7) {
		K_error += 100.0*pow(4.0*(K_min - 0.7), 2);	
	}

	md->K_error = K_error;	
#endif

	md->K_error = K_error;

	// Compute the weighted model error. 

	double weight_spectrum = 80.0;
	double weight_depth    = 15.0;
	double weight_bottom   = 10.0; 
	double weight_K        = 15.0;

	error_model = weight_spectrum*error_spectrum + weight_depth*error_depth + 
		weight_bottom*error_bottom + weight_K*K_error;
	error_model /= weight_spectrum + weight_depth + weight_bottom + weight_K;
	md->model_error = error_model;

#if WRITE_EPGXBHD
	// E
	sprintf(filename, "E_i%d_j%d.csv", md->i, md->j);
	fp = fopen(filename, "a");
	fprintf(fp, "%.4f,", md->model_error);
	fclose(fp);
#endif

#if DEBUG
	printf("\nEXIT samodel_error\n");
#endif

	return error_model;
}


float model_spectral_diff(double ***Rrs_measured, double ***lut_modelled_Rrs, model_data *md) {

	int k_region, k_scene, k_band, n_total_bands;
	double error, Rrs_measured_mean, Rrs_total;

	error = 0.0; 
	Rrs_total = 0.0;
	n_total_bands = 0;

	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {
				if (! approx_equal(Rrs_measured[k_region][k_scene][k_band], 0.0, 1.0e-6) && 
					! approx_equal(lut_modelled_Rrs[k_region][k_scene][k_band], 0.0, 1.0e-6)) {
					error += pow(Rrs_measured[k_region][k_scene][k_band] - lut_modelled_Rrs[k_region][k_scene][k_band], 2.0);
					Rrs_total += Rrs_measured[k_region][k_scene][k_band];
					n_total_bands++;
				}
			}
		}
	}

	Rrs_measured_mean = Rrs_total/((double) n_total_bands);
	error = 100.0*sqrt(error/((double) n_total_bands))/Rrs_measured_mean;

	return error;
}



float model_spectral_diff_mean(double ***Rrs_measured, double ***lut_modelled_Rrs, model_data *md) {

	int k_region, k_scene, k_band, n_total_bands;
	double error, Rrs_measured_mean, Rrs_total, current, previous, nobs;

	error = 0.0; 
	Rrs_total = 0.0;
	n_total_bands = 0;

	for (k_band = 0; k_band < md->max_n_bands; k_band++) {

		if ( ! (k_band == 1 || k_band == 2 || k_band == 3)) {
			continue ;
		}

		current = 0.0;
		previous = 0.0;
		nobs = 0.0;

		// Mean difference across scenes and region. 

		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_region = 0; k_region < md->n_regions; k_region++) {
				if (k_band < md->n_bands[k_scene]) {
					if (! approx_equal(Rrs_measured[k_region][k_scene][k_band], 0.0, 1.0e-6) && 
							! approx_equal(lut_modelled_Rrs[k_region][k_scene][k_band], 0.0, 1.0e-6)) {
						current  += Rrs_measured[k_region][k_scene][k_band];
						previous += lut_modelled_Rrs[k_region][k_scene][k_band];
						nobs += 1.0;
					}
				}
			}
		}

		if (nobs > 0.5) {
			current /= nobs;
			previous /= nobs;
		}

		// error += pow(current - previous, 2.0);
		// Rrs_total += current;
		error += pow((current - previous)/current, 2.0);
		n_total_bands++;
	}

	// Rrs_measured_mean = Rrs_total/((double) n_total_bands);
	// error = 100.0*sqrt(error/((double) n_total_bands))/Rrs_measured_mean;
	error = 100.0*sqrt(error/((double) n_total_bands));

	return error;
}



// Self-consistent backscattering slope (Y_MODEL 1 and 2): Y from the ratio of the model's own
// deep-water reflectance rrs_dp at two wavelengths, which depends on Y, solved by fixed-point
// iteration. The last result is cached, as it is called with the same P, G, X for every band.

static double deep_rrs(int k, double lambda, double P, double G, double X, double Y) {

	// Deep-water rrs at one of two cached wavelengths (k = 0, 1).

	const double S = 0.015;
	static double cached_lambda[2] = {-1.0, -1.0}, a_0[2], a_1[2], a_w[2], b_bw[2];
	double a, bb, u;

	if (cached_lambda[k] != lambda) {
		a_0[k]  = interp_1d(ref_wlens, a0_Lee, n_ref_wlens, lambda);
		a_1[k]  = interp_1d(ref_wlens, a1_Lee, n_ref_wlens, lambda);
		a_w[k]  = interp_1d(ref_wlens, aw_Pope_Fry1997, n_ref_wlens, lambda);
		b_bw[k] = interp_1d(ref_wlens, bbw_Morel_1974, n_ref_wlens, lambda);
		cached_lambda[k] = lambda;
	}

	a  = a_w[k] + (a_0[k] + a_1[k]*log(P))*P + G*exp(-S*(lambda - 440.0));
	bb = b_bw[k] + X*pow(440.0/lambda, Y);
	u  = bb/(a + bb);
	return (LEE_G0 + LEE_G1*u)*u;
}

double samodel_water_Y(double P, double G, double X) {

	static double last_P = -1.0, last_G = -1.0, last_X = -1.0, last_Y = 1.0;
	double Y, Y_new, ratio;
	int k;

	if (P == last_P && G == last_G && X == last_X) {
		return last_Y;
	}

	Y = last_Y; // Start from the previous solution (usually a nearby P, G, X).
	for (k = 0; k < 50; k++) {
#if Y_MODEL == 1
		ratio = deep_rrs(0, 443.0, P, G, X, Y)/deep_rrs(1, 555.0, P, G, X, Y);
		Y_new = 2.0*(1.0 - 1.2*exp(-0.9*ratio));
#else
		ratio = deep_rrs(0, 440.0, P, G, X, Y)/deep_rrs(1, 490.0, P, G, X, Y);
		Y_new = 3.44*(1.0 - 3.17*exp(-2.01*ratio));
#endif
		if (Y_new < Y_MIN) Y_new = Y_MIN;
		if (Y_new > Y_MAX) Y_new = Y_MAX;
		if (fabs(Y_new - Y) < 1.0e-10) {
			Y = Y_new;
			break ;
		}
		Y = Y_new;
	}

	last_P = P;
	last_G = G;
	last_X = X;
	last_Y = Y;
	return Y;
}



void samodel_Rrs(
		 double P, double G, double X, double B[], double q[], double H, double Delta, // parameters
	model_data *md, int k_scene, int k_band, int k_region) {

#if DEBUG
	printf("\nENTER samodel_Rrs\n");
#endif

	const double S = 0.015;
	int i;
	double rho, a_phi, a_g, a, chi, Y, M, b_p, bb, u, K, DuC, DuB, 
	  f, kd, kup_water, kup_bottom, rrs_dp, rrs_C, rrs_B, rrs;

	double a_0     = md->a_0[k_scene][k_band];
	double a_1     = md->a_1[k_scene][k_band];
	double a_w     = md->a_w[k_scene][k_band];
	double b_bw    = md->b_bw[k_scene][k_band];
	double lambda  = md->wavelengths[k_scene][k_band];
	double Rrs440  = md->Rrs440[k_region][k_scene];
	double Rrs490  = md->Rrs490[k_region][k_scene];
	double Rrs750  = md->Rrs750[k_region][k_scene];
	double tide    = md->H_tide[k_scene];
	double sec_theta_view = md->sec_theta_view[k_scene];
	double sec_theta_sun = md->sec_theta_sun[k_scene];
	double H_scene;

	// Compute the water column depth at the time of this scene. H is the depth below
	// the common datum and tide is this scene's tidal offset to that datum.

	H_scene = H + tide;
	if (H_scene < 0.0) H_scene = 0.0; // Exposed at this tide.

	// Compute bottom albedo. 

	md->bottom_albedo = 0.0;
	rho = 0.0;

	for (i = 0; i < md->n_bottoms; i++) {
		md->bottom_albedo += q[i]*B[i]; // q already normalised. 
		rho += q[i]*B[i]*md->bottom_reflectance[md->bottom_type_indexes[i]][k_scene][k_band];
	}

	md->rho[k_region][k_scene][k_band] = rho;

	// Compute a - total absorption coefficient. 

	a_phi = (a_0 + a_1*log(fabs(P)))*fabs(P); // Phytoplankton pigment absorption coefficient.

	a_g = fabs(G)*exp(-S*(lambda - 440.0)); // Gelbstoff absorption coefficient. 

	a = a_w + a_phi + a_g; 

	// Compute bb - total backscattering coefficient.

#if Y_MODEL == 0
	// chi = (Rrs440 - Rrs750)/(Rrs490 - Rrs750);
	chi = Rrs440/Rrs490;
	Y = 3.44*(1.0 - 3.17*exp(-2.01*chi));
	if (Y < 0.0) Y = 0.0;
	if (Y > 2.5) Y = 2.5;  
#elif Y_MODEL == 3
	Y = md->Y_scene[k_scene];
#else
	Y = samodel_water_Y(P, G, X);
#endif
	// Y = 0.5; // Consistent with more turbid waters per Lee, 2001. 
	b_p = X*pow(440.0/lambda, Y);

	bb = b_bw + b_p;

	// Model of Lee et al.

	// Compute u and K. 

	u = bb/(a + bb);
	K = a + bb;
	if (K < 0.0) K = 0.0; 
	if (K > 2.5) K = 2.5;
	md->K[k_scene][k_band] = K;

	// Compute rrs_dp - the remote sensing reflectance for optically deep water.

	rrs_dp = (LEE_G0 + LEE_G1*u)*u;

	md->rrs_dp[k_region][k_scene][k_band] = rrs_dp;

	// Compute the optical path elongation factors for scattered photons.

	DuC = LEE_C0*sqrt(1.0 + LEE_C1*u); // Water column. 
	DuB = LEE_B0*sqrt(1.0 + LEE_B1*u); // Bottom. 

	// Compute rrs. 

	M = sec_theta_sun + DuC*sec_theta_view;
	rrs_C = rrs_dp*(1.0 - exp(-M*K*H_scene)); // Contribution of the water column. 

	M = sec_theta_sun + DuB*sec_theta_view;
	rrs_B = rho/PI*exp(-M*K*H_scene); // Contribution of the bottom. 
#if INTERREFLECTION
	{
		double T = exp(-2.0*D_DIFFUSE*K*H_scene);
		rrs_B /= 1.0 - (R_INTERNAL*T + PI*rrs_dp*(1.0 - T))*rho;
	}
#endif

	md->rrs_bottom[k_region][k_scene][k_band] = rrs_B;
	// printf("\nrrs_bottom = %.7f\n", md->rrs_bottom[k_region][k_scene][k_band]);

	rrs = rrs_C + rrs_B;
	md->rrs_modelled[k_region][k_scene][k_band] = rrs;	

	// Compute Rrs. 

	md->Rrs_modelled[k_region][k_scene][k_band] = 0.5*rrs/(1.0 - 1.5*rrs) + Delta;

#if DEBUG
	printf("\nEXIT samodel_Rrs\n");
#endif
}







void extract_Rrs_data(int i, int j, scene scene_data[], geogrid gridded_data[], int *scene_indexes, 
	int nscenes, int n_spatial, int n_smoothing_radius, int nrows, int ncols, model_data *md, 
	float n_sigma) {

	int k, ii, jj, iii, jjj, k_region, k_scene, k_band; 
	float refl; 
	bool nodata;

	k_region = 0;

	if (n_spatial == 0) {
		n_spatial = 1;
	}

	for (ii = 1 - n_spatial; ii < n_spatial; ii++) {

		if (i + ii < 0) {
			iii = 0;
		} else if (i + ii > nrows - 1) {
			iii = nrows - 1;
		} else { 
			iii = i + ii;
		}

		for (jj = 1 - n_spatial; jj < n_spatial; jj++) {

			if (j + jj < 0) {
				jjj = 0;
			} else if (j + jj > ncols - 1) {
				jjj = ncols - 1;
			} else { 
				jjj = j + jj;
			}

			nodata = false;

			for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
				for (k_band = 0; k_band < md->n_bands[k_scene]; k_band++) {

					k = scene_data[scene_indexes[k_scene]].band_indexes[k_band];
#if LOCALLY_AVERAGE
					refl = smoothed_index_2d(gridded_data[k].array, iii, jjj, nrows, ncols, 
														n_smoothing_radius, gridded_data[k].nodata_value);
#else 
					refl = gridded_data[k].array[iii][jjj];
#endif
					if (! approx_equal(refl, gridded_data[k].nodata_value, 1.0e-6)) {
						md->Rrs_measured[k_region][k_scene][k_band] = refl;
						if (! approx_equal(n_sigma, 0.0, 1.0e-6)) {
							md->Rrs_measured[k_region][k_scene][k_band] += n_sigma*scene_data[scene_indexes[k_scene]].R_sigma[k_band];
						}
					} else {
						nodata = true;
						break ;
					}
				}
				if (nodata) break ;
			}

			if (i == iii && j == jjj) {
				md->origin = k_region;
			}

			if (! nodata) {
				k_region++;
			}
		}
	}

	md->n_regions = k_region;
}



void partial_copy_md(model_data *md, model_data *md_out, bool allocate) {

	// Only copies the essential model data used in the lookup table. 

	int i, j, k;
	int *n_bands;
	double *P;
	double *G;
	double *X;
	double ***Rrs_modelled;
	double **K;
	double *prev;

	md_out->n_bottoms     = md->n_bottoms;
	md_out->n_regions     = md->n_regions;
	md_out->n_scenes      = md->n_scenes; 
	md_out->max_n_regions = md->max_n_regions;
	md_out->max_n_bands   = md->max_n_bands;

	md_out->depth               = md->depth;
	// md_out->depth_prev          = md->depth_prev;
	md_out->depth_sigma         = md->depth_sigma;
	md_out->Rrs_error           = md->Rrs_error;
	md_out->SAM_error           = md->SAM_error;
	md_out->K_error             = md->K_error;
	md_out->depth_error         = md->depth_error;
	md_out->bottom_error        = md->bottom_error;
	md_out->model_error         = md->model_error;
	md_out->K_min               = md->K_min;
	// md_out->K_min_prev          = md->K_min_prev;
	md_out->bottom_albedo       = md->bottom_albedo;
	md_out->index_optical_depth = md->index_optical_depth;
	md_out->bottom_type         = md->bottom_type;
	md_out->converged           = md->converged; 
	md_out->n_iterations        = md->n_iterations;

	for (i = 0; i < 12; i++) {
		md_out->B_type_percent[i] = md->B_type_percent[i];
	}

	if (allocate) {
		n_bands = (int*) malloc(md->n_scenes*sizeof(int));
		md_out->n_bands = n_bands;
	}

	for (i = 0; i < md->n_scenes; i++) {
		md_out->n_bands[i] = md->n_bands[i];
	}

	if (allocate) {
		P = malloc(md->n_scenes*sizeof(double));
		G = malloc(md->n_scenes*sizeof(double));
		X = malloc(md->n_scenes*sizeof(double));	

		md_out->P = P;
		md_out->G = G;
		md_out->X = X;
	}

	for (i = 0; i < md->n_scenes; i++) {
		md_out->P[i] = md->P[i];
		md_out->G[i] = md->G[i];
		md_out->X[i] = md->X[i];
	}

	if (allocate) {
		allocate_double_array_2d(&K, md->n_scenes, md->max_n_bands);
		md_out->K = K;
	}

	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->max_n_bands; j++) {
			md_out->K[i][j] = md->K[i][j];
		}
	}

	if (allocate) {
		allocate_double_array_3d(&Rrs_modelled, md->max_n_regions, md->n_scenes, md->max_n_bands);
		md_out->Rrs_modelled = Rrs_modelled;
	}

	for (i = 0; i < md->max_n_regions; i++) {
		for (j = 0; j < md->n_scenes; j++) {
			for (k = 0; k < md->n_bands[j]; k++) {
				md_out->Rrs_modelled[i][j][k] = md->Rrs_modelled[i][j][k];
			}
		}
	}

#if DELTA
	if (allocate) {
		prev = (double*) malloc(4*md->n_scenes*sizeof(double));
		md_out->prev = prev;
	}

	for (i = 0; i < 4*md->n_scenes; i++) {
		md_out->prev[i] = md->prev[i];
	}
#else
	if (allocate) {
		prev = (double*) malloc(4*md->n_scenes*sizeof(double));
		md_out->prev = prev;
	}

	for (i = 0; i < N_WATER*md->n_scenes; i++) {
		md_out->prev[i] = md->prev[i];
	}
#endif
}



void free_model_data(model_data *md) {

	free(md->P);
	free(md->G);
	free(md->X);
	free(md->n_bands);

	free_double_array_2d(md->K, md->n_scenes);
	free_double_array_3d(md->Rrs_modelled, md->max_n_regions, md->n_scenes);
	free(md);
}



void samodel_display(model_data *md, double *min) {

	int i, j, k, offset;
	double q_sum;

	offset = md->n_regions + (md->n_bottoms + 1)*md->n_regions;

	printf("\nn_scenes  = %d", md->n_scenes);
	printf("\nn_regions = %d", md->n_regions);
	printf("\nn_bottoms = %d\n", md->n_bottoms);

	if (md->start_at_previous) {
		printf("\nhot start = true");
	} else {
		printf("\nhot start = false");
	}

	if (md->converged){
		printf("\nConverged    = true");
	} else {
		printf("\nConverged    = false");
	}
	printf("\nN evals      = %d\n", md->n_iterations);
	printf("\nError        = %5.2f (%%)", md->model_error);
	printf("\nError Rrs    = %5.2f (%%)", md->Rrs_error);
	printf("\nError H      = %5.2f (%%)", md->depth_error);
	printf("\nError B      = %5.2f (%%)", md->bottom_error);
	
	printf("\n\nH            = %5.2f (m)",  md->depth);
	if (md->empirical_depth_present) {
		printf("\nH_start      = %5.2f (m)",  md->h_empirical);
	} else if (md->start_at_previous) {
		printf("\nH_start      = %5.2f (m)",  md->depth_prev);
	} else {
		printf("\nH_start      = %5.2f (m)",  0.0);
	}

	printf("\nIOD          = %5.2f (%%)", md->index_optical_depth);
	printf("\nK_min        = %6.4f (1/m)\n", md->K_min);

	for (i = 0; i < md->n_bottoms; i++) {
		printf("\nB_%s       = %5.2f (%%)", md->B_type_name[md->bottom_type_indexes[i]], md->B_type_percent[md->bottom_type_indexes[i]]);
	}

#if DELTA
	printf("\n\nP = \n\t");
	for (i = 0; i < md->n_scenes; i++) {
		printf("%.8f\t", 0.01*min[offset + 4*i]);
	}

	printf("\nG = \n\t");
	for (i = 0; i < md->n_scenes; i++) {
		printf("%.8f\t", 0.01*min[offset + 1 + 4*i]);
	}

	printf("\nX = \n\t");
	for (i = 0; i < md->n_scenes; i++) {	
		printf("%.8f\t", 0.01*min[offset + 2 + 4*i]);
	}
#else
	printf("\n\nP = \n\t");
	for (i = 0; i < md->n_scenes; i++) {
		printf("%.8f\t", 0.01*min[offset + 3*i]);
	}

	printf("\nG = \n\t");
	for (i = 0; i < md->n_scenes; i++) {
		printf("%.8f\t", 0.01*min[offset + 1 + 3*i]);
	}

	printf("\nX = \n\t");
	for (i = 0; i < md->n_scenes; i++) {	
		printf("%.8f\t", 0.01*min[offset + 2 + 3*i]);
	}
#endif

	printf("\nB = \n");
	for (i = 0; i < md->n_regions; i++) {
	    printf("\t%.4f\n", 0.01*min[md->n_regions + i*(md->n_bottoms + 1)]);
	}

	printf("\nq = \n\t");
	for (i = 0; i < md->n_regions; i++) {

		q_sum = 0.0;
		for (j = 0; j < md->n_bottoms; j++) {
		  q_sum += min[md->n_regions + i*(md->n_bottoms + 1) + j + 1];
		}

		for (j = 0; j < md->n_bottoms; j++) {
		  printf("%.4f\t", min[md->n_regions + i*(md->n_bottoms + 1) + j + 1]/q_sum);
		}
		printf("\n\t");
	}

	printf("\nH = \n\t");
	for (i = 0; i < md->n_regions; i++) {
		printf("%.2f\n\t", min[i]);
	}

#if DELTA
	printf("\nD = \n\t");
	for (i = 0; i < md->n_scenes; i++) {
		printf("%.8f\t", 0.001*min[offset + 3 + 4*i]);
	}
#endif

	printf("\n\nK = \n");
	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->n_bands[i]; j++) {
			printf("%10.0f\t", md->wavelengths[i][j]);
		}
		printf("\n");
		for (j = 0; j < md->n_bands[i]; j++) {
			printf("%10.4f\t", md->K[i][j]);
		}
		printf("\n");
	}

	printf("\n\n1.7/K = \n");
	
	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->n_bands[i]; j++) {
			printf("%10.0f\t", md->wavelengths[i][j]);
		}
		printf("\n");
		for (j = 0; j < md->n_bands[i]; j++) {
			printf("%10.2f\t", 1.7/(md->K[i][j]));
		}
		printf("\n");
	}

	printf("\nRrs_measured vs Rrs_modelled = \n");

	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->n_regions; j++) {
			printf("scene %d, region %d\n", i + 1, j + 1);
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%10.0f\t", md->wavelengths[i][k]);
			}
			printf("\n");
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%.8f\t", md->Rrs_measured[j][i][k]);
			}
			printf("\n");
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%.8f\t", md->Rrs_modelled[j][i][k]);
			}
			printf("\n");
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%10.0f\t", 100.0*fabs(md->Rrs_modelled[j][i][k] - md->Rrs_measured[j][i][k])/md->Rrs_measured[j][i][k]);
			}
			printf("\n\n");
		}
		printf("\n");
	}

	printf("\nRrs_bottom = \n");

	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->n_regions; j++) {
			printf("scene %d, region %d\n", i + 1, j + 1);
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%10.0f\t", md->wavelengths[i][k]);
			}
			printf("\n");
			for (k = 0; k < md->n_bands[i]; k++) {
#if DELTA
				printf("%.8f\t", 0.001*min[offset + 3 + 4*i] + 0.5*md->rrs_bottom[j][i][k]/(1.0 - 1.5*md->rrs_bottom[j][i][k]));
#else 
				printf("%.8f\t", 0.5*md->rrs_bottom[j][i][k]/(1.0 - 1.5*md->rrs_bottom[j][i][k]));
#endif
			}
			printf("\n\n");
		}
	}

	printf("\nrho = \n");

	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->n_regions; j++) {
			printf("scene %d, region %d\n", i + 1, j + 1);
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%10.0f\t", md->wavelengths[i][k]);
			}
			printf("\n");
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%.8f\t", md->rho[j][i][k]);
			}
			printf("\n\n");
		}
	}

	printf("\nRrs_dp = \n");

	for (i = 0; i < md->n_scenes; i++) {
		for (j = 0; j < md->n_regions; j++) {
			printf("scene %d, region %d\n", i + 1, j + 1);
			for (k = 0; k < md->n_bands[i]; k++) {
				printf("%10.0f\t", md->wavelengths[i][k]);
			}
			printf("\n");
			for (k = 0; k < md->n_bands[i]; k++) {
#if DELTA
				printf("%.8f\t", 0.001*min[offset + 3 + 4*i] + 0.5*md->rrs_dp[j][i][k]/(1.0 - 1.5*md->rrs_dp[j][i][k]));
#else
				printf("%.8f\t", 0.5*md->rrs_dp[j][i][k]/(1.0 - 1.5*md->rrs_dp[j][i][k]));
#endif
			}
			printf("\n\n");
		}
	}

}





