#ifndef SPECTRA_H
#define SPECTRA_H

// Spectral data for the semi-analytical model, read from CSV files (model/data/README.md):
// the water tables (a_w, b_bw, a_0, a_1) and the bottom reflectance library.

#define MAX_BOTTOM_TYPES 12

// The water tables.

enum { WATER_AW, WATER_BBW, WATER_A0, WATER_A1, N_WATER_TABLES };

// Directory of the default data files: $PHOTIC_DATA if set, otherwise the build's model/data.

const char *spectra_data_dir(void);

// Read a water file. It may hold any of the columns aw, bbw, a0 and a1; each one read replaces
// the current table. Returns 0 on success.

int spectra_read_water(const char *file);

// Read a bottom library, replacing the current one. Returns 0 on success.

int spectra_read_bottom(const char *file);

// Read a bottom library from a NetCDF file of individual spectra (model/data/README.md,
// bottom_spectra_insitu.nc): one type per class, the median of its spectra's shapes times their
// median reflectance at 550 nm, bounded by the 5th-95th percentiles of that reflectance over the
// class. classes (n_classes names, or 0 for the file's own order) picks and orders the types;
// all_halves uses every spectrum of a class, not just its library half. Returns 0 on success.

int spectra_read_bottom_nc(const char *file, char **classes, int n_classes, int all_halves);

// Load the default files for any water tables not yet read, and the built-in bright sand if no
// bottom library has been read. Returns 0 on success.

int spectra_load_defaults(void);

// Print the tables in use and the files they came from.

void spectra_print(void);

// A water table at a wavelength (nm), by linear interpolation.

double spectra_water(int table, double lambda);

// The bottom library: the number of types, a type's name, its reflectance normalised to 1 at
// 550 nm, and the bounds of its brightness (the reflectance at 550 nm).

int spectra_n_bottoms(void);
const char *spectra_bottom_name(int k);
double spectra_bottom(int k, double lambda);
double spectra_bottom_min(int k);
double spectra_bottom_max(int k);

// Spectral response functions (SCENE ... SRF /file/): a CSV of relative response with one
// column per band. spectra_read_srf returns a handle (the same file read twice gives the same
// handle), or -1. spectra_srf_band gives the band whose centroid is nearest lambda (within
// 20 nm), or -1. The _srf functions average a table over a band's response.

int spectra_read_srf(const char *file);
int spectra_srf_band(int srf, double lambda);
const char *spectra_srf_band_name(int srf, int band);
double spectra_srf_centroid(int srf, int band);
double spectra_water_srf(int table, int srf, int band);
double spectra_bottom_srf(int k, int srf, int band);

#endif
