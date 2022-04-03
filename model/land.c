
// Estimate land/sea threshold in the NIR band

#include "land.h"

float land_sea_threshold(float **array, int nrows, int ncols) {

	int i, j;
	float mean_cutoff = 0.0, max_diff, diff, *cutoffY, *cutoffX;

	cutoffY = malloc(ncols*sizeof(float));
	cutoffX = malloc(nrows*sizeof(float));

	// Along the rows. 

	for (i = 0; i < nrows; i++) {
		max_diff = 0.0;
		cutoff = 0.0;
		for (j = 1; j < ncols - 1; j++) {
			diff = fabs(array[i][j - 1] - 2.0*array[i][j] + array[i][j + 1]);
			if (diff > max_diff) {
				max_diff = diff; 
				cutoff = array[i][j];
			}
		}
		cutoffX[i] = cutoff; 
	}

	// Along the cols. 

	for (j = 0; j < ncols; j++) {
		max_diff = 0.0;
		cutoff = 0.0;
		for (i = 1; i < nrows - 1; i++) {
			diff = fabs(array[i - 1][j] - 2.0*array[i][j] + array[i + 1][j]);
			if (diff > max_diff) {
				max_diff = diff; 
				cutoff = array[i][j];
			}
		}
		cutoffY[j] = cutoff; 
	}

}