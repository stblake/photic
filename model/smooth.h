

#include "common.h"
#include "io.h"

void smooth(float **unsmoothed, float **smoothed, int ncols, int nrows, 
            float weight, int npasses, float spval);
inline float verify(float gp, float spval, float weight, float *total);