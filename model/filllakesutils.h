

#include "common.h"
#include "io.h"

void fill_lakes(float **grid, int ncols, int nrows, int maxwet, 
                bool spvalpos, bool fillspval, float spval);

void fill_lake_onesize(float **grid, int ncols, int nrows, int maxwet, 
        bool spvalpos, bool fillspval, float spval, int xl, int yl, int *nfilled);