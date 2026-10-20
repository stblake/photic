#include "common.h"

// MODEL lyzenga_log_linear (lyzenga.c).
void run_model_lyzenga(void);

// The weighted least-squares solution of A x = y (A is n x p, row-major), by Householder QR:
// x, and the upper-triangular R (p x p, row-major) of sqrt(w) A. Returns 0, or 1 if A is rank
// deficient.
int lyzenga_wls(const double *A, const double *y, const double *w, int n, int p, double *x, double *R);
