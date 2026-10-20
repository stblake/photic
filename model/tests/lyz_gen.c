//
//    lyz_gen.c
//

/* Synthetic test data for MODEL lyzenga_log_linear. */

// Two bands that follow Lyzenga's model exactly, L_m = Lsm_m + A r_m exp(-g_m Z), over a known
// depth field with a varying bottom brightness A, imaged at a tide of TIDE m (Z = depth + TIDE).
// Two bands remove A, so the fit must recover
//   H = a0 + a1 ln(L_1 - Lsm_1) + a2 ln(L_2 - Lsm_2),  a1 = -a2 = 1/(g_2 - g_1),
//   a0 = -a1 ln(r_1/r_2),
// and every sounding's depth. Writes b1.asc, b2.asc, soundings.xyz (depth below the datum,
// positive down, at every third node) and expected.txt (a0 a1 a2).
//
// Usage: lyz_gen /output directory/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NROWS 40
#define NCOLS 40

static const double wlon = 115.0, slat = -21.0, cellsize = 0.0003, tide = 0.5;
static const double Lsm[2] = {0.004, 0.002}, r[2] = {0.05, 0.04}, g[2] = {0.10, 0.25};

static double depth(int i, int j) {
  return 0.5 + 14.0*(i + 0.5*j)/(NROWS - 1 + 0.5*(NCOLS - 1));
}

static double brightness(int i, int j) {
  return 0.7 + 0.3*sin(0.37*i)*cos(0.23*j);
}

static int write_band(const char *dir, int m) {
  char file[512];
  FILE *fp;
  int i, j;
  snprintf(file, sizeof(file), "%s/b%d.asc", dir, m + 1);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  fprintf(fp, "ncols %d\nnrows %d\nxllcorner %.6f\nyllcorner %.6f\ncellsize %.6f\nNODATA_value -9999\n",
    NCOLS, NROWS, wlon, slat, cellsize);
  // ESRI ascii grids are north row first; row 0 is south in photic.
  for (i = NROWS - 1; i >= 0; i--) {
    for (j = 0; j < NCOLS; j++) {
      fprintf(fp, "%.10g ", Lsm[m] + brightness(i, j)*r[m]*exp(-g[m]*(depth(i, j) + tide)));
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

int main(int argc, char **argv) {

  char file[512];
  FILE *fp;
  int i, j;
  double a1 = 1.0/(g[1] - g[0]);

  if (argc != 2) {
    fprintf(stderr, "usage: lyz_gen /output directory/\n");
    return 1;
  }
  if (write_band(argv[1], 0) || write_band(argv[1], 1)) return 1;

  snprintf(file, sizeof(file), "%s/soundings.xyz", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  for (i = 0; i < NROWS; i += 3) {
    for (j = 0; j < NCOLS; j += 3) {
      fprintf(fp, "%.7f %.7f %.4f\n", wlon + j*cellsize, slat + i*cellsize, depth(i, j));
    }
  }
  fclose(fp);

  snprintf(file, sizeof(file), "%s/expected.txt", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  fprintf(fp, "%.10g %.10g %.10g %d\n", -a1*log(r[0]/r[1]), a1, -a1, ((NROWS + 2)/3)*((NCOLS + 2)/3));
  fclose(fp);
  return 0;
}
