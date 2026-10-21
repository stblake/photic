//
//    con_gen.c
//

/* Synthetic test data for MODEL lyzenga_log_linear ... ATTENUATION ... CONSISTENT. */

// Three scenes of three bands that follow the attenuation model exactly,
//   L_sb = Lsm_sb + r_b(x) exp(-K_sb(x) M_s (depth + h)),
// with K_sb bilinear over the grid (2 x 2 nodes represent it exactly) and different in each scene,
// and a bottom r_b(x) whose colour changes from cell to cell, at the soundings too, so that neither
// stage 1 (bottom in the residuals) nor a flat brightness fits it. It is the same bottom in every
// scene. The tides are equal, so the mean attenuation over the scenes, which the images do not
// determine, does not enter the depth: CONSISTENT must recover the scenes' differences of the
// two-way attenuation, g_sb - mean_s g_sb, at the nodes, and the depth at every cell. Writes
// s<s>_b<b>.asc, soundings.xyz (every third node), check.xyz (every interior node, with its depth)
// and expected_gamma.txt (M_s, then g_sb - mean_s g_sb at the four nodes, row-major from the
// south-west, grids scene-major).
//
// Usage: con_gen /output directory/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NROWS 40
#define NCOLS 40
#define NS 3
#define NB 3

static const double wlon = 1.0, slat = 1.0, cellsize = 0.0003;
static const double tide = 0.2, sunz[NS] = {30.0, 50.0, 40.0};
static const double Lsm[NS][NB] = {{0.004, 0.003, 0.002}, {0.005, 0.0035, 0.0025}, {0.0045, 0.0032, 0.0022}};
static const double r0[NB] = {0.05, 0.045, 0.04};
static const double K0[NS][NB] = {{0.08, 0.12, 0.30}, {0.06, 0.09, 0.26}, {0.11, 0.15, 0.34}};
static const double ka[NS][NB] = {{0.3, -0.2, 0.1}, {-0.1, 0.25, 0.2}, {0.15, 0.1, -0.2}};
static const double kb[NS][NB] = {{0.2, 0.1, -0.15}, {0.3, -0.1, 0.05}, {-0.2, 0.2, 0.1}};
static const double kc[NS][NB] = {{0.1, -0.05, 0.2}, {0.0, 0.15, -0.1}, {-0.1, 0.05, 0.15}};

static double depth(int i, int j) {
  return 0.5 + 14.0*(i + 0.5*j)/(NROWS - 1 + 0.5*(NCOLS - 1));
}

// The bottom of band b: a colour that changes across the grid, differently in each band.
static double bottom(int b, int i, int j) {
  return r0[b]*exp(0.4*sin(0.21*i + 1.3*b)*cos(0.17*j - 0.5*b) + 0.2*sin(0.07*(i - j) + b));
}

static double K(int s, int b, double u, double v) {
  return K0[s][b]*(1.0 + ka[s][b]*u + kb[s][b]*v + kc[s][b]*u*v);
}

static double path(int s) {
  return 1.0/cos(asin(sin(sunz[s]*M_PI/180.0)/1.34)) + 1.0;
}

int main(int argc, char **argv) {

  char file[512];
  FILE *fp;
  int s, b, i, j;

  if (argc != 2) {
    fprintf(stderr, "usage: con_gen /output directory/\n");
    return 1;
  }
  for (s = 0; s < NS; s++) {
    for (b = 0; b < NB; b++) {
      snprintf(file, sizeof(file), "%s/s%d_b%d.asc", argv[1], s + 1, b + 1);
      if ((fp = fopen(file, "w")) == NULL) return 1;
      fprintf(fp, "ncols %d\nnrows %d\nxllcorner %.6f\nyllcorner %.6f\ncellsize %.6f\nNODATA_value -9999\n",
        NCOLS, NROWS, wlon, slat, cellsize);
      // ESRI ascii grids are north row first; row 0 is south in photic.
      for (i = NROWS - 1; i >= 0; i--) {
        for (j = 0; j < NCOLS; j++) {
          double k = K(s, b, j/(NCOLS - 1.0), i/(NROWS - 1.0));
          fprintf(fp, "%.10g ", Lsm[s][b] + bottom(b, i, j)*exp(-k*path(s)*(depth(i, j) + tide)));
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  snprintf(file, sizeof(file), "%s/soundings.xyz", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  for (i = 0; i < NROWS; i += 3) {
    for (j = 0; j < NCOLS; j += 3) fprintf(fp, "%.7f %.7f %.4f\n", wlon + j*cellsize, slat + i*cellsize, depth(i, j));
  }
  fclose(fp);

  snprintf(file, sizeof(file), "%s/check.xyz", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  for (i = 1; i < NROWS - 1; i++) {
    for (j = 1; j < NCOLS - 1; j++) fprintf(fp, "%.7f %.7f %.6f\n", wlon + j*cellsize, slat + i*cellsize, depth(i, j));
  }
  fclose(fp);

  snprintf(file, sizeof(file), "%s/expected_gamma.txt", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  for (s = 0; s < NS; s++) fprintf(fp, "%.10g ", path(s));
  fprintf(fp, "\n");
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      for (s = 0; s < NS; s++) {
        for (b = 0; b < NB; b++) {
          double mean = 0.0;
          int t;
          for (t = 0; t < NS; t++) mean += K(t, b, j, i)*path(t)/NS;
          fprintf(fp, "%.10g ", K(s, b, j, i)*path(s) - mean);
        }
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return 0;
}
