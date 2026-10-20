//
//    att_gen.c
//

/* Synthetic test data for MODEL lyzenga_log_linear ... ATTENUATION. */

// Two scenes of three bands that follow the attenuation model exactly,
//   L_sb = Lsm_sb + A r_sb exp(-K_sb(x) M_s (depth + h_s)),
// with K_sb bilinear over the grid (so 2 x 2 nodes represent it exactly), tides h_s, sun zeniths
// theta_s (M_s = sec theta_w + 1, nadir view) and a bottom brightness A that is 1 at the
// soundings and varies elsewhere. Stage 1 must recover K at the nodes exactly, and stage 2 the
// depth at every cell. Writes s<s>_b<b>.asc, soundings.xyz (every third node), check.xyz (every
// interior node, with its depth) and expected_K.txt (K at the four nodes, nodes row-major from
// the south-west, grids scene-major) and index_s1.xyz, scene 1's two depth-invariant bottom indices
// (Lyzenga 1978, INDICES) at every interior node, exact because the depth term cancels.
//
// Usage: att_gen /output directory/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NROWS 40
#define NCOLS 40
#define NS 2
#define NB 3

static const double wlon = 1.0, slat = 1.0, cellsize = 0.0003;  // near 0: photic keeps soundings as floats
static const double tide[NS] = {0.3, -0.4}, sunz[NS] = {30.0, 50.0};
static const double Lsm[NS][NB] = {{0.004, 0.003, 0.002}, {0.005, 0.0035, 0.0025}};
static const double r[NS][NB] = {{0.05, 0.045, 0.04}, {0.055, 0.05, 0.035}};
static const double K0[NS][NB] = {{0.08, 0.12, 0.30}, {0.10, 0.15, 0.35}};
static const double ka[NS][NB] = {{0.3, -0.2, 0.1}, {-0.1, 0.25, 0.2}};
static const double kb[NS][NB] = {{0.2, 0.1, -0.15}, {0.3, -0.1, 0.05}};
static const double kc[NS][NB] = {{0.1, -0.05, 0.2}, {0.0, 0.15, -0.1}};

static double depth(int i, int j) {
  return 0.5 + 14.0*(i + 0.5*j)/(NROWS - 1 + 0.5*(NCOLS - 1));
}

static double brightness(int i, int j) {
  return i%3 == 0 && j%3 == 0 ? 1.0 : 0.6 + 0.5*sin(0.37*i)*cos(0.23*j)*sin(0.11*(i + j));
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
    fprintf(stderr, "usage: att_gen /output directory/\n");
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
          fprintf(fp, "%.10g ", Lsm[s][b] + brightness(i, j)*r[s][b]*exp(-k*path(s)*(depth(i, j) + tide[s])));
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

  // Scene 1's indices: Y_i = sum_j A_ij X_j with b_j = K_j (the path M cancels in the rotation).
  snprintf(file, sizeof(file), "%s/index_s1.xyz", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  for (i = 1; i < NROWS - 1; i++) {
    for (j = 1; j < NCOLS - 1; j++) {
      double bb[NB], x[NB], y[NB - 1], S1, S2;
      int l, ix;
      for (l = 0; l < NB; l++) {
        bb[l] = K(0, l, j/(NCOLS - 1.0), i/(NROWS - 1.0));
        x[l] = log(brightness(i, j)*r[0][l]) - bb[l]*path(0)*(depth(i, j) + tide[0]);
      }
      for (ix = 0; ix < NB - 1; ix++) {
        S1 = 0.0;
        for (l = 0; l <= ix; l++) S1 += bb[l]*bb[l];
        S2 = S1 + bb[ix + 1]*bb[ix + 1];
        y[ix] = -sqrt(S1/S2)*x[ix + 1];
        for (l = 0; l <= ix; l++) y[ix] += bb[ix + 1]*bb[l]/sqrt(S1*S2)*x[l];
      }
      fprintf(fp, "%.7f %.7f %.8f %.8f\n", wlon + j*cellsize, slat + i*cellsize, y[0], y[1]);
    }
  }
  fclose(fp);

  snprintf(file, sizeof(file), "%s/expected_K.txt", argv[1]);
  if ((fp = fopen(file, "w")) == NULL) return 1;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      for (s = 0; s < NS; s++) {
        for (b = 0; b < NB; b++) fprintf(fp, "%.10g ", K(s, b, j, i));
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return 0;
}
