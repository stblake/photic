//
//    con_check.c
//

/* Checks MODEL lyzenga_log_linear ... ATTENUATION ... CONSISTENT on con_gen's data: the scenes'
   differences of the two-way attenuation at the nodes (from K in attenuation.csv) against
   expected_gamma.txt, and the depth at every interior node (depth_points.xyz, from WRITE POINTS)
   against check.xyz. */

// Usage: con_check /run directory/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NS 3
#define NB 3
#define NGRIDS (NS*NB)
#define NNODES 4

int main(int argc, char **argv) {

  char file[512], line[4096], *tok;
  FILE *fp, *fq;
  double M[NS], expected[NNODES][NGRIDS], g[NGRIDS], worst_g = 0.0, gmax = 0.0, worst_h = 0.0, x, y, z, xg, yg, zg;
  int n, m, s, b, n_points = 0, fail = 0;

  if (argc != 2) {
    fprintf(stderr, "usage: con_check /run directory/\n");
    return 1;
  }

  snprintf(file, sizeof(file), "%s/expected_gamma.txt", argv[1]);
  if ((fp = fopen(file, "r")) == NULL) {
    fprintf(stderr, "con_check: cannot read %s\n", file);
    return 1;
  }
  for (s = 0; s < NS; s++) {
    if (fscanf(fp, "%lf", &M[s]) != 1) return 1;
  }
  for (n = 0; n < NNODES; n++) {
    for (m = 0; m < NGRIDS; m++) {
      if (fscanf(fp, "%lf", &expected[n][m]) != 1) {
        fprintf(stderr, "con_check: short %s\n", file);
        return 1;
      }
    }
  }
  fclose(fp);

  // attenuation.csv: LON,LAT,SOUNDINGS, then K of each grid, one node per row; g = K M_s.
  snprintf(file, sizeof(file), "%s/attenuation.csv", argv[1]);
  if ((fp = fopen(file, "r")) == NULL || fgets(line, sizeof(line), fp) == NULL) {
    fprintf(stderr, "con_check: cannot read %s\n", file);
    return 1;
  }
  for (n = 0; n < NNODES; n++) {
    if (fgets(line, sizeof(line), fp) == NULL) {
      fprintf(stderr, "con_check: %s has fewer than %d nodes\n", file, NNODES);
      return 1;
    }
    tok = strtok(line, ",");
    for (m = 0; m < 2 && tok != NULL; m++) tok = strtok(NULL, ",");
    for (m = 0; m < NGRIDS; m++) {
      if ((tok = strtok(NULL, ",\n")) == NULL) {
        fprintf(stderr, "con_check: %s: node %d has no K for grid %d\n", file, n, m);
        return 1;
      }
      g[m] = atof(tok)*M[m/NB];
      gmax = fmax(gmax, g[m]);
    }
    for (b = 0; b < NB; b++) {
      double mean = 0.0;
      for (s = 0; s < NS; s++) mean += g[s*NB + b]/NS;
      for (s = 0; s < NS; s++) worst_g = fmax(worst_g, fabs(g[s*NB + b] - mean - expected[n][s*NB + b]));
    }
  }
  fclose(fp);

  snprintf(file, sizeof(file), "%s/check.xyz", argv[1]);
  fp = fopen(file, "r");
  snprintf(file, sizeof(file), "%s/depth_points.xyz", argv[1]);
  fq = fopen(file, "r");
  if (fp == NULL || fq == NULL) {
    fprintf(stderr, "con_check: cannot read check.xyz or depth_points.xyz\n");
    return 1;
  }
  while (fscanf(fp, "%lf %lf %lf", &x, &y, &z) == 3) {
    if (fscanf(fq, "%lf, %lf, %lf", &xg, &yg, &zg) != 3) {
      fprintf(stderr, "con_check: depth_points.xyz is short\n");
      return 1;
    }
    worst_h = fmax(worst_h, fabs(-zg - z));
    n_points++;
  }
  fclose(fp);
  fclose(fq);

  printf("consistent: largest error in the scenes' differences of attenuation %.2e (of %.3f m^-1); depth at %d cells, "
    "largest error %.2e m\n", worst_g, gmax, n_points, worst_h);
  if (worst_g > 1.0e-5*gmax) {
    printf("consistent: FAIL, the scenes' differences of attenuation are not recovered\n");
    fail = 1;
  }
  // As att_check: float band grids give ~1e-3 m in the red at depth; a bottom colour leaking into
  // the depth would give metres.
  if (n_points < 1000 || worst_h > 5.0e-3) {
    printf("consistent: FAIL, the depth is not recovered\n");
    fail = 1;
  }
  if (! fail) printf("consistent: ok\n");
  return fail;
}
