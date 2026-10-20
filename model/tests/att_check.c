//
//    att_check.c
//

/* Checks MODEL lyzenga_log_linear ... ATTENUATION on att_gen's data: K at the nodes in
   attenuation.csv against expected_K.txt, and the depth from the attenuation at every interior
   node (depth_points.xyz, from WRITE POINTS) against check.xyz, and scene 1's depth-invariant
   indices (index<i>_points.xyz, INDICES) against index_s1.xyz. */

// Usage: att_check /run directory/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NGRIDS 6
#define NNODES 4

int main(int argc, char **argv) {

  char file[512], line[4096], *tok;
  FILE *fp, *fq;
  double expected[NNODES][NGRIDS], got, worst_k = 0.0, worst_h = 0.0, worst_y = 0.0, x, y, z, xg, yg, zg, y1, y2;
  FILE *fr;
  int n_index = 0;
  int n, m, n_points = 0, fail = 0;

  if (argc != 2) {
    fprintf(stderr, "usage: att_check /run directory/\n");
    return 1;
  }

  snprintf(file, sizeof(file), "%s/expected_K.txt", argv[1]);
  if ((fp = fopen(file, "r")) == NULL) {
    fprintf(stderr, "att_check: cannot read %s\n", file);
    return 1;
  }
  for (n = 0; n < NNODES; n++) {
    for (m = 0; m < NGRIDS; m++) {
      if (fscanf(fp, "%lf", &expected[n][m]) != 1) {
        fprintf(stderr, "att_check: short %s\n", file);
        return 1;
      }
    }
  }
  fclose(fp);

  // attenuation.csv: LON,LAT,SOUNDINGS, then K of each grid, one node per row.
  snprintf(file, sizeof(file), "%s/attenuation.csv", argv[1]);
  if ((fp = fopen(file, "r")) == NULL || fgets(line, sizeof(line), fp) == NULL) {
    fprintf(stderr, "att_check: cannot read %s\n", file);
    return 1;
  }
  for (n = 0; n < NNODES; n++) {
    if (fgets(line, sizeof(line), fp) == NULL) {
      fprintf(stderr, "att_check: %s has fewer than %d nodes\n", file, NNODES);
      return 1;
    }
    tok = strtok(line, ",");
    for (m = 0; m < 2 && tok != NULL; m++) tok = strtok(NULL, ",");
    for (m = 0; m < NGRIDS; m++) {
      if ((tok = strtok(NULL, ",\n")) == NULL) {
        fprintf(stderr, "att_check: %s: node %d has no K for grid %d\n", file, n, m);
        return 1;
      }
      got = atof(tok);
      worst_k = fmax(worst_k, fabs(got - expected[n][m])/expected[n][m]);
    }
  }
  fclose(fp);

  // The depth at every interior node, against the truth (photic writes negative down).
  snprintf(file, sizeof(file), "%s/check.xyz", argv[1]);
  fp = fopen(file, "r");
  snprintf(file, sizeof(file), "%s/depth_points.xyz", argv[1]);
  fq = fopen(file, "r");
  if (fp == NULL || fq == NULL) {
    fprintf(stderr, "att_check: cannot read check.xyz or depth_points.xyz\n");
    return 1;
  }
  while (fscanf(fp, "%lf %lf %lf", &x, &y, &z) == 3) {
    if (fscanf(fq, "%lf, %lf, %lf", &xg, &yg, &zg) != 3) {
      fprintf(stderr, "att_check: depth_points.xyz is short\n");
      return 1;
    }
    worst_h = fmax(worst_h, fabs(-zg - z));
    n_points++;
  }
  fclose(fp);
  fclose(fq);

  snprintf(file, sizeof(file), "%s/index_s1.xyz", argv[1]);
  fp = fopen(file, "r");
  snprintf(file, sizeof(file), "%s/index1_points.xyz", argv[1]);
  fq = fopen(file, "r");
  snprintf(file, sizeof(file), "%s/index2_points.xyz", argv[1]);
  fr = fopen(file, "r");
  if (fp == NULL || fq == NULL || fr == NULL) {
    fprintf(stderr, "att_check: cannot read index_s1.xyz or index<i>_points.xyz\n");
    return 1;
  }
  while (fscanf(fp, "%lf %lf %lf %lf", &x, &y, &y1, &y2) == 4) {
    if (fscanf(fq, "%lf, %lf, %lf", &xg, &yg, &zg) != 3) break;
    worst_y = fmax(worst_y, fabs(zg - y1));
    if (fscanf(fr, "%lf, %lf, %lf", &xg, &yg, &zg) != 3) break;
    worst_y = fmax(worst_y, fabs(zg - y2));
    n_index++;
  }
  fclose(fp);
  fclose(fq);
  fclose(fr);

  printf("attenuation: largest relative K error at the nodes %.2e; depth at %d cells, largest error %.2e m; "
    "indices at %d cells, largest error %.2e\n", worst_k, n_points, worst_h, n_index, worst_y);
  // The indices are combinations of ln(L - Lsm), so float rounding gives ~1e-3 in the red at depth.
  if (n_index < 1000 || worst_y > 2.0e-3) {
    printf("attenuation: FAIL, the depth-invariant indices are not recovered\n");
    fail = 1;
  }
  if (worst_k > 1.0e-5) {
    printf("attenuation: FAIL, K is not recovered\n");
    fail = 1;
  }
  // The band grids are floats: at 14 m in the red, L - Lsm is ~1e-7 on 0.0025, so float rounding
  // alone gives ~1e-3 m. A bottom brightness leaking into the depth would give metres.
  if (n_points < 1000 || worst_h > 5.0e-3) {
    printf("attenuation: FAIL, the depth is not recovered\n");
    fail = 1;
  }
  if (! fail) printf("attenuation: ok\n");
  return fail;
}
