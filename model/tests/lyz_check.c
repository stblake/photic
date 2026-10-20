//
//    lyz_check.c
//

/* Checks MODEL lyzenga_log_linear on lyz_gen's data: the coefficients in model.bam against
   expected.txt, and every sounding's residual and status in output.csv. */

// Usage: lyz_check /run directory/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {

  char file[512], line[1024], status[32];
  FILE *fp;
  double expected[3], fitted[3], residual, worst = 0.0;
  int n_expected, n_used = 0, n_other = 0, m, fail = 0;

  if (argc != 2) {
    fprintf(stderr, "usage: lyz_check /run directory/\n");
    return 1;
  }

  snprintf(file, sizeof(file), "%s/expected.txt", argv[1]);
  if ((fp = fopen(file, "r")) == NULL || fscanf(fp, "%lf %lf %lf %d", &expected[0], &expected[1], &expected[2],
      &n_expected) != 4) {
    fprintf(stderr, "lyz_check: cannot read %s\n", file);
    return 1;
  }
  fclose(fp);

  snprintf(file, sizeof(file), "%s/model.bam", argv[1]);
  if ((fp = fopen(file, "r")) == NULL) {
    fprintf(stderr, "lyz_check: cannot read %s\n", file);
    return 1;
  }
  m = 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (strncmp(line, "COEFFICIENTS=", 13) == 0) {
      m = sscanf(line + 13, "%lf,%lf,%lf", &fitted[0], &fitted[1], &fitted[2]);
    }
  }
  fclose(fp);
  if (m != 3) {
    fprintf(stderr, "lyz_check: no COEFFICIENTS in %s\n", file);
    return 1;
  }

  snprintf(file, sizeof(file), "%s/output.csv", argv[1]);
  if ((fp = fopen(file, "r")) == NULL) {
    fprintf(stderr, "lyz_check: cannot read %s\n", file);
    return 1;
  }
  fgets(line, sizeof(line), fp);
  while (fgets(line, sizeof(line), fp) != NULL) {
    // LON,LAT,ROW,COL,DEPTH,MODEL DEPTH,RESIDUAL,WEIGHT,STATUS
    char *field[9], *s = line;
    int k;
    for (k = 0; k < 9; k++) {
      field[k] = s;
      s = strchr(s, ',');
      if (s == NULL) break;
      *s++ = '\0';
    }
    if (k < 8) continue ;
    sscanf(field[8], "%31s", status);
    if (strcmp(status, "used") == 0) {
      n_used++;
      residual = atof(field[6]);
      if (fabs(residual) > worst) worst = fabs(residual);
    } else {
      n_other++;
    }
  }
  fclose(fp);

  printf("lyzenga: coefficients");
  for (m = 0; m < 3; m++) {
    printf(" %.6f (expected %.6f)", fitted[m], expected[m]);
    if (fabs(fitted[m] - expected[m]) > 1.0e-3*(1.0 + fabs(expected[m]))) fail = 1;
  }
  printf("\nlyzenga: %d soundings used (expected %d), %d not; largest residual %.2e m\n", n_used, n_expected,
    n_other, worst);
  if (n_used != n_expected || n_other != 0 || worst > 0.01) fail = 1;
  printf("lyzenga: %s\n", fail ? "FAIL" : "ok");
  return fail;
}
