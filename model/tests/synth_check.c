//
//    synth_check.c
//

/* Checks for the synthetic regression test of the semi-analytical model. */

// Usage:
//
//   synth_check baseline /output.nc/ /baseline.nc/
//       Fails (exit 1) unless the output matches the baseline grid: same size and
//       nodata cells, and |output - baseline| <= 1e-4 + 1e-3*|baseline| everywhere.
//
//   synth_check truth /modelled_H.nc/ /truth_H.asc/
//       Reports the depth error against the synthetic truth. photic writes depths
//       negative down; truth_H.asc is positive down. Informational: always exits 0.

#include "common.h"
#include "io.h"
#include "nc.h"

static int check_baseline(char *out_file, char *base_file) {

  float **out, **base, *lons, *lats, diff, max_diff = 0.0;
  double out_spval, base_spval;
  int nc_o, nr_o, nc_b, nr_b, i, j, n_bad = 0;
  bool out_nodata, base_nodata;

  read_nc(out_file, &out, &nc_o, &nr_o, &lons, &lats, &out_spval);
  read_nc(base_file, &base, &nc_b, &nr_b, &lons, &lats, &base_spval);

  if (nc_o != nc_b || nr_o != nr_b) {
    printf("FAIL  %s: grid is %d x %d, baseline is %d x %d\n", out_file, nc_o, nr_o, nc_b, nr_b);
    return 1;
  }

  for (i = 0; i < nr_b; i++) {
    for (j = 0; j < nc_b; j++) {
      out_nodata = out[i][j] == (float) out_spval;
      base_nodata = base[i][j] == (float) base_spval;
      if (out_nodata != base_nodata) {
        n_bad++;
        continue ;
      }
      if (base_nodata)
        continue ;
      diff = fabs(out[i][j] - base[i][j]);
      if (diff > max_diff)
        max_diff = diff;
      if (diff > 1.0e-4 + 1.0e-3*fabs(base[i][j]))
        n_bad++;
    }
  }

  printf("%s  %s: max |diff| = %.3g, cells outside tolerance = %d\n",
    n_bad ? "FAIL" : "ok  ", out_file, max_diff, n_bad);

  return n_bad ? 1 : 0;
}

static int check_truth(char *model_file, char *truth_file) {

  float **model, **truth, *lons, *lats, h, err, max_err = 0.0;
  double model_spval, truth_spval, wlon, slat, cellsize, sum_sq = 0.0, sum = 0.0;
  int nc_m, nr_m, nc_t, nr_t, i, j, n = 0, n_missing = 0;

  read_nc(model_file, &model, &nc_m, &nr_m, &lons, &lats, &model_spval);
  read_asc(truth_file, &truth, &nc_t, &nr_t, &wlon, &slat, &cellsize, &truth_spval);

  if (nc_m != nc_t || nr_m != nr_t) {
    printf("modelled grid is %d x %d, truth is %d x %d\n", nc_m, nr_m, nc_t, nr_t);
    return 0;
  }

  for (i = 0; i < nr_t; i++) {
    for (j = 0; j < nc_t; j++) {
      h = -model[i][j];
      if (model[i][j] == (float) model_spval || h <= 0.0) {
        n_missing++;
        continue ;
      }
      err = h - truth[i][j];
      sum += err;
      sum_sq += err*err;
      if (fabs(err) > max_err)
        max_err = fabs(err);
      n++;
    }
  }

  printf("\nAccuracy against the synthetic truth (informational):\n");
  printf("  pixels modelled = %d of %d\n", n, n + n_missing);
  if (n > 0) {
    printf("  bias            = %+.4f m\n", sum/n);
    printf("  RMSE            = %.4f m\n", sqrt(sum_sq/n));
    printf("  max abs error   = %.4f m\n", max_err);
  }

  return 0;
}

int main(int argc, char **argv) {

  if (argc == 4 && strcmp(argv[1], "baseline") == 0)
    return check_baseline(argv[2], argv[3]);

  if (argc == 4 && strcmp(argv[1], "truth") == 0)
    return check_truth(argv[2], argv[3]);

  printf("\nUsage: synth_check baseline /output.nc/ /baseline.nc/\n");
  printf("       synth_check truth /modelled_H.nc/ /truth_H.asc/\n\n");
  return 1;
}
