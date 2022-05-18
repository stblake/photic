/*
 *
 *    Smoothing Algorithms
 *
 */

/* Written by Sam Blake. */

/* Started on 28 Feb 2016. */

#include "smoothutils.h"

float verify(float gp, float spval, float weight, float *total) {
  if (gp == spval) {
    return 0.0;
  } else {
    *total += weight; 
    return weight*gp;
  }
}

void smooth(float **input, float **smoothed, int ncols, int nrows, 
            float weight, int npasses, float spval) {

  int n, i, j;
  float **unsmoothed, f, tot, s0, s1, s2, s3, s4, s5, s6, s7, s8;

  f = 1.0/sqrt(2.0);

// Layout of stencil of points:
//
//              N
//         s0  s1  s2 
//      W  s3  s4  s5  E
//         s6  s7  s8
//              S

  allocate_float_array_2d(&unsmoothed, nrows, ncols);
  copy_float_array_2d(input, unsmoothed, nrows, ncols);
  copy_float_array_2d(unsmoothed, smoothed, nrows, ncols);

  for (n = 0; n < npasses; n++) {

  // Interior points.  

    for (i = 1; i < nrows - 1; i++) {
      for (j = 1; j < ncols - 1; j++) {
        if (unsmoothed[i][j] == spval) 
          continue;
        tot = 0.0; 
        s0 = verify(unsmoothed[i - 1][j - 1], spval, f, &tot);
        s1 = verify(unsmoothed[i - 1][j], spval, 1.0, &tot);
        s2 = verify(unsmoothed[i - 1][j + 1], spval, f, &tot);
        s3 = verify(unsmoothed[i][j - 1], spval, 1.0, &tot);
        s4 = verify(unsmoothed[i][j], spval, weight, &tot);
        s5 = verify(unsmoothed[i][j + 1], spval, 1.0, &tot);
        s6 = verify(unsmoothed[i + 1][j - 1], spval, f, &tot);
        s7 = verify(unsmoothed[i + 1][j], spval, 1.0, &tot);
        s8 = verify(unsmoothed[i + 1][j + 1], spval, f, &tot);
        smoothed[i][j] = (s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8)/tot;
      }
    }

  // Western boundary.

    for (i = 1; i < nrows - 1; i++) {
      if (unsmoothed[i][0] == spval) 
        continue;
      tot = 0.0;
      s0 = verify(unsmoothed[i - 1][0], spval, 1.0, &tot);
      s3 = verify(unsmoothed[i][0], spval, weight, &tot);
      s6 = verify(unsmoothed[i + 1][0], spval, 1.0, &tot);
      s1 = verify(unsmoothed[i - 1][1], spval, f, &tot);
      s4 = verify(unsmoothed[i][1], spval, 1.0, &tot);
      s7 = verify(unsmoothed[i + 1][1], spval, f, &tot);
      smoothed[i][0] = (s0 + s3 + s6 + s1 + s4 + s7)/tot;
    } 

  // Eastern boundary. 

    for (i = 1; i < nrows - 1; i++) {
      if (unsmoothed[i][ncols - 1] == spval) 
        continue;
      tot = 0.0;
      s2 = verify(unsmoothed[i - 1][ncols - 1], spval, 1.0, &tot);
      s5 = verify(unsmoothed[i][ncols - 1], spval, weight, &tot);
      s8 = verify(unsmoothed[i + 1][ncols - 1], spval, 1.0, &tot);
      s1 = verify(unsmoothed[i - 1][ncols - 2], spval, f, &tot);
      s4 = verify(unsmoothed[i][ncols - 2], spval, 1.0, &tot);
      s7 = verify(unsmoothed[i + 1][ncols - 2], spval, f, &tot);
      smoothed[i][ncols - 1] = (s2 + s5 + s8 + s1 + s4 + s7)/tot;
    } 

  // Southern boundary. 

    for (j = 0; j < ncols - 1; j++) {
      if (unsmoothed[0][j] == spval) 
        continue; 
      tot = 0.0;
      s6 = verify(unsmoothed[0][j - 1], spval, 1.0, &tot);
      s7 = verify(unsmoothed[0][j], spval, weight, &tot);
      s8 = verify(unsmoothed[0][j + 1], spval, 1.0, &tot);
      s3 = verify(unsmoothed[1][j - 1], spval, f, &tot);
      s4 = verify(unsmoothed[1][j], spval, 1.0, &tot);
      s5 = verify(unsmoothed[1][j + 1], spval, f, &tot);
      smoothed[0][j] = (s6 + s7 + s8 + s3 + s4 + s5)/tot;
    }

  // Northern boundary. 

    for (j = 0; j < ncols - 1; j++) {
      if (unsmoothed[nrows - 1][j] == spval) 
        continue; 
      tot = 0.0;
      s0 = verify(unsmoothed[nrows - 1][j - 1], spval, 1.0, &tot);
      s1 = verify(unsmoothed[nrows - 1][j], spval, weight, &tot);
      s2 = verify(unsmoothed[nrows - 1][j + 1], spval, 1.0, &tot);
      s3 = verify(unsmoothed[nrows - 2][j - 1], spval, f, &tot);
      s4 = verify(unsmoothed[nrows - 2][j], spval, 1.0, &tot);
      s5 = verify(unsmoothed[nrows - 2][j + 1], spval, f, &tot);
      smoothed[nrows - 1][j] = (s0 + s1 + s2 + s3 + s4 + s5)/tot;
    }

  // SW corner.

    if (unsmoothed[0][0] != spval) { 
      tot = 0.0;
      s3 = verify(unsmoothed[1][0], spval, 1.0, &tot);
      s4 = verify(unsmoothed[1][1], spval, f, &tot);
      s6 = verify(unsmoothed[0][0], spval, weight, &tot);
      s7 = verify(unsmoothed[0][1], spval, 1.0, &tot);
      smoothed[0][0] = (s3 + s4 + s6 + s7)/tot;
    }

  // NW corner. 

    if (unsmoothed[nrows - 1][0] != spval) { 
      tot = 0.0;
      s0 = verify(unsmoothed[nrows - 1][0], spval, weight, &tot);
      s1 = verify(unsmoothed[nrows - 1][1], spval, 1.0, &tot);
      s3 = verify(unsmoothed[nrows - 2][0], spval, 1.0, &tot);
      s4 = verify(unsmoothed[nrows - 2][1], spval, f, &tot);
      smoothed[nrows - 1][0] = (s0 + s1 + s3 + s4)/tot;
    }

  // NE corner. 

    if (unsmoothed[nrows - 1][ncols - 1] != spval) { 
      tot = 0.0;
      s1 = verify(unsmoothed[nrows - 1][ncols - 2], spval, 1.0, &tot);
      s2 = verify(unsmoothed[nrows - 1][ncols - 1], spval, weight, &tot);
      s4 = verify(unsmoothed[nrows - 2][ncols - 2], spval, f, &tot);
      s5 = verify(unsmoothed[nrows - 2][ncols - 1], spval, 1.0, &tot);
      smoothed[nrows - 1][ncols - 1] = (s1 + s2 + s4 + s5)/tot;
    }

  // SE corner. 

    if (unsmoothed[0][ncols - 1] != spval) { 
      tot = 0.0;
      s4 = verify(unsmoothed[1][ncols - 2], spval, f, &tot);
      s5 = verify(unsmoothed[1][ncols - 1], spval, 1.0, &tot);
      s7 = verify(unsmoothed[0][ncols - 2], spval, 1.0, &tot);
      s8 = verify(unsmoothed[0][ncols - 1], spval, weight, &tot);
      smoothed[0][ncols - 1] = (s4 + s5 + s7 + s8)/tot;
    }

  // Copy smoothed to unsmoothed. (Except for the final pass.)

    if (n < npasses - 1) 
      copy_float_array_2d(smoothed, unsmoothed, nrows, ncols);
  }

  free_float_array_2d(unsmoothed, nrows);
}




