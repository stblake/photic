
#include "interp.h"
#include "common.h"

// Code for interpolating from one array to another. 

void grid2gridinterp(
  float **grid_cg, int nx_cg, int ny_cg, 
  double wlon_cg, double slat_cg, double cellsize_cg, double spval_cg,
  float **grid_fg, int nx_fg, int ny_fg, 
  double wlon_fg, double slat_fg, double cellsize_fg, double spval_fg) {

  int i, j;
  float v;
  double elon_cg, nlat_cg, lon, lat;

  elon_cg = wlon_cg + ((double) nx_cg - 1)*cellsize_cg;
  nlat_cg = slat_cg + ((double) ny_cg - 1)*cellsize_cg;

  for (i = 0; i < ny_fg; i++){
    for (j = 0; j < nx_fg; j++){
      lon = wlon_fg + ((double) j)*cellsize_fg;
      lat = slat_fg + ((double) i)*cellsize_fg;
      if (lon >= wlon_cg && lon <= elon_cg && lat >= slat_cg && lat <= nlat_cg) {
        v = interp_bilinear(grid_cg, nx_cg, ny_cg, wlon_cg, slat_cg, cellsize_cg, spval_cg, lon, lat);
        if (approx_equal(v, spval_cg, 1.0e-4)) 
          grid_fg[i][j] = (float) spval_fg;
        else
          grid_fg[i][j] = v;
      } else {
        grid_fg[i][j] = (float) spval_fg;
      }
    }
  }
}


// Code for bicubic array interpolation.

float interp_bicubic(float **grid, int nx, int ny, double wlon, double slat, double cellsize, 
      double spval, double lon, double lat) {

  float x, y;

  x = (lon - wlon)/cellsize;
  y = (lat - slat)/cellsize;

  return interp_bicubic_raw(grid, nx, ny, x, y, (float) spval);
}


// Code for bilinear array interpolation. 

float interp_bilinear(float **grid, int nx, int ny, double wlon, double slat, double cellsize, 
      double spval, double lon, double lat) {

  float x, y;

  x = (lon - wlon)/cellsize;
  y = (lat - slat)/cellsize;

  return interp_bilinear_raw(grid, nx, ny, x, y, (float) spval);
}

float interp_bilinear_raw(float **grid, int nx, int ny, float xf, float yf, float spval) {

  int i, j;
  float eps = 1.0e-6, x = xf, y = yf, A, B, C, D, a, b, c, d, interp;

  j = floor(x);
  i = floor(y);

  if (j == nx - 1) j--; 
  if (i == ny - 1) i--;

  // values at the nodes. 
  //
  //  C   |    D
  //      |
  //   -------
  //      |
  //  A   |    B

  A = grid[i][j];
  B = grid[i][j + 1];
  C = grid[i + 1][j];
  D = grid[i + 1][j + 1];

  if (approx_equal(A, spval, 1.0e-4) && approx_equal(B, spval, 1.0e-4) && 
      approx_equal(C, spval, 1.0e-4) && approx_equal(D, spval, 1.0e-4)) {
    return spval;
  }

  // bilinear interpolant

  x -= floor(x);
  y -= floor(y);

  if (approx_equal(A, spval, 1.0e-4) || approx_equal(B, spval, 1.0e-4) || 
      approx_equal(C, spval, 1.0e-4) || approx_equal(D, spval, 1.0e-4)) {

    if (x < 0.5 + eps && y < 0.5 + eps) {
      interp = A;
    } else if (x > 0.5 - eps && y < 0.5 + eps) {
      interp = B;
    } else if (x < 0.5 + eps && y > 0.5 - eps) {
      interp = C;
    } else if (x > 0.5 - eps && y > 0.5 - eps) {
      interp = D;
    }
  } else {
    interp = A*(1.0 - x)*(1.0 - y) + B*x*(1.0 - y) + C*(1.0 - x)*y + D*x*y;
  }

  return interp;
}



// Lagrangian bi-cubic interpolation. 

float interp_bicubic_raw(float **grid, int nx, int ny, float x, float y, float spval)
  {
    int msize = ny - 1;
    int nsize = nx - 1;
    int mmax = msize;
    int nmax = nsize;
    int mmin = 1;
    int nmin = 1;
    float gn = x + 1.0, gm = y + 1.0; 
    double dgm = floor(gm);
    double dgn = floor(gn);

    int igm = (int)dgm;
    int jgn = (int)dgn;
    float fm = gm - (float)igm;
    float fn = gn - (float)jgn;
    if (fm < 1.e-06) {
      fm = 0.0f;
    }
    if (fn < 1.e-06) {
      fn = 0.0f;
    }
    int ms = mmax - 1;
    int ns = nmax - 1;
    int mr = mmin + 1;
    int nr = nmin + 1;

    float e = 0.0f;
    float t1 = 0.0f;
    float t2 = 0.0f;
    float p = 0.0f;
    float h = 0.0f;
    float scinto = 0.0f;

    if (gm >= mmax) {
      if (gn >= nmax) {
        if (grid[mmax - 1][nmax - 1] == spval || 
            grid[ms - 1][nmax - 1] == spval || 
            grid[mmax - 1][ns - 1] == spval) 
        {
          return spval;
        } 
        e = gm - (float)mmax;
        t1 =
          e
            * (grid[mmax - 1][nmax - 1] - grid[ms - 1][nmax - 1]);
        e = gn - (float)nmax;
        t2 =
          e
            * (grid[mmax - 1][nmax - 1] - grid[mmax - 1][ns - 1]);
        scinto = grid[mmax - 1][nmax - 1] + t1 + t2;
        return scinto;
      } else if (gn < nmin) {
        if (grid[mmax - 1][nmin - 1] == spval || 
            grid[ms - 1][nmin - 1] == spval || 
            grid[mmax - 1][nr - 1] == spval)
        {
          return spval;
        }
        e = gm - (float)mmax;
        t1 =
          e
            * (grid[mmax - 1][nmin - 1] - grid[ms - 1][nmin - 1]);
        e = (float)nmin - gn;
        t2 =
          e
            * (grid[mmax - 1][nmin - 1] - grid[mmax - 1][nr - 1]);
        scinto = grid[mmax - 1][nmin - 1] + t1 + t2;
        return scinto;
      } else {
        if (grid[mmax - 1][jgn - 1] == spval || 
            grid[mmax - 1][jgn]  == spval || 
            grid[ms - 1][jgn - 1] == spval || 
            grid[ms - 1][jgn] == spval || 
            grid[ms - 1][jgn - 1] == spval)
        {
          return spval;
        }
        p =
          grid[mmax - 1][jgn - 1]
            + fn
            * (grid[mmax - 1][jgn] - grid[mmax - 1][jgn - 1]);
        h =
          grid[ms - 1][jgn - 1]
            + fn
            * (grid[ms - 1][jgn] - grid[ms - 1][jgn - 1]);
        e = gm - (float)mmax;
        scinto = p + e * (p - h);
        return scinto;
      }
    } else if (gm < mmin) {
      if (gn >= nmax) {
        if (grid[mmin - 1][nmax - 1] == spval || 
            grid[mmin - 1][ns - 1] == spval || 
            grid[mr - 1][nmax - 1] == spval)
        {
          return spval;
        }
        e = gn - (float)nmax;
        t2 =
          e
            * (grid[mmin - 1][nmax - 1] - grid[mmin - 1][ns - 1]);
        e = (float)mmin - gm;
        t1 =
          e
            * (grid[mmin - 1][nmax - 1] - grid[mr - 1][nmax - 1]);
        scinto = grid[mmin - 1][nmax - 1] + t1 + t2;
        return scinto;
      } else if (gn < nmin) {
        if (grid[mmin - 1][nmin - 1] == spval || 
            grid[mmin - 1][nr - 1] == spval || 
            grid[mr - 1][nmin - 1] == spval)
        {
          return spval;
        }
        e = (float)nmin - gn;
        t2 =
          e
            * (grid[mmin - 1][nmin - 1] - grid[mmin - 1][nr - 1]);
        e = (float)mmin - gm;
        t1 =
          e
            * (grid[mmin - 1][nmin - 1] - grid[mr - 1][nmin - 1]);
        scinto = grid[mmin - 1][nmin - 1] + t1 + t2;
        return scinto;
      } else {
        if (grid[mmin - 1][jgn - 1] == spval || 
            grid[mmin - 1][jgn] == spval || 
            grid[mr - 1][jgn - 1] == spval || 
            grid[mr - 1][jgn] == spval)
        {
          return spval;
        }
        e = (float)mmin - gm;
        p =
          grid[mmin - 1][jgn - 1]
            + fn
            * (grid[mmin - 1][jgn] - grid[mmin - 1][jgn - 1]);
        h =
          grid[mr - 1][jgn - 1]
            + fn
            * (grid[mr - 1][jgn] - grid[mr - 1][jgn - 1]);
        scinto = p - e * (h - p);
        return scinto;
      }
    } else if (gn >= nmax) {
      if (grid[igm - 1][nmax - 1] == spval || 
          grid[igm][nmax - 1] == spval || 
          grid[igm - 1][ns - 1] == spval || 
          grid[igm][ns - 1] == spval)
      {
        return spval;
      }
      e = gn - (float)nmax;
      p =
        grid[igm - 1][nmax - 1]
          + fm
          * (grid[igm][nmax - 1] - grid[igm - 1][nmax - 1]);
      h =
        grid[igm - 1][ns - 1]
          + fm
          * (grid[igm][ns - 1] - grid[igm - 1][ns - 1]);
      scinto = p + e * (p - h);
      return scinto;
    } else if (gn < nmin) {
      if (grid[igm - 1][nmin - 1] == spval || 
          grid[igm][nmin - 1] == spval || 
          grid[igm - 1][nr - 1] == spval || 
          grid[igm][nr - 1] == spval)
      {
        return spval;
      }
      e = (float)nmin - gn;
      p =
        grid[igm - 1][nmin - 1]
          + fm
          * (grid[igm][nmin - 1] - grid[igm - 1][nmin - 1]);
      h =
        grid[igm - 1][nr - 1]
          + fm
          * (grid[igm][nr - 1] - grid[igm - 1][nr - 1]);
      scinto = p - e * (h - p);
      return scinto;
    } else if ((gm >= ms) || (gm < mr) || (gn >= ns) || (gn < nr)) {
      if (grid[igm][jgn - 1] == spval || 
          grid[igm][jgn] == spval || 
          grid[igm - 1][jgn - 1] == spval || 
          grid[igm - 1][jgn] == spval)
      {
        return spval;
      }
      p =
        grid[igm][jgn - 1]
          + fn
          * (grid[igm][jgn] - grid[igm][jgn - 1]);
      h =
        grid[igm - 1][jgn - 1]
          + fn
          * (grid[igm - 1][jgn] - grid[igm - 1][jgn - 1]);
      scinto = h + fm * (p - h);
      return scinto;
    } else {
      if (grid[igm - 2][jgn - 2] == spval || 
          grid[igm - 1][jgn - 2] == spval || 
          grid[igm][jgn - 2] == spval || 
          grid[igm + 1][jgn - 2] == spval || 
          grid[igm - 2][jgn - 1] == spval || 
          grid[igm - 1][jgn - 1] == spval || 
          grid[igm][jgn - 1] == spval || 
          grid[igm + 1][jgn - 1] == spval || 
          grid[igm - 2][jgn] == spval || 
          grid[igm - 1][jgn] == spval || 
          grid[igm][jgn] == spval || 
          grid[igm + 1][jgn] == spval || 
          grid[igm - 2][jgn + 1] == spval || 
          grid[igm - 1][jgn + 1] == spval || 
          grid[igm][jgn + 1] == spval || 
          grid[igm + 1][jgn + 1] == spval)
      {
        return spval;
      }
      float s1 = fm + 1.0f;
      float s2 = fm;
      float s3 = fm - 1.0f;
      float s4 = fm - 2.0f;
      float s12 = s1 * s2;
      float s34 = s3 * s4;
      float a = -s2 * s34;
      float b = 3.0f * s1 * s34;
      float c = -3.0f * s12 * s4;
      float d = s12 * s3;
      float x1 =
        a * grid[igm - 2][jgn - 2] + b
          * grid[igm - 1][jgn - 2] + c
          * grid[igm][jgn - 2] + d
          * grid[igm + 1][jgn - 2];
      float x2 =
        a * grid[igm - 2][jgn - 1] + b
          * grid[igm - 1][jgn - 1] + c
          * grid[igm][jgn - 1] + d
          * grid[igm + 1][jgn - 1];
      float x3 =
        a * grid[igm - 2][jgn] + b
          * grid[igm - 1][jgn] + c
          * grid[igm][jgn] + d
          * grid[igm + 1][jgn];
      float x4 =
        a * grid[igm - 2][jgn + 1] + b
          * grid[igm - 1][jgn + 1] + c
          * grid[igm][jgn + 1] + d
          * grid[igm + 1][jgn + 1];
      s1 = fn + 1.0f;
      s2 = fn;
      s3 = fn - 1.0f;
      s4 = fn - 2.0f;
      s12 = s1 * s2;
      s34 = s3 * s4;
      a = -s2 * s34;
      b = 3.0f * s1 * s34;
      c = -3.0f * s12 * s4;
      d = s12 * s3;
      float y = a * x1 + b * x2 + c * x3 + d * x4;
      scinto = y / 36.0f;
      return scinto;
    }
  }




