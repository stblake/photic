
#include "filllakesutils.h"

/* 
  For a given enclosing box of size yl by xl, the number of enclosed 
  points is (xl - 2)(yl - 2).
*/

void fill_lakes(float **grid, int ncols, int nrows, int maxwet, 
                bool spvalpos, bool fillspval, float spval) {

  int xl, yl, nfilled = 0; 

// Iterate over all bounding rectangle sizes. 

  for (xl = 3; xl < maxwet + 3; xl++) {
    for (yl = 3; yl < maxwet + 3; yl++) {
      if ((xl - 2)*(yl - 2) > maxwet) 
        continue ; 
      fill_lake_onesize(grid, ncols, nrows, 
        maxwet, spvalpos, fillspval, spval, xl, yl, &nfilled);
    }
  }

  printf("\nFilled %d points.\n\n", nfilled);
}

void fill_lake_onesize(float **grid, int ncols, int nrows, int maxwet, 
        bool spvalpos, bool fillspval, float spval, int xl, int yl, int *nfilled) {
  
  int i, j, n, m;
  float interp, a, b, c, d, x0, x1, y0, y1;
  bool bounded;

/* The strange indexes (i = -1, i = nrows - yl + 2, j = -1, j = ncols - xl + 2) 
   below fake the presence of ghost points around the edge of the grid. This 
   allows the program to fill lakes on the edge of the grid without the need to 
   allocate a temporary grid of size ny+2 x nx+2. 
*/

  for (i = -1; i < nrows - yl + 2; i++) {
    for (j = -1; j < ncols - xl + 2; j++) {

      bounded = true;

    // Check bounding box is dry. (Do not check corners.)

      for (n = i + 1; n < i + yl - 1; n++) {

      // Western boundary. 

        if (j != -1) {
          if ( !(grid[n][j] > 0.0 || (spvalpos && grid[n][j] == spval)) ) {
            bounded = false;
            break ;
          }
        }

      // Eastern boundary.

        if (j != ncols - xl + 1) { 
          if ( !(grid[n][j + xl - 1] > 0.0 || (spvalpos && grid[n][j + xl - 1] == spval)) ) {
            bounded = false;
            break ;
          }
        }
      }

      if (!bounded) continue ;

      for (m = j + 1; m < j + xl - 1; m++) {

      // Southern boundary. 

        if (i != -1) {
         if ( !(grid[i][m] > 0.0 || (spvalpos && grid[i][m] == spval)) ) {
           bounded = false;
           break ;
         }
        }

      // Northern boundary.
        
        if (i != nrows - yl + 1) {
          if ( !(grid[i + yl - 1][m] > 0.0 || (spvalpos && grid[i + yl - 1][m] == spval)) ) {
            bounded = false;
            break ;
          }
        }
      }

      if (!bounded) continue ;

    // Fill enclosing points with (locally-interpolated) positive values.  
    //printf("\nFilling lake of size %d at: i,j = %d, %d\n", (yl - 2)*(xl - 2), i, j);

      for (n = i + 1; n < i + yl - 1; n++) {
        for (m = j + 1; m < j + xl - 1; m++) {
          if (grid[n][m] > 0.0)
            continue ; // Dry point - no filling required. 
          (*nfilled)++;
          if (fillspval) {
            grid[n][m] = spval;
          } else {
          // Interpolate out (wet) interiour point using (dry) boundary points. 
            a = (j == -1) ? 0.0 : (float) grid[n][j];
            b = (j == ncols - xl + 1) ? 0.0 : (float) grid[n][j + xl - 1];
            c = (i == -1) ? 0.0 : (float) grid[i][m];
            d = (i == nrows - yl + 1) ? 0.0 : (float) grid[i + yl - 1][m];
            x0 = (j == -1 || a == spval) ? 0.0 : expf(-(m - j - 1.0)/2.0f);
            x1 = (j == ncols - xl + 1 || b == spval) ? 0.0 : expf(-(xl - 1.0 - m)/2.0f);
            y0 = (i == -1 || c == spval) ? 0.0 : expf(-(n - i - 1.0)/2.0f);
            y1 = (i == nrows - yl + 1 || d == spval) ? 0.0 : expf(-(yl - 1.0 - n)/2.0f);
            if (x0 + x1 + y0 + y1 == 0.0)
              interp = spval; 
            else
              interp = (a*x0 + b*x1 + c*y0 + d*y1)/(x0 + x1 + y0 + y1);
            grid[n][m] = interp;
          }
        }
      }
    }
  }
}

