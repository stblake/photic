
#include "plot.h"

//
//    PLOT
//

// PLOT GRID /grid name/ [optional] SOUNDINGS ...
// 		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/

// PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ ... 
//    [optional] THRESHOLD /min radiance/ /max radiance/ ... 
//		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/

// PLOT LOGSCATTER /band name 1/ /band name 2/ TYPE /water type/ ... 
//    COORDINATES /x1/ /y1/ /x2/ /y2/ ...

void run_plot() { 

	if (strcmp(parsed[1], "GRID") == 0) {
		run_plot_grid();
		return ;
	} else if (strcmp(parsed[1], "LINE") == 0) {
		run_plot_line();
		return ;
  } else if (strcmp(parsed[1], "LOGSCATTER") == 0) {
    run_plot_logscatter();
	} else {
		printf("\nERROR: unknown command '%s'\n", parsed[1]);
		return ;
	}
}



// PLOT LOGSCATTER /band name 1/ /band name 2/ ... 
//    COORDINATES /x1/ /y1/ /x2/ /y2/ ... [optional] Lsm /Lsm band 1/
//    /Lsm band 2/

void run_plot_logscatter() {

  int n, m, ii, jj, i, j, k, ps1, ps2, segpts, np, npoints, nlines, 
    ncoords, segsizes[LINE_MAX_POINTS], 
    ncolours = 9, colours[] = {2,6,8,4,3,11,12,13,5};
  float coords[LINE_MAX_POINTS], dist, x, y, *radiance1, *radiance2, del, 
    grad, rad, xa, xb, ya, yb, Lsm1, Lsm2, error, *seg1, *seg2,
    xmin, xmax, ymin, ymax, a, b, r;
  bool present, success;

  for (n=0; n < LINE_MAX_POINTS; n++)
    segsizes[n] = 0;

  if (parsed[2][0] == '\0') {
    printf("\nERROR: PLOT LOGSCATTER /band name 1/ /band name 2/ \
COORDINATES /x1/ /y1/ /x2/ /y2/ ... [optional] LSM /Lsm band 1/ /Lsm band 2/\n");
    return ;
  }

// Check first grid is present. 

  present = false;

  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[2]) == 0) {
      present = true;
      ps1 = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
   return ;
  }

// Check second grid is present. 

  present = false;
  
  for (n = 0; n < MAX_GRIDS; n++) {
    if (strcmp(grid_names[n], parsed[3]) == 0) {
      present = true;
      ps2 = n;
      break;
    }
  }

  if (! present) {
    printf("\nERROR: unknown grid '%s'.\n", parsed[3]);
   return ;
  }

// Read coordinates of lines. 

  if (strcmp(parsed[4], "COORDINATES") != 0) {
    printf("\nERROR: PLOT LOGSCATTER /band name 1/ /band name 2/ \
COORDINATES /x1/ /y1/ /x2/ /y2/ ... [optional] LSM /Lsm band 1/ /Lsm band 2/\n");
    return ;
  }

  n = 5;
  nlines = 0;
  ncoords = 0;

  while(parsed[n][0] != '\0' && strcmp(parsed[n], "LSM") != 0) {

    nlines++;

    if ((ncoords + 4) >= LINE_MAX_POINTS) {
      printf("\nERROR: too many lines.\n");
      return ;
    }

    if (parsed[n][0] == '\0' || strcmp(parsed[n], "LSM") == 0) { 
      break ;
    } else if (strcmp(parsed[n], "SAND") == 0) {
      coords[ncoords++] = sand_lon;
      coords[ncoords++] = sand_lat;
      if (strcmp(parsed[n + 1], "DEEP") == 0) {
        coords[ncoords++] = deep_lon;
        coords[ncoords++] = deep_lat;
        n += 2;
      } else {
        coords[ncoords++] = atof(parsed[n + 1]);
        coords[ncoords++] = atof(parsed[n + 2]);
        n += 3;
      }
    } else if (strcmp(parsed[n], "DEEP") == 0) {
      coords[ncoords++] = deep_lon;
      coords[ncoords++] = deep_lat;
      if (strcmp(parsed[n + 1], "SAND") == 0) {
        coords[ncoords++] = sand_lon;
        coords[ncoords++] = sand_lat;
        n += 2;
      } else {
        coords[ncoords++] = atof(parsed[n + 1]);
        coords[ncoords++] = atof(parsed[n + 2]);
        n += 3;  
      }
    } else {
      coords[ncoords++] = atof(parsed[n]);
      coords[ncoords++] = atof(parsed[n + 1]);
      coords[ncoords++] = atof(parsed[n + 2]);
      coords[ncoords++] = atof(parsed[n + 3]);
      n += 4;
    }
  }

// Check grids are on the same extents. 

  if (gridded_data[ps1].nrows != gridded_data[ps2].nrows || 
    gridded_data[ps1].ncols != gridded_data[ps2].ncols || 
    fabs(gridded_data[ps1].slat - gridded_data[ps2].slat) > extent_eps || 
    fabs(gridded_data[ps1].nlat - gridded_data[ps2].nlat) > extent_eps || 
    fabs(gridded_data[ps1].wlon - gridded_data[ps2].wlon) > extent_eps || 
    fabs(gridded_data[ps1].elon - gridded_data[ps2].elon) > extent_eps || 
    fabs(gridded_data[ps1].cellsize - gridded_data[ps2].cellsize) > extent_eps) {
    printf("\nERROR: input grids are not commensurate.\n");
    return ;
  }

// Compute Lsm if not specified by the user.

  if (strcmp(parsed[n], "LSM") == 0) {
    Lsm1 = atof(parsed[++n]);
    Lsm2 = atof(parsed[++n]);
  } else if (Lsm[ps1] != 0.0 && Lsm[ps2] != 0.0) {
    Lsm1 = Lsm[ps1];
    Lsm2 = Lsm[ps2];
  } else {
    Lsm1 = array_min2(
        gridded_data[ps1].array, 
        gridded_data[ps1].nrows, 
        gridded_data[ps1].ncols, 
        gridded_data[ps1].nodata_value);
    Lsm2 = array_min2(
        gridded_data[ps2].array, 
        gridded_data[ps2].nrows, 
        gridded_data[ps2].ncols, 
        gridded_data[ps2].nodata_value);
  }

//  printf("\nLsm1,Lsm2 = %f, %f\n", Lsm1, Lsm2);

// Extract points on each line on each grid. 

  ii = 0;
  npoints = 0;

  for (k = 0; k < nlines; k++) {
    xa = coords[ii++];
    ya = coords[ii++];
    xb = coords[ii++];
    yb = coords[ii++];
    dist = sqrt(pow(xb - xa,2) + pow(yb - ya,2));
    npoints += floor(dist/gridded_data[ps1].cellsize);
  }

  radiance1 = (float*) malloc(npoints*sizeof(float));
  radiance2 = (float*) malloc(npoints*sizeof(float));

  ii = 0;
  m = 0;

  for (k = 0; k < nlines; k++) {
    xa = coords[ii++];
    ya = coords[ii++];
    xb = coords[ii++];
    yb = coords[ii++];

    dist = sqrt(pow(xb - xa,2) + pow(yb - ya,2));
    segpts = floor(dist/gridded_data[ps1].cellsize);
    segsizes[k] = segpts;
    grad = (yb - ya)/(xb - xa);

    for (n = 0; n < segpts; n++) {

      if (m > npoints - 1)
        break ;

      // Compute the point on the line. 

      del = ((float)n)/((float) segpts);
      x = xa + del*(xb - xa);
      y = ya + grad*del*(xb - xa);

      // Compute grid coordinates of point. 

      i = round( ((float) gridded_data[ps1].nrows)*(y - gridded_data[ps1].slat)/
          (gridded_data[ps1].nlat - gridded_data[ps1].slat) );
      j = round( ((float) gridded_data[ps1].ncols)*(x - gridded_data[ps1].wlon)/
          (gridded_data[ps1].elon - gridded_data[ps1].wlon) );

      // Radiance at i,j for band 1. 

      rad = 0.0; 
      np = 0;

      for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
        for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
          if (gridded_data[ps1].array[i + ii][j + jj] != gridded_data[ps1].nodata_value) {
            rad += gridded_data[ps1].array[i + ii][j + jj];
            np++;
          }
        }
      }

      if (np == 0)
        radiance1[m] = 0.0;
      else
        radiance1[m] = rad/((float) np);

      // Radiance at i,j for band 2. 
      
      rad = 0.0; 
      np = 0;

      for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
        for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
          if (gridded_data[ps2].array[i + ii][j + jj] != gridded_data[ps2].nodata_value) {
            rad += gridded_data[ps2].array[i + ii][j + jj];
            np++;
          }
        }
      }

      if (np == 0)
        radiance2[m] = 0.0;
      else
        radiance2[m] = rad/((float) np);

      m++;
    }
  }

// Load PGPLOT. 
#if PGPLOT
  error = cpgopen("?");
  if (error < 1) {
    printf("\nERROR: cannot load graphics library.\n");
    return ;
  }

// Set the page size and aspect ratio to 1.0. 

  cpgpap(pagesize, 1.0);

// Background and text colours. 

  if (background == BLACK) {
    cpgscr(0, 0.0, 0.0, 0.0);
    cpgscr(1, 1.0, 1.0, 1.0);
  } else {
    cpgscr(0, 1.0, 1.0, 1.0);
    cpgscr(1, 0.0, 0.0, 0.0);
  }

// Plot the log-transformed scatter plot of radiances. 

  for (n = 0; n < npoints; n++) {
    if (radiance1[n] < Lsm1)
      radiance1[n] = 0.0; 
    else
      radiance1[n] = log(radiance1[n] - Lsm1);

    if (radiance2[n] < Lsm2)
      radiance2[n] = 0.0; 
    else
      radiance2[n] = log(radiance2[n] - Lsm2);
  }

  int offset = 0;

  xmin = vec_min(radiance1, npoints);
  xmax = vec_max(radiance1, npoints);
  ymin = vec_min(radiance2, npoints);
  ymax = vec_max(radiance2, npoints);


  seg1 = (float*) malloc(npoints*sizeof(float));
  seg2 = (float*) malloc(npoints*sizeof(float));

  bool first = true;
  char xlabel[128], ylabel[128];

  strcpy(xlabel, "log(Ls\\d");
  strcat(xlabel, parsed[2]);
  strcat(xlabel, "\\u - Lsm\\d");
  strcat(xlabel, parsed[2]);
  strcat(xlabel, "\\u)");

  strcpy(ylabel, "log(Ls\\d");
  strcat(ylabel, parsed[3]);
  strcat(ylabel, "\\u - Lsm\\d");
  strcat(ylabel, parsed[3]);
  strcat(ylabel, "\\u)");

  for (k = 0; k < nlines; k++) {

    for (n = 0; n < segsizes[k]; n++) {
      seg1[n] = radiance1[offset + n];
      seg2[n] = radiance2[offset + n];
    }

    plot_scatter("LOG-LOG SCATTER PLOT", xlabel, ylabel, 
      colours[k%ncolours], 16, 
      segsizes[k], seg1, seg2, line_width,
      xmin, xmax, ymin, ymax, first);

    offset += segsizes[k];
    first = false;
  }

  if (nlines == 1) {

    // Compute the lines of best fit.

    success = linear_fit(radiance1, radiance2, segsizes[0], &a, &b, &r);

    printf("\nm, b, r = %.3f, %.3f, %.3f\n", a, b, r);

    // Plot the lines of best fit.
    cpgsls(2);
    cpgsci(3);
    cpgmove(xmin, a*xmin + b);
    cpgdraw(xmax, a*xmax + b);
    cpgsls(1);
  }

// Close the graphics device. 

  cpgclos();

  free(seg1);
  free(seg2);
#endif

  free(radiance1);
  free(radiance2);
}



// PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ ... 
//    [optional] THRESHOLD /min radiance/ /max radiance/ ... 
//		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/

void run_plot_line() {

	int n, ps, error, palette, i, j, np, npoints, ii, jj;
	bool present = false;
	float x[2], y[2], xa, ya, xb, yb, xmin = 0.0, xmax = 0.0, 
	  ymin = 0.0, ymax = 0.0, *xpts, *ypts, *radiance, grad, 
	  rad, del, min_threshold = 0.0, max_threshold = 0.0, 
	  max_slope, slope, radmin, radmax, noise;
	char title[256], snum[32];

// Read the grid. 

  if (parsed[2][0] == '\0') {
  	printf("\nERROR: PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ \
[optional] THRESHOLD /min radiance/ /max radiance/ REGION /xmin/ /ymin/ \
/xmax/ /ymax/ PALETTE /palette name/\n");
  	return ;
  } else {

  	// Check grid is present. 

		for (n = 0; n < MAX_GRIDS; n++) {
			if (strcmp(grid_names[n], parsed[2]) == 0) {
				present = true;
				ps = n;
				break;
			}
		}

		if (! present) {
			printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
			return ;
		}
  }

// Read the line coordinates. 

  if (strcmp(parsed[3], "COORDINATES") != 0) {
  	printf("\nERROR: PLOT LINE /grid name/ COORDINATES /xa/ /ya/ /xb/ /yb/ \
[optional] THRESHOLD /min radiance/ /max radiance/ REGION /xmin/ /ymin/ \
/xmax/ /ymax/ PALETTE /palette name/\n");
  	return ;
  } else {
    if (strcmp(parsed[4], "SAND") == 0) {
      xa = sand_lon;
      ya = sand_lat;
      if (strcmp(parsed[5], "DEEP") == 0) {
        xb = deep_lon;
        yb = deep_lat;
        n = 6;
      } else {
        xb = atof(parsed[5]);
        yb = atof(parsed[6]);
        n = 7;       
      }
    } else if (strcmp(parsed[4], "DEEP") == 0) {
      xa = deep_lon;
      ya = deep_lat;
      if (strcmp(parsed[5], "SAND") == 0) {
        xb = sand_lon;
        yb = sand_lat;
        n = 6;
      } else {
        xb = atof(parsed[5]);
        yb = atof(parsed[6]);
        n = 7;       
      }
    } else {
  	  xa = atof(parsed[4]);
  	  ya = atof(parsed[5]);
  	  xb = atof(parsed[6]);
  	  yb = atof(parsed[7]);
      n = 8;
    }
  }

  // printf("\nxa,ya = %f, %f", xa, ya);
  // printf("\nxb,yb = %f, %f\n", xb, yb);

// Read optional arguments. 

  if (strcmp(parsed[n], "THRESHOLD") == 0) {
  	min_threshold = atof(parsed[++n]);
  	max_threshold = atof(parsed[++n]);
  	n++;
  }

  if (strcmp(parsed[n], "REGION") == 0) {
  	xmin = atof(parsed[++n]);
  	ymin = atof(parsed[++n]);
  	xmax = atof(parsed[++n]);
  	ymax = atof(parsed[++n]);
  	n++;
  }

  palette = GRAYSCALE; // Default.

  if (strcmp(parsed[n], "PALETTE") == 0) {
    if (strcmp(parsed[++n], "jet") == 0)
      palette = JET;
    else if (strcmp(parsed[n], "inverse_jet") == 0)
      palette = INVERSEJET;
    else if (strcmp(parsed[n], "grayscale") == 0)
      palette = GRAYSCALE;
    else if (strcmp(parsed[n], "ocean") == 0)
      palette = OCEAN;
    else if (strcmp(parsed[n], "gmt") == 0)
      palette = GMT;
    else if (strcmp(parsed[n], "rygb") == 0)
      palette = RYGB;
    else if (strcmp(parsed[n], "wor") == 0)
      palette = WOR;
    else if (strcmp(parsed[n], "viridis") == 0)
      palette = VIRIDIS;
    else if (strcmp(parsed[n], "plasma") == 0)
      palette = PLASMA;
    else {
      printf("\nERROR: unknown colour palette '%s'\n", parsed[n]);
    }
  }

// Extract points on the line from the grid. 

  npoints = 1 + floor(sqrt(pow(xb - xa, 2) + pow(yb - ya, 2))/gridded_data[ps].cellsize);

  xpts = (float*) malloc(npoints*sizeof(float));
  ypts = (float*) malloc(npoints*sizeof(float));
  radiance = (float*) malloc(npoints*sizeof(float));

// Gradient of line. 

  grad = (yb - ya)/(xb - xa); // TODO: check for vertical line.

  for (n = 0; n < npoints; n++) {
  	// Compute the point on the line. 
  	del = ((float)n)/((float) npoints);
  	xpts[n] = xa + del*(xb - xa);
  	ypts[n] = ya + grad*del*(xb - xa);
  	// Compute grid coordinates of point. 
  	// i = (ypts[n] - gridded_data[ps].slat)/gridded_data[ps].cellsize;
  	// j = (xpts[n] - gridded_data[ps].wlon)/gridded_data[ps].cellsize;
  	i = round( ((float) gridded_data[ps].nrows)*(ypts[n] - gridded_data[ps].slat)/
          (gridded_data[ps].nlat - gridded_data[ps].slat) );
  	j = round( ((float) gridded_data[ps].ncols)*(xpts[n] - gridded_data[ps].wlon)/
          (gridded_data[ps].elon - gridded_data[ps].wlon) );
  	// Radiance at i,j.
  	rad = 0.0; 
  	np = 0;
  	for (ii = 1 - depth_radius; ii < depth_radius; ii++) {
  		for (jj = 1 - depth_radius; jj < depth_radius; jj++) {
  			if (gridded_data[ps].array[i + ii][j + jj] != gridded_data[ps].nodata_value) {
  				rad += gridded_data[ps].array[i + ii][j + jj];
  				np++;
  			}
  		}
		}
		if (np == 0)
			radiance[n] = 0.0;
		else
			radiance[n] = rad/((float) np);
  }

// Estimate deep water radiance. 

  radmin = vec_min(radiance, npoints);
  radmax = vec_max(radiance, npoints);

  if (min_threshold == 0.0)
	  min_threshold = radmin + noise_threshold*(radmax - radmin);
  
  min_radiance_threshold = min_threshold; // Global.

// Estimate land-sea radiance. 
  
  max_slope = 0.0;

  if (max_threshold == 0.0) {
  	for (n = 1; n < npoints - 1; n++) {
   		slope = radiance[n + 1] - radiance[n]; 
			if (fabs(slope) > max_slope) {
				max_slope = fabs(slope);
				if (radiance[n + 1] > radiance[n])
					max_threshold = radiance[n];
				else 
					max_threshold = radiance[n + 1];
			}
  	}
  }

  max_radiance_threshold = max_threshold; // Global.

#if PGPLOT
// Load PGPLOT. 

  error = cpgopen("?");
  if (error < 1) {
  	printf("\nERROR: cannot load graphics library.\n");
  	return ;
  }

// Background and text colours. 

	if (background == BLACK) {
  	cpgscr(0, 0.0, 0.0, 0.0);
  	cpgscr(1, 1.0, 1.0, 1.0);
  } else {
  	cpgscr(0, 1.0, 1.0, 1.0);
  	cpgscr(1, 0.0, 0.0, 0.0);
  }

// Create two (horizontally-alligned) graphics windows. 

  cpgsubp(2,1);

// Create pixel-array plot.

  x[0] = xa;
  x[1] = xb; 
  y[0] = ya;
  y[1] = yb;

  plot_array(gridded_data[ps], parsed[2], pagesize, textsize, pointsize,
  	xmin, ymin, xmax, ymax, palette, 2, x, y, 0.5, line_width, 
    false, gridded_data[ps], 0.0, 0.0);

// Plot land/sea intersection on pixel-array plot. 

  cpgsci(8);

  for (n = 0; n < npoints - 1; n++) {
  	if ((radiance[n] <= max_threshold && radiance[n + 1] >= max_threshold) || 
  		(radiance[n] >= max_threshold && radiance[n + 1] <= max_threshold) ) {
  		cpgpt1(xpts[n], ypts[n], 20);
  	}
  }

// Plot deep water threshold on pixel-array plot. 

  cpgsci(4);

  for (n = 0; n < npoints - 1; n++) {
  	if ((radiance[n] <= min_threshold && radiance[n + 1] >= min_threshold) || 
  		(radiance[n] >= min_threshold && radiance[n + 1] <= min_threshold) ) {
  		cpgpt1(xpts[n], ypts[n], 20);
  	}
  }


// Reset colour. 

  cpgsci(1);

// Create scatter plot. 

  strcpy(title, "RADIANCE MIN/MAX THRESHOLDS = ");
  sprintf(snum, "%5.2f, %5.2f", min_threshold, max_threshold);
  strcat(title, snum);

  plot_scatter(title, "", "", 3, 20, npoints, xpts, radiance, line_width, 0.0, 0.0, 0.0, 0.0, true);

// Plot the radiance thresholds on the scatter plot.

  cpgsci(4);
  cpgsls(2);
  cpgmove(xpts[0], min_threshold);
  cpgdraw(xpts[npoints - 1], min_threshold);
  cpgsls(1);

  cpgsci(8);
  cpgsls(2);
  cpgmove(xpts[0], max_threshold);
  cpgdraw(xpts[npoints - 1], max_threshold);
  cpgsls(1);

// Close the graphics device. 

  cpgclos();
#endif

// Free memory.

  free(xpts);
  free(ypts);
  free(radiance);

}



// PLOT GRID /grid name/ [optional] SOUNDINGS ...
// 		REGION /xmin/ /ymin/ /xmax/ /ymax/ PALETTE /palette name/ ...
//    LAND /land grid name/ MIN /min/ MAX /max/

void run_plot_grid() {
	
	int i, k, n, ps, ps_land = 0, error;
	bool present = false, plot_soundings = false, plot_land = false;
	float xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, res = 0.0,
		*x, *y, min = 0.0, max = 0.0;

// Read the grid name. 

  if (parsed[2][0] == '\0') {
  	printf("\nERROR: PLOT /grid name/ [optional] REGION /xmin/ /ymin/ /xmax/ /ymax/ \
PALETTE /palette name/ LAND /land grid name/\n");
  	return ;
  } else {

  	// Check grid is present. 

		for (n = 0; n < MAX_GRIDS; n++) {
			if (strcmp(grid_names[n], parsed[2]) == 0) {
				present = true;
				ps = n;
				break;
			}
		}

		if (! present) {
			printf("\nERROR: unknown grid '%s'.\n", parsed[2]);
			return ;
		}
  }

  n = 3;

// Read SOUNDINGS

  if (strcmp(parsed[n], "SOUNDINGS") == 0) {
  	plot_soundings = true;
  	n++;
  	x = (float*) malloc(nsoundings*sizeof(float));
  	y = (float*) malloc(nsoundings*sizeof(float));
  	for (i = 0; i < nsoundings; i++) {
  		x[i] = soundings[i][0];
  		y[i] = soundings[i][1];
  	}
  }

// Read REGION

  if (strcmp(parsed[n], "REGION") == 0) {
  	xmin = atof(parsed[++n]);
  	ymin = atof(parsed[++n]);
  	xmax = atof(parsed[++n]);
  	ymax = atof(parsed[++n]);
  	n++;
  } 

// Read PALETTE 

  palette = JET; // Default.

  if (strcmp(parsed[n], "PALETTE") == 0) {
    if (strcmp(parsed[++n], "jet") == 0)
      palette = JET;
    else if (strcmp(parsed[n], "inverse_jet") == 0)
      palette = INVERSEJET;
    else if (strcmp(parsed[n], "grayscale") == 0)
      palette = GRAYSCALE;
    else if (strcmp(parsed[n], "ocean") == 0)
      palette = OCEAN;
    else if (strcmp(parsed[n], "gmt") == 0)
      palette = GMT;
    else if (strcmp(parsed[n], "rygb") == 0)
      palette = RYGB;
    else if (strcmp(parsed[n], "wor") == 0)
      palette = WOR;
    else if (strcmp(parsed[n], "viridis") == 0)
      palette = VIRIDIS;
    else if (strcmp(parsed[n], "plasma") == 0)
      palette = PLASMA;
    else {
      printf("\nERROR: unknown colour palette '%s'\n", parsed[n]);
    }
    n++;
  }

// Read LAND

  if (strcmp(parsed[n], "LAND") == 0) {
    plot_land = true;

    // Check grid is present. 

    present = false;
    n++;
    for (k = 0; k < MAX_GRIDS; k++) {
      if (strcmp(grid_names[k], parsed[n]) == 0) {
        present = true;
        ps_land = k;
        break;
      }
    }

    if (! present) {
      printf("\nERROR: unknown grid '%s'.\n", parsed[n]);
      return ;
    }
    n++;
  }

// Read MIN and MAX.

  if (strcmp(parsed[n], "MIN") == 0) {
    min = atof(parsed[++n]);
    n++;
  }

  if (strcmp(parsed[n], "MAX") == 0) {
    max = atof(parsed[++n]);
    n++;
  }

#if PGPLOT
// Load PGPLOT. 

  error = cpgopen("?");
  if (error < 1) {
  	printf("\nERROR: cannot load graphics library.\n");
  	return ;
  }

// Background and text colours. 

	if (background == BLACK) {
  		cpgscr(0, 0.0, 0.0, 0.0);
  		cpgscr(1, 1.0, 1.0, 1.0);
  	} else {
  		cpgscr(0, 1.0, 1.0, 1.0);
  		cpgscr(1, 0.0, 0.0, 0.0);
  	}

// Create plot.

  if (plot_soundings) {
  	plot_array(gridded_data[ps], parsed[2], pagesize, textsize, pointsize,
  		xmin, ymin, xmax, ymax, palette, nsoundings, x, y, 1.0, line_width, 
      plot_land, gridded_data[ps_land], min, max);
  	free(x);
  	free(y);
  } else {
  	plot_array(gridded_data[ps], parsed[2], pagesize, textsize, pointsize,
  		xmin, ymin, xmax, ymax, palette, 0, NULL, NULL, 1.0, line_width,
      plot_land, gridded_data[ps_land], min, max);
  }

// Close the graphics device. 

  cpgclos();
#endif

}

