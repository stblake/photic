
#include "write.h"


//
//    WRITE
//

// WRITE GRID /grid name/ /output file name/ REGION /xmin/ /ymin/ /xmax/ /ymax/

// WRITE LINE /grid name/ /output file name/ COORDINATES /xa/ /ya/ /xb/ /yb/

// WRITE POINTS /grid name/ /output file name/

void run_write() {

  if (strcmp(parsed[1], "GRID") == 0) {
    run_write_grid();
  } else if (strcmp(parsed[1], "LINE") == 0) {
    run_write_line();
  } else if (strcmp(parsed[1], "POINTS") == 0) {
    run_write_points();
  } else {
    printf("\nERROR: unknown command '%s'\n", parsed[1]);
  }
}


// WRITE POINTS /grid name/ /output file name/

void run_write_points() {

  int n, ps;
  float *xpts, *ypts, *zpts;
  bool present = false; 

  // Get position of grid name.
      
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

  // Check soundings are already in memory. 

  if (! read_soundings) {
    printf("\nERROR: soundings have not been read into memory.\n");
    return ;
  }

  // Interpolate grid at soundings. 

  xpts = malloc(nsoundings*sizeof(float));
  ypts = malloc(nsoundings*sizeof(float));
  zpts = malloc(nsoundings*sizeof(float));

  for (n = 0; n < nsoundings; n++) {
    xpts[n] = soundings[n][0];
    ypts[n] = soundings[n][1];

    if (xpts[n] > gridded_data[ps].wlon && 
    xpts[n] < gridded_data[ps].elon && 
    ypts[n] > gridded_data[ps].slat && 
    ypts[n] < gridded_data[ps].nlat) {
      zpts[n] = interp_bicubic(gridded_data[ps].array, 
              gridded_data[ps].ncols,
              gridded_data[ps].nrows, 
              gridded_data[ps].wlon, 
              gridded_data[ps].slat, 
              gridded_data[ps].cellsize, 
              gridded_data[ps].nodata_value, 
              xpts[n], 
              ypts[n]);
    } else {
      zpts[n] = gridded_data[ps].nodata_value;
    }
  }

  // Write points to file. 

  write_xyz(parsed[3], xpts, ypts, zpts, nsoundings);

  free(xpts);
  free(ypts);
  free(zpts);
}


// WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/ FILE //

void run_write_line() {

  int n, ps, npoints, i, j, ii, jj, np;
  float xa, ya, xb, yb, rad, grad, del, *xpts, *ypts, *radiance;
  bool present = false;

  if (parsed[2][0] == '\0'){
    printf("\nERROR: WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

  // Get position of grid name.
      
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

  // Read output file name. 

  if (parsed[3][0] == '\0'){
    printf("\nERROR: WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

  // Check output grid name is a netCDF (.nc or .nc4 extension)

  if (strstr(parsed[3], ".csv") == NULL && strstr(parsed[3], ".xyz") == NULL) { // Crude check.
    printf("\nERROR: output file should have the extension .csv or .xyz\n");
    return ;
  }

  // Read LINE.

  if (strcmp(parsed[4], "COORDINATES") != 0) {
    printf("\nERROR: WRITE LINE /grid name/ /output line/ COORDINATES /xa/ /ya/ /xb/ /yb/\n");
    return ;
  }

  if (strcmp(parsed[5], "SAND") == 0) {
    xa = sand_lon;
    ya = sand_lat;
    if (strcmp(parsed[6], "DEEP") == 0) {
      xb = deep_lon;
      yb = deep_lat;
      n = 7;
    } else {
      xb = atof(parsed[6]);
      yb = atof(parsed[7]);
      n = 8;
    }
  } else if (strcmp(parsed[5], "DEEP") == 0) {
    xa = deep_lon;
    ya = deep_lat;
    if (strcmp(parsed[6], "SAND") == 0) {
      xb = sand_lon;
      yb = sand_lat;
      n = 7;
    } else {
      xb = atof(parsed[6]);
      yb = atof(parsed[7]);
      n = 8;
    }
  } else {
    xa = atof(parsed[5]);
    ya = atof(parsed[6]);
    xb = atof(parsed[7]);
    yb = atof(parsed[8]);
    n = 9;
  }

// Extract points on the line from the grid. 

  float dist = sqrt(pow(xb - xa, 2) + pow(yb - ya, 2));
  npoints = 1 + floor(dist/gridded_data[ps].cellsize);

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

//  Write line to xyz file. 

  write_xyz(parsed[3], xpts, ypts, radiance, npoints);

// Free memory.

  free(xpts);
  free(ypts);
  free(radiance);
}



// WRITE GRID /grid name/ /output file name/ REGION /xmin/ /ymin/ /xmax/ /ymax/

void run_write_grid() {

  int n, ps;
  float *lats, *lons, xmin, ymin, xmax, ymax;
  bool present = false, region = false;

  if (parsed[2][0] == '\0'){
    printf("\nERROR: WRITE GRID /grid name/ /output file name/\n");
    return ;
  }
  	
  // Get position of grid name.
  	  
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

  // Read output file name. 

  if (parsed[3][0] == '\0'){
    printf("\nERROR: WRITE GRID /grid name/ /output file name/\n");
    return ;
  }

  // Check output grid name is a netCDF (.nc or .nc4 extension)

  if (strstr(parsed[3], ".nc") == NULL && strstr(parsed[3], ".asc") == NULL) { // Crude check.
    printf("\nERROR: output file should have the extension .nc, .nc4 or .asc\n");
    return ;
  }

  if (strcmp(parsed[4], "REGION") == 0) {
  	region = true;
  	xmin = atof(parsed[5]);
  	ymin = atof(parsed[6]);
  	xmax = atof(parsed[7]);
  	ymax = atof(parsed[8]);
  }

  // Write output grid. 

  if (region) {
  	int ncols = (xmax - xmin)/gridded_data[ps].cellsize;
  	int nrows = (ymax - ymin)/gridded_data[ps].cellsize;
  	lons = (float*) malloc(ncols*sizeof(float));
  	lats = (float*) malloc(nrows*sizeof(float));

  	for (n = 0; n < ncols; n++)
  		lons[n] = xmin + ((float) n)*gridded_data[ps].cellsize;

  	for (n = 0; n < nrows; n++)
  		lats[n] = ymin + ((float) n)*gridded_data[ps].cellsize;
  } else {
  	lons = (float*) malloc(gridded_data[ps].ncols*sizeof(float));
  	lats = (float*) malloc(gridded_data[ps].nrows*sizeof(float));

  	for (n = 0; n < gridded_data[ps].ncols; n++)
  		lons[n] = gridded_data[ps].wlon + ((float) n)*gridded_data[ps].cellsize;

  	for (n = 0; n < gridded_data[ps].nrows; n++)
  		lats[n] = gridded_data[ps].slat + ((float) n)*gridded_data[ps].cellsize;

  	write_nc(parsed[3], 
  		gridded_data[ps].array, 
  		gridded_data[ps].ncols,
  		gridded_data[ps].nrows,
  		lons, lats, 
  		gridded_data[ps].nodata_value);
  }

  free(lons);
  free(lats);
}


