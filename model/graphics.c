
#include "graphics.h"


// 
//  2D histogram. 
// 

void plot_bidimensional_histogram(char *title, char *xlabel, char *ylabel,  
  int npoints, float *xpts, float *ypts, int line_width,
  float xmin, float xmax, float ymin, float ymax, bool first, int palette, 
  int nxbins, int nybins) {

#if PGPLOT

  int i, j, k, cindex;
  float diff, buffer, **density, den_max, dx, dy, lon0, lat0, lon1, lat1, 
    mean, scaled_max;

// Load colour palette. 

  if (palette == JET)
    load_jet_palette();
  else if (palette == INVERSEJET)
    load_inverse_jet_palette();
  else if (palette == OCEAN)
    load_esri_ocean_palette();
  else if (palette == GMT)
    load_gmt_palette();
  else if (palette == RYGB)
    load_rygb_palette();
  else if (palette == WOR)
    load_white_orange_red();
  else if (palette == VIRIDIS)
    load_matplotlib_viridis();
  else if (palette == PLASMA)
    load_matplotlib_plasma();
  else
    load_grayscale_palette();

// Compute plotting region. 

  if (xmin == 0.0 && xmax == 0.0 && ymin == 0.0 && ymax == 0.0) {
    xmin = vec_min(xpts, npoints);
    xmax = vec_max(xpts, npoints);
    ymin = vec_min(ypts, npoints);
    ymax = vec_max(ypts, npoints);
  }

  buffer = 1.05; // 5% buffer

  diff = xmax - xmin;
  xmin -= (buffer - 1.0)*diff;
  xmax += (buffer - 1.0)*diff;

  diff = ymax - ymin;
  ymin -= (buffer - 1.0)*diff;
  ymax += (buffer - 1.0)*diff;

// Set size of plotting window. 

  if (first)
    cpgenv(xmin, xmax, ymin, ymax, 0, 0);

// Set line width. 

  cpgslw(line_width);

// Title. 

  for (i = 0; i < strlen(title); i++)
    if (title[i] == '_') title[i] = ' ';

  cpglab(xlabel, ylabel, trim(title));

// Compute 2D density. 

  allocate_float_array_2d(&density, nybins, nxbins);

  for (i = 0; i < nybins; i++) {
    for (j = 0; j < nxbins; j++) {
      density[i][j] = 0.0;
    }
  }

  for (k = 0; k < npoints; k++) {
    j = floor(((float) nxbins - 1)*(xpts[k] - xmin)/(xmax - xmin));
    i = floor(((float) nybins - 1)*(ypts[k] - ymin)/(ymax - ymin));
    density[i][j] += 1.0;
  }

  den_max = array_max(density, nybins, nxbins);

#if 1
  for (i = 0; i < nybins; i++) {
    for (j = 0; j < nxbins; j++) {
      density[i][j] = sqrt(density[i][j]);
    }
  }
#endif

  mean = array_mean(density, nybins, nxbins, 0.0);
  scaled_max = mean + 3.0*array_stddev(density, nybins, nxbins, 0.0, mean);

// Draw colour ramp. 

  if (palette == GMT)
    cpgscir(17, 36);
  else if (palette == RYGB)
    cpgscir(17, 83);
  else // all other palettes have 64 colours
    cpgscir(17, 80);

  cpgwedg("RI", 0.5, 3.0, 0.0, den_max, " ");

#if 0
  printf("\n\n");
  for (i = 0; i < nybins; i++) {
    for (j = 0; j < nxbins; j++) {
      printf("%6.0f  ", density[i][j]);
    }
    printf("\n");
  }
  printf("\n");
#endif

// Plot points. 

  cpgbbuf();

  dx = (xmax - xmin)/((float) nxbins);
  dy = (ymax - ymin)/((float) nybins);

  for (i = 0; i < nybins; i++) {
    for (j = 0; j < nxbins; j++) {

      if (density[i][j] < 1.0) { // 3.0
        continue ;
      }

      cindex = (int) 63.0*MIN(1.0, density[i][j]/scaled_max);
      cpgsci(cindex + 17);
      lon0 = xmin + ((float) j)*dx;
      lon1 = xmin + ((float) j + 1)*dx;
      lat0 = ymin + ((float) i)*dy;
      lat1 = ymin + ((float) i + 1)*dy;
      cpgrect(lon0, lon1, lat0, lat1);
    }
  }

  cpgebuf();
  cpgsci(1);

// Clean up memory. 

  free_float_array_2d(density, nybins);

#endif
}


//
//  2D log-scaled histogram. 
//

void plot_log_bidimensional_histogram(char *title, char *xlabel, char *ylabel,  
  int npoints, float *xpts, float *ypts, float lsmx, float lsmy, int line_width,
  float xmin, float xmax, float ymin, float ymax, bool first, int palette, 
  int nxbins, int nybins) {

#if PGPLOT

  int i, j, k, cindex;
  float diff, buffer, xlogmin, xlogmax, ylogmin, ylogmax, **log_density, 
    den_max, dx, dy, lat0, lon0, lat1, lon1;

// Load colour palette. 

  if (palette == JET)
    load_jet_palette();
  else if (palette == INVERSEJET)
    load_inverse_jet_palette();
  else if (palette == OCEAN)
    load_esri_ocean_palette();
  else if (palette == GMT)
    load_gmt_palette();
  else if (palette == RYGB)
    load_rygb_palette();
  else if (palette == WOR)
    load_white_orange_red();
  else if (palette == VIRIDIS)
    load_matplotlib_viridis();
  else if (palette == PLASMA)
    load_matplotlib_plasma();
  else
    load_grayscale_palette();

// Compute plotting region. 

  buffer = 1.05; // 5% buffer

  xlogmin = log(xmin - lsmx);
  if (xlogmin < 0.0) {
    xlogmin = 0.0;
  }
  xlogmax = log(xmax - lsmx);
  ylogmin = log(ymin - lsmy);
  if (ylogmin < 0.0) {
    ylogmin = 0.0;
  }
  ylogmax = log(ymax - lsmy);

  diff = xlogmax - xlogmin;
  xlogmin -= (buffer - 1.0)*diff;
  xlogmax += (buffer - 1.0)*diff;

  diff = ylogmax - ylogmin;
  ylogmin -= (buffer - 1.0)*diff;
  ylogmax += (buffer - 1.0)*diff;

// Set size of plotting window. 

  if (first)
    cpgenv(xlogmin, xlogmax, ylogmin, ylogmax, 0, 0);

// Set line width. 

  cpgslw(line_width);

// Title. 

  for (i = 0; i < strlen(title); i++)
    if (title[i] == '_') title[i] = ' ';

  cpglab(xlabel, ylabel, trim(title));

// Create the 2D log-scaled histogram.

  allocate_float_array_2d(&log_density, nybins, nxbins);

  for (i = 0; i < nybins; i++) {
    for (j = 0; j < nxbins; j++) {
      log_density[i][j] = 0.0;
    }
  }

  for (k = 0; k < npoints; k++) {

    j = floor(((float) nxbins - 1)*(xpts[k] - xmin)/(xmax - xmin));
    i = floor(((float) nybins - 1)*(ypts[k] - ymin)/(ymax - ymin));
    log_density[i][j] += 1.0;
  }

  den_max = log(array_max(log_density, nybins, nxbins));

// Plot points. 

  cpgbbuf();

  dx = (xmax - xmin)/((float) nxbins);
  dy = (ymax - ymin)/((float) nybins);

  for (i = 0; i < nybins; i++) {
    for (j = 0; j < nxbins; j++) {
      if (log_density[i][j] == 0.0) {
        continue ;
      }
      cindex = 64 - (int) 63.0*log(log_density[i][j])/den_max;
      cpgsci(cindex + 17);
      lon0 = log(xmin + ((float) j)*dx - lsmx);
      lon1 = log(xmin + ((float) j + 1)*dx - lsmx);
      lat0 = log(ymin + ((float) i)*dy - lsmy);
      lat1 = log(ymin + ((float) i + 1)*dy - lsmy);
      cpgrect(lon0, lon1, lat0, lat1);
    }
  }

  cpgebuf();
  cpgsci(1);

// Clean up memory. 

  free_float_array_2d(log_density, nybins);

#endif
}

//
// Scatter plot
// 

void plot_scatter(char *title, char *xlabel, char *ylabel, int colour, int marker, 
  int npoints, float *xpts, float *ypts, int line_width,
  float xmin, float xmax, float ymin, float ymax, bool first) {

#if PGPLOT

	int i;
	float diff, buffer;

// Compute plotting region. 

  if (xmin == 0.0 && xmax == 0.0 && ymin == 0.0 && ymax == 0.0) {
	  xmin = vec_min(xpts, npoints);
	  xmax = vec_max(xpts, npoints);
	  ymin = vec_min(ypts, npoints);
	  ymax = vec_max(ypts, npoints);
  }

  buffer = 1.05; // 5% buffer

  diff = xmax - xmin;
	xmin -= (buffer - 1.0)*diff;
	xmax += (buffer - 1.0)*diff;

  diff = ymax - ymin;
	ymin -= (buffer - 1.0)*diff;
	ymax += (buffer - 1.0)*diff;

// Set size of plotting window. 

  if (first)
    cpgenv(xmin, xmax, ymin, ymax, 0, 0);

// Set line width. 

  cpgslw(line_width);

// Title. 

  for (i = 0; i < strlen(title); i++)
    if (title[i] == '_') title[i] = ' ';

  cpglab(xlabel, ylabel, trim(title));

// Set colour. 

  cpgsci(colour);

// Plot points. 

  cpgbbuf();

  for (i = 0; i < npoints; i++) {
  	cpgpt1(xpts[i], ypts[i], marker);
  }

  cpgebuf();

  cpgsci(1);

#endif
}


//
// Plotting a pixel-based array
//

void plot_array(geogrid gr, char *title, float pagesize, float textsize, float pointsize,
	float xmin, float ymin, float xmax, float ymax, int palette, 
	int npoints, float *xpts, float *ypts, float aspect_factor, int line_width,
  bool plot_land, geogrid land, float user_min, float user_max) {

#if PGPLOT

  int i, j, cindex;
  float lond, latd, aspect, min, max, lat0, lon0, lat1, lon1;
  bool box = false;

// Load colour palette. 

  load_noaa_land_colour();

  if (palette == JET)
    load_jet_palette();
  else if (palette == INVERSEJET)
    load_inverse_jet_palette();
  else if (palette == OCEAN)
    load_esri_ocean_palette();
  else if (palette == GMT)
    load_gmt_palette();
  else if (palette == RYGB)
    load_rygb_palette();
  else if (palette == WOR)
    load_white_orange_red();
  else if (palette == VIRIDIS)
    load_matplotlib_viridis();
  else if (palette == PLASMA)
    load_matplotlib_plasma();
  else
    load_grayscale_palette();
  

// Compute the aspect ratio. 

  if (xmin == 0.0 && xmax == 0.0 && ymin == 0.0 && ymax == 0.0) {
  	lond = gr.elon - gr.wlon;
  	latd = gr.nlat - gr.slat;
  } else {
  	lond = xmax - xmin;
  	latd = ymax - ymin;
  }

  aspect = latd/lond;

// Compute array min and max. 

  if (user_min == user_max) {
    min = array_min2(gr.array, gr.nrows, gr.ncols, gr.nodata_value);
    max = array_max2(gr.array, gr.nrows, gr.ncols, gr.nodata_value);
  } else {
    min = user_min; 
    max = user_max;
  }

  // printf("\nmin,max = %f, %f\n", min, max);

// Set the page size. 

  cpgpap(pagesize, aspect*aspect_factor);

// Set line width. 

  cpgslw(line_width);

// Set the text size.

  cpgsch(textsize);
  
// Set size of plotting window. 

  if (xmin == 0.0 && xmax == 0.0 && ymin == 0.0 && ymax == 0.0) {
  	cpgenv(gr.wlon, gr.elon, gr.slat, gr.nlat, 0, 0);
  } else {
  	box = true;
  	cpgenv(xmin, xmax, ymin, ymax, 0, 0);
  }

// Label axes. 

  for (i = 0; i < strlen(title); i++)
    if (title[i] == '_') title[i] = ' ';

  cpglab("photic", " ", trim(title));

// Draw colour ramp. 

  if (palette == GMT)
    cpgscir(17, 36);
  else if (palette == RYGB)
    cpgscir(17, 83);
  else // all other palettes have 64 colours
    cpgscir(17, 80);

  cpgwedg("RI", 0.5, 3.0, min, max, " ");

// Draw array plot. 

  cpgbbuf();

  for (i = 0; i < gr.nrows; i++) {
  	lat0 = gr.slat + i*gr.cellsize;
  	lat1 = lat0 + gr.cellsize;
  	for (j = 0; j < gr.ncols; j++) {
  		lon0 = gr.wlon + j*gr.cellsize;
  		lon1 = lon0 + gr.cellsize;

      if (gr.array[i][j] < min || gr.array[i][j] > max)
        continue ; // to get clipped min/max working. 

      if (gr.array[i][j] == gr.nodata_value) {
        cpgsci(0); // white
      } else {
        if (palette == GMT) {
           cindex = (int) 19.0*(gr.array[i][j] - min)/(max - min);
           cpgsci(cindex + 17);
        } else if (palette == RYGB) {
           cindex = (int) 66.0*(gr.array[i][j] - min)/(max - min);
           cpgsci(cindex + 17);
  	    } else { // Grayscale
  		     cindex = (int) 63.0*(gr.array[i][j] - min)/(max - min);
  		     cpgsci(cindex + 17);
  		  }
      }

  		if (box && lon0 > xmin && lon0 < xmax && lat0 > ymin && lat0 < ymax)
  		  cpgrect(lon0, lon1, lat0, lat1);
  		else 
  		  cpgrect(lon0, lon1, lat0, lat1); // what the hell is this?? STB 1016
  	}
  }

  cpgebuf();

// Draw land points. 

  if (plot_land) {

    cpgbbuf();
    cpgsci(15);
    for (i = 0; i < land.nrows; i++) {
      lat0 = land.slat + i*land.cellsize;
      lat1 = lat0 + land.cellsize;      
      for (j = 0; j < land.ncols; j++) {
        lon0 = land.wlon + j*land.cellsize;
        lon1 = lon0 + land.cellsize;
        if (land.array[i][j] == land.nodata_value) {
          cpgrect(lon0, lon1, lat0, lat1);
        }
      }
    }

    cpgebuf();
  }

// Draw points. 

  cpgsci(3);
  for (i = 0; i < npoints; i++) {
  	cpgpt1(xpts[i], ypts[i], 2);
  }

// Draw line between points. 

  if (npoints == 2) {
  	cpgsls(2);
  	for (i = 0; i < npoints - 1; i++) {
  		cpgmove(xpts[i], ypts[i]);
  		cpgdraw(xpts[i + 1], ypts[i + 1]);
  		cpgmove(xpts[i + 1], ypts[i + 1]);
  	}
  	cpgsls(1);
  }

// Redraw frame. 

  cpgsci(1);
  cpgbox("ABCNST", 0.0, 0, "ABCNST", 0.0, 0);

#endif
}











