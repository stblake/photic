
#include "common.h"
#include "math.h"
#if PGPLOT
#include "cpgplot.h"
#include "colours.h"
#endif

#define JET 101
#define GRAYSCALE 102
#define OCEAN 103
#define GMT 104
#define INVERSEJET 105
#define RYGB 106
#define WOR 108
#define VIRIDIS 109
#define PLASMA 110

void plot_scatter(char *title, char *xlabel, char *ylabel, int colour, int marker, int npoints, 
	float *xpts, float *ypts, int line_width, float xmin, float xmax, float ymin, float ymax, bool first);

void plot_array(geogrid gr, char *title, float pagesize, float textsize, float pointsize,
	float xmin, float ymin, float xmax, float ymax, int palette, 
	int npoints, float *xpts, float *ypts, float aspect_factor, int line_width,
  bool plot_land, geogrid land, float user_min, float user_max);

void plot_log_bidimensional_histogram(char *title, char *xlabel, char *ylabel,  
  int npoints, float *xpts, float *ypts, float lsmx, float lsmy, int line_width,
  float xmin, float xmax, float ymin, float ymax, bool first, int palette, 
  int nxbins, int nybins);

void plot_bidimensional_histogram(char *title, char *xlabel, char *ylabel,  
  int npoints, float *xpts, float *ypts, int line_width,
  float xmin, float xmax, float ymin, float ymax, bool first, int palette, 
  int nxbins, int nybins);

