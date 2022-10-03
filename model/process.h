
//	GRID PROCESSING ROUTINES


#include "common.h"
#include "io.h"
#include "interp.h"
#include "filllakesutils.h"
#include "smoothutils.h"
#if PGPLOT
#include "cpgplot.h"
#endif

void run_process_interpolate();
void run_process_smooth();
void run_process_shallow();
void run_process_pansharpen();
void run_process_noise();
void run_process_deglint();
void run_process_image();
void run_process_unite();
void run_process_cloud();
void run_process_mask();
void run_process_land();
void run_process_overlay();
void run_process_polygons();
void run_process_error();
void run_process_relativeerror();
void run_process_mndwi();
void run_process_fill();
void run_process_add();
void run_process_median();
void run_process();