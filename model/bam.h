

#include "common.h"
#include "io.h"
#include "nc.h"
#include "interp.h"
#include "smoothutils.h"
#include "filllakesutils.h"
#include "colours.h"
#include "graphics.h"
#include "linenoise.h"
#include "jerlov.h"
#include "secchi.h"
#include "samodel.h"
#include "process.h"
#include "align.h"
#include "refine.h"
#include "write.h"
#include "plot.h"
#if PGPLOT
#include "cpgplot.h"
#endif

#define DEBUG 1

#define UNDER_CONSTRUCTION printf("\nUNDER CONSTRUCTION...\n\n")


