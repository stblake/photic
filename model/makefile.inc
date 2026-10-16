
# Compiler settings and object lists, shared by makefile and ../validation/makefile.

CC         = cc
CSTD       = -std=gnu17
CFLAGS     = -O3
WARN       =

NETCDF_CFLAGS := $(shell nc-config --cflags)
NETCDF_LIBS   := $(shell nc-config --libs)

PGPLOT_DIR    = /opt/local
PGPLOT_CFLAGS = -I$(PGPLOT_DIR)/include
PGPLOT_LIBS   = -L$(PGPLOT_DIR)/lib -lcpgplot -lpgplot -lX11

CORE   = common.o nc.o io.o

PHOTIC = $(CORE) graphics.o colours.o interp.o filllakesutils.o smoothutils.o \
	jerlov.o secchi.o linenoise.o asa047.o samodelgraphics.o spectra.o samodel.o \
	align.o process.o refine.o write.o plot.o bam.o


# Default spectral data files (model/data), compiled into spectra.o; $PHOTIC_DATA overrides at run time.

PHOTIC_DATA_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))/data
