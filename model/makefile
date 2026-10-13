
# photic makefile.
#
#   make                 builds photic and xyz2nc
#   make utils           builds the standalone grid utilities
#   make WARN=-Wall      builds with warnings
#   make PGPLOT_DIR=...  PGPLOT (or giza) installed somewhere other than MacPorts
#
# Must be compiled as C17 or earlier (see the bool typedef in common.h).

CC         = cc
CSTD       = -std=gnu17
CFLAGS     = -O3
WARN       =

NETCDF_CFLAGS := $(shell nc-config --cflags)
NETCDF_LIBS   := $(shell nc-config --libs)

PGPLOT_DIR    = /opt/local
PGPLOT_CFLAGS = -I$(PGPLOT_DIR)/include
PGPLOT_LIBS   = -L$(PGPLOT_DIR)/lib -lcpgplot -lpgplot -lX11

ALL_CFLAGS = $(CSTD) $(CFLAGS) $(WARN) -MMD -MP $(NETCDF_CFLAGS) $(PGPLOT_CFLAGS)

CORE   = common.o nc.o io.o

PHOTIC = $(CORE) graphics.o colours.o interp.o filllakesutils.o smoothutils.o \
	jerlov.o secchi.o linenoise.o asa047.o samodelgraphics.o samodel.o \
	align.o process.o refine.o write.o plot.o bam.o

UTILS  = asc2nc grid2grid smooth smoothspval emptygrid srtm filllakes

.PHONY: all utils clean

all: photic xyz2nc

utils: $(UTILS)

photic: $(PHOTIC)
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) $(PGPLOT_LIBS) -lm -o $@

xyz2nc asc2nc smoothspval emptygrid: %: $(CORE) %.o
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) -lm -o $@

grid2grid srtm: %: $(CORE) interp.o %.o
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) -lm -o $@

smooth: $(CORE) smoothutils.o smooth.o
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) -lm -o $@

filllakes: $(CORE) filllakesutils.o filllakes.o
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) -lm -o $@

%.o: %.c
	$(CC) $(ALL_CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.d photic xyz2nc $(UTILS)

-include $(wildcard *.d)
