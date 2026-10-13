
# Validation drivers, linked against photic's model objects.
#
#   make            builds photic (in ../model) and the drivers
#   make clean

MODEL = ../model

include $(MODEL)/makefile.inc

MODEL_OBJS = $(addprefix $(MODEL)/, $(filter-out bam.o, $(PHOTIC)))

ALL_CFLAGS = $(CSTD) $(CFLAGS) $(WARN) -I$(MODEL) $(NETCDF_CFLAGS) $(PGPLOT_CFLAGS)

.PHONY: all model clean

all: model
	$(MAKE) drivers/forward_dump

model:
	$(MAKE) -C $(MODEL) photic

drivers/forward_dump: drivers/forward_dump.c $(MODEL_OBJS)
	$(CC) $(ALL_CFLAGS) $^ $(NETCDF_LIBS) $(PGPLOT_LIBS) -lm -o $@

clean:
	rm -f drivers/forward_dump
