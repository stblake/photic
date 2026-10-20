
# photic makefile.
#
#   make                 builds photic and xyz2nc
#   make utils           builds the standalone grid utilities
#   make test            runs the synthetic regression tests (semi-analytical and Lyzenga models)
#   make WARN=-Wall      builds with warnings
#   make PGPLOT_DIR=...  PGPLOT (or giza) installed somewhere other than MacPorts
#
# Must be compiled as C17 or earlier (see the bool typedef in common.h).

include makefile.inc

ALL_CFLAGS = $(CSTD) $(CFLAGS) $(WARN) -MMD -MP -I. $(NETCDF_CFLAGS) $(PGPLOT_CFLAGS)

UTILS  = asc2nc grid2grid smooth smoothspval emptygrid srtm filllakes

.PHONY: all utils test clean

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

# Synthetic regression test. tests/synth_gen writes band grids from known depths with 
# photic's own forward model, photic runs MODEL Lee_Semi_Analytical on them, and the 
# outputs must match tests/baseline/ (make test-baseline regenerates it after an 
# intentional change). The accuracy against the known depths is reported, not checked. 

TEST_OUTPUTS = modelled_H modelled_error modelled_albedo modelled_min_K modelled_bottom_sand \
	s1_P s1_G s1_X s2_P s2_G s2_X s1_K_blue s2_K_blue

tests/synth_gen: $(filter-out bam.o, $(PHOTIC)) tests/synth_gen.o
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) $(PGPLOT_LIBS) -lm -o $@

tests/synth_check: $(CORE) tests/synth_check.o
	$(CC) $(CFLAGS) $^ $(NETCDF_LIBS) -lm -o $@

# MODEL lyzenga_log_linear on two bands that follow its model exactly: the coefficients and
# every sounding's depth must be recovered (tests/lyz_check).
tests/lyz_gen: tests/lyz_gen.o
	$(CC) $(CFLAGS) $^ -lm -o $@

tests/lyz_check: tests/lyz_check.o
	$(CC) $(CFLAGS) $^ -lm -o $@

tests/lyz_run: photic tests/lyz_gen
	rm -rf tests/lyz_run && mkdir -p tests/lyz_run
	tests/lyz_gen tests/lyz_run
	cp tests/lyz.bam tests/lyz_run/
	cd tests/lyz_run && ../../photic lyz.bam < /dev/null > photic.log 2>&1

tests/run: photic tests/synth_gen
	rm -rf tests/run && mkdir -p tests/run
	tests/synth_gen tests/run > /dev/null
	cp tests/synth.bam tests/run/
	cd tests/run && ../../photic synth.bam < /dev/null > photic.log 2>&1

.PHONY: tests/run tests/lyz_run test-baseline

test: tests/run tests/synth_check tests/lyz_run tests/lyz_check
	@fail=0; for f in $(TEST_OUTPUTS); do \
	  tests/synth_check baseline tests/run/$$f.nc tests/baseline/$$f.nc || fail=1; \
	done; \
	tests/synth_check truth tests/run/modelled_H.nc tests/run/truth_H.asc; \
	tests/lyz_check tests/lyz_run || fail=1; \
	if [ $$fail = 0 ]; then echo "\nPASS"; else echo "\nFAIL (see tests/run/photic.log)"; exit 1; fi

test-baseline: tests/run
	mkdir -p tests/baseline
	for f in $(TEST_OUTPUTS); do cp tests/run/$$f.nc tests/baseline/; done

%.o: %.c
	$(CC) $(ALL_CFLAGS) -c $< -o $@

spectra.o: ALL_CFLAGS += -DPHOTIC_DATA_DIR='"$(PHOTIC_DATA_DIR)"'

clean:
	rm -f *.o *.d tests/*.o tests/*.d photic xyz2nc $(UTILS) tests/synth_gen tests/synth_check tests/lyz_gen tests/lyz_check
	rm -rf tests/run

-include $(wildcard *.d tests/*.d)
