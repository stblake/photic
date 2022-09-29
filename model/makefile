

cc= /usr/bin/gcc -O3 -lm 

# cc= gcc -g -O0 -lm -fopenmp 

netcdf=-I/opt/local/include -I/opt/local/include -L/opt/local/lib -lnetcdf 

pgplot=-L/opt/local/lib -I/opt/local/include -lpgplot -lcpgplot -lX11 

all:	photic	xyz2nc

clean:
	rm photic xyz2nc

photic:
	$(cc) $(netcdf) $(pgplot) common.c nc.c io.c graphics.c colours.c interp.c filllakesutils.c smoothutils.c jerlov.c secchi.c linenoise.c asa047.c samodelgraphics.c samodel.c align.c process.c bam.c -o photic

xyz2nc:
	$(cc) $(netcdf) common.c nc.c io.c xyz2nc.c -o xyz2nc
