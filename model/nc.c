

#include "nc.h"

// The following definitions are not externally visible. 
void compress_2d(float **grid, short int *compressed, int ncols, int nrows, 
		double spval, float *add_offset, float *scale_factor, short int *missing_value);

void decompress_2d(short int *compressed, float **grid, int ncols, int nrows,
		float add_offset, float scale_factor, short int compressed_spval, double spval);

// Write a grid to a netCDF file. 

void write_nc(char *file, float **grid, int ncols, int nrows, 
	float *lons, float *lats, double spval) {

	int i, j, n, retval, ncid, lons_varid, lats_varid, grid_varid, 
		grid_dimid[2], lons_dimid, lats_dimid; 
	short int *compressed, missing_value;
	float add_offset, scale_factor, spval_float = spval;
	static char data_type[] = "regular gridded data",
		datum[] = "geographical coordinates, WGS84 projection",
		lon_name[] = "longitude", 
		lat_name[] = "latitude",
		lon_units[] = "degrees_east",
		lat_units[] = "degrees_north";

// Allocate memory for compressed grid. 

	compressed = (short int*) malloc(nrows*ncols*sizeof(short int));
	MEMCHECK(compressed);

// Compress data. 

	compress_2d(grid, compressed, ncols, nrows, spval, 
		&add_offset, &scale_factor, &missing_value);

//	printf("\n add_offset = %f, scale_factor = %f, missing_value = %d\n",
//		add_offset, scale_factor, missing_value);

// Create netCDF file. 

	retval = nc_create(trim(file), NC_CLOBBER, &ncid);
	if (retval) ERR(retval);

// Write global attributes. 

	retval = nc_put_att_text(ncid, NC_GLOBAL, "filename", strlen(file), file);
	if (retval) ERR(retval);
	retval = nc_put_att_text(ncid, NC_GLOBAL, "data_type", strlen(data_type), data_type);
	if (retval) ERR(retval);

// Define dimensions. 

	retval = nc_def_dim(ncid, "lons", ncols, &lons_dimid);
	if (retval) ERR(retval);
	retval = nc_def_dim(ncid, "lats", nrows, &lats_dimid);
	if (retval) ERR(retval);

	grid_dimid[0] = lats_dimid;
	grid_dimid[1] = lons_dimid;

// Define variables. 

	retval = nc_def_var(ncid, "lats", NC_FLOAT, 1, &lats_dimid, &lats_varid);
	if (retval) ERR(retval);

	retval = nc_def_var(ncid, "lons", NC_FLOAT, 1, &lons_dimid, &lons_varid);
	if (retval) ERR(retval);

	retval = nc_def_var(ncid, "bathy", NC_SHORT, 2, grid_dimid, &grid_varid);
	if (retval) ERR(retval);	

// Write attributes. 

	retval = nc_put_att_text(ncid, lons_varid, "standard_name", strlen(lon_name), lon_name);
	if (retval) ERR(retval);
	retval = nc_put_att_text(ncid, lons_varid, "units", strlen(lon_units), lon_units);
	if (retval) ERR(retval);
	retval = nc_put_att_text(ncid, lons_varid, "reference_datum", strlen(datum), datum);
	if (retval) ERR(retval);

	retval = nc_put_att_text(ncid, lats_varid, "standard_name", strlen(lat_name), lat_name);
	if (retval) ERR(retval);
	retval = nc_put_att_text(ncid, lats_varid, "units", strlen(lat_units), lat_units);
	if (retval) ERR(retval);
	retval = nc_put_att_text(ncid, lats_varid, "reference_datum", strlen(datum), datum);
	if (retval) ERR(retval);

// Write add_offset, scale_factor, and missing_value to file. 

	retval = nc_put_att_float(ncid, grid_varid, "add_offset", NC_FLOAT, 1, &add_offset);
	if (retval) ERR(retval);

	retval = nc_put_att_float(ncid, grid_varid, "scale_factor", NC_FLOAT, 1, &scale_factor);
	if (retval) ERR(retval);

	retval = nc_put_att_short(ncid, grid_varid, "missing_value", NC_SHORT, 1, &missing_value);
	if (retval) ERR(retval);

	retval = nc_put_att_float(ncid, grid_varid, "decompressed_missing_value", NC_FLOAT, 1, &spval_float);
	if (retval) ERR(retval);

	retval = nc_enddef(ncid);
	if (retval) ERR(retval);

// Write lats and lons to file. 

	retval = nc_put_var_float(ncid, lons_varid, &lons[0]);
	if (retval) ERR(retval);

	retval = nc_put_var_float(ncid, lats_varid, &lats[0]);
	if (retval) ERR(retval);

// Write compressed grid to file. 

	size_t start[2] = {0, 0};
	size_t count[2] = {nrows, ncols};

	// Note that the two-dimensional grid has been converted to a one-dimensional
	// grid as netCDF requires a contiguous block of memory. 

	retval = nc_put_vara_short(ncid, grid_varid, start, count, &compressed[0]);
	if (retval) ERR(retval);

// Close netCDF file. 

	retval = nc_close(ncid);
	if (retval) ERR(retval);

// Free memory. 

	free(compressed);

}


void read_nc(char *file, float ***grid, int *ncols, int *nrows, 
	float **lons, float **lats, double *spval) {

	int i, j, k, ncid, retval, lons_varid, lats_varid, grid_varid, 
		lons_dimid, lats_dimid, nxy;
	short int *compressed, compressed_spval;
	float add_offset, scale_factor, fspval;
	size_t nx, ny; 

// Open file. 

	retval = nc_open(trim(file), NC_NOWRITE, &ncid);
	if (retval) ERR(retval);

// Get variable ids. 

	retval = nc_inq_varid(ncid, "lons", &lons_varid);
	if (retval) ERR(retval);

	retval = nc_inq_varid(ncid, "lats", &lats_varid);
	if (retval) ERR(retval);	

	retval = nc_inq_varid(ncid, "bathy", &grid_varid);
	if (retval) ERR(retval);

// Read array dimensions. 

	retval = nc_inq_dimid(ncid, "lons", &lons_dimid);
	if (retval) ERR(retval);
	retval = nc_inq_dimlen(ncid, lons_dimid, &nx);
	if (retval) ERR(retval);
	*ncols = ((int) nx);

	retval = nc_inq_dimid(ncid, "lats", &lats_dimid);
	if (retval) ERR(retval);
	retval = nc_inq_dimlen(ncid, lats_dimid, &ny);
	if (retval) ERR(retval);
	*nrows = ((int) ny);

//	printf("\nnrows,ncols = %d, %d", *nrows, *ncols);

// Read add_offset, scale_factor, missing_value, and decompressed_missing_value

	retval = nc_get_att_float(ncid, grid_varid, "add_offset", &add_offset);
	if (retval) ERR(retval);	

	retval = nc_get_att_float(ncid, grid_varid, "scale_factor", &scale_factor);
	if (retval) ERR(retval);

	retval = nc_get_att_float(ncid, grid_varid, "decompressed_missing_value", &fspval);
	if (retval) {
		printf("\nWARNING: 'decompressed_missing_value' not present. Setting to zero.\n");
		fspval = 0.0;
	}
	*spval = ((double) fspval);

	retval = nc_get_att_short(ncid, grid_varid, "missing_value", &compressed_spval);
	if (retval) {
		compressed_spval = SHRT_MIN;
		printf("\nWARNING: 'missing_value' not present. Setting to %d.\n", compressed_spval);
	}

//	printf("\nadd_offset = %f\nscale_factor = %f\ndecompressed_missing_value = %f\n\
//missing_value = %d\n", add_offset, scale_factor, fspval, compressed_spval);

// Allocate memory. 

	(*lons) = (float*) malloc((*ncols)*sizeof(float));
	MEMCHECK(lons);

	(*lats) = (float*) malloc((*nrows)*sizeof(float));
	MEMCHECK(lats);

	nxy = (*nrows)*(*ncols);
	compressed = (short int*) malloc(nxy*sizeof(short int));
	MEMCHECK(compressed);

	allocate_float_array_2d(grid, *nrows, *ncols);

// Read longitudes. 

	retval = nc_get_var(ncid, lons_varid, &(*lons)[0]);
	if (retval) ERR(retval);

// Read latitudes. 

	retval = nc_get_var(ncid, lats_varid, &(*lats)[0]);
	if (retval) ERR(retval);

// Read compressed array. 

	retval = nc_get_var(ncid, grid_varid, &compressed[0]);
	if (retval) ERR(retval);

// Decompress array. 

	decompress_2d(compressed, *grid, *ncols, *nrows, 
		add_offset, scale_factor, compressed_spval, *spval);

// Close netCDF file.

	retval = nc_close(ncid);
	if (retval) ERR(retval);

// Free memory. 

	free(compressed);
}



// Decompress a 1D array of short ints to a 2D array of floats. 

void decompress_2d(short int *compressed, float **grid, int ncols, int nrows,
		float add_offset, float scale_factor, short int compressed_spval, double spval) {

	int i, j, k; 

	k = 0; 
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			if (compressed[k] == compressed_spval)
				grid[i][j] = ((float) spval);
			else
				grid[i][j] = ((float) compressed[k])*scale_factor + add_offset;
			k++;
		}
	}
}


// Compress a 2D array of floats to a 1D array of short ints. 

void compress_2d(float **grid, short int *compressed, int ncols, int nrows, double spval,
				float *add_offset, float *scale_factor, short int *missing_value) {

	int i, j, k; 
	float grmin, grmax, scale, fspv;

	fspv = ((float) spval);

// Assign missing_value. 

	*missing_value = SHRT_MIN;

// Compute min and max of grid. 

	grmin = FLT_MAX; 
	grmax = FLT_MIN;

	for(i = 0; i < nrows; i++) {
		for(j = 0; j < ncols; j++) {
			if (grid[i][j] == fspv)
				continue; 
			if (grid[i][j] < grmin)
				grmin = grid[i][j];
			else if (grid[i][j] > grmax)
				grmax = grid[i][j];
		}
	}

// Assign add_offset. 

	*add_offset = grmin; 

// Assign scale_factor. 

	scale = (grmax - grmin)/((float) SHRT_MAX);
	*scale_factor = scale;

// Compute compressed array. 

	k = 0;
	for(i = 0; i < nrows; i++) {
		for(j = 0; j < ncols; j++) {
			if (grid[i][j] == fspv)
				compressed[k++] = *missing_value;
			else
				compressed[k++] = ((short int) rint((grid[i][j] - grmin)/scale));
		}
	}	

}
















