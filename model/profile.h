
#include "common.h"
#include "nc.h"

void profile_model(
  // inputs
  scene scene_data[], geogrid gridded_data[], int *scene_indexes, int nscenes,
  bool points_present, geogrid points,
  // outputs
  float **depth, float **depth_sigma, float **model_error, float **bottom_albedo,
  float **bottom_sand, float **bottom_seagrass, float **bottom_coral,
  float **K_min, float **bottom_type, float **index_optical_depth);
