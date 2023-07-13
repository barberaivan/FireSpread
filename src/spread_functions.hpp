#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <vector>

// Constants ---------------------------------------------------------------

// Elevation data to standardize distance between pixels
const float elevation_mean = 1163.3;
const float elevation_sd = 399.5;

// Structs ans enums -------------------------------------------------------

// Terrain coefficients names
enum terrain_names {
  northing,
  elev,
  windir,
  slope
};

// Vegetation coefficients names
enum veg_names {
  shrubland,
  subalpine,
  wet,
  dry_a, // araucaria
  dry_b, // cypress
  steppe
};

struct burned_res {
  IntegerMatrix burned_bin;
  IntegerMatrix burned_ids;
  int end;
};

typedef struct _s_burned_compare {
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  NumericVector counts_veg; // need to be numeric to compute divisions later
} burned_compare;

// Functions ---------------------------------------------------------------

burned_res simulate_fire_internal(
  IntegerMatrix vegetation,
  arma::fcube terrain,
  IntegerMatrix ignition_cells,
  arma::frowvec coef,
  int n_veg_types = 6,
  float upper_limit = 1.0,
  double (*prob_fn)(double, double) = R::rbinom
);

burned_compare simulate_fire_compare_cpp(
  IntegerMatrix vegetation,
  arma::fcube terrain,
  IntegerMatrix ignition_cells,
  arma::frowvec coef,
  int n_veg_types = 6,
  float upper_limit = 1.0
);

List simulate_fire_compare(
  IntegerMatrix vegetation,
  arma::fcube terrain,
  IntegerMatrix ignition_cells,
  arma::frowvec coef,
  int n_veg_types = 6,
  float upper_limit = 1.0
);
