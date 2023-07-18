#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <vector>

// Constants ---------------------------------------------------------------

// Elevation data to standardize distance between pixels
constexpr float elevation_mean = 1163.3;
constexpr float elevation_sd = 399.5;

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

struct burned_compare {
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  NumericVector counts_veg; // need to be numeric to compute divisions later
};

// Functions ---------------------------------------------------------------

burned_res simulate_fire_internal(
  const IntegerMatrix& vegetation,
  const arma::fcube& terrain,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  int n_veg_types = 6,
  float upper_limit = 1.0,
  double (*prob_fn)(double, double) = R::rbinom
);

burned_compare simulate_fire_compare_internal(
  const IntegerMatrix& vegetation,
  const arma::fcube& terrain,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  int n_veg_types = 6,
  float upper_limit = 1.0
);

List simulate_fire_compare(
  const IntegerMatrix& vegetation,
  const arma::fcube& terrain,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  int n_veg_types = 6,
  float upper_limit = 1.0
);
