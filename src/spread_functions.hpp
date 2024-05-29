#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <vector>

// Constants ---------------------------------------------------------------

// Elevation data to standardize distance between pixels
constexpr float elevation_mean = 1118.848;
constexpr float elevation_sd = 350.7474;
// from
// <fire_spread/data/NDVI_regional_data/ndvi_optim_and_proportion.rds>

// Structs and enums -------------------------------------------------------

enum veg_names { // not used
  forest,
  shrubland,
  grassland
};

// Coefficients names besides intercepts (they go from 0 to 2)
enum coef_names {
  b_ndvi = 3,
  b_north = 4,
  b_elev = 5,
  b_slope = 6,
  b_wind = 7
};

// Landscape layers names, besides vegetation
enum land_names {
  ndvi,   // pi_ndvi[v] * (ndvi - optim[v]) ^ 2
  north,  // slope-weighted
  elev,   // standardized, but also used to compute slope
  wdir,   // radians
  wspeed  // m/s
};

struct burned_res {
  IntegerMatrix burned_bin;
  IntegerMatrix burned_ids;
  int end;
  int steps_used;
};

struct burned_compare { // used for spatial overlap
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  int steps_used;
};

struct burned_compare_veg { // used for metrics taking vegetation into account
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  NumericVector counts_veg; // numeric to perform division later.
  int steps_used;
};

// Functions ---------------------------------------------------------------

burned_res simulate_fire_internal(
  const arma::fcube& landscape,
  const IntegerMatrix& vegetation,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  float upper_limit = 1.0,
  int steps = 0,
  double (*prob_fn)(double, double) = R::rbinom
);

burned_compare simulate_fire_compare_internal(
  const arma::fcube& landscape,
  const IntegerMatrix& vegetation,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  float upper_limit = 1.0,
  int steps = 0
);

burned_compare_veg simulate_fire_compare_veg_internal(
    const arma::fcube& landscape,
    const IntegerMatrix& vegetation,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 3,
    float upper_limit = 1.0,
    int steps = 0
);

// [[Rcpp::export]]
List simulate_fire_compare(
  const arma::fcube& landscape,
  const IntegerMatrix& vegetation,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  float upper_limit = 1.0,
  int steps = 0
);

// [[Rcpp::export]]
List simulate_fire_compare_veg(
    const arma::fcube& landscape,
    const IntegerMatrix& vegetation,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 3,
    float upper_limit = 1.0,
    int steps = 0
);
// Default arguments are provided here, not at the cpp file.