#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <vector>

// Structs and enums -------------------------------------------------------

// Coefficients names
enum coef_names {
  intercept,
  b_vfi,
  b_tfi,
  b_slope,
  b_wind
};

// Landscape layers names
enum land_names {
  vfi,
  tfi,
  elev,
  wdir,
  wspeed
};

struct burned_res {
  IntegerMatrix burned_bin;
  IntegerMatrix burned_ids;
  int end;
};

struct burned_compare { // used for spatial overlap
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
};

struct burned_compare_veg { // used for metrics taking vegetation into account
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  NumericVector counts_veg; // numeric to perform division later.
};

// Functions ---------------------------------------------------------------

burned_res simulate_fire_internal(
  const arma::fcube& landscape,
  const IntegerMatrix& burnable,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  float upper_limit = 1.0,
  double (*prob_fn)(double, double) = R::rbinom
);

burned_compare simulate_fire_compare_internal(
  const arma::fcube& landscape,
  const IntegerMatrix& burnable,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  float upper_limit = 1.0
);

burned_compare_veg simulate_fire_compare_veg_internal(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    const IntegerMatrix& vegetation,
    int n_veg_types = 6,
    float upper_limit = 1.0
);

List simulate_fire_compare(
  const arma::fcube& landscape,
  const IntegerMatrix& burnable,
  const IntegerMatrix& ignition_cells,
  arma::frowvec coef,
  float upper_limit = 1.0
);

List simulate_fire_compare_veg(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    const IntegerMatrix& vegetation,
    int n_veg_types = 6,
    float upper_limit = 1.0
);

/*
 * Default arguments are provided here, not at the cpp file.
 */