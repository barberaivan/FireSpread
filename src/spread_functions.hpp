#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <vector>

// Structs and enums -------------------------------------------------------

enum b_terrain { // indexes in the coef_terrain vector
  b_slope,
  b_wind
};

// Names of variables used to compute slope and wind effects
enum terrain_names {
  elev,   // m
  wdir,   // radians
  wspeed  // m/s
};

struct burned_res {
  IntegerMatrix burned_layer;
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
    const IntegerMatrix& layer_vegetation,
    const arma::fcube& layer_nd,
    const arma::fcube& layer_terrain,
    arma::frowvec& coef_intercepts,
    arma::fvec& coef_nd,
    arma::frowvec& coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit = 1.0,
    int steps = 0,
    double (*prob_fn)(double, double) = R::rbinom
);

burned_compare simulate_fire_compare_internal(
    const IntegerMatrix& layer_vegetation,
    const arma::fcube& layer_nd,
    const arma::fcube& layer_terrain,
    arma::frowvec& coef_intercepts,
    arma::fvec& coef_nd,
    arma::frowvec& coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit = 1.0,
    int steps = 0
);

burned_compare_veg simulate_fire_compare_veg_internal(
    const IntegerMatrix& layer_vegetation,
    const arma::fcube& layer_nd,
    const arma::fcube& layer_terrain,
    arma::frowvec& coef_intercepts,
    arma::fvec& coef_nd,
    arma::frowvec& coef_terrain,
    const IntegerMatrix& ignition_cells,
    int n_veg = 5,
    float upper_limit = 1.0,
    int steps = 0
);

// [[Rcpp::export]]
List simulate_fire_compare(
    const IntegerMatrix& layer_vegetation,
    const arma::fcube& layer_nd,
    const arma::fcube& layer_terrain,
    arma::frowvec& coef_intercepts,
    arma::fvec& coef_nd,
    arma::frowvec& coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit = 1.0,
    int steps = 0
);

// [[Rcpp::export]]
List simulate_fire_compare_veg(
    const IntegerMatrix& layer_vegetation,
    const arma::fcube& layer_nd,
    const arma::fcube& layer_terrain,
    arma::frowvec& coef_intercepts,
    arma::fvec& coef_nd,
    arma::frowvec& coef_terrain,
    const IntegerMatrix& ignition_cells,
    int n_veg = 5,
    float upper_limit = 1.0,
    int steps = 0
);
// Default arguments are provided here, not at the cpp file.