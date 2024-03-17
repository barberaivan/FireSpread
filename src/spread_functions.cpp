#include "spread_functions.hpp"

// armadillo is used for the matrix representation of the rasters, to use the
// fcube data type.

// Useful armadillo links
//   https://dcgerard.github.io/advancedr/08_cpp_armadillo.html
//   https://arma.sourceforge.net/docs.html#subfcube

using namespace Rcpp;

/*
 * Functions to spread fire.
 *
 * The fire spread model is a cellular automata that spreads towards the 8-
 * neighbouring pixels, with a burn probability defined by a logistic linear
 * model.
 *
 * The landscape is an arma::fcube (a 3D array), where each layer is a matrix:
 *   vfi: vegetation flammability index,
 *   tfi: topographic flammability index,
 *   elev: elevation (m), used to compute directional slope effect,
 *   wdir: direction from where the wind blows (°),
 *   wspeed: wind speed (m/s).
 * As slope and wind effects are directional, they are computed during the
 * simulation.
 *
 * In addition, a burnable binary layer is provided to indicate which cells are
 * burnable (avoiding fire in lakes or high-andean areas).
 *
 * The coef vector holds all the parameters for the logistic regression:
 * {intercept, b_vfi, b_tfi, b_slope, b_wind},
 *
 * The climate is represented by the FWI, and is treated as
 * constant within a fire or landscape. Its value at the ignition point
 * is used to define the mean of a fire-level random effect, with a linear
 * function. This is not included as a pixel-level variable because of its low
 * resolution, which made its values almost constant within landscapes, which
 * would generate high correlation between the fwi and other parameters.
 *
 * The simulations stops when there are no more burning cells or when a number
 * of pre-defined steps are completed. This could represent the duration of
 * fire-prone weather, which is a latent variable.
 *
 * At the bottom there are functions used to simulate and compare fires.
 * Those ending in _veg count the burned cells by vegetation type. To use them
 * it's necessary to provide a vegetation layer (discrete) and the number of
 * vegetation types.
 */

// Constants ---------------------------------------------------------------

// Distance between pixels (m), used to compute slope effect.
// 30 m is the landscape resolution.
const float distances[8] = {
  30 * sqrtf(2), 30, 30 * sqrtf(2),
  30,                30,
  30 * sqrtf(2), 30, 30 * sqrtf(2)
};

constexpr int moves[8][2] = {
  {-1, -1},
  {-1,  0},
  {-1,  1},
  { 0, -1},
  { 0,  1},
  { 1, -1},
  { 1,  0},
  { 1,  1}
};

/* In the case
 * of a 8-pixels neighbourhood, it's a matrix with 2 rows (row and column
 * values) and 8 columns. Its values are {-1, 0, 1}, so that when adding up
 * the row-col ids of a cell and a column of moves, we get the row-col of a
 * neighbour. They are ordered like this:
 * 1 2 3
 * 4   5
 * 6 7 8
 * For example, to get the neighbour 1, we compute
 * focal_row - 1,
 * focal_column - 1.
 * The first column in moves is then
 * -1
 * -1
 * (rows and columns are numbered from the upper-left corner.)
 */

// Angles between cells to compute wind effect. As the wind direction is
// the direction from which the wind comes, these angles must represent where the
// fire would come from if from the neighbours we look at the central pixel.

constexpr float angles[8] = {
  M_PI * 3 / 4, M_PI, M_PI * 5 / 4,
  M_PI / 2,           M_PI * 3 / 2,
  M_PI / 4,      0,   M_PI * 7 / 4
};


// --------------------------------------------------------------------------

//' @title spread_one_cell_prob
//' @description Calculates the probability of a cell spreading fire to another.
//' @return float [0, 1] indicating the probability.
//'
//' @param arma::frowvec landscape_burning: landscape data from burning cell.
//' @param arma::frowvec landscape_neighbour: landscape data from target neighbour.
//' @param arma::frowvec coef: logistic regression parameters.
//' @param int position: relative position of the target in relation to the
//' burning cell. The eight neighbours are labelled from 0 to 7 beginning from
//' the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the slope and wind effects, as they
//'   depend on the angle and distance between burning and target cells.
//' @param float upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires; 0.5 is preferred).

// [[Rcpp::export]]
float spread_one_cell_prob(
  arma::frowvec landscape_burning,
  arma::frowvec landscape_neighbour,
  arma::frowvec coef,
  int position,
  float upper_limit = 1.0
) {

  // wind direction term
  float wdir_term = cosf(angles[position] - landscape_burning(wdir));

  // slope term (from elevation and distance)
  float slope_term = sinf(atanf(
    (landscape_neighbour(elev) - landscape_burning(elev)) / distances[position]
  ));

  // compute linear predictor
  float linpred = coef[intercept];

  linpred += coef[b_vfi] * landscape_neighbour[vfi] +
             coef[b_tfi] * landscape_neighbour[tfi] +
             coef[b_slope] * slope_term +
             coef[b_wind] * wdir_term * landscape_burning[wspeed];

  // burn probability
  float prob = upper_limit / (1 + expf(-linpred));

  return prob;
}

// The same but evaluating the probability and simulating the burn from a
// Bernoulli distribution (here for backwards compatibility)

// [[Rcpp::export]]
int spread_one_cell(
  arma::frowvec landscape_burning,
  arma::frowvec landscape_neighbour,
  arma::frowvec coef,
  int position,
  float upper_limit = 1.0
) {

  float prob = spread_one_cell_prob(
    landscape_burning,
    landscape_neighbour,
    coef,
    position,
    upper_limit
  );

  return (int)R::rbinom(1.0, prob);
}

// -----------------------------------------------------------------------

//' @title simulate_fire_internal
//' @description function to simulate a fire spread given the landscape,
//'   model coefficients and ignition points.
//' @return burned_res: struct containing the following objects:
//'   IntegerMatrix burned_bin, a binary matrix indicating the burned pixels
//'   IntegerMatrix burned_ids, a matrix with a column by burned pixel,
//'     indicating its row (row1) and column (row2) in the landscape,
//'   int end, the number of burned pixels.

//' @param arma::fcube landscape: predictor variables or variables used to
//'   compute the actual predictors.
//'     vfi: vegetation flammability index,
//'     tfi: topographic flammability index,
//'     elev: elevation (m), used to compute directional slope effect,
//'     wdir: direction from where the wind blows (°),
//'     wspeed: wind speed (m/s).
//' @param IntegerMatrix burnable: binary layer indicating cells available to
//'   burn.
//' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
//'   the cell(s) where the fire begun. First row has the row_id, second row has
//'   the col_id.
//' @param arma::frowvec coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param float upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).
//' @param int steps: maximum number of simulation steps allowed. If 0
//'   (the default), a large number is used to avoid limiting fire spread.
//'   The burning of ignition points is considered the first step, so steps = 1
//'   burns only the ignition points. Bear in mind that the simulation may stop
//'   because there are no more burning cells, without reaching the maximum steps
//'   allowed.
//' @param function [unnamed]: function evaluating the burn probability. By
//'   default, it's R::rbinom(), but can be set to a deterministic behaviour to
//'   test whether the probability computation matches R's function. (Just for
//'   testing.)

burned_res simulate_fire_internal(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    float upper_limit,
    int steps,
    double (*prob_fn)(double, double)
  ) {

  // define landscape dimensions
  int n_row = burnable.nrow();
  int n_col = burnable.ncol();
  int n_cell = n_row * n_col;

  // if zero steps are allowed, use a lot to avoid limiting the simulation
  if(steps == 0) steps = n_cell * 10;

  // burned_ids [row-col, cell] will be filled with the row_col ids (rows) of the
  // burning pixels (columns). start and end integers will define the positions
  // limits corresponding to the burning cells in every burn step.
  IntegerMatrix burned_ids(2, n_cell); // check it's filled with 0 // -2147483648 is NA_INTEGER

  int start = 0;
  // end is the last non-empty position in the burned_ids matrix.
  int end = ignition_cells.ncol() - 1;
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.

  // initialize burned_ids and burning_size with ignition_cells
  for(int c = 0; c <= end; c++) {
    for(int r = 0; r < 2; r++) {
      burned_ids(r, c) = ignition_cells(r, c);
    }
  }

  // initialize burning_size
  int burning_size = ignition_cells.ncol(); // == end + 1 - start

  // The burned_bin matrix will indicate whether each pixel is burned or burning
  // (1) or not (0). It's necessary to have this now because it will be used
  // to define burnable neighbours.
  IntegerMatrix burned_bin(n_row, n_col);

  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_bin(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }

  // initialize simulation step, considering the ignition as the first one.
  // The step is updated once the loop is entered.
  int step = 1;

  while(burning_size > 0 && step < steps) {
    // update step
    step++;

    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;

    // Loop over burning cells in the step

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn step
    for(int b = start; b <= end; b++) {

      // Get burning_cell's data
      arma::frowvec landscape_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));

      int neighbours[2][8];
      // get neighbours (adjacent computation here)
      for(int i = 0; i < 8; i++) {
        neighbours[0][i] = burned_ids(0, b) + moves[i][0];
        neighbours[1][i] = burned_ids(1, b) + moves[i][1];
      }

      // Loop over neighbours of the focal burning cell

      for(int n = 0; n < 8; n++) {

        // Is the cell in range?
        bool out_of_range = (
          (neighbours[0][n] < 0) || (neighbours[0][n] >= n_row) || // check rows
          (neighbours[1][n] < 0) || (neighbours[1][n] >= n_col)    // check cols
        );
        if(out_of_range) continue;

        // Is the cell burnable?
        bool burnable_cell =
          (burned_bin(neighbours[0][n], neighbours[1][n]) == 0) && // not burned
          (burnable(neighbours[0][n], neighbours[1][n]) == 1);     // burnable
        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::frowvec landscape_neighbour = landscape.tube(neighbours[0][n], neighbours[1][n]);

        // simulate fire
        float prob = spread_one_cell_prob(
          landscape_burning,
          landscape_neighbour,
          coef,
          n,           // position argument, from 0 to 7
          upper_limit
        );

        int burn = int(prob_fn(1.0, prob));
        if(burn == 0) continue;

        // If burned,
        // store id of recently burned cell and
        // set 1 in burned_bin
        // (but advance end_forward first)
        end_forward++;
        burned_ids(0, end_forward) = neighbours[0][n];
        burned_ids(1, end_forward) = neighbours[1][n];
        burned_bin(neighbours[0][n], neighbours[1][n]) = 1;

      } // end loop over neighbours of burning cell b

    } // end loop over burning cells from this step

    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;

  } // end while

  return {burned_bin, burned_ids, end, step};
}

// a similar function only returning a burned_layer with an integer by burn
// step, used to animate the spread.

// [[Rcpp::export]]
IntegerMatrix simulate_fire_animate(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    float upper_limit = 1.0,
    int steps = 0
) {

  // define landscape dimensions
  int n_row = burnable.nrow();
  int n_col = burnable.ncol();
  int n_cell = n_row * n_col;

  // if zero steps are allowed, use a lot to avoid limiting the simulation
  if(steps == 0) steps = n_cell * 10;

  // burned_ids [row-col, cell] will be filled with the row_col ids (rows) of the
  // burning pixels (columns). start and end integers will define the positions
  // limits corresponding to the burning cells in every burn step.
  IntegerMatrix burned_ids(2, n_cell); // check it's filled with 0 // -2147483648 is NA_INTEGER

  int start = 0;
  // end is the last non-empty position in the burned_ids matrix.
  int end = ignition_cells.ncol() - 1;
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.

  // initialize burned_ids and burning_size with ignition_cells
  for(int c = 0; c <= end; c++) {
    for(int r = 0; r < 2; r++) {
      burned_ids(r, c) = ignition_cells(r, c);
    }
  }

  // initialize burning_size
  int burning_size = ignition_cells.ncol(); // == end + 1 - start

  // The burned_step matrix will indicate the step at which each pixel was
  // burned, starting from the ignition (there will always be at least a "1"
  // pixel). Pixels with 0 are unburned.
  IntegerMatrix burned_step(n_row, n_col);

  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_step(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }

  // initialize simulation step, considering the ignition as the first one.
  // The step is updated once the loop is entered.
  int step = 1;

  while(burning_size > 0 && step < steps) {
    // update step
    step++;
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;

    // Loop over burning cells in the step

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn step
    for(int b = start; b <= end; b++) {

      // Get burning_cell's data
      arma::frowvec landscape_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));

      int neighbours[2][8];
      // get neighbours (adjacent computation here)
      for(int i = 0; i < 8; i++) {
        neighbours[0][i] = burned_ids(0, b) + moves[i][0];
        neighbours[1][i] = burned_ids(1, b) + moves[i][1];
      }

      // Loop over neighbours of the focal burning cell

      for(int n = 0; n < 8; n++) {

        // Is the cell in range?
        bool out_of_range = (
          (neighbours[0][n] < 0) || (neighbours[0][n] >= n_row) || // check rows
          (neighbours[1][n] < 0) || (neighbours[1][n] >= n_col)    // check cols
        );
        if(out_of_range) continue;

        // Is the cell burnable?
        bool burnable_cell =
          (burned_step(neighbours[0][n], neighbours[1][n]) == 0) && // not burned
          (burnable(neighbours[0][n], neighbours[1][n]) == 1);      // burnable
        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::frowvec landscape_neighbour = landscape.tube(neighbours[0][n], neighbours[1][n]);

        // simulate fire
        float prob = spread_one_cell_prob(
          landscape_burning,
          landscape_neighbour,
          coef,
          n,           // position argument, from 0 to 7
          upper_limit
        );

        int burn = R::rbinom(1.0, prob);
        if(burn == 0) continue;

        // If burned,
        // store id of recently burned cell and
        // set step in burned_step
        // (but advance end_forward first)
        end_forward++;
        burned_ids(0, end_forward) = neighbours[0][n];
        burned_ids(1, end_forward) = neighbours[1][n];
        burned_step(neighbours[0][n], neighbours[1][n]) = step;

      } // end loop over neighbours of burning cell b

    } // end loop over burning cells from this step

    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;

  } // end while

  return burned_step;
}

// -----------------------------------------------------------------------

// The same function to be exported to R, only returning the binary burned_bin
// matrix.
// [[Rcpp::export]]
IntegerMatrix simulate_fire(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    float upper_limit = 1.0,
    int steps = 0
) {
  return simulate_fire_internal(
    landscape,
    burnable,
    ignition_cells,
    coef,
    upper_limit,
    steps
  ).burned_bin;
}

// -----------------------------------------------------------------------

// The same function but deterministic, to test if the discrepancy between R and
// cpp is caused by seed problems

// [[Rcpp::export]]
IntegerMatrix simulate_fire_deterministic(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    float upper_limit = 1.0,
    int steps = 0
  ) {

  return simulate_fire_internal(
    landscape,
    burnable,
    ignition_cells,
    coef,
    upper_limit,
    steps,
    [](double _, double x) { return (double)(x >= 0.5); }
  ).burned_bin;
}

// -------------------------------------------------------------------------

// simulate_fire_compare: same as simulate_fire, but returning objects to
// compute discrepancy or similarity metrics (as a List):
//   IntegerMatrix burned_layer: binary matrix indicating the burned pixels,
//   IntegerMatrix burned_ids: ids as [row, col] column-vectors indicating the
//     position of burned pixels.
//
// The _veg versions return the count_veg, a vector with the number of cells
// burned by vegetation type. This is needed to compute similarity metrics
// other than the spatial overlap. They require 2 extra arguments:
//  vegetation, a discrete matrix, and
//  n_veg_types.
// NumericVector counts_veg: burned pixels by veg_type. It has to be numeric
//   and not integer to allow non-integer divisions later.

burned_compare simulate_fire_compare_internal(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    float upper_limit,
    int steps
) {

  burned_res burned = simulate_fire_internal(
    landscape,
    burnable,
    ignition_cells,
    coef,
    upper_limit,
    steps
  );

  IntegerMatrix burned_bin = burned.burned_bin;
  IntegerMatrix burned_ids = burned.burned_ids;
  int end = burned.end;
  int steps_used = burned.steps_used;

  return {burned_bin, burned_ids(_, seq(0, end)), steps_used};
}

burned_compare_veg simulate_fire_compare_veg_internal(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    const IntegerMatrix& vegetation,
    int n_veg_types,
    float upper_limit,
    int steps
) {

  burned_res burned = simulate_fire_internal(
    landscape,
    burnable,
    ignition_cells,
    coef,
    upper_limit,
    steps
  );

  IntegerMatrix burned_bin = burned.burned_bin;
  IntegerMatrix burned_ids = burned.burned_ids;
  int end = burned.end;
  int steps_used = burned.steps_used;

  // Compute burned area by vegetation type
  NumericVector counts_veg(n_veg_types);
  for(int i = 0; i <= end; i++) {
    int veg_i = vegetation(burned_ids(0, i), burned_ids(1, i));
    counts_veg[veg_i]++;
  }

  return {burned_bin, burned_ids(_, seq(0, end)), counts_veg, steps_used};
}

List simulate_fire_compare(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    float upper_limit,
    int steps
) {

  burned_compare burned_com = simulate_fire_compare_internal(
    landscape,
    burnable,
    ignition_cells,
    coef,
    upper_limit,
    steps
  );

  // List to return:
  List L = List::create(Named("burned_layer") = burned_com.burned_layer,
                        Named("burned_ids") = burned_com.burned_ids,
                        Named("steps_used") = burned_com.steps_used);

  return L;
}

List simulate_fire_compare_veg(
    const arma::fcube& landscape,
    const IntegerMatrix& burnable,
    const IntegerMatrix& ignition_cells,
    arma::frowvec coef,
    const IntegerMatrix& vegetation,  // arguments with defaults go last.
    int n_veg_types,
    float upper_limit,
    int steps
) {

  burned_compare_veg burned_com = simulate_fire_compare_veg_internal(
    landscape,
    burnable,
    ignition_cells,
    coef,
    vegetation,
    n_veg_types,
    upper_limit,
    steps
  );

  // List to return:
  List L = List::create(Named("burned_layer") = burned_com.burned_layer,
                        Named("burned_ids") = burned_com.burned_ids,
                        Named("counts_veg") = burned_com.counts_veg,
                        Named("steps_used") = burned_com.steps_used);

  return L;
}