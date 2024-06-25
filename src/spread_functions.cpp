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
 * The landscape is defined by a vegetation (integer) matrix, with classes from
 * 0 to n_veg_types - 1 (99 is non-burnable) and a terrain arma::fcube
 * (a 3D array), with the variables needed to compute slope and wind effects:
 * {elevation, wind direction, wind speed}.
 * Future versions could include an optional landscape with additional,
 * non-directional predictors, as northing, ndvi, elevation, etc. But fitting
 * such a model was complicated, so only veg types remain.
 *
 * The vegetation codes start at zero, and 99 is reserved for non-burnable.
 *
 * Parameters for vegetation (intercepts without reference) are passes by
 * coef_veg, while the effects of slope and wind are in coef_terrain, in
 * this order:
 * {slope, wind}
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

// Distance between pixels, used to compute slope effect.
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
//' @param int vegetation: vegetation class of the target cell, to subset
//'   coef_veg.
//' @param arma::frowvec coef_veg: intercepts for each vegetation type.
//' @param arma::frowvec terrain_burning: terrain data from burning cell.
//' @param arma::frowvec terrain_neighbour: terrain data from target neighbour.
//' @param arma::frowvec coef_terrain: slopes for slope and wind in the
//'   logistic regression
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
  int vegetation,
  arma::frowvec terrain_burning,
  arma::frowvec terrain_neighbour,
  arma::frowvec coef_veg,
  arma::frowvec coef_terrain,
  int position,
  float upper_limit = 1.0
) {

  /* Reminder:
   enum b_terrain {
     b_slope,
     b_wind
   };

   enum terrain_names {
     elev,   // m a.s.l.
     wdir,   // radians
     wspeed  // m/s
   };
  */

  // wind direction term
  float wind_term = cosf(angles[position] - terrain_burning[wdir]) *
                    terrain_burning[wspeed] * coef_terrain[b_wind];

  // slope term (from elevation and distance), only present if uphill
  float slope_term = 0.0f;
  float elev_diff = terrain_neighbour[elev] - terrain_burning[elev];
  if(elev_diff > 0.0f) {
   slope_term += sinf(atanf(elev_diff / distances[position])) *
                 coef_terrain[b_slope];
  }

  // compute linear predictor
  float linpred = coef_veg[vegetation] +
                  slope_term +
                  wind_term;

  // burn probability
  float prob = upper_limit / (1 + expf(-linpred));

  return prob;
}

// The same but evaluating the probability and simulating the burn from a
// Bernoulli distribution (here for backwards compatibility)

// [[Rcpp::export]]
int spread_one_cell(
  int vegetation,
  arma::frowvec terrain_burning,
  arma::frowvec terrain_neighbour,
  arma::frowvec coef_veg,
  arma::frowvec coef_terrain,
  int position,
  float upper_limit = 1.0
) {

  float prob = spread_one_cell_prob(
    vegetation,
    terrain_burning,
    terrain_neighbour,
    coef_veg,
    coef_terrain,
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
//'   int steps_used, the number of simulation steps used.

//' @param arma::fcube terrain: variables used to compute slope and wind
//'   effects: {elevation (m), wind direction (radians), wind speed (m/s)}.
//' @param IntegerMatrix vegetation: from 0 to n_veg_types, with 99 indicating
//'   non-burnable.
//' @param arma::frowvec coef_terrain: slope and wind parameters.
//' @param arma::frowvec coef_veg: intercepts for each vegetation type.
//' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
//'   the cell(s) where the fire begun. First row has the row_id, second row has
//'   the col_id.
//' @param float upper_limit: upper limit for spread probability.
//' @param int steps: maximum number of simulation steps allowed. If 0
//'   (the default), a very large number is set so the simulation is not limited.
//'   The burning of ignition points is considered the first step, so steps = 1
//'   burns only the ignition points. Bear in mind that the simulation may stop
//'   because there are no more burning cells, without reaching the maximum steps
//'   allowed.
//' @param function [unnamed]: function evaluating the burn probability. By
//'   default, it's R::rbinom(), but can be set to a deterministic behaviour to
//'   test whether the probability computation matches R's function. (Just for
//'   testing.)

burned_res simulate_fire_internal(
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit,
    int steps,
    double (*prob_fn)(double, double)
  ) {

  // define landscape dimensions
  int n_row = vegetation.nrow();
  int n_col = vegetation.ncol();
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
      arma::frowvec terrain_burning = terrain.tube(burned_ids(0, b), burned_ids(1, b));

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

        // Extract target vegetation
        int veg_target = vegetation(neighbours[0][n], neighbours[1][n]);

        bool burnable_cell =
          (burned_bin(neighbours[0][n], neighbours[1][n]) == 0) && // not burned
          (veg_target != 99);                                      // burnable
        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::frowvec terrain_neighbour = terrain.tube(neighbours[0][n],
                                                       neighbours[1][n]);

        // simulate fire
        float prob = spread_one_cell_prob(
          veg_target,
          terrain_burning,
          terrain_neighbour,
          coef_veg,
          coef_terrain,
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
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit = 1.0,
    int steps = 0
) {

  // define landscape dimensions
  int n_row = vegetation.nrow();
  int n_col = vegetation.ncol();
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
      arma::frowvec terrain_burning = terrain.tube(burned_ids(0, b), burned_ids(1, b));

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

        // Extract target vegetation
        int veg_target = vegetation(neighbours[0][n], neighbours[1][n]);

        bool burnable_cell =
          (burned_step(neighbours[0][n], neighbours[1][n]) == 0) && // not burned
          (veg_target != 99);                                       // burnable
        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::frowvec terrain_neighbour = terrain.tube(neighbours[0][n], neighbours[1][n]);

        // simulate fire
        float prob = spread_one_cell_prob(
          veg_target,
          terrain_burning,
          terrain_neighbour,
          coef_veg,
          coef_terrain,
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
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit = 1.0,
    int steps = 0
) {
  return simulate_fire_internal(
    vegetation,
    terrain,
    coef_veg,
    coef_terrain,
    ignition_cells,
    upper_limit,
    steps
  ).burned_bin;
}

// -----------------------------------------------------------------------

// The same function but deterministic, to test if the discrepancy between R and
// cpp is caused by seed problems

// [[Rcpp::export]]
IntegerMatrix simulate_fire_deterministic(
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit = 1.0,
    int steps = 0
  ) {

  return simulate_fire_internal(
    vegetation,
    terrain,
    coef_veg,
    coef_terrain,
    ignition_cells,
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
// other than the spatial overlap.
// NumericVector counts_veg: burned pixels by veg_type. It has to be numeric
//   and not integer to allow non-integer divisions later.

burned_compare simulate_fire_compare_internal(
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit,
    int steps
) {

  burned_res burned = simulate_fire_internal(
    vegetation,
    terrain,
    coef_veg,
    coef_terrain,
    ignition_cells,
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
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit,
    int steps
) {

  burned_res burned = simulate_fire_internal(
    vegetation,
    terrain,
    coef_veg,
    coef_terrain,
    ignition_cells,
    upper_limit,
    steps
  );

  IntegerMatrix burned_bin = burned.burned_bin;
  IntegerMatrix burned_ids = burned.burned_ids;
  int end = burned.end;
  int steps_used = burned.steps_used;
  int n_veg_types = coef_veg.size();

  // Compute burned area by vegetation type
  NumericVector counts_veg(n_veg_types);
  for(int i = 0; i <= end; i++) {
    int veg_i = vegetation(burned_ids(0, i), burned_ids(1, i));
    counts_veg[veg_i]++;
  }

  return {burned_bin, burned_ids(_, seq(0, end)), counts_veg, steps_used};
}

List simulate_fire_compare(
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit,
    int steps
) {

  burned_compare burned_com = simulate_fire_compare_internal(
    vegetation,
    terrain,
    coef_veg,
    coef_terrain,
    ignition_cells,
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
    const IntegerMatrix& vegetation,
    const arma::fcube& terrain,
    arma::frowvec coef_veg,
    arma::frowvec coef_terrain,
    const IntegerMatrix& ignition_cells,
    float upper_limit,
    int steps
) {

  burned_compare_veg burned_com = simulate_fire_compare_veg_internal(
    vegetation,
    terrain,
    coef_veg,
    coef_terrain,
    ignition_cells,
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