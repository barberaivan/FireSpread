# Functions to spread fire, written in R to check results with C++.
# For a detailed explanation, see <spread_functions.cpp>.

library(terra)

# Constants ----------------------------------------------------------------

# Elevation data to standardize distance between pixels
elevation_mean = 1118.848
elevation_sd = 350.7474
# from
# <fire_spread/data/NDVI_regional_data/ndvi_optim_and_proportion.rds>

# Distance between pixels, in m / elevation_sd.
# 30 m is the landscape resolution.
distances <- rep(30, 8)
distances[c(1, 3, 6, 8)] <- distances[c(1, 3, 6, 8)] * sqrt(2)
distances <- distances / elevation_sd

# Angles between cells to compute wind effect. As the wind direction is
# the direction from which the wind comes, these angles must represent where the
# fire would come from if from the neighbours we look at the central pixel.
angles <- c(
    135, 180, 225,
    90,       270,
    45,   0,  315
    ) * pi / 180   # in radians

# queen neighbours movements
moves <- matrix(c(-1,-1,-1,  0,0,  1,1,1,
                  -1, 0, 1, -1,1, -1,0,1),
                nrow = 2, byrow = TRUE)
# in the eight-neighbour setting (1 step queen moves), to arrive at neighbours
# we can make eight moves of rows and columns. moves indicates the row movements
# (moves[1, ]) and the column movements (moves[2, ]) to get each neighbour.
# neighbours are ordered row-wise:
# (1, 2, 3,
#  4, NA, 5,
#  6, 7, 8)

# Coefficients names
forest    <- 1
shrubland <- 2
grassland <- 3

b_ndvi    <- 4
b_north   <- 5
b_elev    <- 6
b_slope   <- 7
b_wind    <- 8

# steps (no number)

# Landscape layers names, besides vegetation, must follow this order:
ndvi   <- 1   # pi_ndvi[v] * (ndvi - optim[v]) ^ 2
north  <- 2   # slope-weighted
elev   <- 3   # standardized, but also used to compute slope
wdir   <- 4   # radians
wspeed <- 5   # m/s


# Raster-matrix-array conversion functions --------------------------------

# Make a 3D array from landscape SpatRaster.
land_cube <- function(x) {
  v <- values(x)
  a <- array(NA, dim = c(nrow(x), ncol(x), nlyr(x)),
             dimnames = list(row = NULL, col = NULL, layer = names(x)))
  for(l in 1:nlyr(x)) a[, , l] <- matrix(v[, l], nrow(x), ncol(x),
                                         byrow = TRUE)
  return(a)
}

# Function to turn matrix into SpatRaster (for plotting)
rast_from_mat <- function(m, fill_raster) { # fill_raster is a SpatRaster from terra
  mt <- t(m)
  for(i in 1:nrow(m)) mt[, i] <- m[i, ]
  r <- fill_raster[[1]]
  values(r) <- as.numeric(mt)
  return(r)
}

# Spread functions --------------------------------------------------------

#' @title spread_one_cell_r
#' @description Calculates the probability of a cell spreading fire to another.
#' @return float [0, 1] indicating the probability.
#'
#' @param int vegetation: vegetation class of the target cell, to obtain the
#'   corresponding intercept from coef.
#' @param arma::frowvec landscape_burning: data from burning cell.
#' @param arma::frowvec landscape_neighbour: data from target neighbour.
#' @param arma::frowvec coef: logistic regression parameters.
#' @param int position: relative position of the neighbour in relation to the
#' burning cell. The eight neighbours are labelled from 0 to 7 beginning from
#' the upper-left one (by row):
#'   0 1 2
#'   3   4
#'   5 6 7.
#'   This is necessary to compute the slope and wind effects, as they
#'   depend on the angle and distance between burning and target pixels.
#' @param float upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires; 0.5 is preferred).

spread_one_cell_r <- function(
    vegetation,
    landscape_burning,
    landscape_neighbour,
    coef,
    position,
    upper_limit = 1.0
  ) {

  # wind term
  wdir_term = cos(angles[position] - landscape_burning[wdir])

  # slope term (from elevation and distance)
  slope_term = sin(atan(
    (landscape_neighbour[elev] - landscape_burning[elev]) / distances[position]
  ))

  # compute linear predictor
  linpred <-
    coef[vegetation + 1] + # vegetation classes start at zero
    coef[b_ndvi] * landscape_neighbour[ndvi] +
    coef[b_north] * landscape_neighbour[north] +
    coef[b_elev] * landscape_neighbour[elev] +
    coef[b_slope] * slope_term +
    coef[b_wind] * wdir_term * landscape_burning[wspeed]

  # burn probability
  probs <- plogis(linpred) * upper_limit

  burn <- rbinom(1, size = 1, prob = probs)

  # return both the probability and the burn to check with the cpp function
  result <- c(probs, burn)
  names(result) <- c("probs", "burn")

  return(result)
}

# .........................................................................

#' @title simulate_fire_r
#' @description function to simulate a fire spread given the landscape,
#'   model coefficients and ignition points.
#' @return burned_res: struct containing the following objects:
#'   IntegerMatrix burned_bin, a binary matrix indicating the burned pixels
#'   IntegerMatrix burned_ids, a matrix with a column by burned pixel,
#'     indicating its row (row1) and column (row2) in the landscape,
#'   int end, the number of burned pixels.

#' @param SpatRaster[terra] landscape: predictor variables or variables used to
#'   compute the actual predictors.
#'     veg: vegetation class, in {0, 1, 2}. 99 is non-burnable. This layer is
#'       separated within the function.
#'     ndvi: ndvi-quadratic terms, computed as pi[v] * (ndvi - optim[v]) ^ 2.
#'       pi[forest] = 1, while the others are lower. pi and optim are
#'       estimated at
#'       <fire_spread/data/NDVI_regional_data/ndvi_optim_and_proportion.rds>.
#'     north: cos(aspect), weighted by slope steepness.
#'     elev: elevation (standardized), used to compute directional slope
#'       effect, but also included as a main effect.
#'     wdir: direction from where the wind blows (°),
#'     wspeed: wind speed (m/s).
#' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
#'   the cell(s) where the fire begun. First row has the row_id, second row has
#'   the col_id.
#' @param arma::frowvec coef: parameters in logistic regression to compute the
#'   spread probability as a function of covariates.
#' @param float upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires).
#' @param int steps: maximum number of simulation steps allowed. If 0
#'   (the default), a large number is used to avoid limiting fire spread.
#'   The burning of ignition points is considered the first step, so steps = 1
#'   burns only the ignition points. Bear in mind that the simulation may stop
#'   because there are no more burning cells, without reaching the maximum steps
#'   allowed.
                                                                                                                                    #'   testing.)
#' @param bool plot_animation: whether to plot the fire progress while running
#'   or not (set to FALSE by default).

simulate_fire_r <- function(
    landscape,
    ignition_cells,
    coef,
    upper_limit = 1.0,
    steps = 0,
    plot_animation = FALSE
  ) {

  # define landscape dimensions
  n_row <- nrow(landscape)
  n_col <- ncol(landscape)
  n_cell <- n_row * n_col
  n_layers <- nlyr(landscape)

  # non-limited steps
  if(steps == 0) steps <- n_cell * 10

  # turn landscape into numeric array
  temp <- land_cube(landscape)
  landscape_arr <- temp[, , -1] # without vegetation
  vegetation <- temp[, , 1]

  # Create burn layer, which will be exported.
  burned_bin <- matrix(0, n_row, n_col)

  # Make burning_ids matrix
  burning_ids <- matrix(NA, 2, n_cell)

  # Initialize burning_ids and burned_bin
  for(i in 1:ncol(ignition_cells)) {
    burning_ids[1, i] <- ignition_cells[1, i]
    burning_ids[2, i] <- ignition_cells[2, i]

    burned_bin[ignition_cells[1, i], ignition_cells[2, i]] <- 1
  }

  # positions from where to read the ids of the currently burning cells
  start <- 1
  end <- ncol(ignition_cells)

  # Fire raster for plotting
  burn_raster <- landscape[[1]]
  values(burn_raster) <- 0

  # get burning cells ids
  burning_cells <- cellFromRowCol(burn_raster,
                                  row = burning_ids[1, start:end],
                                  col = burning_ids[2, start:end])

  values(burn_raster)[burning_cells] <- 1
  if (plot_animation) {
    plot(burn_raster, col = c("green", "red"), main = "step 1")
  }

  # spread
  step <- 1
  burning_size <- length(burning_cells)

  while(burning_size > 0 & step < steps) {
    # update step
    step <- step + 1

    # Loop over all the burning cells to burn their neighbours. Use end_forward
    # to update the last position in burning_ids within this loop, without
    # compromising the loop's integrity.
    end_forward <- end

    # Loop over burning cells in the step

    # b is going to keep the position in burning_ids that have to be evaluated
    # in this burn step

    # spread from burning pixels
    for(b in start:end) {
      # Get burning_cells' data
      landscape_burning <- landscape_arr[burning_ids[1, b], burning_ids[2, b], ];

      # get neighbours (adjacent computation here)
      neighbours <- burning_ids[, b] + moves

      # Loop over neighbours of the focal burning cell
      for(n in 1:8) {

        # Is the cell in range?
        out_of_range <- (
          (neighbours[1, n] < 1) | (neighbours[1, n] > n_row) | # check rows
          (neighbours[2, n] < 1) | (neighbours[2, n] > n_col)   # check cols
        )
        if(out_of_range) next # (jumps to next iteration if TRUE)

        # Extract target vegetation
        veg_target <- vegetation[neighbours[1, n], neighbours[2, n]]

        # Is the cell burnable?
        burnable_cell <-
          (burned_bin[neighbours[1, n], neighbours[2, n]] == 0) &
          (veg_target != 99)
        if(!burnable_cell) next

        # obtain data from the neighbour
        landscape_neighbour = landscape_arr[neighbours[1, n], neighbours[2, n], ];

        # simulate fire
        burn <- spread_one_cell_r(
          veg_target,
          landscape_burning,
          landscape_neighbour,
          coef,
          n, # neighbour position (in 1:8)
          upper_limit
        )["burn"] # because it returns also the probability

        if(burn == 0) next

        # If burned,
        # store id of recently burned cell and
        # set 1 to burned_bin
        # (but advance end_forward first)
        end_forward <- end_forward + 1
        burning_ids[1, end_forward] = neighbours[1, n]
        burning_ids[2, end_forward] = neighbours[2, n]
        burned_bin[neighbours[1, n], neighbours[2, n]] <- 1
      } # end loop over neighbours of burning cell b

    } # end loop over burning cells from this step

    # update start and end
    start <- end + 1
    end <- end_forward
    burning_size <- end - start + 1

    # update: burning to burned
    values(burn_raster)[burning_cells] <- 2

    if(burning_size > 0) {
      # update burning_cells (this correspond to the next step)
      burning_cells <- cellFromRowCol(burn_raster,
                                      row = burning_ids[1, start:end],
                                      col = burning_ids[2, start:end])
      values(burn_raster)[burning_cells] <- 1
    }

    if (plot_animation) {
      plot(burn_raster, col = c("green", "red", "black"), main = paste("step", step))
    }

  }

  return(burned_bin)
}

# .........................................................................
# The same function but deterministic, to test if the discrepancy between R and
# cpp is caused by seed problems

simulate_fire_deterministic_r <- function(
    landscape,
    ignition_cells,
    coef,
    upper_limit = 1.0,
    steps = 0,
    plot_animation = FALSE
) {

  # define landscape dimensions
  n_row <- nrow(landscape)
  n_col <- ncol(landscape)
  n_cell <- n_row * n_col
  n_layers <- nlyr(landscape)

  # non-limited steps
  if(steps == 0) steps <- n_cell * 10

  # turn landscape into numeric array
  temp <- land_cube(landscape)
  landscape_arr <- temp[, , -1] # without vegetation
  vegetation <- temp[, , 1]

  # Create burn layer, which will be exported.
  burned_bin = matrix(0, n_row, n_col)

  # Make burning_ids matrix
  burning_ids <- matrix(NA, 2, n_cell)

  # Initialize burning_ids and burned_bin
  for(i in 1:ncol(ignition_cells)) {
    burning_ids[1, i] <- ignition_cells[1, i]
    burning_ids[2, i] <- ignition_cells[2, i]

    burned_bin[ignition_cells[1, i], ignition_cells[2, i]] <- 1
  }

  # positions from where to read the ids of the currently burning cells
  start <- 1
  end <- ncol(ignition_cells)

  # Fire raster for plotting
  burn_raster <- landscape[[1]]
  values(burn_raster) <- 0

  # get burning cells ids
  burning_cells <- cellFromRowCol(burn_raster,
                                  row = burning_ids[1, start:end],
                                  col = burning_ids[2, start:end])

  values(burn_raster)[burning_cells] <- 1
  if (plot_animation) {
    plot(burn_raster, col = c("green", "red"), main = "step 1")
  }

  # spread
  step <- 1
  burning_size <- length(burning_cells)

  while(burning_size > 0 & step < steps) {
    # update burn step
    step <- step + 1

    # Loop over all the burning cells to burn their neighbours. Use end_forward
    # to update the last position in burning_ids within this loop, without
    # compromising the loop's integrity.
    end_forward <- end

    # Loop over burning cells in the step

    # b is going to keep the position in burning_ids that have to be evaluated
    # in this burn step

    # spread from burning pixels
    for(b in start:end) {
      # Get burning_cells' data
      landscape_burning <- landscape_arr[burning_ids[1, b], burning_ids[2, b], ];

      # get neighbours (adjacent computation here)
      neighbours <- burning_ids[, b] + moves

      # Loop over neighbours of the focal burning cell
      for(n in 1:8) {

        # Is the cell in range?
        out_of_range <- (
          (neighbours[1, n] < 1) | (neighbours[1, n] > n_row) | # check rows
          (neighbours[2, n] < 1) | (neighbours[2, n] > n_col)   # check cols
        )
        if(out_of_range) next # (jumps to next iteration if TRUE)

        # Extract target vegetation
        veg_target <- vegetation[neighbours[1, n], neighbours[2, n]]

        # Is the cell burnable?
        burnable_cell <-
          (burned_bin[neighbours[1, n], neighbours[2, n]] == 0) &
          (veg_target != 99)
        if(!burnable_cell) next

        # obtain data from the neighbour
        landscape_neighbour = landscape_arr[neighbours[1, n], neighbours[2, n], ];

        # simulate fire
        burn <- spread_one_cell_r(
          veg_target,
          landscape_burning,
          landscape_neighbour,
          coef,
          n, # neighbour position (in 1:8)
          upper_limit
        )

        # make deterministic!
        if(burn["probs"] < 0.5000000000) next

        # If burned,
        # store id of recently burned cell and
        # set 1 to burned_bin
        # (but advance end_forward first)
        end_forward <- end_forward + 1
        burning_ids[1, end_forward] = neighbours[1, n]
        burning_ids[2, end_forward] = neighbours[2, n]
        burned_bin[neighbours[1, n], neighbours[2, n]] <- 1
      } # end loop over neighbours of burning cell b

    } # end loop over burning cells from this step

    # update start and end
    start <- end + 1
    end <- end_forward
    burning_size <- end - start + 1

    # update: burning to burned
    values(burn_raster)[burning_cells] <- 2

    if(burning_size > 0) {
      # update burning_cells (this correspond to the next step)
      burning_cells <- cellFromRowCol(burn_raster,
                                      row = burning_ids[1, start:end],
                                      col = burning_ids[2, start:end])
      values(burn_raster)[burning_cells] <- 1
    }

    if (plot_animation) {
      plot(burn_raster, col = c("green", "red", "black"), main = paste("step", step))
    }

  }

  return(burned_bin)
}