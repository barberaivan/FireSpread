library(terra)

source("R_spread_functions.R")
source("R_similarity_functions.R")

test_that("Similarity functions", {
  # create data for testing

  n_veg <- 5
  n_nd <- 2      # non-directional terms
  n_terrain <- 3 # it doesn't include steps
  n_b_terrain <- n_terrain - 1 # n for terrain coefficients
  n_layers <- 1 + n_nd + n_terrain # extra 1 is the intercept
  layer_names <- c("vegetation", "vfi", "tfi", "elev", "wdir", "wspeed")

  # landscape raster
  size <- 30
  n_rows <- size
  n_cols <- size
  res <- 30
  n_cells <- n_rows * n_cols

  landscape <- rast(
    ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
    nlyrs = n_layers,
    xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
    names = layer_names
  )

  # fill data
  landscape$vegetation <- sample(0:(n_veg-1), n_rows * n_cols, replace = TRUE)
  landscape$vfi <- rep(0, ncell(landscape))
  landscape$tfi <- rep(0, ncell(landscape))
  landscape$elev <- rnorm(ncell(landscape), 1500, 300)
  landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
  landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

  ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

  set.seed(2343)
  b_int <- rnorm(1)
  b_nd <- rnorm(n_nd)
  b_terrain <- rnorm(n_b_terrain)

  # make arma::cube landscape
  land <- land_cube(landscape)

  # simulate fires
  set.seed(1)
  fire_1 <- simulate_fire_compare_veg(
    layer_vegetation = land[, , 1],
    layer_nd = land[, , 2:3],
    layer_terrain = land[, , 4:6],
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location-1,
    upper_limit = 1.0,
    steps = 0,
    n_veg = 5
  )

  set.seed(1)
  fire_1_ <- simulate_fire_compare_veg(
    layer_vegetation = land[, , 1],
    layer_nd = land[, , 2:3],
    layer_terrain = land[, , 4:6],
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location-1,
    upper_limit = 1.0,
    steps = 0,
    n_veg = 5
  )

  set.seed(2)
  fire_2 <- simulate_fire_compare_veg(
    layer_vegetation = land[, , 1],
    layer_nd = land[, , 2:3],
    layer_terrain = land[, , 4:6],
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location-1,
    upper_limit = 1.0,
    steps = 0,
    n_veg = 5
  )

  # a fire against itself, simulated with the same seed
  expect_equal(fire_1, fire_1_)

  # different fires, r and cpp functions
  similarity_cpp_1_2 <- compare_fires_try(fire_1, fire_2)
  similarity_r_1_2 <- compare_fires_r(fire_1, fire_2)

  expect_equal(similarity_cpp_1_2, similarity_r_1_2, tolerance = 1e-6)

  # a fire against itself, cpp function
  similarity_1_1 <- compare_fires_try(fire_1, fire_1)
  expect_equal(unname(similarity_1_1), rep(1, length(similarity_1_1)))

  # compare overlap sp
  ov_cpp_1_2 <- overlap_spatial(fire_1, fire_2)
  ov_r_1_2 <- overlap_spatial_r(fire_1, fire_2)
  expect_equal(ov_cpp_1_2, ov_r_1_2, tolerance = 1e-6)
  expect_equal(overlap_spatial(fire_1, fire_1),
               overlap_spatial_r(fire_1, fire_1),
               1,
               tolerance = 1e-6)
})