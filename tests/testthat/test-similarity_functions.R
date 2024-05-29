library(terra)

source("R_spread_functions.R")
source("R_similarity_functions.R")

test_that("Similarity functions", {
  # create data for testing

  n_veg <- 3 # vegetation types
  n_coef <- n_veg + 5 # it doesn't include steps
  n_layers <- n_coef - n_veg + 1
  layer_names <- c("vegetation", "ndvi", "north", "elev", "wdir", "wspeed")

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
  landscape$ndvi <- rnorm(ncell(landscape))
  landscape$north <- runif(ncell(landscape))
  landscape$elev <- rnorm(ncell(landscape))
  landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
  landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

  ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

  set.seed(2343)
  coefs <- rnorm(n_coef)

  # simulate fires
  set.seed(1)
  fire_1 <- simulate_fire_compare_veg(
    landscape = land_cube(landscape)[, , -1], # use the array
    vegetation = land_cube(landscape)[, , 1],
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0,
    n_veg_types = n_veg
  )

  set.seed(1)
  fire_1_ <- simulate_fire_compare_veg(
    landscape = land_cube(landscape)[, , -1], # use the array
    vegetation = land_cube(landscape)[, , 1],
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0,
    n_veg_types = n_veg
  )

  set.seed(2)
  fire_2 <- simulate_fire_compare_veg(
    landscape = land_cube(landscape)[, , -1], # use the array
    vegetation = land_cube(landscape)[, , 1],
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0,
    n_veg_types = n_veg
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