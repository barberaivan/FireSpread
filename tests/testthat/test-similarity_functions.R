library(terra)

source("R_spread_functions.R")
source("R_similarity_functions.R")

test_that("Similarity functions", {
  # create data for testing
  n_coef <- 5

  set.seed(2345)
  coefs <- rnorm(n_coef, 0, 3)

  # landscape raster
  size <- 30
  n_rows <- size
  n_cols <- size
  res <- 30
  n_cells <- n_rows * n_cols

  landscape <- rast(
    ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
    nlyrs = n_coef, # no intercept, but wind uses two lawers
    xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
    names = c("vfi", "tfi", "elev", "wdir", "wspeed")
  )

  # fill data
  landscape$vfi <- rnorm(ncell(landscape))
  landscape$tfi <- rnorm(ncell(landscape))
  landscape$elev <- runif(ncell(landscape), 0, 2200)
  landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
  landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

  burnable <- matrix(rbinom(n_rows * n_cols, size = 1, prob = 0.15),
                     n_rows, byrow = T)

  # provide vegetation to compute all metrics
  n_veg_types <- 6
  vegetation <- matrix(sample(0:(n_veg_types-1), n_rows * n_cols, replace = TRUE),
                       n_rows, byrow = T)

  ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

  # simulate fires
  set.seed(1)
  fire_1 <- simulate_fire_compare_veg(
    landscape = land_cube(landscape), # use the array
    burnable = burnable,
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0,
    vegetation = vegetation,
    n_veg_types = n_veg_types
  )

  set.seed(1)
  fire_1_ <- simulate_fire_compare_veg(
    landscape = land_cube(landscape), # use the array
    burnable = burnable,
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0,
    vegetation = vegetation,
    n_veg_types = n_veg_types
  )

  set.seed(2)
  fire_2 <- simulate_fire_compare_veg(
    landscape = land_cube(landscape), # use the array
    burnable = burnable,
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0,
    vegetation = vegetation,
    n_veg_types = n_veg_types
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