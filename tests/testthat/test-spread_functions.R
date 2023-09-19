library(terra)
library(testthat)

source("R_spread_functions.R")

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


test_that("Fire spread functions", {
  set.seed(30)
  fire_r <- simulate_fire_r(
    landscape = landscape, # use the SpatRaster
    burnable = burnable,
    ignition_cells = ig_location,
    coef = coefs,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_cpp <- simulate_fire(
    landscape = land_cube(landscape), # use the array
    burnable = burnable,
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_compare_cpp <- simulate_fire_compare(
    landscape = land_cube(landscape), # use the array
    burnable = burnable,
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
  expect_equal(fire_r, fire_compare_cpp$burned_layer)
})

test_that("Deterministic fire spread functions", {
  fire_r <- simulate_fire_deterministic_r(
    landscape = landscape, # use the SpatRaster
    burnable = burnable,
    ignition_cells = ig_location,
    coef = coefs,
    upper_limit = 1.0
  )

  fire_cpp <- simulate_fire_deterministic(
    landscape = land_cube(landscape), # use the array
    burnable = burnable,
    ignition_cells = ig_location - 1,
    coef = coefs,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
})