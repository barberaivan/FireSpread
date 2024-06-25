library(terra)
library(testthat)

source("R_spread_functions.R")

# create data for testing

n_veg <- 3 # vegetation types
n_terrain <- 3
n_layers <- 1 + n_terrain
layer_names <- c("vegetation", "elev", "wdir", "wspeed")

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
landscape$elev <- rnorm(ncell(landscape), 1500, 300)
landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

set.seed(2343)
cv <- rnorm(n_veg)
ct <- rnorm(n_terrain - 1)

test_that("Fire spread functions", {
  set.seed(30)
  fire_r <- simulate_fire_r(
    landscape = landscape, # use the SpatRaster, with vegetation and terrain
    coef_veg = cv,
    coef_terrain = ct,
    ignition_cells = ig_location,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_cpp <- simulate_fire(
    vegetation = land_cube(landscape)[, , 1],
    terrain = land_cube(landscape)[, , -1], # use the array
    coef_veg = cv,
    coef_terrain = ct,
    ignition_cells = ig_location - 1,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_compare_cpp <- simulate_fire_compare(
    vegetation = land_cube(landscape)[, , 1],
    terrain = land_cube(landscape)[, , -1], # use the array
    coef_veg = cv,
    coef_terrain = ct,
    ignition_cells = ig_location - 1,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
  expect_equal(fire_r, fire_compare_cpp$burned_layer)
})

test_that("Deterministic fire spread functions", {
  fire_r <- simulate_fire_deterministic_r(
    landscape = landscape, # use the SpatRaster
    coef_veg = cv,
    coef_terrain = ct,
    ignition_cells = ig_location,
    upper_limit = 1.0
  )

  fire_cpp <- simulate_fire_deterministic(
    vegetation = land_cube(landscape)[, , 1],
    terrain = land_cube(landscape)[, , -1], # use the array
    coef_veg = cv,
    coef_terrain = ct,
    ignition_cells = ig_location - 1,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
})