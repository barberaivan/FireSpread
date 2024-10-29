library(terra)
library(testthat)

source("R_spread_functions.R")

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

test_that("Fire spread functions", {
  set.seed(30)
  fire_r <- simulate_fire_r(
    landscape = landscape, # SpatRaster
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location,
    upper_limit = 1.0,
    steps = 0,
    plot_animation = F,
    det = F
  )

  set.seed(30)
  fire_cpp <- simulate_fire(
    layer_vegetation = land[, , 1],
    layer_nd = land[, , 2:3],
    layer_terrain = land[, , 4:6],
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location-1,
    upper_limit = 1.0,
    steps = 0
  )

  set.seed(30)
  fire_compare_cpp <- simulate_fire_compare(
    layer_vegetation = land[, , 1],
    layer_nd = land[, , 2:3],
    layer_terrain = land[, , 4:6],
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location-1,
    upper_limit = 1.0,
    steps = 0
  )

  expect_equal(fire_r, fire_cpp)
  expect_equal(fire_r, fire_compare_cpp$burned_layer)
})

test_that("Deterministic fire spread functions", {
  fire_r <- simulate_fire_r(
    landscape = landscape, # SpatRaster
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location,
    upper_limit = 1.0,
    steps = 0,
    plot_animation = F,
    det = T
  )

  fire_cpp <- simulate_fire_deterministic(
    layer_vegetation = land[, , 1],
    layer_nd = land[, , 2:3],
    layer_terrain = land[, , 4:6],
    coef_intercepts = b_int,
    coef_nd = b_nd,
    coef_terrain = b_terrain,
    ignition_cells = ig_location-1,
    upper_limit = 1.0,
    steps = 0
  )

  expect_equal(fire_r, fire_cpp)
})