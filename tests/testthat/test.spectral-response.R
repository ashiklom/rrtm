context("Spectral response functions work")

sensor_rsr_data <- get("sensor.rsr", envir = asNamespace("rrtm"))
p5 <- prospect5(1.4, 40, 10, 0, 0.01, 0.01)
spec <- p5[["reflectance"]]

test_that("identity sensor returns the input spectrum unchanged", {
  idx <- sample(2101, 50)
  expect_equal(spectral_response(spec, "identity"), matrix(spec, nrow = 1))
  expect_equal(
    spectral_response(spec[idx], "identity"),
    matrix(spec[idx], nrow = 1)
  )
})

test_that("sensor names are validated", {
  expect_error(spectral_response(spec, "not-a-sensor"))
})

test_that("sensor names are case-insensitive", {
  expect_equal(
    spectral_response(spec, "landsat8"),
    spectral_response(spec, "LANDSAT8")
  )
})

test_that("spectra convolution for landsat 8 works as expected", {
  l8 <- spectral_response(spec, "landsat8")

  expect_equal(ncol(l8), 7)
  expect_lt(l8[, "L8-1-aerosol"], l8[, "L8-2-blue"])
  expect_lt(l8[, "L8-2-blue"], l8[, "L8-3-green"])
  expect_lt(l8[, "L8-4-red"], l8[, "L8-5-nir"])
})
