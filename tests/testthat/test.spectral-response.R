context("Spectral response functions work")

wl <- seq(400, 2500)
p5 <- prospect5(1.4, 40, 10, 0, 0.01, 0.01)
spec <- p5[["reflectance"]]
specs <- cbind(reflectance = p5[["reflectance"]],
               transmittance = p5[["transmittance"]])
sensors <- setdiff(sensor_list, "identity")

test_that("identity sensor returns the input spectra unchanged, as columns", {
  idx <- c(1, 50, 300, 1000, 2101)
  expect_equal(spectral_response(spec, "identity"), matrix(spec, ncol = 1))
  expect_equal(
    spectral_response(spec[idx], "identity"),
    matrix(spec[idx], ncol = 1)
  )
  expect_equal(spectral_response(specs, "identity"), specs)
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

test_that("real sensors require a full 400:2500 nm spectrum", {
  expect_error(spectral_response(spec[1:100], "landsat8"), "2101 rows")
})

test_that("multiple spectra are returned as columns, bands as rows", {
  for (sensor in sensors) {
    out <- spectral_response(specs, sensor)
    nband <- ncol(sensor.rsr[[sensor]]) - 1
    expect_equal(dim(out), c(nband, 2), info = sensor)
    expect_equal(colnames(out), colnames(specs), info = sensor)
    expect_equal(rownames(out), colnames(sensor.rsr[[sensor]])[-1], info = sensor)
    # Each column is the same as convolving that spectrum alone
    expect_equal(unname(out[, 1]), unname(spectral_response(specs[, 1], sensor)[, 1]),
                 info = sensor)
    expect_equal(unname(out[, 2]), unname(spectral_response(specs[, 2], sensor)[, 1]),
                 info = sensor)
  }
})

test_that("spectra convolution for landsat 8 works as expected", {
  l8 <- spectral_response(spec, "landsat8")

  expect_equal(nrow(l8), 7)
  expect_lt(l8["L8-1-aerosol", ], l8["L8-2-blue", ])
  expect_lt(l8["L8-2-blue", ], l8["L8-3-green", ])
  expect_lt(l8["L8-4-red", ], l8["L8-5-nir", ])
})
