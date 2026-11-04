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

test_that("RSR weights are normalized: band values are weighted averages", {
  flat <- cbind(zero = 0, one = 1, mid = 0.3)[rep(1, 2101), ]
  for (sensor in sensors) {
    out <- spectral_response(cbind(specs, flat), sensor)
    expect_true(all(out >= 0 & out <= 1), info = sensor)
    expect_equal(unname(out[, "zero"]), rep(0, nrow(out)), info = sensor)
    expect_equal(unname(out[, "one"]), rep(1, nrow(out)), info = sensor)
    expect_equal(unname(out[, "mid"]), rep(0.3, nrow(out)), info = sensor)
  }
})

test_that("Landsat 8 bands match naive averages over the OLI band intervals", {
  l8 <- spectral_response(specs, "landsat8")
  # Nominal Landsat 8 OLI band edges (nm), bands 1-7
  edges <- rbind(
    aerosol = c(433, 453), blue = c(450, 515), green = c(525, 600),
    red = c(630, 680), nir = c(845, 885), swir1 = c(1560, 1660),
    swir2 = c(2100, 2300)
  )
  naive <- t(apply(edges, 1, function(e) {
    colMeans(specs[wl >= e[1] & wl <= e[2], , drop = FALSE])
  }))
  expect_true(all(abs(l8 - naive) < 0.01))

  # Leaf reflectance shape: low blue, higher green, low red, high NIR
  refl <- l8[, "reflectance"]
  expect_lt(refl[["L8-2-blue"]], refl[["L8-3-green"]])
  expect_lt(refl[["L8-4-red"]], refl[["L8-3-green"]])
  expect_lt(refl[["L8-3-green"]], refl[["L8-5-nir"]])
  expect_gt(refl[["L8-5-nir"]], 0.3)
})
