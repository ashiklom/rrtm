context("PROSPECT model works")

test_that("PROSPECT 4", {
  x <- do.call(cbind, prospect4(1.4, 40, 0.01, 0.01))
  expect_true(all(x > 0))
  expect_true(all(x < 1))
})

test_that("PROSPECT 5", {
  p5 <- do.call(cbind, prospect5(1.4, 40, 10, 0, 0.01, 0.01))
  expect_true(all(p5 > 0))
  expect_true(all(p5 < 1))
  p5b <- do.call(cbind, prospect5(1.4, 40, 10, 0, 0.01, 0.01))
  expect_true(all(p5b > 0))
  expect_true(all(p5b < 1))
})

test_that("PROSPECT D", {
  pd <- do.call(cbind, prospectd(1.4, 40, 10, 8, 0, 0.01, 0.01))
  expect_true(all(pd > 0))
  expect_true(all(pd < 1))
})

test_that("PROSPECT with extremely high absorptivity works but throws warning", {
  for (Cab in c(seq(100, 1000, 100))) {
    expect_warning(
      p4 <- do.call(cbind, prospect4(1.4, Cab, 0.01, 0.01)),
      "high absorptivity values",
      all = TRUE
    )
    expect_true(all(p4 >= 0))
    expect_true(all(p4 <= 1))
  }
})

test_that("PROSPECT D interface transmissivities use the PROSPECT D refractive index", {
  expect_equal(pd_t21, pd_t12 / refractive_pd^2)
  expect_equal(pd_t12, tav_abs(90, refractive_pd))
  expect_equal(pd_talf, tav_abs(40, refractive_pd))
})

test_that("PROSPECT D matches independent reference values", {
  skip_if_not_installed("expint")
  # Reference values from the `prosail` Python package (J. Gomez-Dans),
  # PROSPECT-D with exact E1, at 450, 550, 680, 800, 1450, 1940, 2200 nm.
  wl <- c(450, 550, 680, 800, 1450, 1940, 2200) - 399
  ref_r <- c(0.0411123436, 0.1115443677, 0.0359676508, 0.4351475168,
             0.1650296676, 0.0373651153, 0.1547468978)
  ref_t <- c(0.0007884939, 0.1060743066, 0.0051708164, 0.4670787360,
             0.2096989879, 0.0489542338, 0.2531362617)
  pd <- prospectd(1.5, 40, 10, 0.1, 2, 0.01, 0.009, e1fun = expint::expint_E1)
  expect_equal(pd$reflectance[wl], ref_r, tolerance = 1e-8)
  expect_equal(pd$transmittance[wl], ref_t, tolerance = 1e-8)
})
