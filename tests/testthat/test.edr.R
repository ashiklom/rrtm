context("EDR R implementation")

pft <- c(1, 1, 2)
lai <- c(1.5, 1, 0.5)
wai <- c(0, 0, 0)
Cab <- c(40, 50, 60)

npft <- length(pft)

cai <- rep(1, 3)
N <- rep(1.4, 3)
Car <- rep(10, 3)
Cw <- rep(0.01, 3)
Cm <- rep(0.01, 3)
orient_factor <- rep(0, 3)
clumping_factor <- rep(0.7, 3)
soil_moisture <- 0.8
direct_sky_frac <- 0.9
czen <- 0.9

wavelengths <- seq(400, 2500)

result <- edr_r(pft, lai, wai, cai,
                N, Cab, Car, Cw, Cm,
                orient_factor, clumping_factor,
                soil_moisture,
                direct_sky_frac, czen)

albedo <- result[["albedo"]]

test_that("EDR returns reasonable output", {
  expect_true(all(is.finite(albedo)))
  expect_true(all(albedo > 0))
  expect_true(all(albedo < 1))
})
