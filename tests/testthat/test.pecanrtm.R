context("Results match PEcAnRTM")

skip_if_not_installed("PEcAnRTM")

vec_equal <- function(x, y) {
  expect_equal(sum((x - y) ^ 2), 0)
}

test_that("tav calculation matches", {
  pecan_90 <- .Fortran("tav_abs", 90, refractive_p45, numeric(2101),
                       PACKAGE = "PEcAnRTM")[[3]]
  rrtm_90 <- tav_abs(90, refractive_p45)
  vec_equal(pecan_90, rrtm_90)
  pecan_40 <- .Fortran("tav_abs", 40, refractive_p45, numeric(2101),
                       PACKAGE = "PEcAnRTM")[[3]]
  rrtm_40 <- tav_abs(40, refractive_p45)
  vec_equal(pecan_40, rrtm_40)
})

test_that("gpm calculation matches", {
  params <- c(1.4, 40, 0.01, 0.01)
  k <- (dataspec_p4 %*% params[-1]) / params[[1]]
  pecan <- .Fortran("gpm", k, refractive_pd, params[[1]], matrix(0, 2100, 2),
                    PACKAGE = "PEcAnRTM")[[4]]
  rrtm <- gpm(k, refractive_pd, params[[1]])
  ## rrtm <- gpm(k, refractive_p45, params[[1]])
  refractive_p45[i]
})

test_that("PROSPECT 4 results match", {
  params <- c(1.4, 40, 0.01, 0.01)
  pecan_4 <- PEcAnRTM::prospect(params, "4")
  rrtm_4 <- do.call(prospect4, as.list(params))
})
