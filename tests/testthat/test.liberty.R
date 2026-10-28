default_args <- list(
  D = 40, xu = 0.045, thick = 1.6, baseline = 0.0006, element = 2,
  c_factor = 200, w_factor = 100, l_factor = 40, p_factor = 1
)

test_that("Liberty model produces sensible results", {
  for (fun in list(liberty, liberty_dawson)) {
    result <- do.call(fun, default_args)
    for (x in c("reflectance", "transmittance")) {
      expect_true(all(result[[x]] >= 0))
      expect_true(all(result[[x]] <= 1))
    }
    expect_true(all(result$reflectance + result$transmittance <= 1))
  }
})

test_that("LIBERTY interface reflectances match Dawson et al. Eqs. 1 and 2", {
  fresnel <- function(theta, n_from, n_to) {
    out <- asin(n_from / n_to * sin(theta))
    0.5 * (sin(theta - out)^2 / sin(theta + out)^2 +
             tan(theta - out)^2 / tan(theta + out)^2)
  }
  weight <- function(theta) 2 * sin(theta) * cos(theta)
  for (n in c(1.4891, 1.4491, 1.3291, 1.2796)) {
    critical <- asin(1 / n)
    me <- integrate(function(t) fresnel(t, 1, n) * weight(t), 0, pi / 2,
                    rel.tol = 1e-12)$value
    mi <- (1 - sin(critical)^2) +
      integrate(function(t) fresnel(t, n, 1) * weight(t), 0, critical,
                rel.tol = 1e-12)$value
    expect_equal(1 - rrtm:::tav_abs(90, n), me, tolerance = 1e-9)
    expect_equal(1 - rrtm:::tav_abs(90, n) / n^2, mi, tolerance = 1e-9)
  }
})

test_that("LIBERTY infinite reflectance is the sum of Dawson et al. Eq. 15", {
  result <- do.call(liberty, default_args)
  i <- c(11, 57, 81, 161, 251, 361)  # 450, 680, 800, 1200, 1650, 2200 nm
  n <- 1.4891 - 0.0005 * (i - 1)
  me <- 1 - rrtm:::tav_abs(90, n)
  mi <- 1 - rrtm:::tav_abs(90, n) / n^2
  k <- with(default_args, D * (
    baseline + unname(as.matrix(rrtm:::dataspec_liberty[i, ])) %*%
      c(c_factor, w_factor, element, l_factor, p_factor)
  ))[, 1]
  M <- 2 * (1 - (k + 1) * exp(-k)) / k^2
  tau <- (1 - mi) * M / (1 - mi * M)
  x <- default_args$xu / (1 - (1 - 2 * default_args$xu) * tau)
  R <- result$RR[i]
  q <- (1 - me) * tau * R / (1 - me * R)
  eq15 <- 2 * x * me
  for (j in 0:5000) {
    eq15 <- eq15 + x * (1 - x)^j * (1 - 2 * x * me) * tau * q^j
  }
  expect_equal(eq15, R, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(liberty_dawson(
    40, 0.045, 1.6, 0.0006, 2, 200, 100, 40, 1
  )$RR[i], R, tolerance = 1e-6)))
})

test_that("LIBERTY is defined for bright leaves where Eq. 16 has no solution", {
  bright <- modifyList(default_args, list(xu = 0.1, thick = 2.5, baseline = 0.0003, element = 1))
  result <- do.call(liberty, bright)
  expect_true(all(is.finite(result$reflectance) & is.finite(result$transmittance)))
  expect_true(all(result$reflectance > 0 & result$transmittance > 0))
  expect_true(all(result$reflectance + result$transmittance <= 1))
  expect_true(all(result$RR <= 1))
  dawson <- suppressWarnings(do.call(liberty_dawson, bright))
  expect_false(all(is.finite(dawson$reflectance)))
})

test_that("LIBERTY loses no light without absorption", {
  # concentrations at 0; the baseline is kept small but nonzero because
  # the slab transmission is evaluated without a series for tiny optical depths
  none <- modifyList(default_args, list(
    baseline = 1e-5, element = 0, c_factor = 0, w_factor = 0, l_factor = 0, p_factor = 0
  ))
  result <- do.call(liberty, none)
  expect_true(all(result$RR > 0.95))
  expect_true(all(result$reflectance + result$transmittance > 0.995))
})
