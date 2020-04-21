context("PROSPECT model works")

test_that("PROSPECT 4", {
  x <- do.call(cbind, prospect4(1.4, 40, 0.01, 0.01))
  expect_true(all(x > 0))
  expect_true(all(x < 1))
})

test_that("PROSPECT 5", {
  p5 <- do.call(cbind, prospect5(1.4, 40, 10, 0.01, 0.01))
  expect_true(all(p5 > 0))
  expect_true(all(p5 < 1))
  p5b <- do.call(cbind, prospect5(1.4, 40, 10, 0.01, 0.01, 8))
  expect_true(all(p5b > 0))
  expect_true(all(p5b < 1))
})

test_that("Multivariate PROSPECT", {
  skip("Temporarily disable broadcasting")
  x <- prospect4(1.4, c(40, 50, 60), 0.01, 0.01)
  expect_equal(dim(x), c(2101, 3, 2))
  expect_true(all(x > 0))
  expect_true(all(x < 1))
})
