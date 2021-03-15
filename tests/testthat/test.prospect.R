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
