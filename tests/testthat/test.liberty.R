test_that("Liberty model produces sensible results", {
  D <- 40
  xu <- 0.045
  thick <- 1.6
  baseline <- 0.0006
  element <- 2
  c_factor <- 200
  w_factor <- 100
  l_factor <- 40
  p_factor <- 1

  result <- liberty(
    D, xu, thick, baseline, element,
    c_factor, w_factor, l_factor, p_factor
  )
  for (x in c("reflectance", "transmittance")) {
    expect_true(all(result[[x]] >= 0))
    expect_true(all(result[[x]] <= 1))
  }
})

