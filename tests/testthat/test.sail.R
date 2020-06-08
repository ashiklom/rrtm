test_that("Volume-scattering function works across angles", {
  angles <- seq(-60, 60, 20)
  l <- angles
  for (sa in angles) {
    for (oa in angles) {
      for (aa in angles) {
        angle_info <- sprintf(paste(
          "solar angle: %f",
          "instrument angle: %f",
          "azimuth angle: %f",
          sep = "\n"
        ), sa, oa, aa)
        s <- volscatt_scalar_v(sa, oa, aa, l)
        v <- volscatt(sa, oa, aa, l)
        expect_identical(unlist(s["chi_s",]), v$chi_s, info = angle_info)
        expect_identical(unlist(s["chi_o",]), v$chi_o, info = angle_info)
        expect_identical(unlist(s["frho",]), v$frho, info = angle_info)
        expect_identical(unlist(s["ftau",]), v$ftau, info = angle_info)
      }
    }
  }
})

test_that("Manual PRO4SAIL works", {
  lrt <- prospect4(1.4, 40, 0.01, 0.01)
  rsoil <- hapke_soil(0.5)
  sail <- foursail(lrt$reflectance, lrt$transmittance, rsoil, 3)
  for (stream in sail) {
    expect_true(all(stream >= 0))
    expect_true(all(stream <= 1))
  }
})

test_that("PRO4SAIL shortcuts work", {
  args <- list(
    N = 1.4, Cab = 40, Car = 8, Canth = 10, Cbrown = 0,
    Cw = 0.01, Cm = 0.01,
    LAI = 3, soil_moisture = 0.5
  )
  for (fun in c(pro4sail_4, pro4sail_5, pro4sail_d)) {
    fun <- pro4sail_5
    cargs <- args[head(names(formals(fun)), -1)]
    sail <- do.call(fun, cargs)
    for (stream in sail) {
      expect_true(all(stream >= 0))
      expect_true(all(stream <= 1))
    }
  }
})
