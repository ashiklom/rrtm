test_that("Volume-scattering function works across angles", {
  angles <- seq(-60, 60, 10)
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
