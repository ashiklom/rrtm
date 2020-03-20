test_that("Volume-scattering function works across angles", {

  angles <- seq(-60, 60, 10)
  for (sa in angles) {
    for (oa in angles) {
      for (aa in angles) {
        s <- volscatt_scalar_v(sa, oa, aa, l)
        v <- volscatt(sa, oa, aa, l)
        expect_identical(unlist(s["chi_s",]), v$chi_s)
        expect_identical(unlist(s["chi_o",]), v$chi_o)
        expect_identical(unlist(s["frho",]), v$frho)
        expect_identical(unlist(s["ftau",]), v$ftau)
      }
    }
  }

})
