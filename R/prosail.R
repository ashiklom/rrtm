#' PRO4SAIL: PROSPECT leaf model coupled to 4SAIL canopy RTM and Hapke soil
#' model
#'
#' @inheritParams prospectd
#' @inheritParams hapke_soil
#' @inheritParams foursail
#' @param ... Additional arguments to [foursail]
#' @inherit foursail return
#' @author Alexey Shiklomanov
#' @export
pro4sail_4 <- function(N, Cab, Cw, Cm, LAI, soil_moisture,
                       ...) {
  lrt <- prospect4(N = N, Cab = Cab, Cw = Cw, Cm = Cm)
  rsoil <- hapke_soil(soil_moisture)
  foursail(lrt$reflectance, lrt$transmittance, rsoil, LAI, ...)
}

#' @rdname pro4sail_4
#' @export
pro4sail_5 <- function(N, Cab, Car, Cbrown, Cw, Cm, LAI, soil_moisture,
                       ...) {
  lrt <- prospect5(N = N, Cab = Cab, Car = Car, Cbrown = Cbrown,
                   Cw = Cw, Cm = Cm)
  rsoil <- hapke_soil(soil_moisture)
  foursail(lrt$reflectance, lrt$transmittance, rsoil, LAI, ...)
}

#' @rdname pro4sail_4
#' @export
pro4sail_d <- function(N, Cab, Car, Canth, Cbrown, Cw, Cm,
                       LAI, soil_moisture,
                       ...) {
  lrt <- prospectd(N = N, Cab = Cab, Car = Car, Canth = Canth,
                   Cbrown = Cbrown, Cw = Cw, Cm = Cm)
  rsoil <- hapke_soil(soil_moisture)
  foursail(lrt$reflectance, lrt$transmittance, rsoil, LAI, ...)
}
