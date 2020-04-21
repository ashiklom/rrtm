#' PROSPECT leaf spectra model
#'
#' @param N Effective number of leaf layers
#' @param Cab Chlorophyll content (ug cm -2)
#' @param Car Carotenoid content (ug cm -2)
#' @param Canth Anthocyanin content (ug cm-2)
#' @param Cbrown "Brown matter" content (unitless)
#' @param Cw Leaf water content (g cm-2), or leaf water thickness (cm)
#' @param Cm Dry matter content (g cm-2)
#' @return Matrix of reflectance and transmittance values
#' @author Alexey Shiklomanov
#' @rdname prospect
#' @export
prospect4 <- function(N, Cab, Cw, Cm) {
  cc <- rbind(Cab, Cw, Cm) / N
  k <-  dataspec_p4 %*% cc
  gpm(k, N, p45_talf, p45_t12, p45_t21)
}

#' @rdname prospect
#' @export
prospect5 <- function(N, Cab, Car, Cw, Cm, Cbrown = 0) {
  cc <- rbind(Cab, Car, Cbrown, Cw, Cm) / N
  k <-  dataspec_p5 %*% cc
  gpm(k, N, p45_talf, p45_t12, p45_t21)
}

#' @rdname prospect
#' @export
prospectd <- function(N, Cab, Car, Canth, Cbrown, Cw, Cm) {
  cc <- rbind(Cab, Car, Canth, Cbrown, Cw, Cm) / N
  k <-  dataspec_pd %*% cc
  gpm(k, N, pd_talf, pd_t12, pd_t21)
}
