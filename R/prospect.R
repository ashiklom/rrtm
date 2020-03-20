#' PROSPECT leaf spectra model
#'
#' @param N Effective number of leaf layers
#' @param Cab Chlorophyll content (ug cm -2)
#' @param Car Carotenoid content (ug cm -2)
#' @param Canth Anthocyanin content (ug cm-2)
#' @param Cbrown "Brown matter" content (unitless)
#' @param Cw Leaf water content (g cm-2), or leaf water thickness (cm)
#' @param Cm Dry matter content (g cm-2)
#' @param kmat Spectral coefficient matrix
#' @return Matrix of reflectance and transmittance values
#' @rdname prospect
#' @author Alexey Shiklomanov
#' @rdname prospect
#' @export
prospect4 <- function(N, Cab, Cw, Cm,
                      kmat = dataspec_p4,
                      refractive = refractive_p45) {
  stopifnot(NCOL(kmat) == 3)
  # TODO: Broadcast division
  cc <- rbind(Cab, Cw, Cm) / N
  k <-  kmat %*% cc 
  gpm(k, refractive, N)
}

#' @rdname prospect
#' @export
prospect5 <- function(N, Cab, Car, Cw, Cm, Cbrown = 0,
                      kmat = dataspec_p5,
                      refractive = refractive_p45) {
  stopifnot(NCOL(kmat) == 5)
  cc <- rbind(Cab, Car, Cbrown, Cw, Cm) / N
  k <-  kmat %*% cc
  gpm(k, refractive, N)
}

#' @rdname prospect
#' @export
prospectd <- function(N, Cab, Car, Canth, Cbrown, Cw, Cm,
                      kmat = dataspec_pd,
                      refractive = refractive_pd) {
  stopifnot(NCOL(kmat) == 6)
  cc <- rbind(Cab, Car, Canth, Cbrown, Cw, Cm) / N
  k <-  kmat %*% cc
  gpm(k, refractive, N)
}
