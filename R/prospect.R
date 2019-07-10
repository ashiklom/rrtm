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
#' @useDynLib rrtm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname prospect
#' @export
prospect4 <- function(N, Cab, Cw, Cm,
                      kmat = dataspec_p4,
                      refractive = refractive_p45) {
  stopifnot(NCOL(kmat) == 3)
  k <-  kmat %*% c(Cab, Cw, Cm) / N
  gpm(k, refractive, N)
}

#' @rdname prospect
#' @export
prospect5 <- function(N, Cab, Car, Cw, Cm, Cbrown = 0,
                      kmat = dataspec_p5,
                      refractive = refractive_p45) {
  stopifnot(NCOL(kmat) == 5)
  k <-  kmat %*% c(Cab, Car, Cbrown, Cw, Cm) / N
  gpm(k, refractive, N)
}

#' @rdname prospect
#' @export
prospectd <- function(N, Cab, Car, Canth, Cbrown, Cw, Cm,
                      kmat = dataspec_pd,
                      refractive = refractive_pd) {
  stopifnot(NCOL(kmat) == 6)
  k <-  kmat %*% c(Cab, Car, Canth, Cbrown, Cw, Cm) / N
  gpm(k, refractive, N)
}

gpm <- function(k, refractive, N) {
  # Global transmittance
  trans <- (1 - k) * exp(-k) + k ^ 2 * gsl::expint_E1(k)

  # reflectivity and transmissivity at the interface
  t12 <- tav_abs(40, refractive)
  tav_90 <- tav_abs(90, refractive)
  t21 <- tav_90 / refractive ^ 2
  r12 <- 1 - t12
  r21 <- 1 - t21
  x <- t12 / tav_90
  y <- x * (tav_90 - 1) + 1 - t12

  # Reflectance and transmittance of the elementary layer (N = 1)
  trans_2 <- trans ^ 2
  t12_t21 <- t12 * t21
  r21_2 <- r21 ^ 2
  ra <- r12 + (t12_t21 * r21 * trans_2) / (1 - r21_2 * trans_2)
  ta <- (t12_t21 * trans) / (1 - r21_2 * trans_2)
  r90 <- (ra - y) / x
  t90 <- ta / x

  t90_2 <- t90 ^ 2
  r90_2 <- r90 ^ 2
  delta <- sqrt((t90_2 - r90_2 - 1)^2 - 4 * r90_2)
  beta <- (1 + r90_2 - t90_2 - delta) / (2 * r90)
  va <- (1 + r90_2 - t90_2 + delta) / (2 * r90)

  va_beta <- pmin(va * (beta - r90), 1e-14)
  vb <- sqrt(beta * (va - r90) / va_beta)

  vbNN <- vb ^ (N - 1)
  vbNNi <- vbNN ^ -1
  vbNN_diff <- vbNN - vbNNi
  vai <- va ^ -1

  s1 <- ta * t90 * vbNN_diff
  s2 <- ta * (va - vai)
  s3 <- va * vbNN - vai * vbNNi - r90 * vbNN_diff
  
  ### Calculate output reflectance and transmittance of the modeled leaf
  RN <- ra + s1 / s3
  TN <- s2 / s3

  cbind("reflectance" = RN, "transmittance" = TN)
}

#' Generalized plate model
#'
#' @references Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R.
#'   (1969), Interaction of isotropic ligth with a compact plant leaf,
#'   Journal of the Optical Society of American, 59:1376-1379.
#' @references Stokes G.G. (1862), On the intensity of the light
#'   reflected from or transmitted through a pile of plates,
#'   Proceedings of the Royal Society of London, 11:545-556.
#'
#' @param k Optical depth (`numeric(nwl)`)
#' @param refractive Refractive index (`numeric(nwl)`)
#' @param N Effective number of leaf layers (`numeric(1)`)
#' @return 
#' @author Alexey Shiklomanov
gpm2 <- function(k, refractive, N) {

  # Reflectance/transmittance of a single layer

  # Can be pre-calculated for inversion
  refractive2 <- refractive * refractive
  t1 <- tav_abs(90, refractive)
  t2 <- tav_abs(40, refractive)
  t12 <- t1 * t1
  t22 <- t2 * t2

  # Depend on k -- cannot be pre-calculated 
  # NOTE: This returns a dim(k) x 1 matrix
  tau <- gsl::expint_E1(k)[,1]
  tau2 <- tau * tau

  x1 <- 1 - t1
  x2 <- t12 * tau2 * (refractive2 - t1)
  x3 <- t12 * tau * refractive2
  x4 <- refractive2 * refractive2 - tau2 * (refractive2 - t1) ^ 2
  x5 <- t2 / t1
  x6 <- x5 * (t1 - 1) + 1 - t2
  r  <- x1 + x2 / x4
  t  <- x3 / x4
  ra <- x5 * r + x6
  ta <- x5 * t

  # Reflectance and transmittance of N layers

  t_2 <- t * t
  r_2 <- r * r
  rt_2 <- r_2 - t_2
  rx2 <- r * 2
  
  delta <- (t_2 - r_2 - 1) ^ 2 - 4 * r_2
  print(min(delta))
  print(sum(delta < 0))
  delta_rt <- sqrt(delta)
  beta <- (1 - rt_2 - delta_rt) / rx2
  va <- (1 + rt_2 + delta_rt) / rx2
  vb <- sqrt(beta * (va - r) / (va * (beta - r)))

  vbnm1 <- vb ^ (N - 1)
  inv_vbnm1 <- 1 / vbnm1
  inv_va <- 1 / va
  
  s1 <- ra * (va * vbnm1 - inv_va * inv_vbnm1) + (ta * t - ra * r) * (vbnm1 - inv_vbnm1)
  s2 <- ta * (va - inv_va)
  s3 <- va * vbnm1 - (1 / va) * inv_vbnm1 - r * (vbnm1 - inv_vbnm1)

  cbind(reflectance = s1 / s3, transmittance = s2 / s3)
}

tav_abs <- function(theta, refractive) {
  if (eq(theta, 0)) {
    result <- 4 * refractive / (refractive + 1) ^ 2
    return(result)
  }

  thetarad <- pi * theta / 180

  refr2 <- refractive * refractive
  refr2p1 <- refr2 + 1
  refr2p12 <- refr2p1 * refr2p1
  refr2p13 <- refr2p12 * refr2p1
  refr2m12 <- (refr2 - 1) ^ 2

  ax <- 0.5 * (refractive + 1) ^ 2
  bx <- -0.25 * refr2m12

  if (eq(thetarad, pi / 2)) {
    b1 <- 0
  } else {
    b1 <- sqrt(bx + (sin(thetarad) ^ 2 - 0.5 * refr2p1) ^ 2)
  }
  b2 <- sin(thetarad) ^ 2 - 0.5 * refr2p1
  b0 <- b1 - b2
  ts <- (bx^2 / (6 * b0^3) + bx / b0 - b0 / 2) -
    (bx^2 / (6 * ax^3) + bx / ax - ax / 2)
  tp1 <- -2 * refr2 * (b0 - ax) / refr2p12
  tp2 <- -2 * refr2 * refr2p1 * log(b0 / ax) / refr2m12
  tp3 <- refr2 * 0.5 * (1 / b0 - 1 / ax)
  tp4 <- 16 * refr2 ^ 2 * (refr2 ^ 2 + 1) *
    log((2 * refr2p1 * b0 - refr2m12) /
          (2 * refr2p1 * ax - refr2m12)) /
    (refr2p13 * refr2m12)
  tp5 <- 16 * refr2 ^ 3 *
    (1 / (2 * refr2p1 * b0 - refr2m12) - 1 /
       (2 * refr2p1 * ax - refr2m12)) / refr2p13
  tp <- tp1 + tp2 + tp3 + tp4 + tp5

  (ts + tp) / (2 * sin(thetarad) ^ 2)
}
