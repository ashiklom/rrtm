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
gpm <- function(k, refractive, N) {
  # global transmittance
  gt0 <- k > 0
  k0 <- k[gt0]
  trans <- matrix(1.0, nrow(k), ncol(k))
  trans[gt0] <- (1 - k) * exp(-k0) + k0 ^ 2 * gsl::expint_E1(k0)

  # Reflectance/transmittance of a single layer
  talf <- tav_abs(40, refractive)
  ralf <- 1 - talf
  t12 <- tav_abs(90, refractive)
  r12 <- 1 - t12
  t21 <- t12 / (refractive ^ 2)
  r21 <- 1 - t21

  denom <- 1 - r21 ^ 2 * trans ^ 2
  Ta <- talf * trans * t21 / denom
  Ra <- ralf + r21 * trans * Ta

  t <- t12 * trans * t21 / denom
  r <- r12 + r21 * trans * t

  # Reflectance/transmittance of N layers
  D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
  r2 <- r ^ 2
  t2 <- t ^ 2
  va <- (1 + r2 - t2 + D) / (2 * r)
  vb <- (1 - r2 + t2 + D) / (2 * t)

  vbNN <- vb ^ (N - 1)
  vbNN2 <- vbNN ^ 2
  va2 <- va ^ 2
  denomx <- va2 * vbNN2 - 1
  Rsub <- va * (vbNN2 - 1) / denomx
  Tsub <- vbNN * (va2 - 1) / denomx

  gt1 <- r + t >= 1
  tgt1 <- t[gt1]
  Tsub[gt1] <- tgt1 / (tgt1 + (1 - tgt1) * (N - 1))
  Rsub[gt1] <- 1 - Tsub[gt1]

  denomy <- 1 - Rsub * r
  RN <- Ra + Ta * Rsub * t / denomy
  TN <- Ta * Tsub / denomy

  rray::rray_bind(reflectance = RN, transmittance = TN, .axis = 3)
}

#' `tav` function
#'
#' @references Stern F. (1964), Transmission of isotropic radiation
#'   across an interface between two dielectrics, Appl. Opt.,
#'   3(1):111-113. Allen W.A. (1973), Transmission of isotropic light
#'   across a dielectric surface in two and three dimensions, J. Opt.
#'   Soc. Am., 63(6):664-666.
#' @param theta Angle of incident radiation (in degrees)
#' @param refractive Refractive index (vector)
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
