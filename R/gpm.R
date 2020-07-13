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
#' @param N Effective number of leaf layers (`numeric(1)`)
#' @param talf,t12,t21 Pre-calculated quantities based on refractive index and angle
#' @param e1fun Function to use for exponential integral. Default is e1_approx,
#'   included here. For greater precision, use `gsl::expint_E1`.
#' @return A length-2 list containing the modeled reflectance and transmittance
#' @author Alexey Shiklomanov
gpm <- function(k, N, talf, t12, t21, e1fun = e1_approx) {
  # global transmittance
  gt0 <- k > 0
  k0 <- k[gt0]
  trans <- rep(1.0, length(k))
  trans[gt0] <- (1 - k0) * exp(-k0) + k0 ^ 2 * e1fun(k0)

  tlt0 <- trans < 0
  if (any(tlt0)) {
    warning(
      sum(tlt0),
      " transmissivity values less than zero ",
      "were set to zero. ",
      "This happens under very high absorptivity values ",
      "(e.g. very high pigment concentrations) ",
      "and when using the fast exponential integral approximation. ",
      "The transmissivity is on the order of <= 1e-7 ",
      "so this warning can probably be ignored.",
      "However, if you want more precise results, you can try setting ",
      "`e1fun = gsl::expint_E1` to use a more precise approximation."
    )
    trans[tlt0] <- 0
  }

  # Reflectance/transmittance of a single layer
  # Commented quantities are precalculated
  ## talf <- tav_abs(40, refractive)
  ralf <- 1 - talf
  ## t12 <- tav_abs(90, refractive)
  r12 <- 1 - t12
  ## t21 <- t12 / (refractive ^ 2)
  r21 <- 1 - t21

  denom <- 1 - r21 ^ 2 * trans ^ 2
  Ta <- talf * trans * t21 / denom
  Ra <- ralf + r21 * trans * Ta

  tt <- t12 * trans * t21 / denom
  rr <- r12 + r21 * trans * tt

  Tsub <- numeric(2101)
  Rsub <- numeric(2101)

  gt1 <- rr + tt >= 1
  tgt1 <- tt[gt1]
  Tsub[gt1] <- tgt1 / (tgt1 + (1 - tgt1) * (N - 1))
  Rsub[gt1] <- 1 - Tsub[gt1]

  # Extremely high absorptivity leads to zero reflectance and transmittance
  inf <- rr == 0 | tt == 0
  Tsub[inf] <- 0
  Rsub[inf] <- 0

  # Reflectance/transmittance of N layers
  r <- rr[!gt1 & !inf]
  t <- tt[!gt1 & !inf]
  D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
  r2 <- r ^ 2
  t2 <- t ^ 2
  va <- (1 + r2 - t2 + D) / (2 * r)
  vb <- (1 - r2 + t2 + D) / (2 * t)

  vbNN <- vb ^ (N - 1)
  vbNN2 <- vbNN ^ 2
  va2 <- va ^ 2
  denomx <- va2 * vbNN2 - 1
  Rsub[!gt1 & !inf] <- va * (vbNN2 - 1) / denomx
  Tsub[!gt1 & !inf] <- vbNN * (va2 - 1) / denomx

  denomy <- 1 - Rsub * rr
  ## result <- matrix(NA_real_, nrow = 2101, ncol = 2,
  ##                  dimnames = list(NULL, c("reflectance", "transmittance")))
  ## result[, "reflectance"] <- Ra + Ta * Rsub * t / denomy
  ## result[, "transmittance"] <- Ta * Tsub / denomy

  # NOTE: A matrix here might be better, but returning a list is much faster.
  result <- list(
    reflectance = Ra + Ta * Rsub * tt / denomy,
    transmittance = Ta * Tsub / denomy
  )
  result
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

#' Swami and Ohija approximation to the Exponential Integral. This is very
#' simple, but differences with more precise (but expensive) approaches are
#' almost completely washed out by PROSPECT math.
e1_approx <- function(x) {
  A <- log((0.56146 / x + 0.65) * (1 + x))
  B <- x^4 * exp(7.7 * x) * (2 + x)^3.7
  (A^-7.7 + B)^-0.13
}
