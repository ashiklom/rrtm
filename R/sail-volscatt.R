#' Volume scattering function
#'
#' @param solar_zenith Solar zenith angle (degrees)
#' @param instrument_zenith Instrument zenith angle (degrees)
#' @param azimuth Sun-instrument azimuth angle (degrees)
#' @param leaf_angle Leaf angle (degrees)
volscatt <- function(solar_zenith, instrument_zenith, azimuth, leaf_angle) {
  d2r <- pi / 180

  # Sun-sensor geoms
  costs <- cos(d2r * solar_zenith)
  sints <- sin(d2r * solar_zenith)
  costo <- cos(d2r * instrument_zenith)
  sinto <- sin(d2r * instrument_zenith)
  cospsi <- cos(d2r * azimuth)
  psirad <- d2r * azimuth

  # Leaf geoms
  costl <- cos(d2r * leaf_angle)
  sintl <- sin(d2r * leaf_angle)

  cs <- costl * costs
  co <- costl * costo
  ss <- sintl * sints
  so <- sintl * sinto

  # Transition angles (beta) for solar (bts) and instrument (bto) directions
  cosbts <- rep_len(5, length(cs))
  i <- abs(ss) > 1e-6
  if (any(i)) cosbts[i] <- -cs[i] / ss[i]

  cosbto <- rep_len(5, length(co))
  i <- abs(so) > 1e-6
  if (any(i)) cosbto[i] <- -co[i] / so[i]

  bts <- rep_len(pi, length(cosbts))
  ds <- rep_len(cs, length(cosbts))
  i <- abs(cosbts) < 1
  if (any(i)) {
    bts[i] <- acos(cosbts[i])
    ds[i] <- ss[i]
  }

  chi_s <- 2 / pi * ((bts - 0.5 * pi) * cs + sin(bts) * ss)

  bto <- rep_len(0, length(cosbto))
  doo <- rep_len(-co, length(cosbto))
  i <- abs(cosbto) < 1
  if (any(i)) {
    bto[i] <- acos(cosbto[i])
    doo[i] <- so[i]
  }
  j <- instrument_zenith < 90
  if (any(j)) {
    bto[j] <- pi
    doo[j] <- co[j]
  }

  chi_o <- 2 / pi * ((bto - 0.5 * pi) * co + sin(bto) * so)

  # Auxiliary azimuth angles bt1, bt2, bt3, for bidirectional scattering
  btran1 <- abs(bts - bto)
  btran2 <- pi - abs(bts + bto - pi)

  bt1 <- rep_len(psirad, length(btran1))
  bt2 <- btran1
  bt3 <- btran2

  i <- psirad >= btran1
  bt1[i] <- btran1[i]
  bt2[i] <- psirad
  bt3[i] <- btran2[i]
  j <- psirad >= btran2
  bt2[j] <- btran2[j]
  bt3[j] <- psirad

  t1 <- 2 * cs * co + ss * so * cospsi
  t2 <- rep_len(0, length(t1))
  i <- bt2 > 0
  t2[i] <- sin(bt2[i]) * (2 * ds[i] * doo[i] +
                            ss[i] * so[i] * cos(bt1[i]) * cos(bt3[i]))

  denom <- 2 * pi^2
  frho <- ((pi - bt2) * t1 + t2) / denom
  ftau <- (-bt2 * t1 + t2) / denom
  frho[frho < 0] <- 0
  ftau[ftau < 0] <- 0

  list(
    chi_s = chi_s,
    chi_o = chi_o,
    frho = frho,
    ftau = ftau
  )
}

#' Original volscatt implementation, based on SAIL Fortran code. Exists for
#' testing of faster vectorized version.
#'
#' @param solar_zenith Solar zenith angle (degrees)
#' @param instrument_zenith Instrument zenith angle (degrees)
#' @param azimuth Sun-instrument azimuth angle (degrees)
#' @param leaf_angle Leaf angle (degrees)
#' @return
#' @author Alexey Shiklomanov
volscatt_scalar <- function(solar_zenith, instrument_zenith, azimuth, leaf_angle) {
  d2r <- pi / 180

  # Sun-sensor geoms
  costs <- cos(d2r * solar_zenith)
  sints <- sin(d2r * solar_zenith)
  costo <- cos(d2r * instrument_zenith)
  sinto <- sin(d2r * instrument_zenith)
  cospsi <- cos(d2r * azimuth)
  psirad <- d2r * azimuth

  # Leaf geoms
  costl <- cos(d2r * leaf_angle)
  sintl <- sin(d2r * leaf_angle)

  cs <- costl * costs
  co <- costl * costo
  ss <- sintl * sints
  so <- sintl * sinto

  # Transition angles (beta) for solar (bts) and instrument (bto) directions
  cosbts <- 5
  if (abs(ss) > 1e-6) cosbts <- -cs / ss
  cosbto <- 5
  if (abs(so) > 1e-6) cosbto <- -co / so

  if (abs(cosbts) < 1) {
    bts <- acos(cosbts)
    ds <- ss
  } else {
    bts <- pi
    ds <- cs
  }

  chi_s <- 2 / pi * ((bts - 0.5 * pi) * cs + sin(bts) * ss)

  if (abs(cosbto) < 1) {
    bto <- acos(cosbto)
    doo <- so
  } else if (instrument_zenith < 90) {
    bto <- pi
    doo <- co
  } else {
    bto <- 0
    doo <- -co
  }

  chi_o <- 2 / pi * ((bto - 0.5 * pi) * co + sin(bto) * so)

  # Auxiliary azimuth angles bt1, bt2, bt3, for bidirectional scattering
  btran1 <- abs(bts - bto)
  btran2 <- pi - abs(bts + bto - pi)

  if (psirad < btran1) {
    bt1 <- psirad
    bt2 <- btran1
    bt3 <- btran2
  } else {
    bt1 <- btran1
    if (psirad < btran2) {
      bt2 <- psirad
      bt3 <- btran2
    } else {
      bt2 <- btran2
      bt3 <- psirad
    }
  }

  t1 <- 2 * cs * co + ss * so * cospsi
  t2 <- 0
  if (bt2 > 0) {
    t2 <- sin(bt2) * (2 * ds * doo + ss * so * cos(bt1) * cos(bt3))
  }

  denom <- 2 * pi^2
  frho <- ((pi - bt2) * t1 + t2) / denom
  ftau <- (-bt2 * t1 + t2) / denom
  if (frho < 0) frho <- 0
  if (ftau < 0) ftau <- 0

  list(
    chi_s = chi_s,
    chi_o = chi_o,
    frho = frho,
    ftau = ftau
  )
}

#' "Vectorized" version of scalar code
volscatt_scalar_v <- Vectorize(volscatt_scalar, "leaf_angle")
