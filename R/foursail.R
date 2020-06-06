#' 4SAIL model of canopy reflectance
#'
#' @param leaf_refl
#' @param leaf_trans
#' @param soil_refl
#' @param LAI
#' @param hot_spot
#' @param LIDFa
#' @param LIDFb
#' @param solar_zenith
#' @param instrument_zenith
#' @param azimuth
#' @return
#' @author Alexey Shiklomanov
foursail <- function(leaf_refl, leaf_trans, soil_refl, LAI,
                     hot_spot = 0,
                     LIDFa = -0.35, LIDFb = -0.15,
                     solar_zenith = 0,
                     instrument_zenith = 0,
                     azimuth = 0) {

  stopifnot(
    LAI >= 0,
    all(leaf_refl >= 0),
    all(leaf_trans >= 0),
    all(soil_refl >= 0),
    length(leaf_trans) %in% c(1, length(leaf_refl)),
    length(soil_refl) %in% c(1, length(leaf_refl))
  )

  if (LAI == 0) {
    return(list(
      bhr = soil_refl,
      dhr = soil_refl,
      hdr = soil_refl,
      bdr = soil_refl
    ))
  }

  # Degrees to radians
  d2r <- pi / 180

  # FLAG 1 -- Geometric quantities (sail_sensgeom)
  cts <- cos(d2r * solar_zenith)
  cto <- cos(d2r * instrument_zenith)
  ctscto <- cts * cto

  tants <- tan(d2r * solar_zenith)
  tanto <- tan(d2r * instrument_zenith)
  cospsi <- cos(d2r * azimuth)

  dso <- sqrt(tants^2 + tanto^2 - 2 * tants * tanto * cospsi)

  # FLAG 2 -- leaf angle distribution (LIDF_fun)
  #
  # Right now, the default (13) values from SAIL, but in theory, this can be
  # extended to continuous number.
  # TODO: Extend to number of angles
  # TODO: Pull out of SAIL calculation (don't have to do it every time)
  litab <- c(seq(5, 75, 10), seq(81, 89, 2))
  dcum_in <- c(seq(10, 80, 10), seq(82, 88, 2))
  F_lidf <- c(dcum(LIDFa, LIDFb, dcum_in), 1)
  lidf <- F_lidf
  for (i in seq(length(litab), 2)) {
    lidf[i] <- lidf[i] - lidf[i - 1]
  }

  # FLAG 3 -- Angular distance; compensation of shadow length (SUITS)

  ##########
  # volscatt
  ##########
  ctl <- cos(d2r * litab)
  vs <- volscatt(solar_zenith, instrument_zenith, azimuth, litab)

  # Extinction coefficients
  ksli <- vs$chi_s / cts
  koli <- vs$chi_o / cto

  # Area scattering coefficient fractions
  sobli <- vs$frho * pi / ctscto
  sofli <- vs$ftau * pi / ctscto
  bfli <- ctl * ctl
  ks <- sum(ksli * lidf)
  ko <- sum(koli * lidf)
  bf <- sum(bfli * lidf)
  sob <- sum(sobli * lidf)
  sof <- sum(sofli * lidf)

  sdb <- 0.5 * (ks + bf)
  sdf <- 0.5 * (ks - bf)
  dob <- 0.5 * (ko + bf)
  dof <- 0.5 * (ko - bf)
  ddb <- 0.5 * (1 + bf)
  ddf <- 0.5 * (1 - bf)

  sigb <- (ddb * leaf_refl) + (ddf * leaf_trans)
  sigf <- (ddf * leaf_refl) + (ddb * leaf_trans)

  att <- 1 - sigf
  m2 <- (att + sigb) * (att - sigb)
  m2[m2 < 0] <- 0
  m <- sqrt(m2)
  sb <- sdb * leaf_refl + sdf * leaf_trans
  sf <- sdf * leaf_refl + sdb * leaf_trans
  vb <- dob * leaf_refl + dof * leaf_trans
  vf <- dof * leaf_refl + dob * leaf_trans
  w <- sob * leaf_refl + sof * leaf_trans

  e1 <- exp(-m * LAI)
  e2 <- e1 * e1
  rinf <- (att - m) / sigb
  rinf2 <- rinf * rinf
  re <- rinf * e1
  denom <- 1 - rinf2 * e2

  J1ks <- jfunc1(ks, m, LAI)
  J2ks <- jfunc2(ks, m, LAI)
  J1ko <- jfunc1(ko, m, LAI)
  J2ko <- jfunc2(ko, m, LAI)

  Ps <- (sf + sb * rinf) * J1ks
  Qs <- (sf * rinf + sb) * J2ks
  Pv <- (vf + vb * rinf) * J1ko
  Qv <- (vf * rinf + vb) * J2ko

  rdd <- rinf * (1. - e2) / denom
  tdd <- (1. - rinf2) * e1 / denom
  tsd <- (Ps - re * Qs) / denom
  rsd <- (Qs - re * Ps) / denom
  tdo <- (Pv - re * Qv) / denom
  rdo <- (Qv - re * Pv) / denom

  tss <- exp(-ks * LAI)
  too <- exp(-ko * LAI)
  z <- jfunc2(ks, ko, LAI)
  g1 <- (z - J1ks * too) / (ko + m)
  g2 <- (z - J1ko * tss) / (ks + m)

  Tv1 <- (vf*rinf+vb)*g1
  Tv2 <- (vf+vb*rinf)*g2
  T1 <- Tv1*(sf+sb*rinf)
  T2 <- Tv2*(sf*rinf+sb)
  T3 <- (rdo*Qs+tdo*Ps)*rinf

  # Multiple scattering contribution to bidirectional canopy reflectance
  rsod <- (T1 + T2 - T3) / (1 - rinf2)

  # Hot spot effect
  alf <- 1e6
  if (hot_spot > 0) {
    alf <- (dso / hot_spot) * 2 / (ks + ko)
  }
  if (alf > 200) alf <- 200
  if (alf == 0) {
    tstoo <- tss
    sumint <- (1 - tss) / (ks * LAI)
  } else {
    # Outside the hotspot
    fhot <- LAI * sqrt(ks * ko)
    x1=0.
    y1=0.
    f1=1.
    fint=(1.-exp(-alf)) * 0.05
    sumint=0.

    for (i in seq(1, 20)) {
      if (i < 20) {
        x2 <- -log(1 - i * fint) / alf
      } else {
        x2 <- 1
      }
      y2 <- -(ko + ks) * LAI * x2 + fhot * (1 - exp(-alf * x2)) / alf
      f2 <- exp(y2)
      sumint <- sumint + (f2 - f1) * (x2 - x1) / (y2 - y1)
      x1 <- x2
      y1 <- y2
      f1 <- f2
    }
    tsstoo <- f1
  }


  rsos <-  w * LAI * sumint

  # Total canopy contribution
  rso <- rsos + rsod

  # Interaction with the soil
  dn <- 1. - soil_refl * rdd

  # bi-hemispherical reflectance factor
  bhr <- rdd + tdd * soil_refl * tdd / dn
  # directional-hemispherical reflectance factor for solar incident flux
  dhr <- rsd + (tsd + tss) * soil_refl * tdd / dn
  # hemispherical-directional reflectance factor in viewing direction
  hdr <- rdo + tdd * soil_refl * (tdo + too) / dn
  # bi-directional reflectance factor
  rsodt <- rsod + ((tss + tsd) * tdo + (tsd + tss * soil_refl * rdd) * too) *
    soil_refl / dn
  rsost <- rsos + tsstoo * soil_refl
  bdr <- rsost + rsodt

  list(bhr = bhr, dhr = dhr, hdr = hdr, bdr = bdr)
}

jfunc1 <- function(k, l, t) {
  del <- (k - l) * t
  out <- numeric(length(l))
  iii <- abs(del) > 1e-3
  out[iii] <- (exp(-l * t) - exp(-k * t)) / (k - l)
  out[!iii] <- 0.5 * t * (exp(-k * t) + exp(-l * t)) * (1 - (del^2) / 12)
  out
}

jfunc2 <- function(k, l, t) {
  (1 - exp(-(k + l) * t)) / (k + l)
}
