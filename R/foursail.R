foursail <- function(leaf_refl, leaf_trans,
                     LIDFa, LIDFb, LIDF_type,
                     LAI,
                     hot_spot,
                     solar_zenith,
                     instrument_zenith,
                     azimuth,
                     soil_refl,
                     # Cumulative LIDF as a function of angle
                     # Default is spherical (integral of f(theta) = sin(theta))
                     cum_LIDF = function(x) 1 - cos(x * pi / 180)) {

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
  litab <- c(seq(5, 75, 10), seq(81, 89, 2))
  F_lidf <- cum_LIDF(litab)

  # FLAG 3 -- Angular distance; compensation of shadow length (SUITS)
  ks <- 0
  ko <- 0
  sob <- 0
  sof <- 0
  sdb <- 0
  sdf <- 0
  dob <- 0
  dof <- 0
  ddb <- 0
  ddf <- 0

  ##########
  # volscatt
  ##########
}
