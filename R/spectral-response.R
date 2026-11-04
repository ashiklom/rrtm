#' List of sensors with available spectral response functions
#'
#' @export
sensor_list <- c(
  "Native" = "identity",
  "AVIRIS NG" = "aviris.ng",
  "AVIRIS Classic" = "aviris.classic",
  "Hyperion" = "hyperion",
  "CHRIS-Proba" = "chris.proba",
  "Landsat 5" = "landsat5",
  "Landsat 7" = "landsat7",
  "Landsat 8" = "landsat8",
  "MODIS" = "modis",
  "VIIRS" = "viirs",
  "AVHRR" = "avhrr",
  "LiCor 6400 chamber" = "licor"
)

#' Spectral response functions for various sensors
#'
#' @name sensor.rsr
NULL

# sensor.rsr is available via lazily-loaded `data`. This silences the relevant
# R CMD check NOTE.
utils::globalVariables("sensor.rsr")

#' Convolution of spectra to sensor RSR
#'
#' Each band value is the average of the spectrum weighted by the band's
#' relative spectral response (RSR). The weights are normalized to sum to 1
#' (after setting the small negative values in some published RSRs to 0), so
#' a band value always lies within the range of the input spectrum over that
#' band; e.g., reflectances between 0 and 1 give band values between 0 and 1.
#'
#' @param spec Full (1 nm, 400:2500 nm) spectrum ((column) vector, or matrix
#'   with each spectrum in its own column). For `sensor = "identity"`, any
#'   length is allowed.
#' @param sensor Sensor name (string). See `sensor_list`.
#' @return A matrix of convolved spectra. Each *column* is a spectrum; each row
#' is a sensor band. For type stability, this *always* returns a matrix (so a
#' single spectrum passed as a vector will return a matrix with nband rows and
#' 1 column).
#' @export
spectral_response <- function(spec, sensor) {
  sensor <- tolower(sensor)
  stopifnot(sensor %in% sensor_list)
  spec <- as.matrix(spec)
  if (sensor == "identity") {
    return(spec)
  }
  if (nrow(spec) != 2101) {
    stop("`spec` must have 2101 rows (400:2500 nm), not ", nrow(spec), ".")
  }
  rsr <- sensor.rsr[[sensor]]
  weights <- pmax(rsr[, -1, drop = FALSE], 0)
  weights <- sweep(weights, 2, colSums(weights), "/")
  x <- spec[rsr[, "index"], , drop = FALSE]
  out <- crossprod(weights, x)
  # A weighted average lies within the range of its inputs; clamp to remove
  # floating point overshoot (e.g. 1 + 1e-15 for a spectrum of all 1s).
  lo <- rep(apply(x, 2, min), each = nrow(out))
  hi <- rep(apply(x, 2, max), each = nrow(out))
  pmin(pmax(out, lo), hi)
} # spectral_response
