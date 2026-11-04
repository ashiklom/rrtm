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
  crossprod(rsr[, -1, drop = FALSE], spec[rsr[, "index"], , drop = FALSE])
} # spectral_response
