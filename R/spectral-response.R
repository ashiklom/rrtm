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
#' @alias sensor.rsr
NULL

# sensor.rsr is available via lazily-loaded `data`. This silences the relevant
# R CMD check NOTE.
utils::globalVariables("sensor.rsr")

#' Convolution of spectra to sensor RSR
#'
#' @param spec Full (1 nm) spectrum ((column) vector or, matrix with each
#'  spectrum in its own column)
#' @param sensor Sensor name (string). See `sensor_list`.
#' @return A matrix of convolved spectra. Each *row* is a spectrum; each column
#' is the sensor band. For type stability, this *always* returns a matrix (so a
#' single spectrum passed as a vector will return a matrix with 1 row and nwl
#' columns).
#' @export
spectral_response <- function(spec, sensor) {
  sensor <- tolower(sensor)
  stopifnot(sensor %in% sensor_list)
  if (sensor == "identity") {
    return(matrix(spec, nrow = nrow(spec) %||% 1))
  }
  rsr <- sensor.rsr[[sensor]]
  spec[rsr[, "index"]] %*% rsr[, -1]
} # spectral_response
