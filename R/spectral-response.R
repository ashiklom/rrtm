#' Sensor spectral response functions
#'
#' @export
sensor_list <- c("identity", "aviris.ng", "aviris.classic",
                 "hyperion", "chris.proba", "landsat5", "landsat7",
                 "landsat8", "modis", "viirs", "avhrr", "licor")

#' Sensor list with proper names
#'
#' @export
sensor_proper <- c("Native", "AVIRIS NG", "AVIRIS Classic",
                   "Hyperion", "CHRIS-Proba", "Landsat 5", "Landsat 7",
                   "Landsat 8", "MODIS", "VIIRS", "AVHRR", "LiCor 6400 chamber")

names(sensor_proper) <- sensor_list

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
