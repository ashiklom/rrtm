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
#' @param spec Full (1 nm) spectrum (vector)
#' @param sensor Sensor name (string). See sensor.list
#' @export
spectral_response <- function(spec, sensor) {
  sensor <- tolower(sensor)
  stopifnot(sensor %in% sensor_list)
  if (sensor == "identity") {
    return(spec)
  }
  rsr <- sensor_rsr[[sensor]]
  m <- spec[rsr[, "index"]] * rsr[, -1]
  .colSums(m, nrow(m), ncol(m))
} # spectral.response
