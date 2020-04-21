#' Simple model of soil reflectance as a function of soil moisture
#'
#' Just the average of [wet_soil] and [dry_soil], weighted by
#' `soil_moisture`. `soil_moisture` of 0 means fully dry soil, and
#'
#' @param soil_moisture Soil moisture index (0-1; 0 = completely dry,
#'   1 = completely wet)
#' @return Spectrum of soil reflectance (400:2500 nm)
#' @author Alexey Shiklomanov
#' @export
hapke_soil <- function(soil_moisture) {
  stopifnot(soil_moisture >= 0, soil_moisture <= 1)
  wet_soil * soil_moisture + dry_soil * (1 - soil_moisture)
}
