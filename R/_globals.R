# sensor.rsr is available via lazily-loaded `data`. This silences the relevant
# R CMD check NOTE.
utils::globalVariables("sensor.rsr")

# This silences the R CMD CHECK note about RcppArmadillo being unused
#' @importFrom RcppArmadillo RcppArmadillo.package.skeleton
NULL
