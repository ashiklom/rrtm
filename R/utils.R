#' Test for equality within machine tolerance
#'
#' @param x,y Numeric values to compare
#' @return `TRUE` if equal, `FALSE` otherwise
#' @author Alexey Shiklomanov
eq <- function(x, y) {
  isTRUE(all.equal(x, y))
}
