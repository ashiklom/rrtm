#' Leaf inclination distribution function from Verhoef dissertation
#'
#' @param a,b LIDF parameters
#' @param theta Angle
#' @return LIDF fraction value
dcum <- function(a, b, theta) {
  d2r <- pi / 180
  thetarad <- theta * d2r
  x0 <- x <- 2 * thetarad
  t <- 1e-6
  dx_abs <- 1
  while (any(dx_abs > t)) {
    y <- a * sin(x) + 0.5 * b * sin(2 * x)
    dx <- 0.5 * (y - x + x0)
    x <- x + dx
    dx_abs <- abs(dx)
  }
  2 * (y + thetarad) / pi
}

if (FALSE) {

  # Some examples
  curve(dcum(0, 0, x), 0, 90)
  curve(dcum(-0.35, -0.15, x), 0, 90, add = TRUE, col = "blue")
  curve(1 - cos(x * pi / 180), 0, 90, add = TRUE, col = "red")
 
}
