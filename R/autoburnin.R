#' Automatically determine burn-in period for MCMC chains
#'
#' Estimates the burn-in period by finding the point where the running mean
#' of the chain stabilizes within a specified tolerance of the overall mean.
#'
#' @param samples A matrix or data.frame of MCMC samples (rows = iterations, columns = parameters)
#' @param tolerance Numeric. Convergence tolerance (default = 0.01)
#' @param min_burnin Integer. Minimum number of iterations to evaluate (default = 10)
#' @return Integer indicating the number of burn-in iterations to discard.
#' @export
autoburnin <- function(samples, tolerance = 0.01, min_burnin = 10) {
  samples <- as.matrix(samples)
  n <- nrow(samples)
  if (n <= min_burnin) return(0L)

  overall_means <- colMeans(samples)
  burn_in <- n

  for (i in (min_burnin + 1):n) {
    running_means <- colMeans(samples[1:i, , drop = FALSE])
    if (max(abs(running_means - overall_means)) <= tolerance) {
      burn_in <- i
      break
    }
  }

  return(burn_in)
}

#' Extract post burn-in samples
#'
#' @param samples MCMC samples
#' @param burn_in Number of iterations to discard
#' @return Matrix of post burn-in samples
#' @export
extract_post_burnin <- function(samples, burn_in) {
  samples[(burn_in + 1):nrow(samples), , drop = FALSE]
}
