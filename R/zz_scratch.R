gpm_serbin <- function(k, refractive, N) {
  # Global transmittance
  trans <- (1 - k) * exp(-k) + k ^ 2 * gsl::expint_E1(k)

  # reflectivity and transmissivity at the interface
  t12 <- tav_abs(40, refractive)
  r12 <- 1 - t12
  tav_90 <- tav_abs(90, refractive)
  t21 <- tav_90 / refractive ^ 2
  r21 <- 1 - t21
  x <- t12 / tav_90
  y <- x * (tav_90 - 1) + 1 - t12

  # Reflectance and transmittance of the elementary layer (N = 1)
  trans_2 <- trans ^ 2
  t12_t21 <- t12 * t21
  r21_2 <- r21 ^ 2
  ra <- r12 + (t12_t21 * r21 * trans_2) / (1 - r21_2 * trans_2)
  ta <- (t12_t21 * trans) / (1 - r21_2 * trans_2)
  r90 <- (ra - y) / x
  t90 <- ta / x

  t90_2 <- t90 ^ 2
  r90_2 <- r90 ^ 2
  delta <- sqrt((t90_2 - r90_2 - 1) ^ 2 - 4 * r90_2)
  beta <- (1 + r90_2 - t90_2 - delta) / (2 * r90)
  va <- (1 + r90_2 - t90_2 + delta) / (2 * r90)

  va_beta <- pmin(va * (beta - r90), 1e-14)
  vb <- sqrt(beta * (va - r90) / va_beta)

  vbNN <- vb ^ (N - 1)
  vbNNi <- 1 / vbNN
  vbNN_diff <- vbNN - vbNNi
  vai <- 1 / va

  s1 <- ta * t90 * vbNN_diff
  s2 <- ta * (va - vai)
  s3 <- va * vbNN - vai * vbNNi - r90 * vbNN_diff
  
  ### Calculate output reflectance and transmittance of the modeled leaf
  RN <- ra + s1 / s3
  TN <- s2 / s3
  print(head(RN))
  print(summary(TN))

  cbind(reflectance = RN, transmittance = TN)
}

gpm_mine <- function(k, refractive, N) {

  # Reflectance/transmittance of a single layer

  # Can be pre-calculated for inversion
  refractive2 <- refractive * refractive
  t1 <- tav_abs(90, refractive)
  t2 <- tav_abs(40, refractive)
  t12 <- t1 * t1
  t22 <- t2 * t2

  # Depend on k -- cannot be pre-calculated 
  # NOTE: This returns a dim(k) x 1 matrix
  tau <- gsl::expint_E1(k)[,1]
  tau2 <- tau * tau

  x1 <- 1 - t1
  x2 <- t12 * tau2 * (refractive2 - t1)
  x3 <- t12 * tau * refractive2
  x4 <- refractive2 * refractive2 - tau2 * (refractive2 - t1) ^ 2
  x5 <- t2 / t1
  x6 <- x5 * (t1 - 1) + 1 - t2
  r  <- x1 + x2 / x4
  t  <- x3 / x4
  ra <- x5 * r + x6
  ta <- x5 * t

  # Reflectance and transmittance of N layers

  t_2 <- t * t
  r_2 <- r * r
  rt_2 <- r_2 - t_2
  rx2 <- r * 2
  
  delta <- (t_2 - r_2 - 1) ^ 2 - 4 * r_2
  print(min(delta))
  print(sum(delta < 0))
  delta_rt <- sqrt(delta)
  beta <- (1 - rt_2 - delta_rt) / rx2
  va <- (1 + rt_2 + delta_rt) / rx2
  vb <- sqrt(beta * (va - r) / (va * (beta - r)))

  vbnm1 <- vb ^ (N - 1)
  inv_vbnm1 <- 1 / vbnm1
  inv_va <- 1 / va
  
  s1 <- ra * (va * vbnm1 - inv_va * inv_vbnm1) + (ta * t - ra * r) * (vbnm1 - inv_vbnm1)
  s2 <- ta * (va - inv_va)
  s3 <- va * vbnm1 - (1 / va) * inv_vbnm1 - r * (vbnm1 - inv_vbnm1)

  cbind(reflectance = s1 / s3, transmittance = s2 / s3)
}
