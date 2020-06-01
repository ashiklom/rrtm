if (FALSE) {
  usethis::use_mit_license("Alexey Shiklomanov")
  usethis::use_roxygen_md()

  usethis::use_rcpp()

  usethis::use_data_raw()

  usethis::use_testthat()

  Rcpp::Rcpp.package.skeleton(path = ".")
}


kmat <- dataspec_p4
cc <- rbind(40, c(0.01, 0.005), 0.01)
N <- 1.4
k <-  kmat %*% cc
refractive <- refractive_p45

out2 <- prospect4(1.4, c(30, 40, 50), c(0.01, 0.005, 0.002), 0.01)
plot(c(400, 2500), c(0, 1), type = 'n',
     xlab = "Wavelength (nm)",
     ylab = "Spectra")
matplot(400:2500, f[,,1], type = "l", add = TRUE)
matplot(400:2500, 1 - f[,,2], type = "l", add = TRUE)

if (FALSE) {
  debug(gpm)
  undebug(gpm)
}

params <- c(1.4, 40, 0.01, 0.01)
k <- (dataspec_p4 %*% params[-1]) / params[[1]]
plot(k, type = 'l')
refractive <- refractive_pd
plot(refractive, type = "l")

fortran <- tibble::tibble(
    tau = gsl::expint_E1(k)[,1],
    t1 = tav_abs(90,refractive),
    t2 = tav_abs(40,refractive),
    x1      = 1-t1,
    x2      = t1^2*tau^2*(refractive^2-t1),
    x3      = t1^2*tau*refractive^2,
    x41 = refractive ^ 4,
    x42 = (refractive ^ 2 - t1),
    x43 = x42 * x42, 
    x44 = x43 * tau * tau,
    x4      = refractive^4-tau^2*(refractive^2-t1)^2,
    x5      = t2/t1,
    x6      = x5*(t1-1)+1-t2,
    r       = x1+x2/x4,
    t       = x3/x4,
    ra      = x5*r+x6,
    ta      = x5*t,
    delta   = (t^2-r^2-1)^2-4*r^2
)

refractive_p45[2001]

sum(fortran$x41 - fortran$x44 < 0)

plot(refractive_p45 - refractive_pd, type = "l", xlim = xlim)
abline(v = i, col = "red", lty = "dashed")

xlim <- c(1200, 1300)
with(fortran, plot(x41 - x44, type = "l", xlim = xlim, ylim = c(0, 1)))
abline(h = 0, col = "blue", lty = "dashed")

i <- which(abs(fortran$x2 / fortran$x4) > 100)
par(mfrow = c(2, 1))
plot(fortran$x2, type = "l", xlim = c(1200, 1300))
abline(v = i, col = "red")
plot(fortran$x4, type = "l", xlim = c(1200, 1300))
abline(v = i, col = "red")
abline(h = 0, col = "blue")


i <- which(fortran$delta < 0)
plot(k, type = "l")
abline(v = i, lty = "dashed", col = "red")

rrtm <- tibble::tibble(
  refractive2 = refractive * refractive,
  t1 = tav_abs(90, refractive),
  t2 = tav_abs(40, refractive),
  t12 = t1 * t1,
  t22 = t2 * t2,
  tau = gsl::expint_E1(k)[,1],
  tau2 = tau * tau,
  x1 = 1 - t1,
  x2 = t12 * tau2 * (refractive2 - t1),
  x3 = t12 * tau * refractive2,
  x4 = refractive2 * refractive2 - tau2 * (refractive2 - t1) ^ 2,
  x5 = t2 / t1,
  x6 = x5 * (t1 - 1) + 1 - t2,
  r  = x1 + x2 / x4,
  t  = x3 / x4,
  ra = x5 * r + x6,
  ta = x5 * t,
  t_2 = t * t,
  r_2 = r * r,
  rt_2 = r_2 - t_2,
  rx2 = r * 2,
  delta = (t_2 - r_2 - 1) ^ 2 - 4 * r_2,
  delta_rt = sqrt(delta)
) %>%
  select(colnames(fortran))

fortran

map2_lgl(rrtm, fortran, ~isTRUE(all.equal(.x, .y)))

##################################################
# Fitting algorithm
##################################################
# Setup
library(rrtm)
library(purrr)

p4l <- function(params) {
  prospect4(params[1], params[2], params[3], params[4])[,,1]
}

true_param <- c(1.4, 40, 0.01, 0.01)
observed <- p4l(true_param)

r_prior <- function() {
  c(
    1 + rexp(1, 1),
    rlnorm(1, log(40), log(5)),
    rlnorm(1, log(0.01), 1),
    rlnorm(1, log(0.01), 1)
  )
}

curve(dexp(x - 1, 1), 1, 10)
curve(dlnorm(x, log(0.01), 1), 0, 0.05)

##################################################

solar_zenith <- 10
instrument_zenith <- 0
azimuth <- 0

p4 <- prospect4(1.4, 40, 0.01, 0.01)
rsoil <- hapke_soil(0.5)

out <- foursail(p4$r, p4$t,
                0, 0, 1,
                3, 0, 10, 0, 0, rsoil)

plot(400:2500, out[[1]], type = 'l')
lines(400:2500, out[[2]], col = 2)
lines(400:2500, out[[3]], col = 3)
lines(400:2500, out[[4]], col = 4)
