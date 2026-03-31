devtools::load_all()

e1_me <- prospect4(1.4, 40, 0.01, 0.01)
e1_gsl <- prospect4(1.4, 40, 0.01, 0.01, e1fun = gsl::expint_E1)
e1_expint <- prospect4(1.4, 40, 0.01, 0.01, e1fun = expint::expint_E1)

if (any(something <- c(TRUE, TRUE, FALSE))) {
  print("true")
} else {
  print("false")
}
print(something)

par(mfrow = c(2, 2))
plot(e1_me$reflectance, e1_gsl$reflectance)
plot(e1_me$transmittance, e1_gsl$transmittance)
plot(e1_me$reflectance - e1_gsl$reflectance)
abline(h = 0)
plot(e1_me$transmittance - e1_gsl$transmittance)
abline(h = 0)

par(mfrow = c(2, 2))
plot(e1_me$reflectance, e1_expint$reflectance)
plot(e1_me$transmittance, e1_expint$transmittance)
plot(e1_me$reflectance - e1_expint$reflectance)
abline(h = 0)
plot(e1_me$transmittance - e1_expint$transmittance)
abline(h = 0)

par(mfrow = c(2, 2))
plot(e1_gsl$reflectance, e1_expint$reflectance)
plot(e1_gsl$transmittance, e1_expint$transmittance)
plot(e1_gsl$reflectance - e1_expint$reflectance)
abline(h = 0)
plot(e1_gsl$transmittance - e1_expint$transmittance)
abline(h = 0)
