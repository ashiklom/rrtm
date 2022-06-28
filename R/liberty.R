#' LIBERTY model
#'
#' @param D Cell diameter. Average leaf cell diameter (m-6). Typical value 40 (20-100).
#' @param xu Intercellular air space. Determinant for the amount of radiative
#' flux passing between cells. Typical value 0.045 (0.01-0.1)
#' @param thick Leaf thickness. Arbitrary value to determine single leaf
#' reflectance and transmittance from infinite reflectance criteria. Typical
#' value 1.6 (1-10).
#' @param baseline Baseline absorption. Wavelength-independent absorption to
#' compensate for changes in absolute reflectance. Typical value 0.0006 for
#' fresh leaves or or 0.0004 for dry leaves.
#' @param element Albino absorption. Absorption in the visible region due to
#' lignin. Typical value 2 (0-4).
#' @param c_factor Chlorophyll content (mg m-2). Typical value 200 (0-600).
#' @param w_factor Water content (g m-2). Typical value 100 (0-500).
#' @param l_factor Lignin and cellulose content. Combined lignin and cellulose
#' content (mg m-2). Typical value 40 (10-80).
#' @param p_factor Nitrogen content (g m-2). Typical value 1  (0.3-2.0).
#' @return Length 4 list containing modeled reflectance, transmittance,
#' wavelengths, and R (?).
#' @author Alexey Shiklomanov
#' @export
liberty <- function(D, xu, thick, baseline, element,
                    c_factor, w_factor, l_factor, p_factor) {

  k_chloro <- dataspec_liberty[, 1]
  k_water <- dataspec_liberty[, 2]
  ke <- dataspec_liberty[, 3]
  k_ligcell <- dataspec_liberty[, 4]
  k_protein <- dataspec_liberty[, 5]

  nwl <- 420
  transR <- numeric(nwl)
  reflR <- numeric(nwl);
  wavelength <- numeric(nwl)
  RR <- numeric(nwl)

  for (i in seq(1, nwl)) {

    coeff <- D * (baseline +
                  (k_chloro[i] * c_factor) +
                  (k_water[i] * w_factor) +
                  (ke[i] * element) +
                  (k_ligcell[i] * l_factor) +
                  (k_protein[i] * p_factor))

    # change of refractive index over wavelength...
    N1 <- 1.4891 - (0.0005 * (i - 1))
    N0 <- 1.0
    in_angle <- 59

    # Index of refraction
    alpha <- in_angle * pi / 180
    beta <- asin((N0 / N1) * sin(alpha))
    para_r <- (tan(alpha - beta)) / (tan(alpha + beta))
    vert_r <- -(sin(alpha - beta)) / (sin(alpha + beta))

    me <- 0
    width <- pi / 180
    pp <- seq(1, 90)
    alpha <- pp * pi / 180
    beta <- asin(N0 / N1 * sin(alpha))
    plus <- alpha + beta
    dif <- alpha - beta
    refl <- 0.5 * ((sin(dif))^2 / (sin(plus))^2 + (tan(dif))^2 / (tan(plus))^2)
    me <- 2 * sum(refl * sin(alpha) * cos(alpha) * width)

    mi <- 0
    mint <- 0
    width <- pi / 180
    critical <- asin(N0 / N1) * 180 / pi
    jj <- seq(1, critical)
    alpha <- jj * pi / 180
    beta <- asin(N0 / N1 * sin(alpha))
    plus <- alpha + beta
    dif <- alpha - beta
    refl <- 0.5 * ((sin(dif))^2 / (sin(plus))^2 + (tan(dif))^2 / (tan(plus))^2)
    mint <- sum(refl * sin(alpha) * cos(alpha) * width)
    mi <- (1 - sin(critical * pi / 180)^2) + 2 * mint

    M <- 2 * (1 - (coeff + 1) * exp(-coeff)) / coeff^2
    T <- ((1 - mi) * M) / (1 - (mi * M))
    x <- xu / (1 - (1 - 2 * xu) * T)
    a <- me*T + x*T - me - T - x*me*T
    b <- 1 + x*me*T - 2*x^2*me^2*T
    c <- 2*me*x^2*T - x*T - 2*x*me
    R <- 0.5
    # Initial guess...
    for (iterations in seq(1, 50)) {
      next_R <- -(a * R^2 + c) / b
      R <- next_R
    }

    # The next bit works out transmittance based upon Benford...
    # setting up unchanging parameters...
    rb <- (2 * x * me) + (x * T) - (x * T * 2 * x * me)
    tb <- sqrt(((R - rb) * (1 - (R * rb))) / R)
    whole <- fix(thick)
    fraction <- thick - whole
    #
    # % 	/* The next bit works out the fractional value */
    # % 	/* for the interval between 1 and 2... */
    top <- tb^(1 + fraction) * ((((1 + tb)^2) - rb^2)^(1 - fraction))
    bot1 <- (1 + tb)^(2 * (1 - fraction)) - rb^2
    bot2 <- 1 + (64 / 3) * fraction * (fraction - 0.5) * (fraction - 1) * 0.001
    tif <- top / (bot1 * bot2)
    rif <- (1 + rb^2 - tb^2 - sqrt((1 + rb^2 - tb^2)^2 - 4 * rb^2 * (1 - tif^2))) / (2 * rb)

    # Now to work out for integral thickness greater than 2...
    if (whole >= 2) {
      prev_t <- 1
      prev_r <- 0
      for (step in seq(1, (whole - 1))) {
        cur_t <- (prev_t * tb) / (1 - (prev_r * rb))
        cur_r <- prev_r + (((prev_t * prev_t) * rb) / (1 - (prev_r * rb)))
        prev_t <- cur_t
        prev_r <- cur_r
      }
    } else {
      cur_t <- 1
      cur_r <- 0
    }
    trans <- (cur_t * tif) / (1 - (rif * cur_r))
    refl <- cur_r + ((cur_t^2 * rif) / (1 - (rif * cur_r)))

    transR[i] <- trans
    reflR[i] <- refl
    wavelength[i] <- 400 + ((i - 1) * 5)
    RR[i] <- R
  }

  LRT <- list(
    reflectance = reflR,
    transmittance = transR,
    wavelength = wavelength,
    RR = RR
  )
  LRT
}

# https://www.mathworks.com/help/matlab/ref/fix.html
fix <- function(x) trunc(x)
