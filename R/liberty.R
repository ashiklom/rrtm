#' LIBERTY model
#'
#' Leaf radiative transfer model designed for conifer needles (Dawson et al.
#' 1998).
#'
#' `liberty()` and `liberty_dawson()` differ only in the equation for the
#' infinite reflectance R (Dawson et al. Eq. 16). `liberty_dawson()` uses Eq. 16
#' as published. Eq. 16 is meant to be the closed form of the series in Eq. 15,
#' but it places the term `2 x me` over the common denominator without
#' multiplying it by that denominator. The correct sum is
#'
#' `R = 2 x me + x (1 - 2 x me) tau (1 - me R) / ((1 - me R) - (1 - x)(1 - me) tau R)`,
#'
#' which `liberty()` uses. The series of Eq. 15 describes a surface layer of
#' cells that loses light only by absorption, so R always lies in \[0, 1\] and
#' tends to 1 as absorption vanishes. Eq. 16 behaves as if the layer returned
#' more light than it received (`me (1 + me) / 2` more at zero absorption): for
#' strong scattering (large `xu`) and weak absorption it has no real solution,
#' and `liberty_dawson()` diverges. Everything else is shared, and both use the
#' internal reflectance of Eq. 2 for light leaving the cell (earlier versions of
#' this package used the reflectance for light entering it).
#'
#' Note that these coefficients are based on the original values from Dawson et
#' al. (1998). Improved specific absorption coefficients are available from
#' later work by Di Vittorio (2009).
#'
#' @references Dawson, T. P., Curran, P. J., & Plummer, S. E. (1998). LIBERTY—Modeling the
#' Effects of Leaf Biochemical Concentration on Reflectance Spectra. Remote
#' Sensing of Environment, 65(1), 50–60.
#' https://doi.org/10.1016/S0034-4257(98)00007-8
#' @references Di Vittorio, A. V. (2009). Enhancing a leaf radiative transfer
#' model to estimate concentrations and in vivo specific absorption
#' coefficients of total carotenoids and chlorophylls a and b from
#' single-needle reflectance and transmittance. Remote Sensing of Environment,
#' 113(9), 1948–1966. https://doi.org/10.1016/j.rse.2009.05.002
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
#' wavelengths, and infinite reflectance (`RR`).
#' @author Alexey Shiklomanov
#' @export
liberty <- function(D, xu, thick, baseline, element,
                    c_factor, w_factor, l_factor, p_factor) {
  liberty_impl(D, xu, thick, baseline, element,
               c_factor, w_factor, l_factor, p_factor, eq16 = FALSE)
}

#' @rdname liberty
#' @export
liberty_dawson <- function(D, xu, thick, baseline, element,
                           c_factor, w_factor, l_factor, p_factor) {
  liberty_impl(D, xu, thick, baseline, element,
               c_factor, w_factor, l_factor, p_factor, eq16 = TRUE)
}

liberty_impl <- function(D, xu, thick, baseline, element,
                         c_factor, w_factor, l_factor, p_factor, eq16) {

  k_chloro <- dataspec_liberty[, 1]
  k_water <- dataspec_liberty[, 2]
  ke <- dataspec_liberty[, 3]
  k_ligcell <- dataspec_liberty[, 4]
  k_protein <- dataspec_liberty[, 5]

  nwl <- 420
  transR <- numeric(nwl)
  reflR <- numeric(nwl)
  wavelength <- numeric(nwl)
  RR <- numeric(nwl)

  # change of refractive index over wavelength...
  N1 <- 1.4891 - (0.0005 * (seq(1, nwl) - 1))
  # Diffuse reflectance of the cell surface for light arriving from the air
  # (me, Eq. 1) and from inside the cell (mi, Eq. 2), in closed form
  t12 <- tav_abs(90, N1)
  me_all <- 1 - t12
  mi_all <- 1 - t12 / N1^2

  for (i in seq(1, nwl)) {

    coeff <- D * (baseline +
                  (k_chloro[i] * c_factor) +
                  (k_water[i] * w_factor) +
                  (ke[i] * element) +
                  (k_ligcell[i] * l_factor) +
                  (k_protein[i] * p_factor))

    me <- me_all[i]
    mi <- mi_all[i]

    M <- 2 * (1 - (coeff + 1) * exp(-coeff)) / coeff^2
    Tt <- ((1 - mi) * M) / (1 - (mi * M))
    x <- xu / (1 - (1 - 2 * xu) * Tt)
    # Reflectance of a single layer (Eq. 17)
    rb <- (2 * x * me) + (x * Tt) - (x * Tt * 2 * x * me)

    if (eq16) {
      a <- me*Tt + x*Tt - me - Tt - x*me*Tt
      b <- 1 + x*me*Tt - 2*x^2*me^2*Tt
      c <- 2*me*x^2*Tt - x*Tt - 2*x*me
      R <- 0.5
      # Initial guess...
      for (iterations in seq(1, 50)) {
        next_R <- -(a * R^2 + c) / b
        R <- next_R
      }
      # Transmittance of a single layer (Eq. 19)
      tb <- sqrt(((R - rb) * (1 - (R * rb))) / R)
    } else {
      # Eq. 15 is R = rb + td * tu * R / (1 - ru * R): adding the surface
      # layer of cells to leaf material of reflectance R, where rb, td are its
      # reflectance and transmittance for light from above, and ru, tu for
      # light from below. Take the root in [0, 1].
      td <- (1 - x) * (1 - 2 * x * me) * Tt
      ru <- me + (1 - x) * (1 - me) * Tt
      tu <- x * (1 - me) * Tt
      bb <- 1 + rb * ru - td * tu
      R <- 2 * rb / (bb + sqrt(max(bb^2 - 4 * ru * rb, 0)))
      # Eq. 19, using R - rb = td * tu * R / (1 - ru * R)
      tb <- sqrt(td * tu * (1 - R * rb) / (1 - R * ru))
    }

    # The next bit works out transmittance based upon Benford...
    # setting up unchanging parameters...
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
