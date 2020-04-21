#' Run ED two-stream radiative transfer model
#'
#' Wrapper around [sw_two_stream] that also calls PROSPECT 5 (see
#' [PEcAnRTM::prospect()]) to generate the leaf reflectance and
#' transmittance spectra, [hapke_soil()] to generate the soil spectrum
#' given `soil_moisture`, and a "flat" incident solar spectra for
#' direct and diffuse light based on `direct_sky_frac`.
#'
#' @inherit sw_two_stream params return
#' @inherit hapke_soil params
#' @param N PROSPECT 5 effective number of leaf mesophyll layers (>1; npft)
#' @param Cab PROSPECT 5 leaf chlorophyll content (ug cm-2; npft)
#' @param Car PROSPECT 5 leaf carotenoid content (ug cm-2; npft)
#' @param Cw PROSPECT 5 leaf water content (g cm-2; npft)
#' @param Cm PROSPECT 5 leaf dry matter content (g cm-2; npft)
#' @param direct_sky_frac Fraction of incident solar radiation that is
#'   direct (0-1; 0 = all diffuse radiation)
#' @param wood_reflect Wood reflectance spectrum. Default [wood_spectrum()]
#' @return
#' @author Alexey Shiklomanov
#' @export
edr_r <- function(pft, lai, wai, cai,
                  N, Cab, Car, Cw, Cm,
                  orient_factor, clumping_factor,
                  soil_moisture,
                  direct_sky_frac,
                  czen,
                  wood_reflect = matrix(rep(wood_spec, length(pft)), 2101),
                  wavelengths = seq(400, 2500)) {
  ncohort <- length(pft)
  npft <- length(N)
  nwl <- length(wavelengths)
  stopifnot(
    length(pft) == ncohort,
    length(lai) == ncohort, all(lai >= 0),
    length(wai) == ncohort, all(wai >= 0),
    length(cai) == ncohort,
    all(cai >= 0), all(cai <= 1),
    length(N) == npft, all(N >= 1),
    length(Cab) == npft, all(Cab > 0),
    length(Car) == npft, all(Car > 0),
    length(Cw) == npft, all(Cw > 0),
    length(Cm) == npft, all(Cm > 0),
    length(orient_factor) == npft,
    all(orient_factor < 1), all(orient_factor > -1),
    length(clumping_factor) == npft,
    all(clumping_factor > 0), all(clumping_factor <= 1),
    length(soil_moisture) == 1,
    soil_moisture <= 1, soil_moisture >= 0,
    length(direct_sky_frac) == 1,
    direct_sky_frac >= 0, direct_sky_frac <= 1,
    length(czen) == 1,
    NROW(wood_reflect) %in% c(2101, nwl)
  )

  # Wavelength indices -- everything relative to 400:2500 (so 400nm is
  # index 1, 2500 is index 2101)
  wli <- wavelengths - 399

  # If using full wood reflectance spectrum, subset to only used
  # wavelengths
  if (nwl != NROW(wood_reflect)) wood_reflect <- wood_reflect[wli, ]

  leaf_spectra <- mapply(prospect5, N, Cab, Car, Cw, Cm, SIMPLIFY = FALSE)
  leaf_reflect <- Reduce(
    cbind,
    Map(function(x) x[["reflectance"]], leaf_spectra)
  )[wli,]
  leaf_trans <- Reduce(
    cbind,
    Map(function(x) x[["transmittance"]], leaf_spectra)
  )[wli,]

  # Soil reflectance as a function of soil moisture
  soil_reflect <- hapke_soil(soil_moisture)[wli]

  # "Flat" spectra of incident solar radiation
  down0_sky <- rep(direct_sky_frac, nwl)
  down_sky <- rep(1 - direct_sky_frac, nwl)

  # Wood does not transmit in the VIS or NIR
  wood_trans <- wood_reflect
  wood_trans[] <- 0

  sw_two_stream(
    czen = czen,
    iota_g = soil_reflect,
    pft = pft,
    lai = lai,
    wai = wai,
    cai = cai,
    orient_factor = orient_factor,
    clumping_factor = clumping_factor,
    leaf_reflect = leaf_reflect,
    leaf_trans = leaf_trans,
    wood_reflect = wood_reflect,
    wood_trans = wood_trans,
    down_sky = down_sky,
    down0_sky = down0_sky,
    wavelengths = wavelengths
  )
}

#' R implementation of ED2 two-stream radiation model
#'
#' @param czen Cosine of angle of incidence (`cosaio`)
#' @param iota_g Ground (soil + snow) albedo (nwl)
#' @param pft PFT identities of each cohort, as integer (ncoh)
#' @param lai Leaf area index of each cohort (ncoh)
#' @param wai Wood area index of each cohort (ncoh)
#' @param cai Crown area of each cohort (ncoh)
#' @param orient_factor Orient factor (npft)
#' @param clumping_factor Clumping factor (npft)
#' @param leaf_reflect Leaf reflectance spectra (nwl * npft)
#' @param leaf_trans Leaf transmittance spectra (nwl * npft)
#' @param wood_reflect Wood reflectance spectra (nwl * npft)
#' @param wood_trans Wood transmittance spectra (nwl * npft)
#' @param down_sky Normalized diffuse solar spectrum (nwl)
#' @param down0_sky Normalized direct solar spectrum (nwl)
#' @param wavelengths Numeric vector of wavelengths to use, in nm
#'   (nwl). Default is 400:2500.
#' @return
#' @author Alexey Shiklomanov
#' @useDynLib rrtm
#' @export
sw_two_stream <- function(czen,
                          iota_g,
                          pft,
                          lai, wai, cai,
                          orient_factor, clumping_factor,
                          leaf_reflect, leaf_trans,
                          wood_reflect, wood_trans,
                          down_sky, down0_sky,
                          wavelengths = seq(400, 2500)
                          ) {

  # Sanity checks
  nwl <- length(wavelengths)
  stopifnot(
    NROW(leaf_reflect) == nwl,
    NROW(leaf_trans) == nwl,
    NROW(wood_reflect) == nwl,
    NROW(wood_trans) == nwl,
    NROW(iota_g) == nwl
  )

  ncoh <- length(pft)
  stopifnot(
    length(lai) == ncoh,
    length(wai) == ncoh,
    length(cai) == ncoh
  )

  stopifnot(abs(czen) <= 1)

  # All spectra are between 0 and 1
  stopifnot(
    !any(iota_g > 1 | iota_g < 0),
    !any(leaf_reflect > 1 | leaf_reflect < 0),
    !any(leaf_trans > 1 | leaf_trans < 0),
    !any(wood_reflect > 1 | wood_reflect < 0),
    !any(wood_trans > 1 | wood_trans < 0)
  )

  # Incident radiation has to sum to 1 across all wavelengths
  stopifnot(
    !any(down_sky > 1 | down_sky < 0),
    !any(down0_sky > 1 | down0_sky < 0),
    all(down_sky + down0_sky == 1)
  )

  ##########
  # Unpack PFT-specific parameters across cohorts for easier vectorization
  leaf_reflect <- leaf_reflect[, pft]
  leaf_trans <- leaf_trans[, pft]
  wood_reflect <- wood_reflect[, pft]
  wood_trans <- wood_trans[, pft]
  orient_factor <- orient_factor[pft]
  clumping_factor <- clumping_factor[pft]

  # Calculations from sfc_rad
  leaf_scatter <- leaf_reflect + leaf_trans
  wood_scatter <- wood_reflect + wood_trans
  leaf_backscatter <- (leaf_scatter +
                         0.25 * (leaf_reflect - leaf_trans) *
                         (1 + orient_factor) ^ 2) / (2 * leaf_scatter)
  wood_backscatter <- (wood_scatter +
                         0.25 * (wood_reflect - wood_trans) *
                         (1 + orient_factor) ^ 2) / (2 * wood_scatter)
  phi1 <- 0.5 - orient_factor * (0.633 + 0.33 * orient_factor)
  phi2 <- 0.877 * (1 - 2 * phi1)
  mu_bar <- (1 - phi1 * log(1 + phi2 / phi1) / phi2) / phi2
  mu_bar[orient_factor == 0] <- 1
  stopifnot(all(is.finite(mu_bar)))
  ##########

  # Size of solution matrix
  nsiz <- 2 * ncoh + 2

  # Loop over bands -- skipping because this is already vectorized

  # Loop over cohorts -- vectorizing here
  elai <- clumping_factor * lai
  etai <- elai + wai

  leaf_weight <- elai / etai
  wood_weight <- 1 - leaf_weight

  # Inverse optical depth of direct radiation
  proj_area <- phi1 + phi2 * czen
  mu0 <- -etai / log((1 - cai) + cai * exp(-proj_area * etai / (cai * czen)))

  # Inverse optical depth of diffuse radiation
  mu <- -etai / log((1 - cai) + cai * exp(-etai / mu_bar))

  # Backscatter coefficients for diffuse radiation
  iota_ratio <- 1 / (2 * (1 + phi2 * mu0)) *
    (1 - phi1 * mu0 / (1 + phi2 * mu0) *
       log((1 + (phi1 + phi2) * mu0) / (phi1 * mu0)))
  stopifnot(all(is.finite(iota_ratio)))
  beta0 <- iota_ratio * (1 + mu0 / mu)
  stopifnot(all(is.finite(beta0)))
  epsil0 <- 1 - 2 * beta0

  # Transmissivity of direct radiation
  expm0_minus <- exp(-etai / mu0)

  #################
  # Define boundary conditions
  #################
  z <- ncoh + 1
  elai[z] <- 0
  etai[z] <- 0
  leaf_weight[z] <- 0.5
  wood_weight[z] <- 0.5
  proj_area[z] <- 0.5
  mu0[z] <- czen / proj_area[z]
  mu[z] <- 1
  iota_ratio[z] <- 0.5 * (1 - 0.5 * mu0[z] *
                            log(1 / (0.5 * mu0[z]) + 1))
  beta0[z] <- iota_ratio[z] * (mu0[z] + mu[z]) / mu[z]
  epsil0[z] <- 1 - 2 * beta0[z]
  expm0_minus[z] <- 1

  # Direct radiation profile via exponential attentuation
  # TODO: This has to be a matrix because Down0 is a spectrum.
  down0 <- matrix(0, nwl, z)
  down0[, z] <- down0_sky
  for (j in seq(ncoh, 1)) {
    down0[, j] <- down0[, j + 1] * expm0_minus[j]
  }

  # Convert scalar quantities to matrices, so I can do elementwise multiplication
  vec2mat <- function(x) matrix(rep(x, nwl), nrow = nwl, byrow = TRUE)
  etai <- vec2mat(etai)
  leaf_weight <- vec2mat(leaf_weight)
  wood_weight <- vec2mat(wood_weight)
  mu <- vec2mat(mu)
  mu0 <- vec2mat(mu0)
  epsil0 <- vec2mat(epsil0)
  expm0_minus <- vec2mat(expm0_minus)

  # Diffuse radiation properties
  # All of these are wavelength-dependent quantities (nwl x ncoh)
  iota <- leaf_weight[, -z] * leaf_scatter + wood_weight[, -z] * wood_scatter
  beta <- leaf_weight[, -z] * leaf_backscatter + wood_weight[, -z] * wood_backscatter
  epsil <- 1 - 2 * beta
  lambda <- sqrt((1 - epsil * iota) * (1 - iota)) / mu[, -z]

  # Ancillary variables for right-hand side
  iota_mu <- iota / mu[, -z]
  iota_mu0 <- iota / mu0[, -z]
  down0_mu0 <- down0[, -1] / mu0[, -z]
  mu02 <- mu0[, -z] ^ 2
  lambda2 <- lambda ^ 2

  a_aux <- -((1 - epsil * iota) * iota_mu + epsil0[, -z] * iota_mu0) * down0_mu0
  s_aux <- -((1 - iota) * epsil0[, -z] * iota_mu + iota_mu0) * down0_mu0
  delta <- (a_aux + s_aux) * mu02 / (2 * (1 - lambda2 * mu02))
  upsilon <- (a_aux - s_aux) * mu02 / (2 * (1 - lambda2 * mu02))

  # Upwelling and downwelling radiation
  iez <- sqrt((1 - iota) / (1 - epsil * iota))
  gamm_plus <- 0.5 * (1 + iez)
  gamm_minus <- 0.5 * (1 - iez)

  # Transmissivity of diffuse light
  expl_plus <- exp(lambda * etai[, -z])
  expl_minus <- exp(-lambda * etai[, -z])

  # Define boundary conditions for above
  iota <- cbind(iota, rep(1, nwl))
  beta <- cbind(beta, rep(0, nwl))
  epsil <- cbind(epsil, 1 - 2 * beta[, z])
  lambda <- cbind(lambda, rep(0, nwl))
  a_aux <- cbind(a_aux, -epsil0[, z] * down0_sky / (mu0[, z] ^ 2))
  s_aux <- cbind(s_aux, -iota[, z] * down0_sky / (mu0[, z] ^ 2))
  delta <- cbind(delta, 0.5 * (a_aux[, z] + s_aux[, z]) * mu0[, z] ^ 2)
  upsilon <- cbind(upsilon, 0.5 * (a_aux[, z] - s_aux[, z]) * mu0[, z] ^ 2)
  gamm_plus <- cbind(gamm_plus, rep(1, nwl))
  gamm_minus <- cbind(gamm_minus, rep(0, nwl))
  expl_plus <- cbind(expl_plus, rep(1, nwl))
  expl_minus <- cbind(expl_minus, rep(1, nwl))

  mmat <- array(0, c(nsiz, nsiz, nwl))
  yvec <- matrix(0, nsiz, nwl)

  # Bottom (1) and top boundary conditions
  mmat[1, 1, ] <- (gamm_minus[, 1] - iota_g * gamm_plus[, 1]) * expl_minus[, 1]
  mmat[1, 2, ] <- (gamm_plus[, 1] - iota_g * gamm_minus[, 1]) * expl_plus[, 1]
  mmat[nsiz, nsiz-1, ] <- gamm_plus[, z]
  mmat[nsiz, nsiz, ] <- gamm_minus[, z]
  yvec[1, ] <- iota_g * down0[, 1] - (upsilon[, 1] - iota_g * delta[, 1]) * expm0_minus[1]
  yvec[nsiz, ] <- down_sky - delta[, z]

  for (k in seq_len(ncoh)) {
    kp1 <- k + 1
    k2 <- 2 * k
    k2m1 <- k2 - 1
    k2p1 <- k2 + 1
    k2p2 <- k2 + 2

    yvec[k2, ] <- delta[, kp1] * expm0_minus[, kp1] - delta[, k]
    yvec[k2p1, ] <- upsilon[, kp1] * expm0_minus[, kp1] - upsilon[, k]

    mmat[k2, k2m1, ] <- gamm_plus[, k]
    mmat[k2, k2, ] <- gamm_minus[, k]
    mmat[k2, k2p1, ] <- -gamm_plus[, kp1] * expl_minus[, kp1]
    mmat[k2, k2p2, ] <- -gamm_minus[, kp1] * expl_plus[, kp1]
    mmat[k2p1, k2m1, ] <- gamm_minus[, k]
    mmat[k2p1, k2, ] <- gamm_plus[, k]
    mmat[k2p1, k2p1, ] <- -gamm_minus[, kp1] * expl_minus[, kp1]
    mmat[k2p1, k2p2, ] <- -gamm_plus[, kp1] * expl_plus[, kp1]
  }

  stopifnot(is.finite(sum(mmat)))

  # Solve the radiation balance at each wavelength
  xvec_t <- solvearray(mmat, yvec)
  xvec <- t(xvec_t)

  # Store the solution in matrices (nwl x (ncoh + 1))
  down <- matrix(0, nwl, ncoh + 1)
  up <- matrix(0, nwl, ncoh + 1)
  for (k in seq_len(ncoh + 1)) {
    k2 <- 2 * k
    k2m1 <- k2 - 1
    down[, k] <- xvec[, k2m1] * gamm_plus[, k] * expl_minus[, k] +
      xvec[, k2] * gamm_minus[, k] * expl_plus[, k] +
      delta[, k] * expm0_minus[, k]
    up[, k] <- xvec[, k2m1] * gamm_minus[, k] * expl_minus[, k] +
      xvec[, k2] * gamm_plus[, k] * expl_plus[, k] +
      upsilon[, k] * expm0_minus[, k]
  }

  # Integrate light levels
  light_level <- matrix(0, nwl, ncoh)
  light_beam_level <- matrix(0, nwl, ncoh)
  light_diff_level <- matrix(0, nwl, ncoh)

  for (k in seq_len(ncoh)) {
    kp1 <- k + 1
    light_level[, k] <- light_level[, k] + 0.5 * (down[, k] + down[, kp1] + down0[, k] + down0[, kp1])
    light_beam_level[, k] <- light_beam_level[, k] + 0.5 * (down0[, k] + down0[, kp1])
    light_diff_level[, k] <- light_diff_level[, k] + 0.5 * (down[, k] + down[, kp1])
  }

  # Albedo is "up" at top of canopy
  list(
    albedo = up[, z],                    # Albedo is upwelling radiation profile at top of canopy (nwl)
    up = up,                             # Upwelling radiation profile, by cohort + top of canopy (nwl x (ncoh + 1))
    down = down,                         # Downwelling radiation, by cohort + top of canopy (nwl x (ncoh + 1))
    light_level = light_level,           # Light level, by cohort (nwl x ncoh)
    light_beam_level = light_beam_level, # Direct light level, by cohort (nwl x ncoh)
    light_diff_level = light_diff_level  # Diffuse light level, by cohort (nwl x ncoh)
  )
}
