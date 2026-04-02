#!/usr/bin/env Rscript

check_requirements <- function() {
  required <- c("nimble")
  missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "))
  }
  library(nimble)
}

get_script_path <- function() {
  file_arg <- "--file="
  matches <- grep(paste0("^", file_arg), commandArgs(FALSE), value = TRUE)
  if (length(matches) == 0) {
    return(NULL)
  }
  normalizePath(sub(file_arg, "", matches[[1]]), mustWork = FALSE)
}

get_project_root <- function() {
  script_path <- get_script_path()
  if (is.null(script_path)) {
    return(normalizePath(getwd()))
  }
  normalizePath(file.path(dirname(script_path), "..", ".."))
}

load_rrtm_namespace <- function(project_root) {
  if (
    requireNamespace("pkgload", quietly = TRUE) &&
      file.exists(file.path(project_root, "DESCRIPTION"))
  ) {
    pkgload::load_all(
      project_root,
      export_all = FALSE,
      helpers = FALSE,
      quiet = TRUE
    )
  }
  if (!"rrtm" %in% loadedNamespaces()) {
    loadNamespace("rrtm")
  }
  asNamespace("rrtm")
}

get_rrtm_object <- function(rrtm_ns, name) {
  get(name, envir = rrtm_ns, inherits = FALSE)
}

read_observation <- function() {
  obs <- utils::read.csv(
    "https://raw.githubusercontent.com/ashiklom/spectra_db/master/nasa_fft/spectra/nasa_fft%7CAK06_ABBA_M2%7C2009.csvy",
    comment.char = "#",
    col.names = c("wavelength", "observed")
  )
  prospect_wl <- seq(400, 2500, 1)
  obs_prospect <- stats::approx(
    obs[["wavelength"]],
    obs[["observed"]],
    prospect_wl
  )[["y"]]
  list(obs = obs, prospect_wl = prospect_wl, obs_prospect = obs_prospect)
}

default_inits <- function() {
  list(
    N_minus_1 = 0.4,
    Cab = 40,
    Car = 10,
    Cbrown = 0.1,
    Canth = 10,
    Cw = 0.01,
    Cm = 0.01,
    rsd = 0.05
  )
}

parameter_names <- function() {
  c("N", "Cab", "Car", "Cbrown", "Canth", "Cw", "Cm", "rsd")
}

invalid_reflectance <- function(n_wl) {
  rep(-1, n_wl)
}

PROSPECT_WL <- seq(400, 2500, 1)
PROSPECTD_FUN <- NULL

prospectd_reflectance_direct <- function(N, Cab, Car, Cbrown, Canth, Cw, Cm) {
  tryCatch(
    {
      spec <- withCallingHandlers(
        PROSPECTD_FUN(N, Cab, Car, Cbrown, Canth, Cw, Cm),
        warning = function(w) {
          invokeRestart("muffleWarning")
        }
      )
      refl <- spec[["reflectance"]]
      if (length(refl) != length(PROSPECT_WL)) {
        return(invalid_reflectance(length(PROSPECT_WL)))
      }
      if (any(!is.finite(refl)) || any(refl < 0)) {
        return(invalid_reflectance(length(PROSPECT_WL)))
      }
      refl
    },
    error = function(e) {
      invalid_reflectance(length(PROSPECT_WL))
    }
  )
}

prospectd_rcall <- nimble::nimbleRcall(
  prototype = function(
    N = double(0),
    Cab = double(0),
    Car = double(0),
    Cbrown = double(0),
    Canth = double(0),
    Cw = double(0),
    Cm = double(0)
  ) {},
  returnType = double(1),
  Rfun = "prospectd_reflectance_direct"
)

prospectd_rewrite_nf <- nimble::nimbleFunction(
  run = function(
    N = double(0),
    Cab = double(0),
    Car = double(0),
    Cbrown = double(0),
    Canth = double(0),
    Cw = double(0),
    Cm = double(0),
    dataspec = double(2),
    talf = double(1),
    t12 = double(1),
    t21 = double(1)
  ) {
    returnType(double(1))

    n_wl <- dim(dataspec)[1]
    invalid <- numeric(length = n_wl, init = FALSE)
    refl <- numeric(length = n_wl, init = FALSE)
    conc <- numeric(length = 6, init = FALSE)

    for (i in 1:n_wl) {
      invalid[i] <- -1.0
    }
    if (N <= 1.0) {
      return(invalid)
    }

    conc[1] <- Cab / N
    conc[2] <- Car / N
    conc[3] <- Canth / N
    conc[4] <- Cbrown / N
    conc[5] <- Cw / N
    conc[6] <- Cm / N

    for (i in 1:n_wl) {
      k_i <- 0.0
      for (j in 1:6) {
        k_i <- k_i + dataspec[i, j] * conc[j]
      }

      trans_i <- 1.0
      if (k_i > 0.0) {
        A <- log((0.56146 / k_i + 0.65) * (1.0 + k_i))
        B <- k_i^4 * exp(7.7 * k_i) * (2.0 + k_i)^3.7
        e1_i <- (A^(-7.7) + B)^(-0.13)
        trans_i <- (1.0 - k_i) * exp(-k_i) + k_i^2 * e1_i
      }
      if (trans_i < 0.0) {
        trans_i <- 0.0
      }

      ralf <- 1.0 - talf[i]
      r12 <- 1.0 - t12[i]
      r21 <- 1.0 - t21[i]

      denom <- 1.0 - (r21 * r21) * (trans_i * trans_i)
      if (abs(denom) < 1e-12) {
        return(invalid)
      }

      Ta <- talf[i] * trans_i * t21[i] / denom
      Ra <- ralf + r21 * trans_i * Ta
      tt_i <- t12[i] * trans_i * t21[i] / denom
      rr_i <- r12 + r21 * trans_i * tt_i

      Tsub <- 0.0
      Rsub <- 0.0
      if (rr_i == 0.0) {
        Tsub <- 0.0
        Rsub <- 0.0
      } else if (tt_i == 0.0) {
        Tsub <- 0.0
        Rsub <- 0.0
      } else if (rr_i + tt_i >= 1.0) {
        Tsub <- tt_i / (tt_i + (1.0 - tt_i) * (N - 1.0))
        Rsub <- 1.0 - Tsub
      } else {
        disc <- (1.0 + rr_i + tt_i) *
          (1.0 + rr_i - tt_i) *
          (1.0 - rr_i + tt_i) *
          (1.0 - rr_i - tt_i)
        if (disc < 0.0) {
          return(invalid)
        }

        D <- sqrt(disc)
        if (abs(rr_i) < 1e-12) {
          return(invalid)
        }
        if (abs(tt_i) < 1e-12) {
          return(invalid)
        }

        va <- (1.0 + rr_i * rr_i - tt_i * tt_i + D) / (2.0 * rr_i)
        vb <- (1.0 - rr_i * rr_i + tt_i * tt_i + D) / (2.0 * tt_i)
        vbNN <- vb^(N - 1.0)
        vbNN2 <- vbNN * vbNN
        va2 <- va * va
        denomx <- va2 * vbNN2 - 1.0
        if (abs(denomx) < 1e-12) {
          return(invalid)
        }

        Rsub <- va * (vbNN2 - 1.0) / denomx
        Tsub <- vbNN * (va2 - 1.0) / denomx
      }

      denomy <- 1.0 - Rsub * rr_i
      if (abs(denomy) < 1e-12) {
        return(invalid)
      }

      refl_i <- Ra + Ta * Rsub * tt_i / denomy
      if (refl_i < 0.0) {
        return(invalid)
      }
      if (refl_i > 1.0) {
        return(invalid)
      }
      refl[i] <- refl_i
    }

    return(refl)
  }
)

direct_code <- nimble::nimbleCode({
  N_minus_1 ~ dexp(1)
  N <- N_minus_1 + 1
  Cab ~ dlnorm(log(40), 0.8)
  Car ~ dlnorm(log(10), 0.8)
  Cbrown ~ dexp(0.2)
  Canth ~ dlnorm(log(10), 0.8)
  Cw ~ dlnorm(-4.456, 1.216)
  Cm ~ dlnorm(-5.15, 1.328)
  rsd ~ dexp(10)
  mu[1:n_wl] <- prospectd_rcall(N, Cab, Car, Cbrown, Canth, Cw, Cm)
  for (i in 1:n_wl) {
    obs[i] ~ dnorm(mu[i], sd = rsd)
  }
})

rewrite_code <- nimble::nimbleCode({
  N_minus_1 ~ dexp(1)
  N <- N_minus_1 + 1
  Cab ~ dlnorm(log(40), 0.8)
  Car ~ dlnorm(log(10), 0.8)
  Cbrown ~ dexp(0.2)
  Canth ~ dlnorm(log(10), 0.8)
  Cw ~ dlnorm(-4.456, 1.216)
  Cm ~ dlnorm(-5.15, 1.328)
  rsd ~ dexp(10)
  mu[1:n_wl] <- prospectd_rewrite_nf(
    N,
    Cab,
    Car,
    Cbrown,
    Canth,
    Cw,
    Cm,
    dataspec_pd[1:n_wl, 1:6],
    pd_talf[1:n_wl],
    pd_t12[1:n_wl],
    pd_t21[1:n_wl]
  )
  for (i in 1:n_wl) {
    obs[i] ~ dnorm(mu[i], sd = rsd)
  }
})

build_constants <- function(obs_prospect, lookup, variant) {
  constants <- list(n_wl = length(obs_prospect))
  if (identical(variant, "rewrite")) {
    constants$dataspec_pd <- lookup$dataspec_pd
    constants$pd_talf <- lookup$pd_talf
    constants$pd_t12 <- lookup$pd_t12
    constants$pd_t21 <- lookup$pd_t21
  }
  constants
}

build_model <- function(code, constants, obs_prospect, inits) {
  nimble::nimbleModel(
    code,
    constants = constants,
    data = list(obs = obs_prospect),
    inits = inits,
    check = FALSE,
    calculate = FALSE
  )
}

run_variant_mcmc <- function(
  code,
  constants,
  obs_prospect,
  inits,
  iterations,
  burnin
) {
  model <- build_model(code, constants, obs_prospect, inits)
  model$calculate()
  conf <- nimble::configureMCMC(
    model,
    monitors = c("N", "Cab", "Car", "Cbrown", "Canth", "Cw", "Cm", "rsd"),
    useConjugacy = FALSE
  )
  mcmc <- nimble::buildMCMC(conf)
  c_model <- nimble::compileNimble(model)
  c_mcmc <- nimble::compileNimble(mcmc, project = model)
  samples <- nimble::runMCMC(
    c_mcmc,
    niter = iterations,
    nburnin = burnin,
    samplesAsCodaMCMC = FALSE,
    summary = FALSE,
    setSeed = TRUE
  )
  list(
    model = model,
    compiled_model = c_model,
    compiled_mcmc = c_mcmc,
    samples = samples
  )
}

normalize_samples <- function(samples) {
  samples <- as.matrix(samples)
  if (!"N" %in% colnames(samples) && "N_minus_1" %in% colnames(samples)) {
    samples <- cbind(samples, N = samples[, "N_minus_1"] + 1)
  }
  samples[, parameter_names(), drop = FALSE]
}

call_rewrite_forward <- function(params, lookup) {
  refl <- prospectd_rewrite_nf(
    params[["N"]],
    params[["Cab"]],
    params[["Car"]],
    params[["Cbrown"]],
    params[["Canth"]],
    params[["Cw"]],
    params[["Cm"]],
    lookup$dataspec_pd,
    lookup$pd_talf,
    lookup$pd_t12,
    lookup$pd_t21
  )
  as.numeric(refl)
}

call_direct_forward <- function(params) {
  prospectd_reflectance_direct(
    params[["N"]],
    params[["Cab"]],
    params[["Car"]],
    params[["Cbrown"]],
    params[["Canth"]],
    params[["Cw"]],
    params[["Cm"]]
  )
}

summarize_predictions <- function(samples, forward_fun, summary_samples) {
  n_keep <- min(nrow(samples), summary_samples)
  keep <- samples[seq_len(n_keep), , drop = FALSE]
  predmat <- matrix(NA_real_, n_keep, length(PROSPECT_WL))

  for (i in seq_len(n_keep)) {
    mu <- forward_fun(stats::setNames(
      as.list(keep[i, parameter_names()[1:7]]),
      parameter_names()[1:7]
    ))
    predmat[i, ] <- stats::rnorm(length(mu), mean = mu, sd = keep[i, "rsd"])
  }

  list(
    pred_mean = colMeans(predmat, na.rm = TRUE),
    pred_lo = apply(predmat, 2, stats::quantile, probs = 0.05, na.rm = TRUE),
    pred_hi = apply(predmat, 2, stats::quantile, probs = 0.95, na.rm = TRUE)
  )
}

compare_forward_models <- function(params, lookup) {
  direct <- call_direct_forward(params)
  rewrite <- call_rewrite_forward(params, lookup)
  diff <- rewrite - direct
  list(
    max_abs_diff = max(abs(diff), na.rm = TRUE),
    rmse = sqrt(mean(diff^2, na.rm = TRUE)),
    direct = direct,
    rewrite = rewrite
  )
}

print_variant_summary <- function(name, result, comparison = NULL) {
  samples <- normalize_samples(result$samples)
  posterior_mean <- colMeans(samples)
  message(sprintf("[%s] posterior means:", name))
  print(round(posterior_mean, 4))
  if (!is.null(comparison)) {
    message(sprintf(
      "[%s] rewrite-vs-direct forward comparison: max_abs_diff=%.6f rmse=%.6f",
      name,
      comparison$max_abs_diff,
      comparison$rmse
    ))
  }
}

run_variant <- function(
  name,
  obs_prospect,
  lookup,
  iterations,
  burnin,
  summary_samples
) {
  constants <- build_constants(obs_prospect, lookup, name)
  inits <- default_inits()

  message(sprintf("Running %s variant", name))
  code <- if (identical(name, "direct")) direct_code else rewrite_code
  result <- run_variant_mcmc(
    code,
    constants,
    obs_prospect,
    inits,
    iterations,
    burnin
  )
  samples <- normalize_samples(result$samples)

  comparison <- compare_forward_models(
    stats::setNames(
      as.list(colMeans(samples)[parameter_names()[1:7]]),
      parameter_names()[1:7]
    ),
    lookup
  )

  summary <- if (identical(name, "direct")) {
    summarize_predictions(samples, call_direct_forward, summary_samples)
  } else {
    summarize_predictions(
      samples,
      function(params) call_rewrite_forward(params, lookup),
      summary_samples
    )
  }

  print_variant_summary(name, result, comparison)
  list(samples = samples, comparison = comparison, prediction_summary = summary)
}

parse_args <- function(args) {
  options <- list(
    variant = "both",
    iterations = 10000,
    burnin = 2000,
    summary_samples = 25
  )
  if (length(args) == 0) {
    return(options)
  }

  for (arg in args) {
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- pieces[[1]]
    value <- if (length(pieces) > 1) pieces[[2]] else NULL
    if (identical(key, "--variant") && !is.null(value)) {
      options$variant <- value
    }
    if (identical(key, "--iterations") && !is.null(value)) {
      options$iterations <- as.integer(value)
    }
    if (identical(key, "--burnin") && !is.null(value)) {
      options$burnin <- as.integer(value)
    }
    if (identical(key, "--summary-samples") && !is.null(value)) {
      options$summary_samples <- as.integer(value)
    }
  }
  options
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  check_requirements()
  options <- parse_args(args)
  project_root <- get_project_root()
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  rrtm_ns <- load_rrtm_namespace(project_root)
  PROSPECTD_FUN <<- get_rrtm_object(rrtm_ns, "prospectd")
  lookup <- list(
    dataspec_pd = get_rrtm_object(rrtm_ns, "dataspec_pd"),
    pd_talf = get_rrtm_object(rrtm_ns, "pd_talf"),
    pd_t12 = get_rrtm_object(rrtm_ns, "pd_t12"),
    pd_t21 = get_rrtm_object(rrtm_ns, "pd_t21")
  )
  obs_data <- read_observation()

  variants <- if (identical(options$variant, "both")) {
    c("direct", "rewrite")
  } else {
    options$variant
  }

  results <- list()
  for (variant in variants) {
    tstart <- Sys.time()
    results[[variant]] <- run_variant(
      variant,
      obs_data$obs_prospect,
      lookup,
      options$iterations,
      options$burnin,
      options$summary_samples
    )
    tend <- Sys.time()
    message("Variant: ", variant)
    message("Timing: ", tend - tstart)
  }

  invisible(results)
}

if (sys.nframe() == 0) {
  main()
}
