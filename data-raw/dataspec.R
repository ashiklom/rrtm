devtools::load_all()

read_raw_data <- function(file) {
  raw_lines <- readLines(file, warn = FALSE)
  vals <- strsplit(raw_lines, ",|[[:space:]]")
  varnames <- lapply(vals, "[[", 1)
  values_all <- lapply(vals, function(x) as.numeric(trimws(c(x[-1]))))
  values <- lapply(values_all, na.omit)
  out <- matrix(0, 2101, length(values))
  colnames(out) <- varnames
  for (i in seq_along(values)) {
    v <- values[[i]]
    nv <- length(v)
    out[seq(1, nv), i] <- v
  }
  out
}

refractive_p45_raw <- readLines(file.path("data-raw", "refractive_4.dat"))
refractive_p45 <- as.numeric(trimws(strsplit(refractive_p45_raw, ",")[[1]]))

dataspec_p4 <- read_raw_data(file.path("data-raw", "prospect4_raw.dat"))
dataspec_p5 <- read_raw_data(file.path("data-raw", "prospect5_raw.dat"))
dataspec_pd_full <- read_raw_data(file.path("data-raw", "prospectd_raw.dat"))

dataspec_pd <- dataspec_pd_full[, !grepl("refractive", colnames(dataspec_pd_full))]
refractive_pd <- dataspec_pd_full[, "refractive"]

# Pre-calculate tav -- this is a relatively expensive operation
p45_talf <- tav_abs(40, refractive_p45)
p45_t12 <- tav_abs(90, refractive_p45)
p45_t21 <- p45_t12 / (refractive_p45 ^ 2)
pd_talf <- tav_abs(40, refractive_pd)
pd_t12 <- tav_abs(90, refractive_pd)
pd_t21 <- p45_t12 / (refractive_pd ^ 2)

# Read wood reflectance
wood_spec_file <- file.path("extdata", "wood-reflect.dat")
wood_spec <- scan(wood_spec_file)
wood_spec <- c(wood_spec, tail(wood_spec, 1))

# Read soil reflectance for Hapke model
soil_refl_file <- file.path("extdata", "hapke-soil.dat")
soil_refl <- read.table(soil_refl_file, header = FALSE)
dry_soil <- soil_refl[, 2]
wet_soil <- soil_refl[, 3]

# Read liberty input data
dataspec_liberty <- read.table("data-raw/liberty_raw.dat")

usethis::use_data(
  dataspec_p4,
  dataspec_p5,
  dataspec_pd,
  refractive_p45,
  refractive_pd,
  p45_talf, p45_t12, p45_t21,
  pd_talf, pd_t12, pd_t21,
  wood_spec,
  dry_soil,
  wet_soil,
  dataspec_liberty,
  overwrite = TRUE,
  internal = TRUE
)
