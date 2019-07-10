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

dataspec_pd <- dataspec_pd_full[, !grep("refractive", colnames(dataspec_pd_full))]
refractive_pd <- dataspec_pd_full[, "refractive"]

usethis::use_data(
  dataspec_p4,
  dataspec_p5,
  dataspec_pd,
  refractive_p45,
  refractive_pd
)
