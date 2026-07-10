set.seed(2026)
options(survey.lonely.psu = "adjust")
options(na.action = "na.omit")

required_packages <- c(
  "haven", "survey", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "tibble", "broom", "data.table", "splines"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install required R packages before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
})

dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

root_dir <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(root_dir, "data", "raw")
processed_dir <- file.path(root_dir, "data", "processed")
results_dir <- file.path(root_dir, "results")
docs_dir <- file.path(root_dir, "docs")

dir_create(raw_dir)
dir_create(processed_dir)
dir_create(results_dir)

safe_log <- function(x) {
  ifelse(is.na(x) | x <= 0, NA_real_, log(x))
}

safe_iqr <- function(x) {
  as.numeric(stats::IQR(x, na.rm = TRUE))
}

pct_effect <- function(beta, scale = 1) {
  (exp(beta * scale) - 1) * 100
}

pct_ci <- function(beta, se, scale = 1, z = 1.96) {
  tibble(
    pct_low = pct_effect(beta - z * se, scale),
    pct_high = pct_effect(beta + z * se, scale)
  )
}

svy_term_row <- function(fit, term, design) {
  coefs <- coef(summary(fit))
  if (!term %in% rownames(coefs)) {
    stop("Model term not found: ", term, call. = FALSE)
  }
  beta <- unname(coefs[term, "Estimate"])
  se <- unname(coefs[term, "Std. Error"])
  df <- survey::degf(design)
  t_value <- beta / se
  p_value <- 2 * stats::pt(abs(t_value), df = df, lower.tail = FALSE)
  tibble(term = term, estimate = beta, std.error = se, statistic = t_value, p.value = p_value, df = df)
}

nhanes_design <- function(dat, weight = "WTSA2YR") {
  survey::svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = stats::as.formula(paste0("~", weight)),
    nest = TRUE,
    data = dat
  )
}

pah_dictionary <- tibble::tribble(
  ~var, ~lod_var, ~abbrev, ~parent_pah, ~ring_count, ~role,
  "URXP01", "URDP01LC", "1-OHNa", "Naphthalene", 2L, "context",
  "URXP03", "URDP03LC", "3-OHFlu", "Fluorene", 3L, "context",
  "URXP04", "URDP04LC", "2-OHFlu", "Fluorene", 3L, "context",
  "URXP06", "URDP06LC", "1-OHPhe", "Phenanthrene", 3L, "core",
  "URXP10", "URDP10LC", "1-OHPyr", "Pyrene", 4L, "core",
  "URXP25", "URDP25LC", "2/3-OHPhe", "Phenanthrene", 3L, "core"
)

steroid_outcomes <- tibble::tribble(
  ~outcome, ~label, ~role,
  "ln_SSTES1", "Estrone sulfate", "primary",
  "ln_SSTESO", "Estrone", "active_estrogen_comparator",
  "ln_SSTEST", "Estradiol", "active_estrogen_comparator"
)

primary_covariates <- c(
  "RIDAGEYR", "race_eth", "education", "INDFMPIR", "BMXBMI",
  "ln_LBXCOT", "ln_URXUCR"
)
