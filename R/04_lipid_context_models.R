source("R/00_setup.R")

dat <- readRDS(file.path(processed_dir, "analysis_dataset.rds"))

lipid_covariates <- c("RIDAGEYR", "race_eth", "education", "INDFMPIR", "BMXBMI", "ln_LBXCOT")

fit_lipid_linear <- function(outcome, predictor = "ln_SSTES1") {
  variables <- c(outcome, predictor, lipid_covariates, "SDMVPSU", "SDMVSTRA", "WTSSTS2Y")
  model_data <- dat %>%
    filter(if_all(all_of(variables), ~ !is.na(.x))) %>%
    filter(WTSSTS2Y > 0)
  design <- nhanes_design(model_data, "WTSSTS2Y")
  form <- stats::as.formula(paste(outcome, "~", paste(c(predictor, lipid_covariates), collapse = " + ")))
  fit <- survey::svyglm(form, design = design)
  est <- svy_term_row(fit, predictor, design)
  scale <- safe_iqr(model_data[[predictor]])
  z <- stats::qt(0.975, df = est$df)
  tibble(
    outcome = outcome,
    predictor = predictor,
    n = nrow(model_data),
    beta = est$estimate,
    se = est$std.error,
    p_value = est$p.value,
    scale = scale,
    change_per_iqr = est$estimate * scale,
    ci_low = (est$estimate - z * est$std.error) * scale,
    ci_high = (est$estimate + z * est$std.error) * scale,
    model = "survey_weighted_linear"
  )
}

fit_lipid_binary <- function(outcome, predictor = "ln_SSTES1") {
  variables <- c(outcome, predictor, lipid_covariates, "SDMVPSU", "SDMVSTRA", "WTSSTS2Y")
  model_data <- dat %>%
    filter(if_all(all_of(variables), ~ !is.na(.x))) %>%
    filter(WTSSTS2Y > 0)
  design <- nhanes_design(model_data, "WTSSTS2Y")
  form <- stats::as.formula(paste(outcome, "~", paste(c(predictor, lipid_covariates), collapse = " + ")))
  fit <- survey::svyglm(form, design = design, family = quasipoisson(link = "log"))
  est <- svy_term_row(fit, predictor, design)
  scale <- safe_iqr(model_data[[predictor]])
  z <- stats::qt(0.975, df = est$df)
  tibble(
    outcome = outcome,
    predictor = predictor,
    n = nrow(model_data),
    beta = est$estimate,
    se = est$std.error,
    p_value = est$p.value,
    scale = scale,
    prevalence_ratio = exp(est$estimate * scale),
    ci_low = exp((est$estimate - z * est$std.error) * scale),
    ci_high = exp((est$estimate + z * est$std.error) * scale),
    model = "survey_weighted_modified_poisson"
  )
}

continuous_results <- purrr::map_dfr(
  c("LBXTC", "nonHDL_C", "LBDHDD"),
  fit_lipid_linear
) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

binary_results <- purrr::map_dfr(
  c("high_total_cholesterol", "high_nonHDL_C"),
  fit_lipid_binary
) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

readr::write_csv(continuous_results, file.path(results_dir, "e1s_lipid_continuous_models.csv"))
readr::write_csv(binary_results, file.path(results_dir, "e1s_lipid_binary_models.csv"))

same_sample_checks <- purrr::map_dfr(c("ln_URXP25", "ln_URXP10", "ln_URXP06"), function(exposure) {
  variables <- c("ln_SSTES1", "LBXTC", exposure, primary_covariates, "SDMVPSU", "SDMVSTRA", "WTSA2YR")
  model_data <- dat %>%
    filter(if_all(all_of(variables), ~ !is.na(.x))) %>%
    filter(WTSA2YR > 0)
  design <- nhanes_design(model_data, "WTSA2YR")
  pah_e1s <- survey::svyglm(
    stats::as.formula(paste("ln_SSTES1 ~", paste(c(exposure, primary_covariates), collapse = " + "))),
    design = design
  )
  pah_lipid <- survey::svyglm(
    stats::as.formula(paste("LBXTC ~", paste(c(exposure, primary_covariates), collapse = " + "))),
    design = design
  )
  bind_rows(
    svy_term_row(pah_e1s, exposure, design) %>% mutate(check = "PAH_to_E1S"),
    svy_term_row(pah_lipid, exposure, design) %>% mutate(check = "PAH_to_total_cholesterol")
  ) %>%
    transmute(
      exposure = exposure,
      check = check,
      n = nrow(model_data),
      beta = estimate,
      se = std.error,
      p_value = p.value
    )
})

readr::write_csv(same_sample_checks, file.path(results_dir, "same_sample_adjacent_checks.csv"))
