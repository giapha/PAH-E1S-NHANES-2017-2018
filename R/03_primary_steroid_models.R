source("R/00_setup.R")

dat <- readRDS(file.path(processed_dir, "analysis_dataset.rds"))

fit_svy_linear <- function(data, outcome, exposure, covariates, weight = "WTSA2YR", scale_iqr = TRUE) {
  variables <- c(outcome, exposure, covariates, "SDMVPSU", "SDMVSTRA", weight)
  model_data <- data %>%
    filter(if_all(all_of(variables), ~ !is.na(.x))) %>%
    filter(.data[[weight]] > 0)
  if (nrow(model_data) < 50) {
    return(tibble(outcome = outcome, exposure = exposure, n = nrow(model_data), status = "too_few_observations"))
  }
  design <- nhanes_design(model_data, weight)
  form <- stats::as.formula(paste(outcome, "~", paste(c(exposure, covariates), collapse = " + ")))
  fit <- survey::svyglm(form, design = design)
  est <- svy_term_row(fit, exposure, design)
  scale <- if (scale_iqr) safe_iqr(model_data[[exposure]]) else 1
  ci <- pct_ci(est$estimate, est$std.error, scale, z = stats::qt(0.975, df = est$df))
  tibble(
    outcome = outcome,
    exposure = exposure,
    n = nrow(model_data),
    beta = est$estimate,
    se = est$std.error,
    p_value = est$p.value,
    scale = scale,
    percent_change = pct_effect(est$estimate, scale),
    percent_low = ci$pct_low,
    percent_high = ci$pct_high,
    status = "modeled"
  )
}

exposure_terms <- paste0("ln_", pah_dictionary$var)

screen_results <- tidyr::crossing(
  outcome = steroid_outcomes$outcome,
  exposure = exposure_terms
) %>%
  mutate(result = purrr::map2(outcome, exposure, ~ fit_svy_linear(dat, .x, .y, primary_covariates))) %>%
  select(result) %>%
  tidyr::unnest(result) %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(steroid_outcomes, by = "outcome") %>%
  left_join(pah_dictionary %>% mutate(exposure = paste0("ln_", var)), by = "exposure")

readr::write_csv(screen_results, file.path(results_dir, "primary_pah_steroid_screen.csv"))

reservoir_outcomes <- c("ln_E1S_E1", "ln_E1S_E2", "ln_E1S_active_pool")
core_terms <- pah_dictionary %>%
  filter(role == "core") %>%
  mutate(exposure = paste0("ln_", var)) %>%
  pull(exposure)

reservoir_results <- tidyr::crossing(
  outcome = reservoir_outcomes,
  exposure = core_terms
) %>%
  mutate(result = purrr::map2(outcome, exposure, ~ fit_svy_linear(dat, .x, .y, primary_covariates))) %>%
  select(result) %>%
  tidyr::unnest(result) %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(pah_dictionary %>% mutate(exposure = paste0("ln_", var)), by = "exposure")

readr::write_csv(reservoir_results, file.path(results_dir, "reservoir_ratio_models.csv"))

weight_sensitivity <- purrr::map_dfr(core_terms, function(exposure) {
  fit_svy_linear(dat, "ln_SSTES1", exposure, primary_covariates, weight = "WTSSTS2Y")
}) %>%
  left_join(pah_dictionary %>% mutate(exposure = paste0("ln_", var)), by = "exposure")

readr::write_csv(weight_sensitivity, file.path(results_dir, "alternative_weight_e1s_models.csv"))
