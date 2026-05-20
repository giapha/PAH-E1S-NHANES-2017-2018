source("R/00_setup.R")

dat <- readRDS(file.path(processed_dir, "analysis_dataset.rds"))

exposure_terms <- paste0("ln_", pah_dictionary$var)

complete_exposure_data <- dat %>%
  select(all_of(exposure_terms)) %>%
  filter(if_all(everything(), ~ !is.na(.x)))

spearman_correlation <- stats::cor(
  complete_exposure_data,
  method = "spearman",
  use = "pairwise.complete.obs"
)

readr::write_csv(
  as.data.frame(spearman_correlation) %>% tibble::rownames_to_column("exposure"),
  file.path(results_dir, "pah_spearman_correlation.csv")
)

mixture_note <- tibble(
  item = c("mixture_scope", "interpretation"),
  value = c(
    "The repository provides exposure-correlation and single-component reproducibility code. Full BKMR/WQS/qgcomp objects can be regenerated from the documented source files and analysis definitions.",
    "Mixture diagnostics are interpreted as exposure-structure support, not as causal source apportionment."
  )
)

readr::write_csv(mixture_note, file.path(results_dir, "mixture_sensitivity_note.csv"))
