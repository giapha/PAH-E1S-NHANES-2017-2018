source("R/00_setup.R")

discovery <- readr::read_csv(
  file.path(results_dir, "PAH_discovery_E1S_forest_canonical.csv"),
  show_col_types = FALSE
)
lipid <- readr::read_csv(
  file.path(results_dir, "PAH_E1S_lipid_scan_canonical.csv"),
  show_col_types = FALSE
)

expect_close <- function(actual, expected, tolerance, label) {
  if (length(actual) != 1 || is.na(actual) || abs(actual - expected) > tolerance) {
    stop(
      label, " did not reproduce: observed ", actual,
      ", expected ", expected, " +/- ", tolerance,
      call. = FALSE
    )
  }
}

pyrene <- discovery %>% filter(analyte == "1-OH-pyrene")
phenanthrene <- discovery %>% filter(analyte == "2/3-OH-phenanthrene")
total_cholesterol <- lipid %>%
  filter(block == "E1S->lipid", outcome == "Total cholesterol")

expect_close(pyrene$n, 524, 0, "1-OH-pyrene N")
expect_close(pyrene$df, 15, 0, "1-OH-pyrene design df")
expect_close(pyrene$pct, 28.3, 0.05, "1-OH-pyrene percent change")
expect_close(pyrene$lo, 12.3, 0.05, "1-OH-pyrene lower CI")
expect_close(pyrene$hi, 46.6, 0.05, "1-OH-pyrene upper CI")
expect_close(pyrene$fdr, 0.00722, 0.000005, "1-OH-pyrene FDR")

expect_close(phenanthrene$pct, 22.8, 0.05, "2/3-OH-phenanthrene percent change")
expect_close(phenanthrene$lo, 7.8, 0.05, "2/3-OH-phenanthrene lower CI")
expect_close(phenanthrene$hi, 39.8, 0.05, "2/3-OH-phenanthrene upper CI")
expect_close(phenanthrene$fdr, 0.0128, 0.00005, "2/3-OH-phenanthrene FDR")

expect_close(total_cholesterol$effect, 9.4, 0.05, "E1S-total cholesterol effect")
expect_close(total_cholesterol$fdr, 0.00264, 0.000005, "E1S-total cholesterol FDR")

message("Manuscript-value verification passed.")
