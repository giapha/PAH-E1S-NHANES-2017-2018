#!/usr/bin/env Rscript

# Final NHANES epidemiology recalibration for the PAH-E1S manuscript.
#
# The script applies one inference policy to every journal-facing survey model:
#   1. Taylor-linearized NHANES design with the component-specific weight.
#   2. Domain analysis from all positive-weight records.
#   3. Denominator df = number of represented PSUs - number of strata.
#   4. Design-t confidence intervals and P values.
#   5. Benjamini-Hochberg correction within explicitly named families.
#   6. Parsimonious three-knot natural splines (10th, 50th, 90th percentiles)
#      with a design-based F test of the nonlinear component.

suppressMessages({
  library(data.table)
  library(haven)
  library(splines)
  library(survey)
})

set.seed(2026)
options(survey.lonely.psu = "adjust")

ROOT <- normalizePath(".", mustWork = TRUE)
TAB <- file.path(ROOT, "results")
AUDIT <- file.path(ROOT, "results/audit/design_based_inference")
RAW <- file.path(ROOT, "data/raw")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT, recursive = TRUE, showWarnings = FALSE)

canonical_files <- c(
  "PAH_discovery_E1S_forest_canonical.csv",
  "PAH_endocrine_axis_scan_canonical.csv",
  "PAH_steroid_specificity_scan_canonical.csv",
  "PAH_E1S_lipid_scan_canonical.csv",
  "PAH_hepatic_substrate_selectivity.csv",
  "PAH_replication_2013_2016_active_hormones.csv",
  "PAH_reviewer_credibility_sensitivity_20260711.csv"
)

df <- readRDS(file.path(ROOT, "data/processed/analysis_dataset.rds"))
df$race_eth <- factor(df$race_eth)
df$education <- factor(df$education)

base_covariates <- c(
  "ln_URXUCR", "RIDAGEYR", "BMXBMI", "INDFMPIR",
  "race_eth", "education", "ln_LBXCOT"
)
base_no_bmi <- setdiff(base_covariates, "BMXBMI")
design_vars <- c("SDMVPSU", "SDMVSTRA")

iqr_value <- function(x) {
  as.numeric(diff(quantile(x, c(0.25, 0.75), na.rm = TRUE)))
}

make_domain_design <- function(data, variables, weight, extra_domain = NULL) {
  data <- as.data.frame(data)
  positive <- !is.na(data[[weight]]) & data[[weight]] > 0
  full <- data[positive, , drop = FALSE]
  domain <- complete.cases(full[, unique(c(variables, weight, design_vars)), drop = FALSE])
  if (!is.null(extra_domain)) domain <- domain & extra_domain[positive]
  full$.analysis_domain <- domain
  all_design <- svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = as.formula(paste0("~", weight)),
    data = full,
    nest = TRUE
  )
  design <- subset(all_design, .analysis_domain)
  list(design = design, data = design$variables, df = degf(design))
}

fit_term <- function(data, outcome, exposure, covariates = base_covariates,
                     weight = "WTSSTS2Y", family = gaussian(),
                     extra_domain = NULL) {
  variables <- unique(c(outcome, exposure, covariates))
  obj <- make_domain_design(data, variables, weight, extra_domain)
  model <- svyglm(
    as.formula(paste(outcome, "~", exposure, "+", paste(covariates, collapse = " + "))),
    design = obj$design,
    family = family
  )
  beta <- unname(coef(model)[exposure])
  se <- unname(sqrt(diag(vcov(model)))[exposure])
  statistic <- beta / se
  critical <- qt(0.975, obj$df)
  list(
    n = nrow(obj$data),
    df = obj$df,
    beta = beta,
    se = se,
    statistic = statistic,
    p = 2 * pt(-abs(statistic), df = obj$df),
    ci_beta_low = beta - critical * se,
    ci_beta_high = beta + critical * se,
    iqr = iqr_value(obj$data[[exposure]]),
    design = obj$design,
    model = model,
    data = obj$data
  )
}

as_percent_iqr <- function(fit) {
  data.table(
    n = fit$n,
    df = fit$df,
    beta = fit$beta,
    se = fit$se,
    iqr = fit$iqr,
    pct = (exp(fit$beta * fit$iqr) - 1) * 100,
    lo = (exp(fit$ci_beta_low * fit$iqr) - 1) * 100,
    hi = (exp(fit$ci_beta_high * fit$iqr) - 1) * 100,
    t = fit$statistic,
    p = fit$p
  )
}

map <- data.table(
  code = c("ln_URXP01", "ln_URXP03", "ln_URXP04", "ln_URXP06", "ln_URXP25", "ln_URXP10"),
  analyte = c(
    "1-OH-naphthalene", "3-OH-fluorene", "2-OH-fluorene",
    "1-OH-phenanthrene", "2/3-OH-phenanthrene", "1-OH-pyrene"
  ),
  ring = c(2L, 3L, 3L, 3L, 3L, 4L)
)

# -------------------------------------------------------------------------
# Discovery family: six OH-PAH biomarkers versus E1S.
# -------------------------------------------------------------------------
discovery <- rbindlist(lapply(seq_len(nrow(map)), function(i) {
  fit <- fit_term(df, "ln_SSTES1", map$code[i])
  cbind(data.table(analyte = map$analyte[i], ring = map$ring[i]), as_percent_iqr(fit))
}))
discovery[, fdr := p.adjust(p, method = "BH")]
fwrite(discovery, file.path(TAB, "PAH_discovery_E1S_forest_canonical.csv"))

# -------------------------------------------------------------------------
# Endocrine selectivity family: two leads across 13 serum endpoints.
# -------------------------------------------------------------------------
df$ln_ratio_E1S_E1 <- df$ln_SSTES1 - df$ln_SSTESO
df$ln_E1_E2 <- df$ln_SSTESO - df$ln_SSTEST
df$ln_E1S_E2 <- df$ln_SSTES1 - df$ln_SSTEST

endpoints <- data.table(
  variable = c(
    "ln_SSTES1", "ln_SSTESO", "ln_ratio_E1S_E1", "ln_E1_E2", "ln_E1S_E2",
    "ln_SSTEST", "ln_SSTAND", "ln_SST17H", "ln_SSTPG4", "ln_SSTSHBG",
    "ln_SSTFSH", "ln_SSTLUH", "ln_SSTAMH"
  ),
  endpoint = c(
    "Estrone sulfate (E1S)", "Estrone (E1)", "E1S/E1 ratio (reservoir index)",
    "E1/E2 (parent : active estrogen)", "E1S/E2 (reservoir : active estrogen)",
    "Estradiol (E2, active)", "Androstenedione", "17-OH-progesterone",
    "Progesterone", "SHBG", "FSH", "LH", "AMH"
  ),
  axis = c(
    rep("estrogen reservoir", 5), "active estrogen", "androgen",
    "steroidogenesis", "steroidogenesis", "binding",
    "HPG gonadotropin", "HPG gonadotropin", "gonadal"
  )
)

lead_map <- c(ln_URXP10 = "1-OH-pyrene", ln_URXP25 = "2/3-OH-phenanthrene")
endocrine <- rbindlist(lapply(names(lead_map), function(exposure) {
  rbindlist(lapply(seq_len(nrow(endpoints)), function(i) {
    fit <- fit_term(df, endpoints$variable[i], exposure)
    cbind(
      data.table(
        pah = unname(lead_map[exposure]),
        endpoint = endpoints$endpoint[i],
        axis = endpoints$axis[i]
      ),
      as_percent_iqr(fit)
    )
  }))
}))
endocrine[, fdr_family := fifelse(
  grepl("ratio|E1/E2|E1S/E2", endpoint),
  "reservoir-ratio",
  "measured endocrine analytes"
)]
endocrine[, fdr := p.adjust(p, method = "BH"), by = .(pah, fdr_family)]
fwrite(endocrine, file.path(TAB, "PAH_endocrine_axis_scan_canonical.csv"))

steroid_rows <- data.table(
  endpoint = c(
    "Estrone sulfate (E1S)", "Estrone (E1)", "Estradiol (E2, active)",
    "Androstenedione", "17-OH-progesterone", "Progesterone", "SHBG", "AMH"
  ),
  steroid = c(
    "E1S (estrone SULFATE)", "E1 (estrone)", "Estradiol", "Androstenedione",
    "17-OH-progesterone", "Progesterone", "SHBG", "AMH"
  ),
  type = c("SULFATED", rep("free", 5), "binding-protein", "free")
)
steroid_specificity <- merge(
  endocrine,
  steroid_rows,
  by = "endpoint",
  all = FALSE,
  allow.cartesian = TRUE
)
steroid_specificity[, pah := fifelse(
  pah == "1-OH-pyrene", "ln_URXP10", "ln_URXP25"
)]
steroid_specificity <- steroid_specificity[, .(pah, steroid, type, n, df, pct, lo, hi, t, p, fdr)]
fwrite(steroid_specificity, file.path(TAB, "PAH_steroid_specificity_scan_canonical.csv"))

# -------------------------------------------------------------------------
# Lipid families: continuous E1S-lipid, clinical E1S-lipid, ratios, direct PAH.
# -------------------------------------------------------------------------
read_xpt_df <- function(name) as.data.frame(read_xpt(file.path(RAW, name)))
tc <- read_xpt_df("TCHOL_J.XPT")[, c("SEQN", "LBXTC")]
hdl <- read_xpt_df("HDL_J.XPT")[, c("SEQN", "LBDHDD")]
lipid_df <- df
if (!"LBXTC" %in% names(lipid_df)) {
  lipid_df <- merge(lipid_df, tc, by = "SEQN", all.x = TRUE)
}
if (!"LBDHDD" %in% names(lipid_df)) {
  lipid_df <- merge(lipid_df, hdl, by = "SEQN", all.x = TRUE)
}
lipid_df$nonHDL <- lipid_df$LBXTC - lipid_df$LBDHDD
lipid_df$high_TC <- ifelse(is.na(lipid_df$LBXTC), NA_integer_, as.integer(lipid_df$LBXTC >= 200))
lipid_df$high_nonHDL <- ifelse(is.na(lipid_df$nonHDL), NA_integer_, as.integer(lipid_df$nonHDL >= 130))

lipid_results <- list()
for (outcome in c("LBXTC", "nonHDL", "LBDHDD")) {
  fit <- fit_term(lipid_df, outcome, "ln_SSTES1")
  lipid_results[[length(lipid_results) + 1L]] <- data.table(
    block = "E1S->lipid",
    predictor = "E1S",
    outcome = c(LBXTC = "Total cholesterol", nonHDL = "Non-HDL-C", LBDHDD = "HDL-C")[[outcome]],
    unit = "mg/dL per E1S IQR",
    n = fit$n,
    df = fit$df,
    effect = fit$beta * fit$iqr,
    lo = fit$ci_beta_low * fit$iqr,
    hi = fit$ci_beta_high * fit$iqr,
    t = fit$statistic,
    p = fit$p
  )
}
for (outcome in c("high_TC", "high_nonHDL")) {
  fit <- fit_term(lipid_df, outcome, "ln_SSTES1", family = quasipoisson())
  lipid_results[[length(lipid_results) + 1L]] <- data.table(
    block = "E1S->clinical",
    predictor = "E1S",
    outcome = c(
      high_TC = "High total cholesterol (>=200)",
      high_nonHDL = "High non-HDL-C (>=130)"
    )[[outcome]],
    unit = "prevalence ratio per E1S IQR",
    n = fit$n,
    df = fit$df,
    effect = exp(fit$beta * fit$iqr),
    lo = exp(fit$ci_beta_low * fit$iqr),
    hi = exp(fit$ci_beta_high * fit$iqr),
    t = fit$statistic,
    p = fit$p
  )
}
for (outcome in c("LBXTC", "nonHDL")) {
  fit <- fit_term(lipid_df, outcome, "ln_ratio_E1S_E1")
  lipid_results[[length(lipid_results) + 1L]] <- data.table(
    block = "ratio->lipid",
    predictor = "E1S_E1",
    outcome = c(LBXTC = "Total cholesterol", nonHDL = "Non-HDL-C")[[outcome]],
    unit = "mg/dL per ratio IQR",
    n = fit$n,
    df = fit$df,
    effect = fit$beta * fit$iqr,
    lo = fit$ci_beta_low * fit$iqr,
    hi = fit$ci_beta_high * fit$iqr,
    t = fit$statistic,
    p = fit$p
  )
}
for (exposure in names(lead_map)) {
  for (outcome in c("LBXTC", "nonHDL")) {
    fit <- fit_term(lipid_df, outcome, exposure)
    lipid_results[[length(lipid_results) + 1L]] <- data.table(
      block = "PAH->lipid (boundary)",
      predictor = unname(lead_map[exposure]),
      outcome = c(LBXTC = "Total cholesterol", nonHDL = "Non-HDL-C")[[outcome]],
      unit = "mg/dL per PAH IQR",
      n = fit$n,
      df = fit$df,
      effect = fit$beta * fit$iqr,
      lo = fit$ci_beta_low * fit$iqr,
      hi = fit$ci_beta_high * fit$iqr,
      t = fit$statistic,
      p = fit$p
    )
  }
}
lipid <- rbindlist(lipid_results, fill = TRUE)
lipid[, fdr_family := fifelse(
  block == "E1S->lipid", "continuous E1S-lipid",
  fifelse(
    block == "E1S->clinical", "clinical E1S-lipid",
    fifelse(block == "ratio->lipid", "reservoir-ratio lipid", "direct PAH-lipid")
  )
)]
lipid[, fdr := p.adjust(p, method = "BH"), by = fdr_family]
fwrite(lipid, file.path(TAB, "PAH_E1S_lipid_scan_canonical.csv"))

# -------------------------------------------------------------------------
# Hepatic/liver comparator family under the environmental subsample weight.
# -------------------------------------------------------------------------
biopro <- read_xpt_df("BIOPRO_J.XPT")[, c(
  "SEQN", "LBXSTB", "LBXSGTSI", "LBXSATSI", "LBXSASSI", "LBXSAPSI"
)]
hepatic_df <- df
missing_biopro <- setdiff(names(biopro), names(hepatic_df))
if (length(missing_biopro) > 0) {
  hepatic_df <- merge(
    hepatic_df,
    biopro[, c("SEQN", missing_biopro), drop = FALSE],
    by = "SEQN",
    all.x = TRUE
  )
}
for (v in c("LBXSTB", "LBXSGTSI", "LBXSATSI", "LBXSASSI", "LBXSAPSI")) {
  hepatic_df[[paste0("ln_", v)]] <- log(hepatic_df[[v]])
}
hepatic_panel <- data.table(
  outcome = c(
    "ln_SSTES1", "ln_SSTESO", "ln_LBXSTB", "ln_LBXSGTSI",
    "ln_LBXSAPSI", "ln_LBXSATSI", "ln_LBXSASSI"
  ),
  label = c(
    "E1S (sulfate conjugate)", "Estrone E1 (unconjugated)",
    "Total bilirubin (OATP1B1 substrate)", "GGT", "Alkaline phosphatase", "ALT", "AST"
  ),
  cls = c(
    "sulfate conjugate", "unconjugated estrogen",
    rep("OATP1B1 substrate / liver", 5)
  )
)
hepatic <- rbindlist(lapply(seq_len(nrow(hepatic_panel)), function(i) {
  fit <- fit_term(
    hepatic_df,
    hepatic_panel$outcome[i],
    "ln_URXP10",
    weight = "WTSA2YR"
  )
  cbind(
    hepatic_panel[i, .(label, cls)],
    as_percent_iqr(fit)
  )
}))
hepatic[, fdr := p.adjust(p, method = "BH")]
fwrite(hepatic, file.path(TAB, "PAH_hepatic_substrate_selectivity.csv"))

# -------------------------------------------------------------------------
# Detection-limit and BMI sensitivity analyses.
# -------------------------------------------------------------------------
df$bmi_group <- factor(
  ifelse(df$BMXBMI < 25, "BMI <25", "BMI >=25"),
  levels = c("BMI <25", "BMI >=25")
)
df$detected_1OHPyr <- as.integer(df$URDP10LC == 0)
df$URXP10_lod2 <- ifelse(df$URDP10LC == 1, df$URXP10 / sqrt(2), df$URXP10)
df$ln_URXP10_lod2 <- log(df$URXP10_lod2)

sensitivity <- list()
add_continuous_sensitivity <- function(analysis, exposure, population,
                                       extra_domain = NULL, covariates = base_covariates,
                                       detection_percent = NA_real_, p_heterogeneity = NA_real_) {
  fit <- fit_term(
    df, "ln_SSTES1", exposure, covariates = covariates,
    extra_domain = extra_domain
  )
  row <- as_percent_iqr(fit)
  row[, `:=`(
    analysis = analysis,
    analyte = fifelse(grepl("URXP25", exposure), "2/3-hydroxyphenanthrene", "1-hydroxypyrene"),
    exposure = exposure,
    population = population,
    detection_percent = detection_percent,
    percent_change = pct,
    ci_low = lo,
    ci_high = hi,
    p_heterogeneity = p_heterogeneity
  )]
  row[, .(
    analysis, analyte, exposure, population, n, detection_percent,
    df, beta, se, iqr, percent_change, ci_low, ci_high, t, p, p_heterogeneity
  )]
}

sensitivity[[length(sensitivity) + 1L]] <- add_continuous_sensitivity(
  "Primary CDC substitution", "ln_URXP10", "Full analytic sample",
  detection_percent = mean(df$URDP10LC == 0, na.rm = TRUE) * 100
)
sensitivity[[length(sensitivity) + 1L]] <- add_continuous_sensitivity(
  "Alternative LOD/2 substitution", "ln_URXP10_lod2", "Full analytic sample",
  detection_percent = mean(df$URDP10LC == 0, na.rm = TRUE) * 100
)
sensitivity[[length(sensitivity) + 1L]] <- add_continuous_sensitivity(
  "Detectable-only continuous", "ln_URXP10", "1-OHPyr at or above the LOD",
  extra_domain = !is.na(df$URDP10LC) & df$URDP10LC == 0,
  detection_percent = 100
)

binary_fit <- fit_term(df, "ln_SSTES1", "detected_1OHPyr")
critical_binary <- qt(0.975, binary_fit$df)
sensitivity[[length(sensitivity) + 1L]] <- data.table(
  analysis = "Binary detectability",
  analyte = "1-hydroxypyrene",
  exposure = "detected_1OHPyr",
  population = "Full analytic sample; detected versus below LOD",
  n = binary_fit$n,
  detection_percent = mean(binary_fit$data$detected_1OHPyr == 1) * 100,
  df = binary_fit$df,
  beta = binary_fit$beta,
  se = binary_fit$se,
  iqr = NA_real_,
  percent_change = (exp(binary_fit$beta) - 1) * 100,
  ci_low = (exp(binary_fit$beta - critical_binary * binary_fit$se) - 1) * 100,
  ci_high = (exp(binary_fit$beta + critical_binary * binary_fit$se) - 1) * 100,
  t = binary_fit$statistic,
  p = binary_fit$p,
  p_heterogeneity = NA_real_
)

interaction_p <- function(exposure) {
  variables <- unique(c("ln_SSTES1", exposure, "bmi_group", base_no_bmi))
  obj <- make_domain_design(df, variables, "WTSSTS2Y")
  model <- svyglm(
    as.formula(paste(
      "ln_SSTES1 ~", exposure, "* bmi_group +",
      paste(base_no_bmi, collapse = " + ")
    )),
    design = obj$design
  )
  term <- grep(
    paste0(exposure, ":bmi_group|bmi_group.*:", exposure),
    names(coef(model)),
    value = TRUE
  )
  beta <- unname(coef(model)[term])
  se <- unname(sqrt(diag(vcov(model)))[term])
  2 * pt(-abs(beta / se), df = obj$df)
}

for (exposure in c("ln_URXP25", "ln_URXP10")) {
  p_heterogeneity <- interaction_p(exposure)
  for (level in levels(df$bmi_group)) {
    sensitivity[[length(sensitivity) + 1L]] <- add_continuous_sensitivity(
      "BMI-stratified association",
      exposure,
      level,
      extra_domain = !is.na(df$bmi_group) & df$bmi_group == level,
      covariates = base_no_bmi,
      detection_percent = if (exposure == "ln_URXP10") {
        mean(df$URDP10LC[df$bmi_group == level] == 0, na.rm = TRUE) * 100
      } else {
        mean(df$URDP25LC[df$bmi_group == level] == 0, na.rm = TRUE) * 100
      },
      p_heterogeneity = p_heterogeneity
    )
  }
}
sensitivity <- rbindlist(sensitivity, fill = TRUE)
fwrite(sensitivity, file.path(TAB, "PAH_reviewer_credibility_sensitivity_20260711.csv"))

# -------------------------------------------------------------------------
# Independent-cycle active-hormone comparison (2013-2016).
# -------------------------------------------------------------------------
rep_raw <- file.path(RAW, "replication_2013_2016")
read_rep <- function(stem) as.data.table(read_xpt(file.path(rep_raw, paste0(stem, ".XPT"))))
build_cycle <- function(suffix) {
  pah <- read_rep(paste0("PAH_", suffix))
  tst <- read_rep(paste0("TST_", suffix))
  dem <- read_rep(paste0("DEMO_", suffix))
  cr <- read_rep(paste0("ALB_CR_", suffix))
  bmx <- read_rep(paste0("BMX_", suffix))
  cot <- read_rep(paste0("COT_", suffix))
  pcols <- intersect(
    c("SEQN", "WTSA2YR", "URXP01", "URXP03", "URXP04", "URXP06", "URXP10", "URXP25"),
    names(pah)
  )
  merged <- Reduce(function(a, b) merge(a, b, by = "SEQN", all.x = TRUE), list(
    pah[, ..pcols],
    tst[, intersect(c("SEQN", "LBXTST", "LBXEST", "LBXSHBG"), names(tst)), with = FALSE],
    dem[, intersect(c(
      "SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "INDFMPIR",
      "SDMVPSU", "SDMVSTRA"
    ), names(dem)), with = FALSE],
    cr[, intersect(c("SEQN", "URXUCR"), names(cr)), with = FALSE],
    bmx[, intersect(c("SEQN", "BMXBMI"), names(bmx)), with = FALSE],
    cot[, intersect(c("SEQN", "LBXCOT"), names(cot)), with = FALSE]
  ))
  merged[, cycle := suffix]
  merged
}
replication_df <- rbindlist(list(build_cycle("H"), build_cycle("I")), fill = TRUE)
replication_df <- replication_df[
  RIAGENDR == 1 & RIDAGEYR >= 20 & !is.na(WTSA2YR) & WTSA2YR > 0
]
replication_df[, w2 := WTSA2YR / 2]
replication_df[, race_eth := factor(RIDRETH1)]
replication_df[, education := factor(fifelse(DMDEDUC2 %in% c(7, 9), NA_real_, DMDEDUC2))]
for (v in c("URXP10", "URXP25", "URXUCR", "LBXCOT", "LBXTST", "LBXEST", "LBXSHBG")) {
  replication_df[[paste0("ln_", v)]] <- log(replication_df[[v]])
}
rep_outcomes <- c(
  ln_LBXEST = "Estradiol (active E2)",
  ln_LBXTST = "Testosterone",
  ln_LBXSHBG = "SHBG"
)
replication <- rbindlist(lapply(names(lead_map), function(exposure) {
  rbindlist(lapply(names(rep_outcomes), function(outcome) {
    fit <- fit_term(
      replication_df,
      outcome,
      exposure,
      weight = "w2"
    )
    cbind(
      data.table(lead = unname(lead_map[exposure]), outcome = unname(rep_outcomes[outcome])),
      as_percent_iqr(fit)
    )
  }))
}))
replication[, fdr := p.adjust(p, method = "BH")]
fwrite(replication, file.path(TAB, "PAH_replication_2013_2016_active_hormones.csv"))

# -------------------------------------------------------------------------
# Three-knot natural-spline diagnostics and design-t prediction bands.
# -------------------------------------------------------------------------
spline_stats <- list()
spline_curves <- list()
for (i in seq_len(nrow(map))) {
  exposure <- map$code[i]
  obj <- make_domain_design(
    df,
    unique(c("ln_SSTES1", exposure, base_covariates)),
    "WTSSTS2Y"
  )
  x <- obj$data[[exposure]]
  boundary <- as.numeric(quantile(x, c(0.10, 0.90), na.rm = TRUE))
  knot <- as.numeric(quantile(x, 0.50, na.rm = TRUE))
  basis <- ns(x, knots = knot, Boundary.knots = boundary)
  linear_space <- cbind(1, x)
  projection <- qr.solve(linear_space, basis)
  nonlinear_residual <- basis - linear_space %*% projection
  sv <- svd(nonlinear_residual)
  nonlinear_vector <- sv$v[, 1, drop = FALSE]
  nonlinear_score <- as.numeric(nonlinear_residual %*% nonlinear_vector)
  obj$design$variables$.nonlinear <- nonlinear_score

  linear_model <- svyglm(
    as.formula(paste("ln_SSTES1 ~", exposure, "+", paste(base_covariates, collapse = " + "))),
    design = obj$design
  )
  spline_model <- svyglm(
    as.formula(paste(
      "ln_SSTES1 ~", exposure, "+ .nonlinear +",
      paste(base_covariates, collapse = " + ")
    )),
    design = obj$design
  )
  nonlinear_test <- regTermTest(
    spline_model,
    ~.nonlinear,
    df = obj$df,
    method = "Wald"
  )

  mode_factor <- function(v) names(sort(table(v), decreasing = TRUE))[1]
  grid <- seq(
    quantile(x, 0.05, na.rm = TRUE),
    quantile(x, 0.95, na.rm = TRUE),
    length.out = 100
  )
  grid_basis <- predict(basis, grid)
  grid_residual <- grid_basis - cbind(1, grid) %*% projection
  grid_nonlinear <- as.numeric(grid_residual %*% nonlinear_vector)
  newdata <- data.frame(
    ln_URXUCR = median(obj$data$ln_URXUCR, na.rm = TRUE),
    RIDAGEYR = median(obj$data$RIDAGEYR, na.rm = TRUE),
    BMXBMI = median(obj$data$BMXBMI, na.rm = TRUE),
    INDFMPIR = median(obj$data$INDFMPIR, na.rm = TRUE),
    ln_LBXCOT = median(obj$data$ln_LBXCOT, na.rm = TRUE),
    race_eth = factor(mode_factor(obj$data$race_eth), levels = levels(df$race_eth)),
    education = factor(mode_factor(obj$data$education), levels = levels(df$education))
  )
  newdata <- newdata[rep(1, length(grid)), , drop = FALSE]
  newdata[[exposure]] <- grid
  newdata$.nonlinear <- grid_nonlinear

  critical <- qt(0.975, obj$df)
  for (model_name in c("Linear", "Natural spline")) {
    model <- if (model_name == "Linear") linear_model else spline_model
    pred <- predict(model, newdata = newdata, se.fit = TRUE)
    pred_se <- sqrt(attr(pred, "var"))
    spline_curves[[length(spline_curves) + 1L]] <- data.table(
      short = map$analyte[i],
      code = exposure,
      x = grid,
      fit = as.numeric(pred),
      lo = as.numeric(pred) - critical * pred_se,
      hi = as.numeric(pred) + critical * pred_se,
      model = model_name,
      df = obj$df
    )
  }
  spline_stats[[length(spline_stats) + 1L]] <- data.table(
    short = map$analyte[i],
    code = exposure,
    ring = map$ring[i],
    n = nrow(obj$data),
    df1 = 1L,
    df2 = obj$df,
    F_nonlinear = as.numeric(nonlinear_test$Ftest),
    p_nonlinear = as.numeric(nonlinear_test$p),
    knot_p10 = boundary[1],
    knot_p50 = knot,
    knot_p90 = boundary[2],
    AIC_linear = AIC(linear_model)[["AIC"]],
    AIC_spline = AIC(spline_model)[["AIC"]]
  )
}
spline_stats <- rbindlist(spline_stats)
spline_stats[, delta_AIC := AIC_linear - AIC_spline]
spline_stats[, shape := fifelse(
  p_nonlinear < 0.05,
  "evidence of curvature",
  "no design-F evidence of curvature"
)]
spline_curves <- rbindlist(spline_curves)
fwrite(spline_stats, file.path(TAB, "PAH_dose_response_designF_20260711.csv"))
fwrite(spline_curves, file.path(TAB, "PAH_dose_response_designT_curves_20260711.csv"))

# -------------------------------------------------------------------------
# Machine-readable run summary.
# -------------------------------------------------------------------------
summary_table <- rbindlist(list(
  discovery[analyte %in% c("1-OH-pyrene", "2/3-OH-phenanthrene"), .(
    family = "discovery", item = analyte, n, df, effect = pct, lo, hi, p, fdr
  )],
  lipid[block %in% c("E1S->lipid", "E1S->clinical"), .(
    family = fdr_family, item = outcome, n, df, effect, lo, hi, p, fdr
  )]
), fill = TRUE)
fwrite(summary_table, file.path(AUDIT, "journal_facing_recalibrated_headlines.csv"))

cat("\n== Discovery (design t) ==\n")
print(discovery[, .(
  analyte, n, df, pct = round(pct, 1), lo = round(lo, 1), hi = round(hi, 1),
  p = signif(p, 4), fdr = signif(fdr, 4)
)])
cat("\n== Dose-response shape (three-knot natural spline; design F) ==\n")
print(spline_stats[, .(
  short, n, df2, F_nonlinear = round(F_nonlinear, 3),
  p_nonlinear = signif(p_nonlinear, 4), delta_AIC = round(delta_AIC, 2), shape
)])
cat("\n== Lipid context (family-specific FDR) ==\n")
print(lipid[block %in% c("E1S->lipid", "E1S->clinical"), .(
  block, outcome, n, df, effect = round(effect, 3), lo = round(lo, 3),
  hi = round(hi, 3), p = signif(p, 4), fdr = signif(fdr, 4)
)])
cat("\nWROTE recalibrated canonical tables to:", TAB, "\n")
cat("AUDIT backup and summaries:", AUDIT, "\n")
