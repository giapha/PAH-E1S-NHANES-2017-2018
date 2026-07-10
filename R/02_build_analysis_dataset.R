source("R/00_setup.R")

read_xpt_local <- function(file_name) {
  path <- file.path(raw_dir, file_name)
  if (!file.exists(path)) {
    stop("Missing raw NHANES file: ", file_name, call. = FALSE)
  }
  haven::read_xpt(path) %>%
    mutate(SEQN = as.integer(SEQN))
}

demo <- read_xpt_local("DEMO_J.XPT")
pah <- read_xpt_local("PAH_J.XPT")
sstst <- read_xpt_local("SSTST_J.XPT")
bmx <- read_xpt_local("BMX_J.XPT")
albcr <- read_xpt_local("ALB_CR_J.XPT")
cot <- read_xpt_local("COT_J.XPT")
tchol <- read_xpt_local("TCHOL_J.XPT")
hdl <- read_xpt_local("HDL_J.XPT")
trigly <- read_xpt_local("TRIGLY_J.XPT")
biopro <- read_xpt_local("BIOPRO_J.XPT")
bpxo <- read_xpt_local("BPXO_J.XPT")
bpq <- read_xpt_local("BPQ_J.XPT")

source_inventory <- tibble(
  file_name = c(
    "DEMO_J.XPT", "PAH_J.XPT", "SSTST_J.XPT", "BMX_J.XPT",
    "ALB_CR_J.XPT", "COT_J.XPT", "TCHOL_J.XPT", "HDL_J.XPT",
    "TRIGLY_J.XPT", "BIOPRO_J.XPT", "BPXO_J.XPT", "BPQ_J.XPT"
  ),
  n_rows = c(
    nrow(demo), nrow(pah), nrow(sstst), nrow(bmx), nrow(albcr), nrow(cot),
    nrow(tchol), nrow(hdl), nrow(trigly), nrow(biopro), nrow(bpxo), nrow(bpq)
  ),
  n_cols = c(
    ncol(demo), ncol(pah), ncol(sstst), ncol(bmx), ncol(albcr), ncol(cot),
    ncol(tchol), ncol(hdl), ncol(trigly), ncol(biopro), ncol(bpxo), ncol(bpq)
  )
)
readr::write_csv(source_inventory, file.path(results_dir, "source_inventory.csv"))

bpxo_derived <- bpxo %>%
  mutate(
    SBP_mean = rowMeans(across(c(BPXOSY1, BPXOSY2, BPXOSY3)), na.rm = TRUE),
    DBP_mean = rowMeans(across(c(BPXODI1, BPXODI2, BPXODI3)), na.rm = TRUE),
    SBP_mean = if_else(is.nan(SBP_mean), NA_real_, SBP_mean),
    DBP_mean = if_else(is.nan(DBP_mean), NA_real_, DBP_mean)
  ) %>%
  select(SEQN, SBP_mean, DBP_mean)

dat <- pah %>%
  inner_join(sstst, by = "SEQN") %>%
  inner_join(demo, by = "SEQN") %>%
  left_join(bmx %>% select(SEQN, BMXBMI, BMXWAIST), by = "SEQN") %>%
  left_join(albcr %>% select(SEQN, URXUCR, URDACT), by = "SEQN") %>%
  left_join(cot %>% select(SEQN, LBXCOT), by = "SEQN") %>%
  left_join(tchol %>% select(SEQN, LBXTC), by = "SEQN") %>%
  left_join(hdl %>% select(SEQN, LBDHDD), by = "SEQN") %>%
  left_join(trigly %>% select(SEQN, WTSAF2YR, LBXTR, LBDLDL, LBDLDLM, LBDLDLN), by = "SEQN") %>%
  left_join(biopro %>% select(SEQN, LBXSATSI, LBXSASSI, LBXSGTSI, LBXSCR, LBXSBU), by = "SEQN") %>%
  left_join(bpxo_derived, by = "SEQN") %>%
  left_join(bpq %>% select(SEQN, BPQ100D), by = "SEQN") %>%
  mutate(
    sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
    race_eth = factor(
      RIDRETH3,
      levels = c(3, 1, 2, 4, 6, 7),
      labels = c("NH White", "Mexican American", "Other Hispanic", "NH Black", "NH Asian", "Other/Multi")
    ),
    education = factor(
      DMDEDUC2,
      levels = c(5, 4, 3, 2, 1),
      labels = c("College+", "Some college", "HS/GED", "9-11th", "<9th grade")
    ),
    adult_male = RIAGENDR == 1 & RIDAGEYR >= 20,
    ln_URXUCR = safe_log(URXUCR),
    ln_LBXCOT = safe_log(LBXCOT + 0.01),
    ln_SSTES1 = safe_log(SSTES1),
    ln_SSTESO = safe_log(SSTESO),
    ln_SSTEST = safe_log(SSTEST),
    ln_SSTAND = safe_log(SSTAND),
    ln_SST17H = safe_log(SST17H),
    ln_SSTPG4 = safe_log(SSTPG4),
    ln_SSTSHBG = safe_log(SSTSHBG),
    ln_SSTFSH = safe_log(SSTFSH),
    ln_SSTLUH = safe_log(SSTLUH),
    ln_SSTAMH = safe_log(SSTAMH),
    ln_E1S_E1 = safe_log(SSTES1SI / SSTESOSI),
    ln_E1S_E2 = safe_log(SSTES1SI / SSTESTSI),
    ln_E1S_active_pool = safe_log(SSTES1SI / (SSTESOSI + SSTESTSI)),
    nonHDL_C = LBXTC - LBDHDD,
    high_total_cholesterol = as.integer(LBXTC >= 200),
    high_nonHDL_C = as.integer(nonHDL_C >= 130),
    cholesterol_med = if_else(BPQ100D == 1, 1L, 0L, missing = NA_integer_)
  )

for (exposure in pah_dictionary$var) {
  dat[[paste0("ln_", exposure)]] <- safe_log(dat[[exposure]])
}

dat <- dat %>%
  filter(adult_male, WTSA2YR > 0)

cohort_variables <- c(
  paste0("ln_", pah_dictionary$var), "ln_SSTEST", primary_covariates,
  "SDMVPSU", "SDMVSTRA", "WTSSTS2Y"
)
primary_cohort <- dat %>%
  filter(if_all(all_of(cohort_variables), ~ !is.na(.x))) %>%
  filter(WTSSTS2Y > 0)

analysis_flow <- tibble(
  step = c(
    "PAH_J x SSTST_J",
    "PAH_J x SSTST_J x DEMO_J",
    "Adult male with positive WTSA2YR",
    "Prespecified six-biomarker cohort with active-estrogen comparator"
  ),
  n = c(
    nrow(pah %>% inner_join(sstst, by = "SEQN")),
    nrow(pah %>% inner_join(sstst, by = "SEQN") %>% inner_join(demo, by = "SEQN")),
    nrow(dat),
    nrow(primary_cohort)
  )
)

readr::write_csv(analysis_flow, file.path(results_dir, "analysis_flow.csv"))
saveRDS(dat, file.path(processed_dir, "analysis_dataset_full_adult_males.rds"))
saveRDS(primary_cohort, file.path(processed_dir, "analysis_dataset.rds"))
readr::write_csv(pah_dictionary, file.path(results_dir, "pah_dictionary.csv"))
readr::write_csv(steroid_outcomes, file.path(results_dir, "steroid_outcomes.csv"))
