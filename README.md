# PAH-E1S NHANES 2017-2018 Analysis Code

This repository contains reproducibility code for the manuscript:

`Urinary phenanthrene/pyrene PAH biomarkers are associated with an estrone sulfate-centered steroid-reservoir phenotype in U.S. adult men`

Repository URL: <https://github.com/giapha/PAH-E1S-NHANES-2017-2018>

The analysis uses public NHANES 2017-2018 files to evaluate urinary monohydroxylated polycyclic aromatic hydrocarbon biomarkers, serum estrone sulfate, active estrogen comparators, lipid context, and survey design variables in adult men.

## Data

No participant-level data are stored in this repository. The scripts download or read public NHANES 2017-2018 XPT files listed in `docs/data_manifest.csv`.

Place raw XPT files in `data/raw/`, or run:

```r
source("R/01_download_nhanes_files.R")
```

## Reproduce

From the repository root:

```r
source("run_all.R")
```

The scripts write generated aggregate outputs to `results/` and intermediate RDS files to `data/processed/`. Aggregate result tables are included in this repository for transparent review; raw NHANES XPT files and derived participant-level RDS files are intentionally not tracked.

## Main Scripts

| Script | Purpose |
|---|---|
| `R/00_setup.R` | Shared package loading, file paths, variables, and helper functions. |
| `R/01_download_nhanes_files.R` | Download public NHANES 2017-2018 source files. |
| `R/02_build_analysis_dataset.R` | Merge source files and derive analysis variables. |
| `R/03_primary_steroid_models.R` | PAH-steroid screen, retained PAH-E1S models, and reservoir-ratio analyses. |
| `R/04_lipid_context_models.R` | E1S-lipid and same-sample adjacent association checks. |
| `R/05_mixture_sensitivity.R` | Optional mixture and robustness sensitivity summaries. |
| `R/06_session_info.R` | Save R session information. |

## Interpretation Boundary

This code estimates cross-sectional associations. It does not establish temporal sequence, causal mediation, tissue-level mechanism, receptor activation, enzyme flux, or disease causality.

## License

Code is released under the MIT License. NHANES data remain governed by NCHS/CDC public data terms.
