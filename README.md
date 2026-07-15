# PAH-E1S NHANES 2017-2018 Analysis Code

This repository contains reproducibility code for the manuscript:

`Urinary phenanthrene and pyrene biomarkers are associated with selective serum estrone sulfate accumulation in adult men`

Repository URL: <https://github.com/giapha/PAH-E1S-NHANES-2017-2018>

The analysis uses public NHANES 2017-2018 files to evaluate urinary monohydroxylated polycyclic aromatic hydrocarbon biomarkers, serum estrone sulfate, active estrogen comparators, lipid context, and survey design variables in adult men. Public NHANES 2013-2016 files are used only for the independent active-hormone comparison because estrone sulfate was unavailable in those cycles.

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
| `R/03_journal_facing_models.R` | Design-based discovery, endocrine-selectivity, lipid, hepatic-comparator, sensitivity, replication, and spline analyses used in the manuscript. |
| `R/04_verify_manuscript_values.R` | Fail-fast checks for the submitted headline values and survey design degrees of freedom. |
| `R/06_session_info.R` | Save R session information. |

## Inference Policy

All journal-facing models use Taylor-linearized NHANES survey designs with the component-specific examination weight. Analytic domains are defined within all positive-weight records. Confidence intervals and two-sided P values use the survey design degrees of freedom (represented primary sampling units minus strata), and Benjamini-Hochberg correction is applied within the prespecified outcome families documented in the script. Dose-response shape is assessed with three-knot natural splines and a design-based F test of the nonlinear component.

Unweighted BKMR, WQS, and qgcomp analyses are not part of the submitted manuscript or the reproducibility route.

## Locked Headline Check

The July 2026 clean run reproduces the manuscript estimates:

- 1-hydroxypyrene: +28.3% E1S per IQR (95% CI 12.3 to 46.6; FDR = 0.0072)
- 2/3-hydroxyphenanthrene: +22.8% E1S per IQR (95% CI 7.8 to 39.8; FDR = 0.0128)
- pyrene-associated estrone: +5.3%, nominal only (endocrine-family FDR = 0.055)

## Interpretation Boundary

This code estimates cross-sectional associations. It does not establish temporal sequence, causal mediation, tissue-level mechanism, receptor activation, enzyme flux, or disease causality.

## License

Code is released under the MIT License. NHANES data remain governed by NCHS/CDC public data terms.
