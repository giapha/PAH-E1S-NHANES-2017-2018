# Reproducibility Check

Date: 2026-07-11

## Scope

This repository contains clean reproducibility code for the NHANES PAH, steroid-reservoir, and lipid-context analyses used in the Environmental Research submission package.

It includes:

- public NHANES source-file manifest;
- data download script;
- analysis dataset construction script;
- survey-weighted PAH-steroid models with design-t inference;
- reservoir-ratio models;
- E1S-lipid context models;
- hepatic-substrate comparator models;
- LOD and BMI sensitivity analyses;
- independent-cycle active-hormone comparisons;
- survey-weighted spline diagnostics with design-based nonlinearity tests;
- aggregate generated result tables.

It intentionally excludes:

- raw NHANES XPT files;
- derived participant-level RDS files;
- exploratory manuscript drafting scripts;
- internal audit notes;
- private local paths.

## Verification

The repository was run from a clean repository root with:

```r
source("run_all.R")
```

All scripts were run from a clean repository root. The generated primary outputs reproduce the submitted Table 2 and Table 3 rounded values, including N = 524, 15 design degrees of freedom, +28.3% (95% CI 12.3, 46.6) for 1-hydroxypyrene, and +22.8% (95% CI 7.8, 39.8) for 2/3-hydroxyphenanthrene.

The submitted analysis does not use unweighted mixture models.

## Environment

The run environment is recorded in `results/sessionInfo.txt`.
