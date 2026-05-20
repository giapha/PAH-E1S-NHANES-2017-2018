# Reproducibility Check

Date: 2026-05-20

## Scope

This repository contains clean reproducibility code for the NHANES 2017-2018 PAH, steroid-reservoir, and lipid-context analyses used in the EHP submission package.

It includes:

- public NHANES source-file manifest;
- data download script;
- analysis dataset construction script;
- survey-weighted PAH-steroid models;
- reservoir-ratio models;
- E1S-lipid context models;
- same-sample adjacent association checks;
- exposure-correlation sensitivity summary;
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

All R scripts parsed successfully. The generated primary model outputs reproduce the EHP package Table 2 and Table 3 rounded values for the retained PAH-E1S and reservoir-ratio findings.

## Environment

The run environment is recorded in `results/sessionInfo.txt`.

