source("R/00_setup.R")

capture.output(
  sessionInfo(),
  file = file.path(results_dir, "sessionInfo.txt")
)

