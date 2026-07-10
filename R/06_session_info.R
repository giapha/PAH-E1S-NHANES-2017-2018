source("R/00_setup.R")

session_lines <- capture.output(sessionInfo())
writeLines(
  sub("[[:space:]]+$", "", session_lines),
  file.path(results_dir, "sessionInfo.txt")
)
