source("R/00_setup.R")

manifest <- readr::read_csv(file.path(docs_dir, "data_manifest.csv"), show_col_types = FALSE)

download_one <- function(file_name, source_url) {
  destination <- file.path(raw_dir, file_name)
  relative_destination <- file.path("data", "raw", file_name)
  if (file.exists(destination)) {
    return(tibble(file_name = file_name, status = "already_present", path = relative_destination))
  }
  utils::download.file(source_url, destination, mode = "wb", quiet = TRUE)
  tibble(file_name = file_name, status = "downloaded", path = relative_destination)
}

download_log <- purrr::pmap_dfr(
  manifest[, c("file_name", "source_url")],
  download_one
)

readr::write_csv(download_log, file.path(results_dir, "download_log.csv"))
