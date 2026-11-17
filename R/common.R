# Shared, privacy-aware helpers. No external API calls or personal paths.
project_root <- function() {
  path <- getOption("chikv.project.root", getwd())
  normalizePath(path, mustWork = TRUE)
}
private_data_path <- function(name) file.path(project_root(), "data", "private", name)
private_result_path <- function(name) {
  out <- file.path(project_root(), "results", "private")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  file.path(out, name)
}
require_columns <- function(dat, columns, name = "input") {
  absent <- setdiff(columns, names(dat))
  if (length(absent)) stop(name, " is missing columns: ", paste(absent, collapse = ", "))
  invisible(TRUE)
}
read_private_csv <- function(name) {
  f <- private_data_path(name)
  if (!file.exists(f)) stop("Required restricted input is absent: ", name)
  readr::read_csv(f, show_col_types = FALSE)
}
safe_write_csv <- function(x, name) {
  # ALL generated rows, including apparently anonymous case rows, stay private.
  readr::write_csv(x, private_result_path(name))
}
make_onset_date <- function(x) {
  d <- as.Date(x)
  if (anyNA(d)) stop("Missing or invalid onset dates. Expect ISO 8601 YYYY-MM-DD.")
  d
}
check_case_ids <- function(x) {
  if (anyNA(x) || any(!nzchar(trimws(as.character(x)))) || anyDuplicated(x)) {
    stop("Case IDs must be present and unique (IDs stay in private files).")
  }
  invisible(TRUE)
}
