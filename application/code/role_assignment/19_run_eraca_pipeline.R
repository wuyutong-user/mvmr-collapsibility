#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(prefix, default = "") {
  hit <- args[grepl(paste0("^", prefix, "="), args)]
  if (length(hit)) sub(paste0("^", prefix, "="), "", hit[[1]]) else default
}

step <- value_arg("--step", "all")
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
skip_tests <- "--skip-tests" %in% args
rscript <- file.path(R.home("bin"), "Rscript")

run_script <- function(file, extra_args = character()) {
  path <- file.path(script_dir, file)
  if (!file.exists(path)) stop("Pipeline script not found: ", path, call. = FALSE)
  cat("\n=== Running ", file, " ===\n", sep = "")
  status <- system2(rscript, c(shQuote(path), shQuote(extra_args)))
  if (!identical(status, 0L)) stop(file, " exited with status ", status, call. = FALSE)
  invisible(status)
}

run_tests <- function() {
  run_script(file.path("tests", "test_eraca_synthetic_triads.R"))
}

common <- if (overwrite) "--overwrite" else character()

if (step == "test") {
  run_tests()
} else if (step == "14") {
  run_script("14_prepare_all_trait_instruments.R", c(common, if (dry_run) "--dry-run" else character()))
} else if (step == "15") {
  run_script("15_run_bidirectional_uvmr.R", common)
} else if (step == "15-aggregate") {
  run_script("15_run_bidirectional_uvmr.R", "--aggregate-only")
} else if (step == "16") {
  run_script("16_build_local_marginal_graphs.R")
} else if (step == "17") {
  run_script("17_run_local_conditional_learning.R", common)
} else if (step == "17-aggregate") {
  run_script("17_run_local_conditional_learning.R", "--aggregate-only")
} else if (step == "18") {
  run_script("18_assign_roles_and_compare_response.R")
} else if (step == "all") {
  if (!skip_tests) run_tests()
  run_script("14_prepare_all_trait_instruments.R", common)
  run_script("15_run_bidirectional_uvmr.R", common)
  run_script("16_build_local_marginal_graphs.R")
  run_script("17_run_local_conditional_learning.R", common)
  run_script("18_assign_roles_and_compare_response.R")
} else {
  stop(
    "Unknown --step=", step,
    ". Use test, 14, 15, 15-aggregate, 16, 17, 17-aggregate, 18, or all.",
    call. = FALSE
  )
}

cat("\nERACA pipeline request completed: step=", step, "\n", sep = "")
