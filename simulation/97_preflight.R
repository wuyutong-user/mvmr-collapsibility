#!/usr/bin/env Rscript
# 97_preflight.R
# Mandatory server-side preflight before production Monte Carlo.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
project_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
} else getwd()

options(mvmr.sim.project_dir = project_dir)
setwd(project_dir)

# Verify the frozen source manifest before parsing/running code.
manifest_file <- file.path(project_dir, "SOURCE_MANIFEST_MD5.tsv")
if (!file.exists(manifest_file)) {
  stop("Frozen source manifest is missing: SOURCE_MANIFEST_MD5.tsv")
}
manifest <- read.delim(
  manifest_file, stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(c("file", "md5") %in% names(manifest))) {
  stop("Malformed SOURCE_MANIFEST_MD5.tsv")
}
manifest_paths <- file.path(project_dir, manifest$file)
if (any(!file.exists(manifest_paths))) {
  stop("One or more source-manifest files are missing.")
}
observed_md5 <- unname(tools::md5sum(manifest_paths))
if (!identical(tolower(observed_md5), tolower(manifest$md5))) {
  bad <- manifest$file[tolower(observed_md5) != tolower(manifest$md5)]
  stop("Frozen source manifest mismatch: ", paste(bad, collapse = ", "))
}
cat(sprintf("SOURCE MANIFEST PASS: %d files\n", nrow(manifest)))

# Parse every project R file before sourcing anything.
r_files <- list.files(project_dir, pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  parse(file = f)
}
cat(sprintf("PARSE PASS: %d R files\n", length(r_files)))

source("bootstrap.R")
config <- default_sim_config(file.path(project_dir, "simulation_output"))
ensure_output_dirs(config)
config <- apply_frozen_strength_calibration(config, require = FALSE)

cat("Running formal validation tests...\n")
validation <- run_validation_tests(config, project_dir = project_dir, verbose = TRUE)
require_mvmr <- identical(Sys.getenv("REQUIRE_MVMR", "0"), "1")
if (require_mvmr && !isTRUE(validation$mvmr_package_crosscheck)) {
  stop(
    "REQUIRE_MVMR=1 but the MVMR::strength_mvmr package cross-check did not run. ",
    "Install MVMR and rerun preflight."
  )
}

cat("Running exact theorem modules...\n")
th1 <- run_theorem1_calibration(config, write_output = TRUE)
th23 <- run_binary_scale_calibration(config, write_output = TRUE)

# Small continuous smoke test
sc <- build_primary_scenarios(config)
scc <- sc[
  sc$role == "independent_cause" &
  sc$outcome_type == "continuous" &
  sc$beta_x_target == 0.10 &
  sc$instrument_config == "all_direct_classes",
  , drop = FALSE
]
cat("Running 10-replicate continuous smoke test...\n")
rawc <- run_primary_scenario(
  scc, config, nrep = 10L, n_cores = 1L,
  project_dir = project_dir, write_output = FALSE,
  resume = FALSE
)
stopifnot(all(rawc$converged))

# Small binary smoke test
scb <- sc[
  sc$role == "independent_cause" &
  sc$outcome_type == "binary" &
  sc$beta_x_target == 0.10 &
  sc$instrument_config == "all_direct_classes",
  , drop = FALSE
]
cat("Running 5-replicate binary smoke test...\n")
rawb <- run_primary_scenario(
  scb, config, nrep = 5L, n_cores = 1L,
  project_dir = project_dir, write_output = FALSE,
  resume = FALSE
)
stopifnot(all(rawb$converged))

audit <- list(
  timestamp = Sys.time(),
  protocol_version = config$protocol_version,
  source_manifest_pass = TRUE,
  parse_pass = TRUE,
  validation_pass = isTRUE(validation$passed),
  mvmr_package_crosscheck = isTRUE(validation$mvmr_package_crosscheck),
  theorem1_max_abs_identity_error = max(abs(th1$identity_error)),
  theorem23_max_abs_decomposition_error = max(abs(th23$decomposition_error)),
  continuous_smoke_all_converged = all(rawc$converged),
  binary_smoke_all_converged = all(rawb$converged),
  continuous_smoke_median_F_min = median(rawc$F_min, na.rm = TRUE),
  binary_smoke_median_F_min = median(rawb$F_min, na.rm = TRUE),
  binary_smoke_mean_prevalence = mean(rawb$prevalence, na.rm = TRUE),
  frozen_strength_calibration_present =
    file.exists(strength_calibration_file(config)),
  effect_rms_primary_used = config$effect_rms_primary
)
saveRDS(audit, file.path(config$output_dir, "validation", "preflight_audit.rds"))
capture.output(
  print(audit),
  file = file.path(config$output_dir, "validation", "preflight_audit.txt")
)

cat("\nPREFLIGHT PASS.\n")
cat("Review simulation_output/validation/preflight_audit.txt before production.\n")
