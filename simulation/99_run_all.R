#!/usr/bin/env Rscript
# 99_run_all.R
#
# Production runner.
# A single persistent PSOCK cluster is created and reused across all scenarios.
# Deterministic replicate-specific seeds make load-balanced scheduling reproducible.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
project_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
} else getwd()

options(mvmr.sim.project_dir = project_dir)
setwd(project_dir)
source("bootstrap.R")

config <- default_sim_config(file.path(project_dir, "simulation_output"))
ensure_output_dirs(config)
config <- apply_frozen_strength_calibration(config, require = TRUE)

# Mandatory preflight marker.
preflight_file <- file.path(config$output_dir, "validation", "preflight_audit.rds")
if (!file.exists(preflight_file)) {
  stop("Preflight audit not found. Run: Rscript 97_preflight.R")
}
preflight <- readRDS(preflight_file)
if (!isTRUE(preflight$parse_pass) || !isTRUE(preflight$validation_pass)) {
  stop("Preflight audit is not valid.")
}
if (!isTRUE(preflight$frozen_strength_calibration_present) ||
    !isTRUE(all.equal(
      unname(preflight$effect_rms_primary_used),
      unname(config$effect_rms_primary),
      tolerance = 0
    ))) {
  stop(
    "Preflight was not run after the frozen strength calibration. ",
    "Run Rscript 97_preflight.R again."
  )
}

run_full <- identical(Sys.getenv("RUN_FULL_SIM", "0"), "1")
if (!run_full) {
  cat("RUN_FULL_SIM is not 1; production grid not launched.\n")
  quit(save = "no", status = 0)
}

n_cores <- as.integer(Sys.getenv(
  "SIM_CORES",
  as.character(config$hpc_cores_default)
))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L

resume <- !identical(Sys.getenv("SIM_RESUME", "1"), "0")
chunk_size <- as.integer(Sys.getenv(
  "SIM_CHUNK_SIZE",
  as.character(config$replicate_chunk_size)
))
if (!is.finite(chunk_size) || chunk_size < 1L) chunk_size <- config$replicate_chunk_size

cat(sprintf(
  "Protocol %s | cores=%d | chunk_size=%d | resume=%s\n",
  config$protocol_version, n_cores, chunk_size, resume
))

cl <- NULL
if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterCall(cl, function(pd) {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    setwd(pd)
    options(mvmr.sim.project_dir = pd)
    source(file.path(pd, "bootstrap.R"))
    NULL
  }, project_dir)
}

scenarios <- build_primary_scenarios(config)

for (i in seq_len(nrow(scenarios))) {
  cat(sprintf(
    "[%d/%d] %s | %s\n",
    i, nrow(scenarios), format(Sys.time()), scenarios$scenario_id[i]
  ))
  run_primary_scenario(
    scenarios[i, , drop = FALSE],
    config = config,
    n_cores = n_cores,
    project_dir = project_dir,
    write_output = TRUE,
    cluster = cl,
    chunk_size = chunk_size,
    resume = resume,
    output_group = "primary"
  )
}

cat("Summarising primary-grid raw results...\n")
summarise_raw_directory(config, output_group = "primary", summary_label = "primary")
cat("PRIMARY PRODUCTION RUN COMPLETE.\n")
