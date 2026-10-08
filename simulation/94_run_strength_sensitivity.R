#!/usr/bin/env Rscript
# 94_run_strength_sensitivity.R
# Runs strong/borderline/weak conditional-strength sensitivity after pilot freeze.

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

rms <- config$strength_rms_sensitivity
if (any(!is.finite(rms))) {
  stop(
    "strength_rms_sensitivity contains NA. Run 98_pilot_strength.R, ",
    "then rerun 97_preflight.R."
  )
}

n_cores <- as.integer(Sys.getenv("SIM_CORES", as.character(config$hpc_cores_default)))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
resume <- !identical(Sys.getenv("SIM_RESUME", "1"), "0")
chunk_size <- as.integer(Sys.getenv("SIM_CHUNK_SIZE", as.character(config$replicate_chunk_size)))

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

grid <- build_strength_sensitivity_scenarios(
  config = config,
  calibrated_rms = rms
)

for (i in seq_len(nrow(grid))) {
  cat(sprintf("[strength %d/%d] %s\n", i, nrow(grid), grid$scenario_id[i]))
  run_primary_scenario(
    grid[i, , drop = FALSE],
    config = config,
    nrep = config$nrep_sensitivity,
    n_cores = n_cores,
    project_dir = project_dir,
    write_output = TRUE,
    cluster = cl,
    chunk_size = chunk_size,
    resume = resume,
    output_group = "sensitivity/strength"
  )
}

summarise_raw_directory(
  config,
  output_group = "sensitivity/strength",
  summary_label = "sensitivity_strength"
)
cat("Strength sensitivity complete.\n")
