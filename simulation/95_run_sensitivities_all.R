#!/usr/bin/env Rscript
# 95_run_sensitivities_all.R
# Runs prespecified sensitivity grids after the primary production run.

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

preflight_file <- file.path(config$output_dir, "validation", "preflight_audit.rds")
if (!file.exists(preflight_file)) stop("Run Rscript 97_preflight.R first.")

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

run_grid <- function(grid, label) {
  cat("Sensitivity grid: ", label, " | scenarios=", nrow(grid), "\n", sep = "")
  for (i in seq_len(nrow(grid))) {
    cat(sprintf("[%s %d/%d] %s\n", label, i, nrow(grid), grid$scenario_id[i]))
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
      output_group = file.path("sensitivity", label)
    )
  }
  summarise_raw_directory(
    config,
    output_group = file.path("sensitivity", label),
    summary_label = paste0("sensitivity_", label)
  )
}

run_grid(build_same_sample_scenarios(config), "same_sample")
run_grid(build_alignment_sensitivity_scenarios(config), "alignment")
run_grid(build_n_sensitivity_scenarios(config), "sample_size")
run_grid(build_binary_prevalence_sensitivity_scenarios(config), "binary_prevalence")
run_grid(build_effect_size_sensitivity_scenarios(config), "effect_size_005")

cat("Sensitivity run complete.\n")
cat("Strength-level sensitivity is intentionally excluded until pilot-calibrated RMS values are frozen.\n")
