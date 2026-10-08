#!/usr/bin/env Rscript
# 98_pilot_strength.R
# Prespecified strength pilot. Run BEFORE production and freeze selected RMS values.
#
# Reference scenario:
#   independent cause, continuous outcome, beta_X=0.10,
#   all direct genetic classes, non-overlapping two-sample.
#
# The same candidate RMS grid is evaluated once. The RMS values closest to
# conditional F targets 25, 10, and 5 are reported. These are recommendations;
# they must be frozen in 00_config.R before strength sensitivity analyses.

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

preflight_file <- file.path(config$output_dir, "validation", "preflight_audit.rds")
if (!file.exists(preflight_file)) {
  stop("Run Rscript 97_preflight.R before the strength pilot.")
}
preflight <- readRDS(preflight_file)
if (!isTRUE(preflight$parse_pass) || !isTRUE(preflight$validation_pass)) {
  stop("Existing preflight audit is not valid.")
}

nrep <- as.integer(Sys.getenv("PILOT_NREP", "100"))
n_cores <- as.integer(Sys.getenv("SIM_CORES", as.character(config$hpc_cores_default)))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
candidate_rms <- seq(0.02, 0.14, by = 0.01)

cl <- NULL
if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(pd) {
    Sys.setenv(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    setwd(pd)
    options(mvmr.sim.project_dir = pd)
    source(file.path(pd, "bootstrap.R"))
    NULL
  }, project_dir)
}

# Run one grid and retain the full F curve.
pilot <- pilot_calibrate_effect_rms(
  target_F = config$conditional_F_target_strong,
  config = config,
  candidate_rms = candidate_rms,
  nrep_pilot = nrep,
  n_cores = n_cores,
  cluster = cl,
  project_dir = project_dir
)
tab <- pilot$grid

targets <- c(
  strong = config$conditional_F_target_strong,
  borderline = config$conditional_F_target_borderline,
  weak = config$conditional_F_target_weak
)

pick <- lapply(names(targets), function(nm) {
  target <- targets[[nm]]
  j <- which.min(abs(tab$median_F_min - target))
  data.frame(
    strength_level = nm,
    target_F = target,
    recommended_effect_rms = tab$effect_rms[j],
    achieved_median_F_min = tab$median_F_min[j],
    achieved_mean_F_min = tab$mean_F_min[j],
    stringsAsFactors = FALSE
  )
})
recommend <- do.call(rbind, pick)

saveRDS(
  list(grid = tab, recommendation = recommend, nrep = nrep),
  file.path(config$output_dir, "validation", "strength_pilot.rds")
)
write.csv(
  tab,
  file.path(config$output_dir, "validation", "strength_pilot_grid.csv"),
  row.names = FALSE
)
write.csv(
  recommend,
  file.path(config$output_dir, "validation", "strength_rms_recommendation.csv"),
  row.names = FALSE
)

frozen <- list(
  protocol_version = config$protocol_version,
  pilot_timestamp = Sys.time(),
  pilot_nrep = nrep,
  reference_scenario = paste(
    "independent_cause", "continuous", "beta_X=0.10",
    "all_direct_classes", "two_sample", sep = " | "
  ),
  target_F = targets,
  effect_rms_primary =
    recommend$recommended_effect_rms[recommend$strength_level == "strong"],
  strength_rms_sensitivity = setNames(
    recommend$recommended_effect_rms,
    recommend$strength_level
  ),
  achieved_median_F_min = setNames(
    recommend$achieved_median_F_min,
    recommend$strength_level
  )
)
atomic_save_rds(
  frozen,
  file.path(config$output_dir, "validation", "strength_calibration_frozen.rds")
)

cat("Strength pilot complete.\n\n")
print(recommend)
cat("\nThe deterministic pilot recommendation has been frozen to:\n")
cat("simulation_output/validation/strength_calibration_frozen.rds\n")
cat("Production scripts will load this frozen calibration automatically.\n")
cat("Now rerun: Rscript 97_preflight.R\n")
