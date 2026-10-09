# 00_config.R
# Master configuration and deterministic seed helpers.
# Simulation protocol: MVMR collapsibility / adjustment-response framework, v1.1

`%||%` <- function(x, y) if (is.null(x)) y else x

stable_seed <- function(...) {
  txt <- paste(..., collapse = "|")
  ints <- utf8ToInt(txt)
  if (length(ints) == 0L) return(1L)
  # Deterministic, platform-independent integer hash in the valid R seed range.
  h <- 104729
  mod <- 2147483646
  for (i in seq_along(ints)) {
    h <- (h * 1009 + ints[i] * (i + 37)) %% mod
  }
  as.integer(h + 1)
}




scenario_rng_key <- function(sc) {
  # Common-random-number key: deliberately excludes beta_X, strength scale,
  # alignment target, N, prevalence and sample overlap. Thus controlled
  # sensitivity contrasts reuse the same underlying genotype/U/residual streams.
  paste(sc$role, sc$outcome_type, sc$instrument_config, sep = "|")
}

architecture_rng_seed <- function(config, instrument_config) {
  # Deliberately independent of effect magnitude and alignment target so
  # controlled sensitivity analyses reuse the same latent architecture draws.
  stable_seed(config$master_seed, "architecture", instrument_config)
}

binary_intercept_rng_seed <- function(config, sc) {
  # Same reference random stream across beta_X / strength / alignment / N /
  # prevalence / sample-design contrasts; only the linear predictor changes.
  stable_seed(config$master_seed, "binary_intercept", sc$role, sc$instrument_config)
}

bind_rows_fill <- function(lst) {
  if (length(lst) == 0L) return(data.frame())
  all_names <- unique(unlist(lapply(lst, names), use.names = FALSE))
  out <- lapply(lst, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[all_names]
  })
  do.call(rbind, out)
}


atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Atomic save failed: ", path)
  }
  invisible(path)
}

default_sim_config <- function(output_dir = file.path(getwd(), "simulation_output")) {
  list(
    protocol_version = "1.1",
    master_seed = 20260905L,
    alpha = 0.05,
    ci_level = 0.95,
    theorem_tolerance = 1e-10,

    # Primary inference convention.
    # "weighted_lm" reproduces the weighted-regression SE convention used by the
    # legacy simulation and TwoSampleMR::mv_multiple; fixed-effect SEs are also
    # retained as diagnostics.
    se_method_primary = "weighted_lm",
    p_reference_primary = "normal",

    # Six-DAG finite-sample primary design
    n_exposure = 10000L,
    n_outcome = 10000L,
    n_same_sample = 10000L,
    beta_x_primary = c(0, 0.10),
    path_coefficient = 0.30,
    u_coefficient = 0.30,
    residual_sd = 1.0,

    # Genotypes / instruments
    K = 100L,
    maf = 0.30,
    effect_rms_primary = 0.08,
    direct_shared_rho_primary = 0.80,
    direct_shared_rho_sensitivity = c(0, -0.80, 0.95),

    # Direct genetic-source configurations for individual-level six-DAG simulations.
    # These are NOT the theorem's summary-association classes.
    direct_instrument_config = data.frame(
      instrument_config = c(
        "all_direct_classes",
        "no_direct_x",
        "no_direct_z",
        "no_direct_shared"
      ),
      n_direct_x = c(30L, 0L, 60L, 50L),
      n_direct_z = c(30L, 60L, 0L, 50L),
      n_direct_shared = c(40L, 40L, 40L, 0L),
      stringsAsFactors = FALSE
    ),

    # Summary-association configurations for exact Theorem-1 calibration.
    theorem_instrument_config = data.frame(
      instrument_config = c(
        "all_summary_classes",
        "no_x_only",
        "no_z_only",
        "no_shared"
      ),
      n_x_only = c(30L, 0L, 60L, 50L),
      n_z_only = c(30L, 60L, 0L, 50L),
      n_shared = c(40L, 40L, 40L, 0L),
      stringsAsFactors = FALSE
    ),

    # Theorem calibration geometry
    theorem_rho_primary = 0.80,
    theorem_rho_sensitivity = c(0, -0.80, 0.95),
    theorem_effect_rms = 0.08,

    # Binary outcome
    binary_prevalence_primary = 0.10,
    binary_prevalence_sensitivity = 0.50,
    binary_reference_n = 100000L,

    # Monte Carlo replicates
    nrep_continuous_primary = 2000L,
    nrep_binary_primary = 1000L,
    nrep_sensitivity = 1000L,

    # HPC execution
    hpc_cores_default = 100L,
    replicate_chunk_size = 100L,
    resume_default = TRUE,

    # Practical instrument-selection sensitivity
    practical_selection_p = 5e-8,

    # Strength calibration targets (pilot calibration is optional and prespecified)
    conditional_F_target_strong = 25,
    conditional_F_target_borderline = 10,
    conditional_F_target_weak = 5,
    # Fill these ONLY after running 98_pilot_strength.R and freezing the pilot result.
    strength_rms_sensitivity = c(
      strong = NA_real_,
      borderline = NA_real_,
      weak = NA_real_
    ),

    # Numerics
    min_variance = 1e-14,
    min_instruments_univariable = 2L,
    min_instruments_multivariable = 3L,

    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
}

ensure_output_dirs <- function(config) {
  dirs <- c(
    config$output_dir,
    file.path(config$output_dir, "raw"),
    file.path(config$output_dir, "raw", "primary"),
    file.path(config$output_dir, "raw", "sensitivity"),
    file.path(config$output_dir, "summary"),
    file.path(config$output_dir, "theorem"),
    file.path(config$output_dir, "validation"),
    file.path(config$output_dir, "logs")
  )
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}


strength_calibration_file <- function(config) {
  file.path(config$output_dir, "validation", "strength_calibration_frozen.rds")
}

apply_frozen_strength_calibration <- function(config, require = FALSE) {
  f <- strength_calibration_file(config)
  if (!file.exists(f)) {
    if (isTRUE(require)) {
      stop(
        "Frozen strength calibration not found. Run 98_pilot_strength.R first: ", f
      )
    }
    return(config)
  }

  cal <- readRDS(f)
  if (!identical(cal$protocol_version, config$protocol_version)) {
    stop("Frozen strength calibration protocol version does not match code protocol.")
  }
  stopifnot(
    is.finite(cal$effect_rms_primary),
    length(cal$strength_rms_sensitivity) == 3L,
    all(is.finite(cal$strength_rms_sensitivity))
  )
  config$effect_rms_primary <- unname(cal$effect_rms_primary)
  config$strength_rms_sensitivity <- cal$strength_rms_sensitivity
  config$strength_calibration <- cal
  config
}
