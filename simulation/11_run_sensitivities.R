# 11_run_sensitivities.R
# Prespecified targeted sensitivities. These do not alter the theorem definitions.

build_same_sample_scenarios <- function(config = default_sim_config()) {
  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == "all_direct_classes", ,
    drop = FALSE
  ]
  grid <- expand.grid(
    role = roles,
    beta_x_target = config$beta_x_primary,
    outcome_type = c("continuous", "binary"),
    stringsAsFactors = FALSE
  )
  grid$instrument_config <- dc$instrument_config
  grid$n_direct_x <- dc$n_direct_x
  grid$n_direct_z <- dc$n_direct_z
  grid$n_direct_shared <- dc$n_direct_shared
  grid$sample_design <- "same_sample"
  grid$sample_overlap <- 1
  grid$direct_shared_rho <- config$direct_shared_rho_primary
  grid$effect_rms <- config$effect_rms_primary
  grid$n_exposure <- config$n_same_sample
  grid$n_outcome <- config$n_same_sample
  grid$target_prevalence <- ifelse(
    grid$outcome_type == "binary",
    config$binary_prevalence_primary, NA_real_
  )
  grid$architecture_seed <- architecture_rng_seed(config, "all_direct_classes")
  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i], grid$outcome_type[i],
      sprintf("bx%.2f", grid$beta_x_target[i]),
      "all_direct_classes", "same_sample", sep = "__"
    ),
    character(1)
  )
  grid
}

build_alignment_sensitivity_scenarios <- function(config = default_sim_config()) {
  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == "all_direct_classes", ,
    drop = FALSE
  ]
  grid <- expand.grid(
    role = roles,
    beta_x_target = config$beta_x_primary,
    outcome_type = c("continuous", "binary"),
    direct_shared_rho = config$direct_shared_rho_sensitivity,
    stringsAsFactors = FALSE
  )
  grid$instrument_config <- dc$instrument_config
  grid$n_direct_x <- dc$n_direct_x
  grid$n_direct_z <- dc$n_direct_z
  grid$n_direct_shared <- dc$n_direct_shared
  grid$sample_design <- "two_sample"
  grid$sample_overlap <- 0
  grid$effect_rms <- config$effect_rms_primary
  grid$n_exposure <- config$n_exposure
  grid$n_outcome <- config$n_outcome
  grid$target_prevalence <- ifelse(
    grid$outcome_type == "binary",
    config$binary_prevalence_primary, NA_real_
  )
  grid$architecture_seed <- architecture_rng_seed(config, "all_direct_classes")
  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i], grid$outcome_type[i],
      sprintf("bx%.2f", grid$beta_x_target[i]),
      sprintf("rho%.2f", grid$direct_shared_rho[i]),
      "alignment_sensitivity", sep = "__"
    ),
    character(1)
  )
  grid
}

# Pilot effect-size calibration against a reference scenario.
# The resulting multiplier must be frozen before the final simulation is run.
pilot_calibrate_effect_rms <- function(
    target_F = 25,
    role = "independent_cause",
    outcome_type = "continuous",
    beta_x_target = 0.10,
    instrument_config = "all_direct_classes",
    config = default_sim_config(),
    candidate_rms = seq(0.02, 0.12, by = 0.01),
    nrep_pilot = 100L,
    n_cores = 1L,
    cluster = NULL,
    project_dir = getwd()) {

  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == instrument_config, ,
    drop = FALSE
  ]
  if (nrow(dc) != 1L) stop("Unknown instrument_config.")

  rows <- vector("list", length(candidate_rms))
  for (i in seq_along(candidate_rms)) {
    rms <- candidate_rms[i]
    sc <- data.frame(
      role = role,
      beta_x_target = beta_x_target,
      instrument_config = instrument_config,
      outcome_type = outcome_type,
      n_direct_x = dc$n_direct_x,
      n_direct_z = dc$n_direct_z,
      n_direct_shared = dc$n_direct_shared,
      sample_design = "two_sample",
      sample_overlap = 0,
      direct_shared_rho = config$direct_shared_rho_primary,
      effect_rms = rms,
      n_exposure = config$n_exposure,
      n_outcome = config$n_outcome,
      target_prevalence = if (outcome_type == "binary")
        config$binary_prevalence_primary else NA_real_,
      architecture_seed = architecture_rng_seed(config, instrument_config),
      scenario_id = paste("pilot", role, outcome_type, rms, sep = "__"),
      stringsAsFactors = FALSE
    )
    raw <- run_primary_scenario(
      sc, config, nrep = nrep_pilot, n_cores = n_cores,
      project_dir = project_dir, write_output = FALSE,
      cluster = cluster, resume = FALSE,
      output_group = "pilot"
    )
    rows[[i]] <- data.frame(
      effect_rms = rms,
      median_F_min = median(raw$F_min[raw$converged], na.rm = TRUE),
      mean_F_min = mean(raw$F_min[raw$converged], na.rm = TRUE)
    )
  }

  tab <- do.call(rbind, rows)
  tab$distance_to_target <- abs(tab$median_F_min - target_F)
  best <- tab[which.min(tab$distance_to_target), , drop = FALSE]
  list(target_F = target_F, best = best, grid = tab)
}


build_n_sensitivity_scenarios <- function(
    config = default_sim_config(),
    N_values = c(5000L, 50000L)) {

  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == "all_direct_classes", ,
    drop = FALSE
  ]

  grid <- expand.grid(
    role = roles,
    beta_x_target = config$beta_x_primary,
    outcome_type = c("continuous", "binary"),
    N_value = N_values,
    stringsAsFactors = FALSE
  )
  grid$instrument_config <- dc$instrument_config
  grid$n_direct_x <- dc$n_direct_x
  grid$n_direct_z <- dc$n_direct_z
  grid$n_direct_shared <- dc$n_direct_shared
  grid$sample_design <- "two_sample"
  grid$sample_overlap <- 0
  grid$direct_shared_rho <- config$direct_shared_rho_primary
  grid$effect_rms <- config$effect_rms_primary
  grid$n_exposure <- grid$N_value
  grid$n_outcome <- grid$N_value
  grid$target_prevalence <- ifelse(
    grid$outcome_type == "binary",
    config$binary_prevalence_primary, NA_real_
  )
  grid$architecture_seed <- architecture_rng_seed(config, "all_direct_classes")
  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i], grid$outcome_type[i],
      sprintf("bx%.2f", grid$beta_x_target[i]),
      paste0("N", grid$N_value[i]), "N_sensitivity", sep = "__"
    ),
    character(1)
  )
  grid
}

build_binary_prevalence_sensitivity_scenarios <- function(
    config = default_sim_config()) {

  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == "all_direct_classes", ,
    drop = FALSE
  ]
  grid <- expand.grid(
    role = roles,
    beta_x_target = config$beta_x_primary,
    target_prevalence = config$binary_prevalence_sensitivity,
    stringsAsFactors = FALSE
  )
  grid$outcome_type <- "binary"
  grid$instrument_config <- dc$instrument_config
  grid$n_direct_x <- dc$n_direct_x
  grid$n_direct_z <- dc$n_direct_z
  grid$n_direct_shared <- dc$n_direct_shared
  grid$sample_design <- "two_sample"
  grid$sample_overlap <- 0
  grid$direct_shared_rho <- config$direct_shared_rho_primary
  grid$effect_rms <- config$effect_rms_primary
  grid$n_exposure <- config$n_exposure
  grid$n_outcome <- config$n_outcome
  grid$architecture_seed <- architecture_rng_seed(config, "all_direct_classes")
  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i], "binary",
      sprintf("bx%.2f", grid$beta_x_target[i]),
      sprintf("prev%.2f", grid$target_prevalence[i]),
      "prevalence_sensitivity", sep = "__"
    ),
    character(1)
  )
  grid
}

build_strength_sensitivity_scenarios <- function(
    config = default_sim_config(),
    calibrated_rms = c(strong = 0.08, borderline = 0.05, weak = 0.03)) {

  # IMPORTANT: calibrated_rms must be frozen after the prespecified pilot.
  # The labels are not accepted as evidence of actual strength; realised conditional F
  # remains the reported strength diagnostic.
  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == "all_direct_classes", ,
    drop = FALSE
  ]
  grid <- expand.grid(
    role = roles,
    beta_x_target = config$beta_x_primary,
    outcome_type = c("continuous", "binary"),
    strength_level = names(calibrated_rms),
    stringsAsFactors = FALSE
  )
  grid$effect_rms <- unname(calibrated_rms[grid$strength_level])
  grid$instrument_config <- dc$instrument_config
  grid$n_direct_x <- dc$n_direct_x
  grid$n_direct_z <- dc$n_direct_z
  grid$n_direct_shared <- dc$n_direct_shared
  grid$sample_design <- "two_sample"
  grid$sample_overlap <- 0
  grid$direct_shared_rho <- config$direct_shared_rho_primary
  grid$n_exposure <- config$n_exposure
  grid$n_outcome <- config$n_outcome
  grid$target_prevalence <- ifelse(
    grid$outcome_type == "binary",
    config$binary_prevalence_primary, NA_real_
  )
  grid$architecture_seed <- architecture_rng_seed(config, "all_direct_classes")
  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i], grid$outcome_type[i],
      sprintf("bx%.2f", grid$beta_x_target[i]),
      grid$strength_level[i], "strength_sensitivity", sep = "__"
    ),
    character(1)
  )
  grid
}


build_effect_size_sensitivity_scenarios <- function(
    config = default_sim_config(),
    beta_x_values = 0.05) {

  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config[
    config$direct_instrument_config$instrument_config == "all_direct_classes", ,
    drop = FALSE
  ]
  grid <- expand.grid(
    role = roles,
    beta_x_target = beta_x_values,
    outcome_type = c("continuous", "binary"),
    stringsAsFactors = FALSE
  )
  grid$instrument_config <- dc$instrument_config
  grid$n_direct_x <- dc$n_direct_x
  grid$n_direct_z <- dc$n_direct_z
  grid$n_direct_shared <- dc$n_direct_shared
  grid$sample_design <- "two_sample"
  grid$sample_overlap <- 0
  grid$direct_shared_rho <- config$direct_shared_rho_primary
  grid$effect_rms <- config$effect_rms_primary
  grid$n_exposure <- config$n_exposure
  grid$n_outcome <- config$n_outcome
  grid$target_prevalence <- ifelse(
    grid$outcome_type == "binary",
    config$binary_prevalence_primary, NA_real_
  )
  grid$architecture_seed <- architecture_rng_seed(config, "all_direct_classes")
  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i], grid$outcome_type[i],
      sprintf("bx%.2f", grid$beta_x_target[i]),
      "effect_size_sensitivity", sep = "__"
    ),
    character(1)
  )
  grid
}
