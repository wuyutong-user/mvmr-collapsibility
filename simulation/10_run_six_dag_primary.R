# 10_run_six_dag_primary.R
# Primary six-DAG finite-sample simulation:
#   - strict same-set IVW/MVMR-IVW is primary
#   - non-overlapping two-sample design is primary
#   - practical different-set estimates are stored as a sensitivity in the same replicate

make_architecture_for_scenario <- function(scenario, config) {
  make_direct_architecture(
    K = config$K,
    n_direct_x = scenario$n_direct_x,
    n_direct_z = scenario$n_direct_z,
    n_direct_shared = scenario$n_direct_shared,
    shared_rho = scenario$direct_shared_rho,
    effect_rms = scenario$effect_rms,
    maf = config$maf,
    seed = scenario$architecture_seed
  )
}

build_primary_scenarios <- function(config = default_sim_config()) {
  roles <- get_role_registry(config$path_coefficient)$role
  dc <- config$direct_instrument_config

  grid <- expand.grid(
    role = roles,
    beta_x_target = config$beta_x_primary,
    instrument_config = dc$instrument_config,
    outcome_type = c("continuous", "binary"),
    stringsAsFactors = FALSE
  )
  grid <- merge(grid, dc, by = "instrument_config", sort = FALSE)
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

  grid$architecture_seed <- vapply(
    grid$instrument_config,
    function(x) architecture_rng_seed(config, x),
    integer(1)
  )

  grid$scenario_id <- vapply(
    seq_len(nrow(grid)),
    function(i) paste(
      grid$role[i],
      grid$outcome_type[i],
      sprintf("bx%.2f", grid$beta_x_target[i]),
      grid$instrument_config[i],
      "two_sample",
      sep = "__"
    ),
    character(1)
  )
  grid
}

prepare_scenario_truth <- function(scenario_row, config = default_sim_config()) {
  scenario <- as.list(scenario_row)
  arch <- make_architecture_for_scenario(scenario, config)

  binary_cal <- NULL
  if (identical(scenario$outcome_type, "binary")) {
    binary_cal <- calibrate_binary_intercept(
      role = scenario$role,
      architecture = arch,
      beta_x = scenario$beta_x_target,
      target_prevalence = scenario$target_prevalence,
      role_registry = get_role_registry(config$path_coefficient),
      u_coefficient = config$u_coefficient,
      residual_sd = config$residual_sd,
      n_reference = config$binary_reference_n,
      seed = binary_intercept_rng_seed(config, scenario)
    )
  }

  list(
    scenario = scenario,
    architecture = arch,
    binary_calibration = binary_cal
  )
}

generate_samples_for_replicate <- function(scenario_truth, config, replicate) {
  sc <- scenario_truth$scenario
  arch <- scenario_truth$architecture

  rng_key <- scenario_rng_key(sc)
  seed_e <- stable_seed(config$master_seed, rng_key, replicate, "exposure")
  seed_y <- stable_seed(config$master_seed, rng_key, replicate, "outcome")

  if (sc$outcome_type == "continuous") {
    if (sc$sample_design == "two_sample") {
      exposure_pop <- generate_continuous_population(
        sc$role, sc$n_exposure, arch, sc$beta_x_target,
        role_registry = get_role_registry(config$path_coefficient),
        u_coefficient = config$u_coefficient,
        residual_sd = config$residual_sd,
        seed = seed_e, id_prefix = "E"
      )
      outcome_pop <- generate_continuous_population(
        sc$role, sc$n_outcome, arch, sc$beta_x_target,
        role_registry = get_role_registry(config$path_coefficient),
        u_coefficient = config$u_coefficient,
        residual_sd = config$residual_sd,
        seed = seed_y, id_prefix = "Y"
      )
    } else if (sc$sample_design == "same_sample") {
      one <- generate_continuous_population(
        sc$role, sc$n_exposure, arch, sc$beta_x_target,
        role_registry = get_role_registry(config$path_coefficient),
        u_coefficient = config$u_coefficient,
        residual_sd = config$residual_sd,
        seed = seed_e, id_prefix = "S"
      )
      exposure_pop <- outcome_pop <- one
    } else stop("Unknown sample_design.")
  } else {
    eta0 <- scenario_truth$binary_calibration$intercept
    if (sc$sample_design == "two_sample") {
      exposure_pop <- generate_binary_population(
        sc$role, sc$n_exposure, arch, sc$beta_x_target, eta0,
        role_registry = get_role_registry(config$path_coefficient),
        u_coefficient = config$u_coefficient,
        residual_sd = config$residual_sd,
        seed = seed_e, id_prefix = "E"
      )
      outcome_pop <- generate_binary_population(
        sc$role, sc$n_outcome, arch, sc$beta_x_target, eta0,
        role_registry = get_role_registry(config$path_coefficient),
        u_coefficient = config$u_coefficient,
        residual_sd = config$residual_sd,
        seed = seed_y, id_prefix = "Y"
      )
    } else if (sc$sample_design == "same_sample") {
      one <- generate_binary_population(
        sc$role, sc$n_exposure, arch, sc$beta_x_target, eta0,
        role_registry = get_role_registry(config$path_coefficient),
        u_coefficient = config$u_coefficient,
        residual_sd = config$residual_sd,
        seed = seed_e, id_prefix = "S"
      )
      exposure_pop <- outcome_pop <- one
    } else stop("Unknown sample_design.")
  }

  list(exposure = exposure_pop, outcome = outcome_pop)
}

run_one_replicate <- function(scenario_truth, config, replicate) {
  sc <- scenario_truth$scenario
  role_spec <- get_role_spec(sc$role, get_role_registry(config$path_coefficient))
  # Exact population X/Z summary support is algebraically available for the
  # continuous DGP and for binary roles in which Z is not downstream of binary Y.
  pop_truth_available <- !(
    sc$outcome_type == "binary" &&
    sc$role %in% c("collider", "downstream_surrogate_outcome")
  )
  if (pop_truth_available) {
    pop_xz <- population_xz_summary_truth(
      role = sc$role,
      architecture = scenario_truth$architecture,
      beta_x = sc$beta_x_target,
      role_registry = get_role_registry(config$path_coefficient)
    )
    pop_class <- classify_summary_support(pop_xz$gamma_x_true, pop_xz$gamma_z_true)
  } else {
    pop_class <- rep(NA_character_, nrow(scenario_truth$architecture))
  }

  ans <- tryCatch({
    samp <- generate_samples_for_replicate(scenario_truth, config, replicate)

    sd <- make_summary_data(
      exposure_population = samp$exposure,
      outcome_population = samp$outcome,
      outcome_type = sc$outcome_type,
      min_variance = config$min_variance
    )

    # Strict primary calibration requires all prespecified K SNP rows to survive
    # GWAS estimation. A failed binary SNP fit is a replicate failure, not a
    # silently reduced instrument set.
    if (nrow(sd) != config$K || any(!sd$outcome_gwas_converged)) {
      stop("strict_same_set_requires_all_K_outcome_GWAS_fits")
    }
    required_finite <- with(sd,
      is.finite(beta_x) & is.finite(se_x) &
      is.finite(beta_z) & is.finite(se_z) &
      is.finite(cov_xz) &
      is.finite(beta_y) & is.finite(se_y) & se_y > 0
    )
    if (!all(required_finite)) {
      stop("strict_same_set_contains_nonfinite_summary_associations")
    }

    same <- fit_strict_same_set(sd, se_method = config$se_method_primary)
    if (!isTRUE(same$converged)) stop("same_set_fit_failed: ", same$reason)

    # Role-implied prediction is separate from the algebraic sample identity.
    # For continuous outcomes the DGP supplies population gamma_X/gamma_Z and
    # alpha_Z exactly, so we can evaluate the theorem-implied structural response
    # using the realised outcome weighting matrix. For realistic binary outcomes
    # this numeric structural component is deliberately left unidentified.
    if (sc$outcome_type == "continuous") {
      w_realised <- 1 / sd$se_y^2
      delta_role_population <- weighted_projection(
        pop_xz$gamma_x_true, pop_xz$gamma_z_true, w_realised
      )
      role_predicted_response <-
        role_spec$alpha_z_theorem * delta_role_population
      role_response_error <- same$response - role_predicted_response
    } else {
      delta_role_population <- NA_real_
      role_predicted_response <- NA_real_
      role_response_error <- NA_real_
    }

    strength <- conditional_strength(
      sd, idx = seq_len(nrow(sd)), min_variance = config$min_variance
    )

    practical <- fit_practical_different_set(
      sd,
      p_threshold = config$practical_selection_p,
      se_method = config$se_method_primary
    )

    zcrit <- qnorm(1 - (1 - config$ci_level) / 2)
    target <- sc$beta_x_target

    cover_u <- as.integer(
      target >= same$beta_unadj - zcrit * same$se_unadj &&
      target <= same$beta_unadj + zcrit * same$se_unadj
    )
    cover_a <- as.integer(
      target >= same$beta_adj - zcrit * same$se_adj &&
      target <= same$beta_adj + zcrit * same$se_adj
    )

    data.frame(
      simulation_module = "six_dag_finite_sample",
      scenario_id = sc$scenario_id,
      replicate = replicate,
      role = sc$role,
      outcome_type = sc$outcome_type,
      sample_design = sc$sample_design,
      sample_overlap = sc$sample_overlap,
      beta_x_target = target,
      model_alpha_z = role_spec$alpha_z_theorem,
      nominal_adjustment_class = role_spec$nominal_adjustment_class,
      instrument_config = sc$instrument_config,
      K = config$K,
      n_same_set = same$n_instruments,
      n_direct_x = sc$n_direct_x,
      n_direct_z = sc$n_direct_z,
      n_direct_shared = sc$n_direct_shared,
      direct_shared_rho_target = sc$direct_shared_rho,
      effect_rms = sc$effect_rms,
      population_summary_support_available = pop_truth_available,
      n_population_x_only = if (pop_truth_available) sum(pop_class == "X_only") else NA_integer_,
      n_population_z_only = if (pop_truth_available) sum(pop_class == "Z_only") else NA_integer_,
      n_population_shared = if (pop_truth_available) sum(pop_class == "shared") else NA_integer_,
      n_population_null_both = if (pop_truth_available) sum(pop_class == "null_both") else NA_integer_,

      rho_w = same$rho_w,
      delta_xz = same$delta_xz,
      condition_number = same$condition_number,
      F_x_given_z = strength$F_x_given_z,
      F_z_given_x = strength$F_z_given_x,
      F_min = strength$F_min,
      F_x_given_z_package = strength$F_x_given_z_package,
      F_z_given_x_package = strength$F_z_given_x_package,
      F_min_package = strength$F_min_package,

      beta_unadj = same$beta_unadj,
      se_unadj = same$se_unadj,
      se_unadj_fixed = same$se_unadj_fixed,
      se_unadj_weighted_lm = same$se_unadj_weighted_lm,
      p_unadj = same$p_unadj,
      p_unadj_t = same$p_unadj_t,
      sigma_unadj = same$sigma_unadj,
      beta_adj = same$beta_adj,
      se_adj = same$se_adj,
      se_adj_fixed = same$se_adj_fixed,
      se_adj_weighted_lm = same$se_adj_weighted_lm,
      p_adj = same$p_adj,
      p_adj_t = same$p_adj_t,
      sigma_adj = same$sigma_adj,
      alpha_z = same$alpha_z,
      se_alpha_z = same$se_alpha_z,
      response = same$response,
      algebraic_predicted_response = same$predicted_response,
      identity_error = same$identity_error,
      delta_role_population = delta_role_population,
      role_predicted_response = role_predicted_response,
      role_response_error = role_response_error,

      cover_unadj = cover_u,
      cover_adj = cover_a,
      reject_unadj = as.integer(same$p_unadj < config$alpha),
      reject_adj = as.integer(same$p_adj < config$alpha),

      practical_converged = practical$converged,
      practical_beta_unadj = if (isTRUE(practical$converged)) practical$beta_unadj else NA_real_,
      practical_se_unadj = if (isTRUE(practical$converged)) practical$se_unadj else NA_real_,
      practical_se_unadj_fixed = if (isTRUE(practical$converged)) practical$se_unadj_fixed else NA_real_,
      practical_se_unadj_weighted_lm = if (isTRUE(practical$converged)) practical$se_unadj_weighted_lm else NA_real_,
      practical_p_unadj = if (isTRUE(practical$converged)) practical$p_unadj else NA_real_,
      practical_beta_adj = if (isTRUE(practical$converged)) practical$beta_adj else NA_real_,
      practical_se_adj = if (isTRUE(practical$converged)) practical$se_adj else NA_real_,
      practical_se_adj_fixed = if (isTRUE(practical$converged)) practical$se_adj_fixed else NA_real_,
      practical_se_adj_weighted_lm = if (isTRUE(practical$converged)) practical$se_adj_weighted_lm else NA_real_,
      practical_p_adj = if (isTRUE(practical$converged)) practical$p_adj else NA_real_,
      practical_response = if (isTRUE(practical$converged)) practical$response else NA_real_,
      practical_n_unadj = practical$n_unadj %||% NA_integer_,
      practical_n_adj = practical$n_adj %||% NA_integer_,

      prevalence = if (sc$outcome_type == "binary") samp$outcome$prevalence else NA_real_,
      calibrated_prevalence = if (sc$outcome_type == "binary")
        scenario_truth$binary_calibration$calibrated_prevalence else NA_real_,

      converged = TRUE,
      failure_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      simulation_module = "six_dag_finite_sample",
      scenario_id = sc$scenario_id,
      replicate = replicate,
      role = sc$role,
      outcome_type = sc$outcome_type,
      sample_design = sc$sample_design,
      sample_overlap = sc$sample_overlap,
      beta_x_target = sc$beta_x_target,
      instrument_config = sc$instrument_config,
      K = config$K,
      converged = FALSE,
      failure_reason = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  ans
}

run_replicates <- function(
    scenario_truth,
    config,
    replicates,
    n_cores = 1L,
    project_dir = getwd(),
    cluster = NULL) {

  reps <- as.integer(replicates)
  if (length(reps) == 0L) return(data.frame())

  if (!is.null(cluster)) {
    out <- parallel::parLapplyLB(
      cluster, reps,
      function(r, scenario_truth, config) {
        run_one_replicate(scenario_truth, config, r)
      },
      scenario_truth = scenario_truth,
      config = config
    )
  } else if (n_cores <= 1L) {
    out <- lapply(reps, function(r) run_one_replicate(scenario_truth, config, r))
  } else {
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
    out <- parallel::parLapplyLB(
      cl, reps,
      function(r, scenario_truth, config) {
        run_one_replicate(scenario_truth, config, r)
      },
      scenario_truth = scenario_truth,
      config = config
    )
  }
  bind_rows_fill(out)
}

scenario_chunk_dir <- function(config, scenario_id, output_group = "primary") {
  fn <- gsub("[^A-Za-z0-9_.-]", "_", scenario_id)
  d <- file.path(config$output_dir, "raw", output_group, "chunks", fn)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

run_primary_scenario <- function(
    scenario_row,
    config = default_sim_config(),
    nrep = NULL,
    n_cores = 1L,
    project_dir = getwd(),
    write_output = TRUE,
    cluster = NULL,
    chunk_size = config$replicate_chunk_size,
    resume = config$resume_default,
    output_group = "primary") {

  ensure_output_dirs(config)
  truth <- prepare_scenario_truth(scenario_row, config)

  if (is.null(nrep)) {
    nrep <- if (scenario_row$outcome_type == "continuous")
      config$nrep_continuous_primary else config$nrep_binary_primary
  }

  reps_all <- seq_len(as.integer(nrep))
  chunk_size <- max(1L, as.integer(chunk_size))
  chunk_id <- ceiling(reps_all / chunk_size)
  chunks <- split(reps_all, chunk_id)

  cdir <- scenario_chunk_dir(config, scenario_row$scenario_id, output_group)
  parts <- vector("list", length(chunks))

  for (j in seq_along(chunks)) {
    reps <- chunks[[j]]
    chunk_file <- file.path(
      cdir,
      sprintf("rep_%06d_%06d.rds", min(reps), max(reps))
    )

    use_checkpoint <- FALSE
    if (isTRUE(resume) && file.exists(chunk_file)) {
      part <- tryCatch(readRDS(chunk_file), error = function(e) NULL)
      if (!is.null(part)) {
        expected <- reps
        observed <- sort(unique(part$replicate))
        use_checkpoint <- identical(observed, expected)
      }
      if (!use_checkpoint) {
        warning("Invalid checkpoint will be recomputed: ", chunk_file)
        unlink(chunk_file)
      }
    }

    if (!use_checkpoint) {
      part <- run_replicates(
        scenario_truth = truth,
        config = config,
        replicates = reps,
        n_cores = n_cores,
        project_dir = project_dir,
        cluster = cluster
      )
      if (write_output) {
        atomic_save_rds(part, chunk_file)
      }
    }
    parts[[j]] <- part
  }

  raw <- bind_rows_fill(parts)
  raw <- raw[order(raw$replicate), , drop = FALSE]
  if (!identical(sort(unique(raw$replicate)), reps_all)) {
    stop("Combined replicate set is incomplete or duplicated.")
  }

  failure_rate <- mean(!(raw$converged %in% TRUE))
  if (is.finite(failure_rate) && failure_rate > 0.01) {
    warning(sprintf(
      "Scenario %s has %.2f%% failed replicates; investigate before interpretation.",
      scenario_row$scenario_id, 100 * failure_rate
    ))
  }

  if (write_output) {
    fn <- gsub("[^A-Za-z0-9_.-]", "_", scenario_row$scenario_id)
    final_file <- file.path(config$output_dir, "raw", output_group, paste0(fn, ".rds"))
    atomic_save_rds(raw, final_file)
  }

  raw
}
