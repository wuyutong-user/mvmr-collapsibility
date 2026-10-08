# 08_theorem_calibration.R
# Exact same-set, summary-association-level calibration of Theorem 1.

make_theorem_summary_data <- function(geometry, beta_x_target, alpha_z) {
  x <- geometry$gamma_x_true
  z <- geometry$gamma_z_true
  w <- geometry$weight
  y <- beta_x_target * x + alpha_z * z

  data.frame(
    snp_id = geometry$snp_id,
    beta_x = x,
    se_x = rep(1e-6, length(x)),
    beta_z = z,
    se_z = rep(1e-6, length(x)),
    cov_xz = rep(0, length(x)),
    beta_y = y,
    se_y = 1 / sqrt(w),
    outcome_gwas_converged = TRUE,
    stringsAsFactors = FALSE
  )
}

run_theorem1_calibration <- function(config = default_sim_config(), write_output = TRUE) {
  ensure_output_dirs(config)
  registry <- get_role_registry(config$path_coefficient)

  cfg <- config$theorem_instrument_config
  results <- list()
  rr <- 0L

  for (i in seq_len(nrow(cfg))) {
    cc <- cfg[i, , drop = FALSE]
    rhos <- if (cc$n_shared == 0L) {
      0
    } else {
      unique(c(
        config$theorem_rho_primary,
        config$theorem_rho_sensitivity
      ))
    }

    for (rho in rhos) {
      for (role in registry$role) {
        rr <- rr + 1L
        alpha <- registry$alpha_z_theorem[registry$role == role]

        w <- seq(0.75, 1.25, length.out = config$K)
        geom <- make_summary_geometry(
          K = config$K,
          n_x_only = cc$n_x_only,
          n_z_only = cc$n_z_only,
          n_shared = cc$n_shared,
          rho_target = rho,
          w = w,
          effect_rms = config$theorem_effect_rms,
          seed = stable_seed(config$master_seed, "theorem1", cc$instrument_config, rho)
        )
        sd <- make_theorem_summary_data(
          geom,
          beta_x_target = 0.10,
          alpha_z = alpha
        )
        fit <- fit_strict_same_set(sd, se_method = "fixed")

        if (!isTRUE(fit$converged)) stop("Theorem calibration fit failed.")
        if (abs(fit$identity_error) >= config$theorem_tolerance) {
          stop("Theorem 1 identity test failed.")
        }

        results[[rr]] <- data.frame(
          simulation_module = "theorem1_exact",
          role = role,
          instrument_config = cc$instrument_config,
          n_x_only = cc$n_x_only,
          n_z_only = cc$n_z_only,
          n_shared = cc$n_shared,
          rho_target = rho,
          rho_w = fit$rho_w,
          delta_xz = fit$delta_xz,
          alpha_z_true = alpha,
          alpha_z_fit = fit$alpha_z,
          beta_x_target = 0.10,
          beta_unadj = fit$beta_unadj,
          beta_adj = fit$beta_adj,
          response = fit$response,
          predicted_response = fit$predicted_response,
          identity_error = fit$identity_error,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, results)
  if (write_output) {
    saveRDS(out, file.path(config$output_dir, "theorem", "theorem1_exact_calibration.rds"))
    write.csv(out,
              file.path(config$output_dir, "theorem", "theorem1_exact_calibration.csv"),
              row.names = FALSE)
  }
  out
}
