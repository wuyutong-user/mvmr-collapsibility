# 09_binary_scale_calibration.R
# Exact summary-level calibration of Theorems 2 and 3.

construct_scale_residual <- function(x, z, w, target_delta_scale) {
  l <- linear_response_functional(x, z, w)
  den <- sum(l^2)
  if (den <= 1e-20) {
    if (abs(target_delta_scale) <= 1e-14) return(rep(0, length(x)))
    stop("Cannot create non-zero scale residual: response functional is zero.")
  }
  target_delta_scale * l / den
}

run_binary_scale_calibration <- function(config = default_sim_config(), write_output = TRUE) {
  ensure_output_dirs(config)

  # Use an all-class positive-alignment geometry to make structural response non-zero.
  cc <- config$theorem_instrument_config[
    config$theorem_instrument_config$instrument_config == "all_summary_classes", ,
    drop = FALSE
  ]
  w <- seq(0.75, 1.25, length.out = config$K)
  geom <- make_summary_geometry(
    K = config$K,
    n_x_only = cc$n_x_only,
    n_z_only = cc$n_z_only,
    n_shared = cc$n_shared,
    rho_target = config$theorem_rho_primary,
    w = w,
    effect_rms = config$theorem_effect_rms,
    seed = stable_seed(config$master_seed, "binary_scale_geometry")
  )

  x <- geom$gamma_x_true
  z <- geom$gamma_z_true
  delta <- weighted_projection(x, z, w)
  S_grid <- c(0.75, 1.00, 1.50)
  modes <- c("none", "reinforce", "partial_offset", "exact_offset", "scale_only")
  rows <- list()
  k <- 0L

  for (S in S_grid) {
    for (mode in modes) {
      k <- k + 1L

      if (mode == "scale_only") {
        alpha_star <- 0
        beta_star <- 0.10
        delta_struct <- 0
        target_delta_scale <- 0.05
      } else {
        alpha_star <- 0.30
        beta_star <- 0.10
        delta_struct <- (alpha_star / S) * delta
        target_delta_scale <- switch(
          mode,
          none = 0,
          reinforce = 0.50 * delta_struct,
          partial_offset = -0.50 * delta_struct,
          exact_offset = -1.00 * delta_struct
        )
      }

      r_scale <- construct_scale_residual(x, z, w, target_delta_scale)
      ell <- (beta_star * x + alpha_star * z) / S + r_scale

      sd <- data.frame(
        snp_id = geom$snp_id,
        beta_x = x,
        se_x = rep(1e-6, length(x)),
        beta_z = z,
        se_z = rep(1e-6, length(x)),
        cov_xz = 0,
        beta_y = ell,
        se_y = 1 / sqrt(w),
        outcome_gwas_converged = TRUE,
        stringsAsFactors = FALSE
      )

      fit <- fit_strict_same_set(sd, se_method = "fixed")
      lfun <- linear_response_functional(x, z, w)
      delta_scale_observed <- sum(lfun * r_scale)
      delta_total <- fit$response
      decomposition_error <- delta_total - (delta_struct + delta_scale_observed)

      if (abs(decomposition_error) >= config$theorem_tolerance) {
        stop("Binary structural+scale decomposition test failed.")
      }

      rows[[k]] <- data.frame(
        simulation_module = "binary_common_scale_exact",
        S = S,
        residual_mode = mode,
        beta_x_star = beta_star,
        alpha_z_star = alpha_star,
        delta_xz = delta,
        delta_struct = delta_struct,
        delta_scale_target = target_delta_scale,
        delta_scale_observed = delta_scale_observed,
        delta_total = delta_total,
        decomposition_error = decomposition_error,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)

  if (write_output) {
    saveRDS(out, file.path(config$output_dir, "theorem", "binary_scale_calibration.rds"))
    write.csv(out,
              file.path(config$output_dir, "theorem", "binary_scale_calibration.csv"),
              row.names = FALSE)
  }
  out
}
