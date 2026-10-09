# 12_summarise_results.R
# Monte Carlo operating-characteristic summaries.

wilson_ci <- function(x, n, level = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- qnorm(1 - (1 - level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt((p * (1 - p) / n) + z^2 / (4 * n^2)) / den
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

rmse <- function(x, target) sqrt(mean((x - target)^2, na.rm = TRUE))

summarise_scenario <- function(raw, config = default_sim_config()) {
  if (nrow(raw) == 0L) stop("Empty raw data.")
  ok <- raw$converged %in% TRUE
  d <- raw[ok, , drop = FALSE]
  if (nrow(d) == 0L) stop("No converged replicates.")

  target <- unique(d$beta_x_target)
  if (length(target) != 1L) stop("summarise_scenario expects one target.")
  target <- target[1]

  var_u <- var(d$beta_unadj, na.rm = TRUE)
  var_a <- var(d$beta_adj, na.rm = TRUE)
  sd_u <- sd(d$beta_unadj, na.rm = TRUE)
  sd_a <- sd(d$beta_adj, na.rm = TRUE)

  rej_u <- sum(d$reject_unadj == 1, na.rm = TRUE)
  rej_a <- sum(d$reject_adj == 1, na.rm = TRUE)
  n_rej_u <- sum(is.finite(d$reject_unadj))
  n_rej_a <- sum(is.finite(d$reject_adj))
  ci_u <- wilson_ci(rej_u, n_rej_u, config$ci_level)
  ci_a <- wilson_ci(rej_a, n_rej_a, config$ci_level)

  is_null <- abs(target) < 1e-15

  data.frame(
    scenario_id = unique(d$scenario_id)[1],
    role = unique(d$role)[1],
    outcome_type = unique(d$outcome_type)[1],
    sample_design = unique(d$sample_design)[1],
    beta_x_target = target,
    instrument_config = unique(d$instrument_config)[1],
    error_interpretation = if (unique(d$outcome_type)[1] == "binary")
      "target-relative total-logOR error (may include scale-dependent variation)" else
      "continuous-outcome bias",
    n_requested = nrow(raw),
    n_converged = nrow(d),
    failure_rate = 1 - nrow(d) / nrow(raw),
    failure_gt_1pct = (1 - nrow(d) / nrow(raw)) > 0.01,

    mean_beta_unadj = mean(d$beta_unadj, na.rm = TRUE),
    mean_beta_adj = mean(d$beta_adj, na.rm = TRUE),
    mcse_mean_beta_unadj = sd_u / sqrt(nrow(d)),
    mcse_mean_beta_adj = sd_a / sqrt(nrow(d)),

    bias_unadj = mean(d$beta_unadj - target, na.rm = TRUE),
    bias_adj = mean(d$beta_adj - target, na.rm = TRUE),
    mcse_bias_unadj = sd_u / sqrt(nrow(d)),
    mcse_bias_adj = sd_a / sqrt(nrow(d)),

    empirical_var_unadj = var_u,
    empirical_var_adj = var_a,
    empirical_sd_unadj = sd_u,
    empirical_sd_adj = sd_a,

    primary_se_method = unique(config$se_method_primary),
    mean_model_se_unadj = mean(d$se_unadj, na.rm = TRUE),
    mean_model_se_adj = mean(d$se_adj, na.rm = TRUE),
    median_model_se_unadj = median(d$se_unadj, na.rm = TRUE),
    median_model_se_adj = median(d$se_adj, na.rm = TRUE),
    mean_fixed_se_unadj = mean(d$se_unadj_fixed, na.rm = TRUE),
    mean_fixed_se_adj = mean(d$se_adj_fixed, na.rm = TRUE),
    mean_weighted_lm_se_unadj = mean(d$se_unadj_weighted_lm, na.rm = TRUE),
    mean_weighted_lm_se_adj = mean(d$se_adj_weighted_lm, na.rm = TRUE),

    variance_ratio_adj_over_unadj = var_a / var_u,
    relative_efficiency_unadj_over_adj = var_u / var_a,

    rmse_unadj = rmse(d$beta_unadj, target),
    rmse_adj = rmse(d$beta_adj, target),

    coverage_unadj = mean(d$cover_unadj, na.rm = TRUE),
    coverage_adj = mean(d$cover_adj, na.rm = TRUE),

    type1_error_unadj = if (is_null) mean(d$reject_unadj, na.rm = TRUE) else NA_real_,
    type1_error_adj = if (is_null) mean(d$reject_adj, na.rm = TRUE) else NA_real_,
    power_unadj = if (!is_null) mean(d$reject_unadj, na.rm = TRUE) else NA_real_,
    power_adj = if (!is_null) mean(d$reject_adj, na.rm = TRUE) else NA_real_,

    rejection_ci_lower_unadj = ci_u["lower"],
    rejection_ci_upper_unadj = ci_u["upper"],
    rejection_ci_lower_adj = ci_a["lower"],
    rejection_ci_upper_adj = ci_a["upper"],

    mean_F_x_given_z = mean(d$F_x_given_z, na.rm = TRUE),
    mean_F_z_given_x = mean(d$F_z_given_x, na.rm = TRUE),
    median_F_min = median(d$F_min, na.rm = TRUE),
    median_F_min_package = median(d$F_min_package, na.rm = TRUE),

    mean_response = mean(d$response, na.rm = TRUE),
    mean_role_predicted_response = if (all(is.na(d$role_predicted_response)))
      NA_real_ else mean(d$role_predicted_response, na.rm = TRUE),
    mean_role_response_error = if (all(is.na(d$role_response_error)))
      NA_real_ else mean(d$role_response_error, na.rm = TRUE),
    rmse_role_response_error = if (all(is.na(d$role_response_error)))
      NA_real_ else sqrt(mean(d$role_response_error^2, na.rm = TRUE)),

    mean_rho_w = mean(d$rho_w, na.rm = TRUE),
    mean_delta_xz = mean(d$delta_xz, na.rm = TRUE),
    median_condition_number = median(d$condition_number, na.rm = TRUE),
    max_abs_identity_error = max(abs(d$identity_error), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summarise_raw_directory <- function(
    config = default_sim_config(),
    output_group = "primary",
    summary_label = gsub("[^A-Za-z0-9_.-]", "_", output_group)) {

  ensure_output_dirs(config)
  fs <- list.files(
    file.path(config$output_dir, "raw", output_group),
    pattern = "\\.rds$", full.names = TRUE, recursive = FALSE
  )
  if (length(fs) == 0L) stop("No raw RDS files found for group: ", output_group)
  out <- lapply(fs, function(f) summarise_scenario(readRDS(f), config))
  tab <- do.call(rbind, out)
  saveRDS(
    tab,
    file.path(config$output_dir, "summary", paste0(summary_label, "_scenario_summary.rds"))
  )
  write.csv(
    tab,
    file.path(config$output_dir, "summary", paste0(summary_label, "_scenario_summary.csv")),
    row.names = FALSE
  )
  tab
}
