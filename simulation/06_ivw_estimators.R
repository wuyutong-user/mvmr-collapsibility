# 06_ivw_estimators.R
# No-intercept weighted summary regressions.
#
# POINT ESTIMATES exactly match the theorem setup.
#
# For finite-sample inference we retain two uncertainty conventions:
#   1) fixed-effect SE: sqrt(diag((B'WB)^-1))
#   2) weighted-lm SE: fixed-effect SE * residual sigma
#      This is the SE returned by summary(lm(..., weights=...)) and matches
#      the legacy simulation / TwoSampleMR::mv_multiple convention.
#
# The primary finite-sample simulation uses config$se_method_primary
# (default "weighted_lm"), while exact theorem modules can use "fixed".
# Both SEs are saved so empirical precision can be audited independently.

normal_p_two_sided <- function(z) 2 * pnorm(-abs(z))

t_p_two_sided <- function(t, df) {
  if (!is.finite(df) || df <= 0) return(NA_real_)
  2 * pt(-abs(t), df = df)
}

safe_solve <- function(M) {
  tryCatch(solve(M), error = function(e) NULL)
}

select_primary_se <- function(se_fixed, se_weighted_lm,
                              method = c("weighted_lm", "fixed")) {
  method <- match.arg(method)
  if (method == "weighted_lm") se_weighted_lm else se_fixed
}

wls_diagnostics <- function(y, X, w, coef, p) {
  resid <- as.vector(y - X %*% coef)
  Q <- sum(w * resid^2)
  df <- length(y) - p
  sigma <- if (df > 0) sqrt(Q / df) else NA_real_
  list(residual = resid, Q = Q, df = df, sigma = sigma)
}

ivw_unadjusted <- function(summary_data,
                           idx = seq_len(nrow(summary_data)),
                           se_method = c("weighted_lm", "fixed")) {
  se_method <- match.arg(se_method)
  d <- summary_data[idx, , drop = FALSE]
  ok <- with(d, is.finite(beta_x) & is.finite(beta_y) &
                 is.finite(se_y) & se_y > 0)
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 2L) {
    return(list(converged = FALSE, reason = "too_few_instruments"))
  }

  w <- 1 / d$se_y^2
  x <- d$beta_x
  y <- d$beta_y
  den <- sum(w * x^2)
  if (!is.finite(den) || den <= 0) {
    return(list(converged = FALSE, reason = "degenerate_x_direction"))
  }

  beta <- sum(w * x * y) / den
  se_fixed <- sqrt(1 / den)
  diag <- wls_diagnostics(y, matrix(x, ncol = 1L), w, beta, p = 1L)
  se_weighted_lm <- se_fixed * diag$sigma
  se_primary <- select_primary_se(se_fixed, se_weighted_lm, se_method)

  z_normal <- beta / se_primary
  p_normal <- normal_p_two_sided(z_normal)
  p_t <- t_p_two_sided(z_normal, diag$df)

  list(
    converged = TRUE,
    reason = NA_character_,
    beta_x = beta,
    se_x = se_primary,
    se_fixed = se_fixed,
    se_weighted_lm = se_weighted_lm,
    residual_sigma = diag$sigma,
    Q = diag$Q,
    Q_df = diag$df,
    z_x = z_normal,
    p_x = p_normal,
    p_x_normal = p_normal,
    p_x_t = p_t,
    n_instruments = nrow(d),
    snp_id = d$snp_id,
    weight = w
  )
}

ivw_adjusted <- function(summary_data,
                         idx = seq_len(nrow(summary_data)),
                         se_method = c("weighted_lm", "fixed")) {
  se_method <- match.arg(se_method)
  d <- summary_data[idx, , drop = FALSE]
  ok <- with(d, is.finite(beta_x) & is.finite(beta_z) &
                 is.finite(beta_y) & is.finite(se_y) & se_y > 0)
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 3L) {
    return(list(converged = FALSE, reason = "too_few_instruments"))
  }

  w <- 1 / d$se_y^2
  B <- cbind(d$beta_x, d$beta_z)
  M <- crossprod(B, w * B)
  Minv <- safe_solve(M)
  if (is.null(Minv)) {
    return(list(converged = FALSE, reason = "singular_design"))
  }

  coef <- as.vector(Minv %*% crossprod(B, w * d$beta_y))
  se_fixed <- sqrt(diag(Minv))
  diagw <- wls_diagnostics(d$beta_y, B, w, coef, p = 2L)
  se_weighted_lm <- se_fixed * diagw$sigma
  se_primary <- select_primary_se(se_fixed, se_weighted_lm, se_method)

  z_normal <- coef / se_primary
  p_normal <- normal_p_two_sided(z_normal)
  p_t <- vapply(z_normal, t_p_two_sided, numeric(1), df = diagw$df)

  list(
    converged = TRUE,
    reason = NA_character_,
    beta_x = coef[1],
    se_x = se_primary[1],
    se_fixed_x = se_fixed[1],
    se_weighted_lm_x = se_weighted_lm[1],
    z_x = z_normal[1],
    p_x = p_normal[1],
    p_x_normal = p_normal[1],
    p_x_t = p_t[1],
    alpha_z = coef[2],
    se_alpha_z = se_primary[2],
    se_fixed_alpha_z = se_fixed[2],
    se_weighted_lm_alpha_z = se_weighted_lm[2],
    z_alpha_z = z_normal[2],
    p_alpha_z = p_normal[2],
    p_alpha_z_normal = p_normal[2],
    p_alpha_z_t = p_t[2],
    residual_sigma = diagw$sigma,
    Q = diagw$Q,
    Q_df = diagw$df,
    n_instruments = nrow(d),
    snp_id = d$snp_id,
    weight = w,
    condition_number = kappa(M, exact = TRUE)
  )
}

fit_strict_same_set <- function(summary_data,
                                idx = seq_len(nrow(summary_data)),
                                se_method = c("weighted_lm", "fixed")) {
  se_method <- match.arg(se_method)
  u <- ivw_unadjusted(summary_data, idx, se_method = se_method)
  a <- ivw_adjusted(summary_data, idx, se_method = se_method)
  if (!isTRUE(u$converged) || !isTRUE(a$converged)) {
    return(list(
      converged = FALSE,
      reason = paste(na.omit(c(u$reason, a$reason)), collapse = ";")
    ))
  }

  if (!identical(u$snp_id, a$snp_id)) {
    stop("Strict same-set violation: SNP IDs differ.")
  }

  d <- summary_data[match(u$snp_id, summary_data$snp_id), , drop = FALSE]
  w <- 1 / d$se_y^2
  delta <- weighted_projection(d$beta_x, d$beta_z, w)
  rho <- weighted_correlation(d$beta_x, d$beta_z, w)
  response <- u$beta_x - a$beta_x
  predicted <- a$alpha_z * delta
  identity_error <- response - predicted

  list(
    converged = TRUE,
    reason = NA_character_,
    se_method = se_method,

    beta_unadj = u$beta_x,
    se_unadj = u$se_x,
    se_unadj_fixed = u$se_fixed,
    se_unadj_weighted_lm = u$se_weighted_lm,
    p_unadj = u$p_x,
    p_unadj_normal = u$p_x_normal,
    p_unadj_t = u$p_x_t,
    sigma_unadj = u$residual_sigma,
    Q_unadj = u$Q,
    Q_df_unadj = u$Q_df,

    beta_adj = a$beta_x,
    se_adj = a$se_x,
    se_adj_fixed = a$se_fixed_x,
    se_adj_weighted_lm = a$se_weighted_lm_x,
    p_adj = a$p_x,
    p_adj_normal = a$p_x_normal,
    p_adj_t = a$p_x_t,
    sigma_adj = a$residual_sigma,
    Q_adj = a$Q,
    Q_df_adj = a$Q_df,

    alpha_z = a$alpha_z,
    se_alpha_z = a$se_alpha_z,
    se_alpha_z_fixed = a$se_fixed_alpha_z,
    se_alpha_z_weighted_lm = a$se_weighted_lm_alpha_z,

    response = response,
    predicted_response = predicted,
    identity_error = identity_error,
    delta_xz = delta,
    rho_w = rho,
    condition_number = a$condition_number,
    n_instruments = length(u$snp_id),
    snp_id = u$snp_id
  )
}

practical_instrument_sets <- function(summary_data, p_threshold = 5e-8) {
  zx <- summary_data$beta_x / summary_data$se_x
  zz <- summary_data$beta_z / summary_data$se_z
  px <- normal_p_two_sided(zx)
  pz <- normal_p_two_sided(zz)
  list(
    unadjusted = which(is.finite(px) & px < p_threshold),
    adjusted = which((is.finite(px) & px < p_threshold) |
                     (is.finite(pz) & pz < p_threshold))
  )
}

fit_practical_different_set <- function(summary_data,
                                        p_threshold = 5e-8,
                                        se_method = c("weighted_lm", "fixed")) {
  se_method <- match.arg(se_method)
  sets <- practical_instrument_sets(summary_data, p_threshold)
  u <- ivw_unadjusted(summary_data, sets$unadjusted, se_method = se_method)
  a <- ivw_adjusted(summary_data, sets$adjusted, se_method = se_method)
  if (!isTRUE(u$converged) || !isTRUE(a$converged)) {
    return(list(
      converged = FALSE,
      reason = paste(na.omit(c(u$reason, a$reason)), collapse = ";"),
      n_unadj = length(sets$unadjusted),
      n_adj = length(sets$adjusted)
    ))
  }

  list(
    converged = TRUE,
    reason = NA_character_,
    beta_unadj = u$beta_x,
    se_unadj = u$se_x,
    se_unadj_fixed = u$se_fixed,
    se_unadj_weighted_lm = u$se_weighted_lm,
    p_unadj = u$p_x,
    beta_adj = a$beta_x,
    se_adj = a$se_x,
    se_adj_fixed = a$se_fixed_x,
    se_adj_weighted_lm = a$se_weighted_lm_x,
    p_adj = a$p_x,
    alpha_z = a$alpha_z,
    se_alpha_z = a$se_alpha_z,
    response = u$beta_x - a$beta_x,
    n_unadj = length(sets$unadjusted),
    n_adj = length(sets$adjusted)
  )
}

linear_response_functional <- function(x, z, w) {
  B <- cbind(x, z)
  M <- crossprod(B, w * B)
  Minv <- safe_solve(M)
  if (is.null(Minv)) stop("Singular design in response functional.")
  e1 <- c(1, 0)
  a_un <- w * x / sum(w * x^2)
  a_adj <- as.vector((w * B) %*% (Minv %*% e1))
  a_un - a_adj
}
