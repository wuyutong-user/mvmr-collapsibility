# 14_egger_sensitivity.R
# Optional intercept-model sensitivity.
# This module is deliberately excluded from Theorem-1 exact calibration because
# the theorem is derived for no-intercept weighted summary regression.

egger_unadjusted <- function(summary_data, idx = seq_len(nrow(summary_data))) {
  d <- summary_data[idx, , drop = FALSE]
  ok <- with(d, is.finite(beta_x) & is.finite(beta_y) & is.finite(se_y) & se_y > 0)
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 3L) return(list(converged = FALSE, reason = "too_few_instruments"))

  w <- 1 / d$se_y^2
  B <- cbind(1, d$beta_x)
  M <- crossprod(B, w * B)
  Minv <- safe_solve(M)
  if (is.null(Minv)) return(list(converged = FALSE, reason = "singular_design"))
  coef <- as.vector(Minv %*% crossprod(B, w * d$beta_y))
  se <- sqrt(diag(Minv))
  z <- coef / se

  list(
    converged = TRUE,
    intercept = coef[1],
    se_intercept = se[1],
    p_intercept = normal_p_two_sided(z[1]),
    beta_x = coef[2],
    se_x = se[2],
    p_x = normal_p_two_sided(z[2])
  )
}

egger_adjusted <- function(summary_data, idx = seq_len(nrow(summary_data))) {
  d <- summary_data[idx, , drop = FALSE]
  ok <- with(d, is.finite(beta_x) & is.finite(beta_z) & is.finite(beta_y) &
                 is.finite(se_y) & se_y > 0)
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 4L) return(list(converged = FALSE, reason = "too_few_instruments"))

  w <- 1 / d$se_y^2
  B <- cbind(1, d$beta_x, d$beta_z)
  M <- crossprod(B, w * B)
  Minv <- safe_solve(M)
  if (is.null(Minv)) return(list(converged = FALSE, reason = "singular_design"))
  coef <- as.vector(Minv %*% crossprod(B, w * d$beta_y))
  se <- sqrt(diag(Minv))
  z <- coef / se

  list(
    converged = TRUE,
    intercept = coef[1],
    se_intercept = se[1],
    p_intercept = normal_p_two_sided(z[1]),
    beta_x = coef[2],
    se_x = se[2],
    p_x = normal_p_two_sided(z[2]),
    alpha_z = coef[3],
    se_alpha_z = se[3],
    p_alpha_z = normal_p_two_sided(z[3])
  )
}
