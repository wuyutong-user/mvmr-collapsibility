# 07_conditional_strength.R
# Conditional instrument-strength diagnostics for two-exposure summary-data MVMR.
#
# PRIMARY:
#   F_TS = Q_x / [L - (K - 1)] = Q_x / (L - 1) for K=2 exposures,
#   following Sanderson, Spiller & Bowden (Statistics in Medicine, 2021,
#   doi:10.1002/sim.9133, Eq. 7).
#
# The nuisance coefficient delta is estimated by a no-intercept IVW regression
# of the target SNP-exposure association on the other exposure association,
# which is one of the consistent estimators explicitly allowed in the paper.
#
# SECONDARY:
#   A package-compatible diagnostic is also returned to mirror the public
#   MVMR::strength_mvmr implementation: no-intercept OLS delta and Q/L.
#
# The covariance between SNP-X and SNP-Z association estimates is retained
# because X and Z are estimated in the same exposure GWAS sample.

conditional_q_two_exposure <- function(
    target_beta,
    other_beta,
    target_se,
    other_se,
    cov_target_other,
    min_variance = 1e-14) {

  ok <- is.finite(target_beta) & is.finite(other_beta) &
    is.finite(target_se) & is.finite(other_se) &
    is.finite(cov_target_other)

  y <- target_beta[ok]
  x <- other_beta[ok]
  sy <- target_se[ok]
  sx <- other_se[ok]
  cv <- cov_target_other[ok]

  L <- length(y)
  if (L < 3L) {
    return(list(
      Q = NA_real_, F = NA_real_, delta = NA_real_, L = L,
      Q_package = NA_real_, F_package = NA_real_, delta_package = NA_real_
    ))
  }

  # Paper-compatible nuisance coefficient: no-intercept IVW regression.
  w_delta <- 1 / pmax(sy^2, min_variance)
  den_ivw <- sum(w_delta * x^2)
  delta_ivw <- if (den_ivw > 0) sum(w_delta * x * y) / den_ivw else 0

  v <- sy^2 + delta_ivw^2 * sx^2 - 2 * delta_ivw * cv
  v <- pmax(v, min_variance)
  Q <- sum((y - delta_ivw * x)^2 / v)

  # K=2, so Eq. 7 denominator is L-(K-1)=L-1.
  F <- Q / (L - 1L)

  # Public MVMR::strength_mvmr-compatible calculation for audit/crosswalk:
  # delta from unweighted no-intercept OLS, then Q divided by L.
  den_ols <- sum(x^2)
  delta_ols <- if (den_ols > 0) sum(x * y) / den_ols else 0
  v_pkg <- sy^2 + delta_ols^2 * sx^2 - 2 * delta_ols * cv
  v_pkg <- pmax(v_pkg, min_variance)
  Q_pkg <- sum((y - delta_ols * x)^2 / v_pkg)
  F_pkg <- Q_pkg / L

  list(
    Q = Q,
    F = F,
    delta = delta_ivw,
    L = L,
    Q_package = Q_pkg,
    F_package = F_pkg,
    delta_package = delta_ols
  )
}

conditional_strength <- function(summary_data, idx = seq_len(nrow(summary_data)),
                                 min_variance = 1e-14) {
  d <- summary_data[idx, , drop = FALSE]

  fx <- conditional_q_two_exposure(
    target_beta = d$beta_x,
    other_beta = d$beta_z,
    target_se = d$se_x,
    other_se = d$se_z,
    cov_target_other = d$cov_xz,
    min_variance = min_variance
  )

  fz <- conditional_q_two_exposure(
    target_beta = d$beta_z,
    other_beta = d$beta_x,
    target_se = d$se_z,
    other_se = d$se_x,
    cov_target_other = d$cov_xz,
    min_variance = min_variance
  )

  fvals <- c(fx$F, fz$F)
  fmin <- if (any(is.finite(fvals))) min(fvals[is.finite(fvals)]) else NA_real_

  fvals_pkg <- c(fx$F_package, fz$F_package)
  fmin_pkg <- if (any(is.finite(fvals_pkg))) {
    min(fvals_pkg[is.finite(fvals_pkg)])
  } else NA_real_

  list(
    # Primary Eq. 7 quantities
    F_x_given_z = fx$F,
    F_z_given_x = fz$F,
    F_min = fmin,
    Q_x = fx$Q,
    Q_z = fz$Q,
    delta_x_on_z = fx$delta,
    delta_z_on_x = fz$delta,

    # Package-compatible crosswalk
    F_x_given_z_package = fx$F_package,
    F_z_given_x_package = fz$F_package,
    F_min_package = fmin_pkg,
    Q_x_package = fx$Q_package,
    Q_z_package = fz$Q_package,
    delta_x_on_z_package = fx$delta_package,
    delta_z_on_x_package = fz$delta_package,

    L = min(fx$L, fz$L)
  )
}
