# 05_gwas_summary.R
# Vectorised SNP-trait summary association estimation.

center_genotypes <- function(G) {
  sweep(G, 2L, colMeans(G), FUN = "-")
}

gwas_continuous <- function(G, y, min_variance = 1e-14) {
  N <- nrow(G)
  Gc <- center_genotypes(G)
  yc <- y - mean(y)
  sxx <- colSums(Gc^2)
  beta <- as.vector(crossprod(Gc, yc)) / sxx
  syy <- sum(yc^2)
  rss <- pmax(syy - beta^2 * sxx, 0)
  sigma2 <- rss / pmax(N - 2L, 1L)
  se <- sqrt(pmax(sigma2 / sxx, min_variance))
  list(beta = beta, se = se)
}

gwas_continuous_pair <- function(G, x, z, min_variance = 1e-14) {
  N <- nrow(G)
  Gc <- center_genotypes(G)
  xc <- x - mean(x)
  zc <- z - mean(z)
  sxx_g <- colSums(Gc^2)

  bx <- as.vector(crossprod(Gc, xc)) / sxx_g
  bz <- as.vector(crossprod(Gc, zc)) / sxx_g

  sx <- sum(xc^2)
  sz <- sum(zc^2)
  rss_x <- pmax(sx - bx^2 * sxx_g, 0)
  rss_z <- pmax(sz - bz^2 * sxx_g, 0)

  var_x <- pmax((rss_x / pmax(N - 2L, 1L)) / sxx_g, min_variance)
  var_z <- pmax((rss_z / pmax(N - 2L, 1L)) / sxx_g, min_variance)

  # Sampling covariance of the two SNP association estimates when X and Z
  # are measured in the same exposure GWAS sample.
  xz <- sum(xc * zc)
  resid_cross <- xz - bx * bz * sxx_g
  cov_xz <- (resid_cross / pmax(N - 2L, 1L)) / sxx_g

  list(
    beta_x = bx,
    se_x = sqrt(var_x),
    beta_z = bz,
    se_z = sqrt(var_z),
    cov_xz = cov_xz
  )
}

gwas_binary <- function(G, y, min_variance = 1e-14) {
  K <- ncol(G)
  beta <- se <- rep(NA_real_, K)
  converged <- rep(FALSE, K)

  for (j in seq_len(K)) {
    X <- cbind(1, G[, j])
    fit <- tryCatch(
      suppressWarnings(glm.fit(x = X, y = y, family = binomial())),
      error = function(e) NULL
    )
    if (is.null(fit) || !isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) next

    R <- tryCatch(qr.R(fit$qr), error = function(e) NULL)
    if (is.null(R) || nrow(R) < 2L || ncol(R) < 2L) next
    V <- tryCatch(chol2inv(R), error = function(e) NULL)
    if (is.null(V) || !is.finite(V[2,2]) || V[2,2] <= 0) next

    beta[j] <- fit$coefficients[2]
    se[j] <- sqrt(max(V[2,2], min_variance))
    converged[j] <- TRUE
  }

  list(beta = beta, se = se, converged = converged)
}

make_summary_data <- function(
    exposure_population,
    outcome_population,
    outcome_type = c("continuous", "binary"),
    min_variance = 1e-14) {

  outcome_type <- match.arg(outcome_type)
  if (ncol(exposure_population$G) != ncol(outcome_population$G)) {
    stop("Exposure and outcome populations use different K.")
  }

  ex <- gwas_continuous_pair(
    exposure_population$G,
    exposure_population$X,
    exposure_population$Z,
    min_variance = min_variance
  )

  if (outcome_type == "continuous") {
    oy <- gwas_continuous(
      outcome_population$G,
      outcome_population$Y,
      min_variance = min_variance
    )
    cy <- rep(TRUE, length(oy$beta))
  } else {
    oy <- gwas_binary(
      outcome_population$G,
      outcome_population$Y,
      min_variance = min_variance
    )
    cy <- oy$converged
  }

  K <- length(ex$beta_x)
  data.frame(
    snp_id = sprintf("snp_%03d", seq_len(K)),
    beta_x = ex$beta_x,
    se_x = ex$se_x,
    beta_z = ex$beta_z,
    se_z = ex$se_z,
    cov_xz = ex$cov_xz,
    beta_y = oy$beta,
    se_y = oy$se,
    outcome_gwas_converged = cy,
    stringsAsFactors = FALSE
  )
}
