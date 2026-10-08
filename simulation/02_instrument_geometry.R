# 02_instrument_geometry.R
# Genetic architecture and theorem-space weighted geometry.
#
# IMPORTANT TERMINOLOGY
#   direct_class  : biological/direct genetic source in individual-level DGP
#   summary_class : population summary-association class used by the theorems
#                   X_only / Z_only / shared
# These labels are intentionally separate.

scale_rms <- function(x, target_rms) {
  nz <- abs(x) > 0
  if (!any(nz)) return(x)
  current <- sqrt(mean(x[nz]^2))
  if (!is.finite(current) || current <= 0) stop("Cannot scale zero vector.")
  x * (target_rms / current)
}

weighted_projection <- function(x, z, w) {
  stopifnot(length(x) == length(z), length(x) == length(w))
  den <- sum(w * x * x)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(w * x * z) / den
}

weighted_correlation <- function(x, z, w) {
  nx <- sum(w * x * x)
  nz <- sum(w * z * z)
  if (nx <= 0 || nz <= 0) return(NA_real_)
  sum(w * x * z) / sqrt(nx * nz)
}

make_genotype_matrix <- function(N, maf, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(maf) == 1L) maf <- rep(maf, K)
  stopifnot(length(maf) == K, all(maf > 0), all(maf < 1))
  G <- matrix(
    rbinom(N * K, size = 2L, prob = rep(maf, each = N)),
    nrow = N, ncol = K
  )
  colnames(G) <- sprintf("snp_%03d", seq_len(K))
  G
}

make_summary_geometry <- function(
    K,
    n_x_only,
    n_z_only,
    n_shared,
    rho_target,
    w = rep(1, K),
    effect_rms = 0.08,
    seed = 1L,
    tol = 1e-10) {

  stopifnot(n_x_only + n_z_only + n_shared == K)
  if (n_shared == 0L && abs(rho_target) > tol) {
    stop("Non-zero rho_target is impossible when n_shared == 0.")
  }
  if (abs(rho_target) >= 1) stop("|rho_target| must be < 1.")
  if (any(w <= 0) || length(w) != K) stop("w must be positive and length K.")

  set.seed(seed)
  perm <- sample.int(K)
  idx_x <- sort(perm[seq_len(n_x_only)])
  pos <- n_x_only
  idx_z <- if (n_z_only > 0L) sort(perm[pos + seq_len(n_z_only)]) else integer(0)
  pos <- pos + n_z_only
  idx_s <- if (n_shared > 0L) sort(perm[pos + seq_len(n_shared)]) else integer(0)

  # Construct in transformed coordinates xt = sqrt(W) x, zt = sqrt(W) z.
  sw <- sqrt(w)
  xt <- numeric(K)

  if (n_shared > 0L) {
    xt[idx_s] <- abs(rnorm(n_shared)) + 0.5
  }
  if (n_x_only > 0L) {
    # Keep X-only nonzero while allowing high full-vector weighted correlation.
    xt[idx_x] <- 0.05 * (abs(rnorm(n_x_only)) + 0.5)
  }

  zt <- numeric(K)

  if (n_shared == 0L) {
    if (n_z_only == 0L) stop("Z vector is identically zero.")
    zt[idx_z] <- rnorm(n_z_only)
  } else if (abs(rho_target) <= tol) {
    xs <- xt[idx_s]
    q <- rnorm(n_shared)
    q <- q - xs * (sum(q * xs) / sum(xs^2))
    if (sqrt(sum(q^2)) < 1e-12) {
      q <- seq_len(n_shared)
      q <- q - xs * (sum(q * xs) / sum(xs^2))
    }
    q <- q / sqrt(sum(q^2))
    zt[idx_s] <- q
    if (n_z_only > 0L) zt[idx_z] <- rnorm(n_z_only)
  } else {
    xs <- xt[idx_s]
    c_shared <- sum(xs^2)
    nx <- sum(xt^2)
    max_abs_rho <- sqrt(c_shared / nx)
    if (abs(rho_target) >= max_abs_rho - 1e-12) {
      stop(sprintf(
        "Requested |rho|=%.4f exceeds support-constrained maximum %.4f.",
        abs(rho_target), max_abs_rho
      ))
    }

    sgn <- sign(rho_target)
    zt[idx_s] <- sgn * xs
    # Required transformed Z norm for the exact full-vector rho.
    nz_target <- c_shared^2 / (rho_target^2 * nx)
    extra_norm2 <- nz_target - c_shared
    if (extra_norm2 < -1e-10) stop("Internal geometry feasibility error.")
    extra_norm2 <- max(extra_norm2, 0)

    q <- rnorm(n_shared)
    q <- q - xs * (sum(q * xs) / c_shared)
    nq <- sqrt(sum(q^2))
    can_q <- is.finite(nq) && nq > 1e-12

    if (n_z_only > 0L && can_q) {
      frac_zonly <- 0.50
    } else if (n_z_only > 0L) {
      frac_zonly <- 1.00
    } else {
      frac_zonly <- 0.00
    }

    if (can_q && extra_norm2 > 0) {
      zt[idx_s] <- zt[idx_s] +
        q / nq * sqrt(extra_norm2 * (1 - frac_zonly))
    }
    if (n_z_only > 0L && extra_norm2 > 0) {
      rz <- rnorm(n_z_only)
      rz <- rz / sqrt(sum(rz^2))
      zt[idx_z] <- rz * sqrt(extra_norm2 * frac_zonly)
    }
  }

  x <- xt / sw
  z <- zt / sw
  x <- scale_rms(x, effect_rms)
  z <- scale_rms(z, effect_rms)

  rho_achieved <- weighted_correlation(x, z, w)
  delta_achieved <- weighted_projection(x, z, w)

  if (n_shared > 0L && abs(rho_achieved - rho_target) > 1e-8) {
    stop(sprintf("Failed exact rho construction: target %.8f, achieved %.8f",
                 rho_target, rho_achieved))
  }
  if (n_shared == 0L && abs(rho_achieved) > 1e-12) {
    stop("No-shared geometry did not produce weighted orthogonality.")
  }

  cls <- rep(NA_character_, K)
  cls[idx_x] <- "X_only"
  cls[idx_z] <- "Z_only"
  cls[idx_s] <- "shared"

  d <- data.frame(
    snp_id = sprintf("snp_%03d", seq_len(K)),
    summary_class = cls,
    gamma_x_true = x,
    gamma_z_true = z,
    weight = w,
    stringsAsFactors = FALSE
  )
  attr(d, "rho_w_true") <- rho_achieved
  attr(d, "delta_xz_true") <- delta_achieved
  d
}

make_direct_architecture <- function(
    K,
    n_direct_x,
    n_direct_z,
    n_direct_shared,
    shared_rho = 0.80,
    effect_rms = 0.08,
    maf = 0.30,
    seed = 1L) {

  stopifnot(n_direct_x + n_direct_z + n_direct_shared == K)
  if (n_direct_shared == 0L && abs(shared_rho) > 1e-12) {
    shared_rho <- 0
  }
  if (abs(shared_rho) >= 1) stop("|shared_rho| must be < 1.")

  set.seed(seed)
  perm <- sample.int(K)
  ix <- sort(perm[seq_len(n_direct_x)])
  pos <- n_direct_x
  iz <- if (n_direct_z > 0L) sort(perm[pos + seq_len(n_direct_z)]) else integer(0)
  pos <- pos + n_direct_z
  is <- if (n_direct_shared > 0L) sort(perm[pos + seq_len(n_direct_shared)]) else integer(0)

  bx <- numeric(K)
  bz <- numeric(K)

  if (n_direct_x > 0L) bx[ix] <- rnorm(n_direct_x)
  if (n_direct_z > 0L) bz[iz] <- rnorm(n_direct_z)

  if (n_direct_shared > 0L) {
    a <- rnorm(n_direct_shared)
    a <- a / sqrt(sum(a^2))
    if (abs(shared_rho) <= 1e-12) {
      b <- rnorm(n_direct_shared)
      b <- b - a * sum(a * b)
      b <- b / sqrt(sum(b^2))
    } else {
      q <- rnorm(n_direct_shared)
      q <- q - a * sum(a * q)
      q <- q / sqrt(sum(q^2))
      b <- shared_rho * a + sqrt(1 - shared_rho^2) * q
    }
    bx[is] <- a
    bz[is] <- b
  }

  bx <- scale_rms(bx, effect_rms)
  bz <- scale_rms(bz, effect_rms)

  cls <- rep(NA_character_, K)
  cls[ix] <- "direct_X"
  cls[iz] <- "direct_Z"
  cls[is] <- "direct_shared"

  if (length(maf) == 1L) maf <- rep(maf, K)

  data.frame(
    snp_id = sprintf("snp_%03d", seq_len(K)),
    direct_class = cls,
    beta_x_direct = bx,
    beta_z_direct = bz,
    maf = maf,
    stringsAsFactors = FALSE
  )
}


# Population SNP-association truth for X and Z under the continuous structural
# equations. This is used to keep "direct genetic source" distinct from the
# theorem's population summary-association classes.
population_xz_summary_truth <- function(role, architecture, beta_x,
                                        role_registry = get_role_registry()) {
  spec <- get_role_spec(role, role_registry)
  bx <- architecture$beta_x_direct
  bz <- architecture$beta_z_direct

  if (role %in% c("confounder", "upstream_surrogate_exposure")) {
    gamma_z <- bz
    gamma_x <- bx + spec$edge_z_to_x * gamma_z

  } else if (role == "independent_cause") {
    gamma_x <- bx
    gamma_z <- bz

  } else if (role == "collider") {
    gamma_x <- bx
    gamma_y <- beta_x * gamma_x
    gamma_z <- bz + spec$edge_x_to_z * gamma_x +
      spec$edge_y_to_z * gamma_y

  } else if (role == "downstream_surrogate_outcome") {
    gamma_x <- bx
    gamma_y <- beta_x * gamma_x
    gamma_z <- bz + spec$edge_y_to_z * gamma_y

  } else if (role == "downstream_surrogate_exposure") {
    gamma_x <- bx
    gamma_z <- bz + spec$edge_x_to_z * gamma_x

  } else {
    stop("Unknown role: ", role)
  }

  data.frame(
    snp_id = architecture$snp_id,
    gamma_x_true = gamma_x,
    gamma_z_true = gamma_z,
    stringsAsFactors = FALSE
  )
}

classify_summary_support <- function(gamma_x, gamma_z, tol = 1e-12) {
  xnz <- abs(gamma_x) > tol
  znz <- abs(gamma_z) > tol
  out <- rep("null_both", length(gamma_x))
  out[xnz & !znz] <- "X_only"
  out[!xnz & znz] <- "Z_only"
  out[xnz & znz] <- "shared"
  out
}
