# 04_dgp_binary.R
# Binary-outcome DGP. The common U term explicitly enters the logistic
# outcome model, correcting the omission in the legacy binary code.

binary_pre_y_variables <- function(
    role, N, architecture, beta_x, role_registry,
    u_coefficient, residual_sd, seed) {

  spec <- get_role_spec(role, role_registry)
  set.seed(seed)
  K <- nrow(architecture)
  G <- make_genotype_matrix(N, architecture$maf, K)
  gs <- genetic_scores(G, architecture)

  U <- rnorm(N)
  eX <- rnorm(N, sd = residual_sd)
  eZ <- rnorm(N, sd = residual_sd)

  bzx <- spec$edge_z_to_x
  bzy <- spec$edge_z_to_y

  if (role %in% c("confounder", "upstream_surrogate_exposure")) {
    Z <- gs$gz + u_coefficient * U + eZ
    X <- gs$gx + bzx * Z + u_coefficient * U + eX
  } else if (role == "independent_cause") {
    X <- gs$gx + u_coefficient * U + eX
    Z <- gs$gz + u_coefficient * U + eZ
  } else {
    X <- gs$gx + u_coefficient * U + eX
    Z <- NULL
  }

  lp_no_intercept <- beta_x * X + u_coefficient * U
  if (role %in% c("confounder", "independent_cause")) {
    lp_no_intercept <- lp_no_intercept + bzy * Z
  }

  list(lp = lp_no_intercept)
}

calibrate_binary_intercept <- function(
    role,
    architecture,
    beta_x,
    target_prevalence = 0.10,
    role_registry = get_role_registry(),
    u_coefficient = 0.30,
    residual_sd = 1,
    n_reference = 100000L,
    seed = 1L) {

  ref <- binary_pre_y_variables(
    role = role, N = n_reference, architecture = architecture,
    beta_x = beta_x, role_registry = role_registry,
    u_coefficient = u_coefficient, residual_sd = residual_sd, seed = seed
  )

  f <- function(a) mean(plogis(a + ref$lp)) - target_prevalence
  eta0 <- uniroot(f, interval = c(-20, 20), tol = 1e-10)$root
  achieved <- mean(plogis(eta0 + ref$lp))
  list(intercept = eta0, calibrated_prevalence = achieved)
}

generate_binary_population <- function(
    role,
    N,
    architecture,
    beta_x,
    intercept,
    role_registry = get_role_registry(),
    u_coefficient = 0.30,
    residual_sd = 1,
    seed = 1L,
    id_prefix = "S") {

  spec <- get_role_spec(role, role_registry)
  set.seed(seed)
  K <- nrow(architecture)
  G <- make_genotype_matrix(N, architecture$maf, K)
  gs <- genetic_scores(G, architecture)

  U <- rnorm(N)
  eX <- rnorm(N, sd = residual_sd)
  eZ <- rnorm(N, sd = residual_sd)

  bzx <- spec$edge_z_to_x
  bzy <- spec$edge_z_to_y
  bxz <- spec$edge_x_to_z
  byz <- spec$edge_y_to_z

  if (role == "confounder") {
    Z <- gs$gz + u_coefficient * U + eZ
    X <- gs$gx + bzx * Z + u_coefficient * U + eX
    lp <- intercept + beta_x * X + bzy * Z + u_coefficient * U
    Y <- rbinom(N, 1L, plogis(lp))

  } else if (role == "collider") {
    X <- gs$gx + u_coefficient * U + eX
    lp <- intercept + beta_x * X + u_coefficient * U
    Y <- rbinom(N, 1L, plogis(lp))
    Z <- gs$gz + bxz * X + byz * Y + u_coefficient * U + eZ

  } else if (role == "independent_cause") {
    X <- gs$gx + u_coefficient * U + eX
    Z <- gs$gz + u_coefficient * U + eZ
    lp <- intercept + beta_x * X + bzy * Z + u_coefficient * U
    Y <- rbinom(N, 1L, plogis(lp))

  } else if (role == "upstream_surrogate_exposure") {
    Z <- gs$gz + u_coefficient * U + eZ
    X <- gs$gx + bzx * Z + u_coefficient * U + eX
    lp <- intercept + beta_x * X + u_coefficient * U
    Y <- rbinom(N, 1L, plogis(lp))

  } else if (role == "downstream_surrogate_outcome") {
    X <- gs$gx + u_coefficient * U + eX
    lp <- intercept + beta_x * X + u_coefficient * U
    Y <- rbinom(N, 1L, plogis(lp))
    Z <- gs$gz + byz * Y + u_coefficient * U + eZ

  } else if (role == "downstream_surrogate_exposure") {
    X <- gs$gx + u_coefficient * U + eX
    Z <- gs$gz + bxz * X + u_coefficient * U + eZ
    lp <- intercept + beta_x * X + u_coefficient * U
    Y <- rbinom(N, 1L, plogis(lp))

  } else {
    stop("Unknown role: ", role)
  }

  list(
    id = sprintf("%s%07d", id_prefix, seq_len(N)),
    G = G,
    X = X,
    Z = Z,
    Y = Y,
    U = U,
    outcome_type = "binary",
    role = role,
    beta_x_target = beta_x,
    intercept = intercept,
    prevalence = mean(Y)
  )
}
