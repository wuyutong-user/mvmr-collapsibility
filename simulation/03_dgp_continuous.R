# 03_dgp_continuous.R
# Continuous-outcome individual-level DGP for the six Fig.1 structures.

genetic_scores <- function(G, architecture) {
  maf <- architecture$maf
  Gc <- sweep(G, 2L, 2 * maf, FUN = "-")
  list(
    gx = as.vector(Gc %*% architecture$beta_x_direct),
    gz = as.vector(Gc %*% architecture$beta_z_direct)
  )
}

generate_continuous_population <- function(
    role,
    N,
    architecture,
    beta_x,
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
  eY <- rnorm(N, sd = residual_sd)

  bzx <- spec$edge_z_to_x
  bzy <- spec$edge_z_to_y
  bxz <- spec$edge_x_to_z
  byz <- spec$edge_y_to_z

  if (role == "confounder") {
    Z <- gs$gz + u_coefficient * U + eZ
    X <- gs$gx + bzx * Z + u_coefficient * U + eX
    Y <- beta_x * X + bzy * Z + u_coefficient * U + eY

  } else if (role == "collider") {
    X <- gs$gx + u_coefficient * U + eX
    Y <- beta_x * X + u_coefficient * U + eY
    Z <- gs$gz + bxz * X + byz * Y + u_coefficient * U + eZ

  } else if (role == "independent_cause") {
    X <- gs$gx + u_coefficient * U + eX
    Z <- gs$gz + u_coefficient * U + eZ
    Y <- beta_x * X + bzy * Z + u_coefficient * U + eY

  } else if (role == "upstream_surrogate_exposure") {
    Z <- gs$gz + u_coefficient * U + eZ
    X <- gs$gx + bzx * Z + u_coefficient * U + eX
    Y <- beta_x * X + u_coefficient * U + eY

  } else if (role == "downstream_surrogate_outcome") {
    X <- gs$gx + u_coefficient * U + eX
    Y <- beta_x * X + u_coefficient * U + eY
    Z <- gs$gz + byz * Y + u_coefficient * U + eZ

  } else if (role == "downstream_surrogate_exposure") {
    X <- gs$gx + u_coefficient * U + eX
    Y <- beta_x * X + u_coefficient * U + eY
    Z <- gs$gz + bxz * X + u_coefficient * U + eZ

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
    outcome_type = "continuous",
    role = role,
    beta_x_target = beta_x
  )
}
