# 01_role_registry.R
# Frozen Fig.1 six-role registry.
# No simulation function is allowed to infer a role from coefficient movement.

get_role_registry <- function(path_coefficient = 0.30) {
  p <- path_coefficient
  data.frame(
    role = c(
      "confounder",
      "collider",
      "independent_cause",
      "upstream_surrogate_exposure",
      "downstream_surrogate_outcome",
      "downstream_surrogate_exposure"
    ),
    edge_z_to_x = c(p, 0, 0, p, 0, 0),
    edge_z_to_y = c(p, 0, p, 0, 0, 0),
    edge_x_to_z = c(0, p, 0, 0, 0, p),
    edge_y_to_z = c(0, p, 0, 0, p, 0),
    alpha_z_theorem = c(p, 0, p, 0, 0, 0),
    nominal_adjustment_class = c(
      "Necessary",
      "Overadjustment",
      "Target/instrument-set dependent",
      "Overadjustment",
      "Overadjustment",
      "Overadjustment"
    ),
    stringsAsFactors = FALSE
  )
}

get_role_spec <- function(role, registry = get_role_registry()) {
  if (length(role) != 1L || !(role %in% registry$role)) {
    stop("Unknown role: ", paste(role, collapse = ", "))
  }
  as.list(registry[match(role, registry$role), , drop = FALSE])
}

validate_role_registry <- function(registry = get_role_registry()) {
  expected <- c(
    "confounder", "collider", "independent_cause",
    "upstream_surrogate_exposure",
    "downstream_surrogate_outcome",
    "downstream_surrogate_exposure"
  )
  stopifnot(identical(registry$role, expected))
  stopifnot(all(registry$alpha_z_theorem[c(2,4,5,6)] == 0))
  stopifnot(all(registry$alpha_z_theorem[c(1,3)] != 0))
  invisible(TRUE)
}
