#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else normalizePath(sys.frames()[[1]]$ofile)
test_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
code_dir <- normalizePath(file.path(test_dir, ".."), mustWork = FALSE)

source(file.path(code_dir, "13_eraca_config_controlled_data.R"))
source(file.path(code_dir, "helpers_eraca.R"))

assert_true <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS:", label, "\n")
}

assert_equal <- function(value, expected, label) {
  if (!identical(value, expected)) {
    stop(
      "FAIL: ", label, "\nExpected: ", paste(expected, collapse = ", "),
      "\nObserved: ", paste(value, collapse = ", "),
      call. = FALSE
    )
  }
  cat("PASS:", label, "\n")
}

cfg <- eraca_config()
validate_eraca_config(cfg, check_files = FALSE)

assert_equal(length(cfg$candidates$name), 12L, "configuration contains 12 candidate traits")
assert_equal(length(cfg$outcomes$name), 14L, "configuration contains 14 outcomes")

tasks <- build_direction_tasks(cfg)
assert_equal(nrow(tasks), 388L, "direction registry contains 388 tasks")
assert_equal(length(unique(tasks$task_id)), 388L, "direction task IDs are unique")
assert_equal(sum(tasks$family == "bmi_candidate"), 24L, "BMI-candidate directions total 24")
assert_equal(sum(tasks$family == "candidate_outcome"), 336L, "candidate-outcome directions total 336")
assert_equal(sum(tasks$family == "bmi_outcome"), 28L, "BMI-outcome directions total 28")

triads <- build_triad_registry(cfg)
assert_equal(nrow(triads), 168L, "triad registry contains 168 combinations")
assert_equal(length(unique(triads$triad_id)), 168L, "triad IDs are unique")

N <- "not_supported"
S <- "robust_supported"
edges <- function(X_to_Z = N, Z_to_X = N, Z_to_Y = N, Y_to_Z = N) {
  c(X_to_Z = X_to_Z, Z_to_X = Z_to_X, Z_to_Y = Z_to_Y, Y_to_Z = Y_to_Z)
}

assert_equal(assign_role_from_edges(edges(Z_to_X = S, Z_to_Y = S))$role,
             "confounder", "confounder template")
assert_equal(assign_role_from_edges(edges(X_to_Z = S, Y_to_Z = S))$role,
             "collider", "collider template")
assert_equal(assign_role_from_edges(edges(Z_to_Y = S))$role,
             "independent_cause", "independent-cause template")
assert_equal(assign_role_from_edges(edges(Z_to_X = S))$role,
             "upstream_surrogate_of_exposure", "upstream-surrogate template")
assert_equal(assign_role_from_edges(edges(Y_to_Z = S))$role,
             "downstream_surrogate_of_outcome", "outcome-surrogate template")
assert_equal(assign_role_from_edges(edges(X_to_Z = S))$role,
             "downstream_surrogate_of_exposure", "exposure-surrogate template")
assert_equal(assign_role_from_edges(edges(X_to_Z = S, Z_to_Y = S))$role,
             "mediator_out_of_library", "mediator exits the six-role library")
assert_equal(assign_role_from_edges(edges(X_to_Z = S, Z_to_X = S))$role,
             "feedback_or_cycle", "bidirectional pair exits as cycle")
assert_equal(assign_role_from_edges(edges(), overlapping_measure = TRUE)$role,
             "overlapping_measure", "overlapping adiposity measure exit")
assert_equal(assign_role_from_edges(edges())$role,
             "insufficient_evidence", "empty graph exits as insufficient evidence")

robust <- classify_edge_evidence(
  pval = 1e-6, p_fdr = 0.001, b = 0.2, nsnp = 12,
  mean_f = 30, secondary_concordant = TRUE, unstable = FALSE
)
assert_equal(robust, "robust_supported", "robust edge evidence tier")

null_edge <- classify_edge_evidence(
  pval = 0.6, p_fdr = 0.8, b = 0.01, nsnp = 12,
  mean_f = 30, secondary_concordant = FALSE, unstable = FALSE
)
assert_equal(null_edge, "not_supported", "null MR result is recorded as not supported")

x <- c(0.1, 0.2, 0.3)
z <- c(0.2, 0.1, 0.4)
beta_x <- 0.8
alpha_z <- 0.5
y <- beta_x * x + alpha_z * z
identity <- compute_weighted_response_identity(x, z, y, rep(1, 3), beta_x, alpha_z)
assert_true(abs(identity$identity_error) < 1e-12, "weighted coefficient-response identity")
assert_true(abs(identity$observed_response - identity$predicted_response) < 1e-12,
            "observed response equals alpha_Z times delta_XZ")

cat("All ERACA synthetic-triad tests passed.\n")
