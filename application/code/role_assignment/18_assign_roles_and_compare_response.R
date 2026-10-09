#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "13_eraca_config_controlled_data.R"))
source(file.path(script_dir, "helpers_eraca.R"))

cfg <- eraca_config()
ensure_eraca_directories(cfg)

consensus_path <- file.path(cfg$paths$graphs, "conditional_graph_consensus.csv")
marginal_path <- file.path(cfg$paths$graphs, "marginal_edge_ledger.csv")
response_path <- file.path(cfg$paths$graphs, "same_set_response_decomposition.csv")
for (path in c(consensus_path, marginal_path, response_path)) {
  if (!file.exists(path)) eraca_stop("Required role-assignment input not found: ", path)
}

conditional <- utils::read.csv(consensus_path, stringsAsFactors = FALSE, check.names = FALSE)
marginal <- utils::read.csv(marginal_path, stringsAsFactors = FALSE, check.names = FALSE)
response <- utils::read.csv(response_path, stringsAsFactors = FALSE, check.names = FALSE)
triads <- build_triad_registry(cfg)
allow_partial <- eraca_env_flag("ERACA_ALLOW_PARTIAL", FALSE)

available_triads <- unique(conditional$triad_id)
if (!allow_partial && length(available_triads) != 168L) {
  eraca_stop(
    "Role assignment requires 168 completed triads; found ", length(available_triads),
    ". Finish conditional-learning checkpoints or set ERACA_ALLOW_PARTIAL=true for a diagnostic run."
  )
}
triads <- triads[triads$triad_id %in% available_triads, , drop = FALSE]

edge_label <- function(triad_id, from, to) {
  c_row <- conditional[
    conditional$triad_id == triad_id & conditional$from == from & conditional$to == to,
    , drop = FALSE
  ]
  m_row <- marginal[
    marginal$triad_id == triad_id & marginal$from == from & marginal$to == to,
    , drop = FALSE
  ]
  if (!nrow(c_row)) return("not_estimable")
  if (identical(c_row$graph_resolution[[1]], "conditional")) {
    if (isTRUE(c_row$final_edge_present[[1]])) "conditional_supported" else "not_supported"
  } else if (nrow(m_row)) {
    m_row$consensus_status[[1]]
  } else {
    "not_estimable"
  }
}

original_edge_label <- function(triad_id, from, to) {
  row <- marginal[
    marginal$triad_id == triad_id & marginal$from == from & marginal$to == to,
    , drop = FALSE
  ]
  if (nrow(row)) row$consensus_status[[1]] else "not_estimable"
}

required_role_edges <- function(role) {
  switch(
    role,
    confounder = c("Z_to_X", "Z_to_Y"),
    collider = c("X_to_Z", "Y_to_Z"),
    independent_cause = "Z_to_Y",
    upstream_surrogate_of_exposure = "Z_to_X",
    downstream_surrogate_of_outcome = "Y_to_Z",
    downstream_surrogate_of_exposure = "X_to_Z",
    character()
  )
}

theoretical_response_class <- function(role) {
  if (role %in% c("confounder", "independent_cause")) {
    return("projection_dependent_alpha_Z_may_be_nonzero")
  }
  if (role %in% c(
    "collider", "upstream_surrogate_of_exposure",
    "downstream_surrogate_of_outcome", "downstream_surrogate_of_exposure"
  )) return("structural_alpha_Z_zero_under_prespecified_model")
  "outside_or_unresolved"
}

adjustment_class <- function(role) {
  switch(
    role,
    confounder = "necessary_adjustment_under_valid_MVMR_assumptions",
    independent_cause = "usually_unnecessary_instrument_set_dependent",
    collider = "overadjustment",
    upstream_surrogate_of_exposure = "overadjustment",
    downstream_surrogate_of_outcome = "overadjustment",
    downstream_surrogate_of_exposure = "overadjustment",
    mediator_out_of_library = "estimand_changes_from_total_toward_direct_effect",
    "not_assigned"
  )
}

role_rows <- list()
for (i in seq_len(nrow(triads))) {
  triad <- triads[i, , drop = FALSE]
  X <- "BMI"
  Z <- triad$candidate
  Y <- triad$outcome
  states <- c(
    X_to_Z = edge_label(triad$triad_id, X, Z),
    Z_to_X = edge_label(triad$triad_id, Z, X),
    Z_to_Y = edge_label(triad$triad_id, Z, Y),
    Y_to_Z = edge_label(triad$triad_id, Y, Z)
  )
  original <- c(
    X_to_Z = original_edge_label(triad$triad_id, X, Z),
    Z_to_X = original_edge_label(triad$triad_id, Z, X),
    Z_to_Y = original_edge_label(triad$triad_id, Z, Y),
    Y_to_Z = original_edge_label(triad$triad_id, Y, Z)
  )
  assigned <- assign_role_from_edges(states, overlapping_measure = triad$overlapping_measure)
  c_rows <- conditional[conditional$triad_id == triad$triad_id, , drop = FALSE]
  graph_resolution <- if (all(c_rows$graph_resolution == "conditional")) "conditional" else "marginal_unresolved"
  mrsl_status <- paste(sort(unique(c_rows$mrsl_status)), collapse = ";")

  required <- required_role_edges(assigned$role)
  original_required <- original[required]
  strong_required <- length(required) > 0L && all(original_required %in% c(
    "robust_supported", "literature_supported"
  ))
  if (assigned$confidence == "unclassified") {
    confidence <- "unclassified"
  } else if (graph_resolution == "conditional" && strong_required) {
    confidence <- "high"
  } else if (graph_resolution == "conditional") {
    confidence <- "moderate"
  } else {
    confidence <- "low_exploratory"
  }

  downgrade <- character()
  if (graph_resolution != "conditional") downgrade <- c(downgrade, "conditional_graph_unresolved")
  if (any(original %in% c("unstable", "weak_instrument", "not_estimable"))) {
    downgrade <- c(downgrade, "edge_evidence_unstable_or_weak")
  }
  if (assigned$role %in% c(
    "independent_cause", "upstream_surrogate_of_exposure",
    "downstream_surrogate_of_outcome", "downstream_surrogate_of_exposure"
  )) downgrade <- c(downgrade, "role_depends_on_unsupported_edge_not_proven_absent")

  role_rows[[triad$triad_id]] <- data.frame(
    triad_id = triad$triad_id,
    candidate = Z,
    outcome = Y,
    role = assigned$role,
    role_confidence = confidence,
    role_reason = assigned$reason,
    X_to_Z = states[["X_to_Z"]],
    Z_to_X = states[["Z_to_X"]],
    Z_to_Y = states[["Z_to_Y"]],
    Y_to_Z = states[["Y_to_Z"]],
    graph_resolution = graph_resolution,
    mrsl_status = mrsl_status,
    theoretical_response = theoretical_response_class(assigned$role),
    adjustment_interpretation = adjustment_class(assigned$role),
    downgrade_reasons = paste(unique(downgrade), collapse = ";"),
    stringsAsFactors = FALSE
  )
}
roles <- do.call(rbind, role_rows)

combined <- merge(roles, response, by = c("triad_id", "candidate", "outcome"), all.x = TRUE)
outcome_type_map <- setNames(cfg$outcomes$analysis_type, cfg$outcomes$name)
combined$outcome_type <- unname(outcome_type_map[combined$outcome])
combined$response_concordance <- vapply(seq_len(nrow(combined)), function(i) {
  row <- combined[i, , drop = FALSE]
  if (!identical(row$status[[1]], "estimated") || !isTRUE(row$identity_verified[[1]])) {
    return("indeterminate_estimation_or_identity")
  }
  if (identical(row$outcome_type[[1]], "binary")) {
    return("descriptive_binary_total_response_scale_sensitive")
  }
  if (row$theoretical_response[[1]] == "outside_or_unresolved") return("not_applicable")
  if (row$theoretical_response[[1]] == "structural_alpha_Z_zero_under_prespecified_model") {
    if (is.finite(row$alpha_z_pval[[1]]) && row$alpha_z_pval[[1]] >= 0.05) {
      return("consistent_with_zero_alpha_Z")
    }
    return("inconsistent_or_binary_scale_sensitive")
  }
  if (is.finite(row$alpha_z_pval[[1]]) && row$alpha_z_pval[[1]] < 0.05) {
    return("consistent_with_nonzero_alpha_Z")
  }
  "compatible_but_nondiagnostic"
}, character(1))

if (file.exists(cfg$paths$application_response)) {
  legacy <- utils::read.csv(cfg$paths$application_response, stringsAsFactors = FALSE, check.names = FALSE)
  required_legacy <- c("outcome_abbr", "candidate_trait", "OR_unadjusted", "OR_adjusted")
  if (all(required_legacy %in% names(legacy))) {
    legacy$legacy_logOR_response <- log(legacy$OR_unadjusted) - log(legacy$OR_adjusted)
    legacy <- legacy[, c(
      "outcome_abbr", "candidate_trait", "OR_unadjusted", "OR_adjusted",
      "Delta_OR_percent", "legacy_logOR_response"
    )]
    names(legacy)[1:2] <- c("outcome", "candidate")
    combined <- merge(combined, legacy, by = c("candidate", "outcome"), all.x = TRUE)
    combined$legacy_minus_same_set_response <-
      combined$legacy_logOR_response - combined$observed_response_same_set
  }
}

main_candidates <- combined[
  combined$role %in% c(
    "confounder", "collider", "independent_cause",
    "upstream_surrogate_of_exposure", "downstream_surrogate_of_outcome",
    "downstream_surrogate_of_exposure"
  ) & combined$role_confidence %in% c("high", "moderate") &
    combined$response_concordance %in% c(
      "consistent_with_zero_alpha_Z", "consistent_with_nonzero_alpha_Z",
      "compatible_but_nondiagnostic",
      "descriptive_binary_total_response_scale_sensitive"
    ),
  , drop = FALSE
]
confidence_order <- match(main_candidates$role_confidence, c("high", "moderate"))
main_candidates <- main_candidates[order(confidence_order, main_candidates$role, main_candidates$candidate), ]

write_csv_atomic(roles, file.path(cfg$paths$output_dir, "role_assignments_168.csv"))
write_csv_atomic(combined, file.path(cfg$paths$output_dir, "role_response_concordance_168.csv"))
write_csv_atomic(main_candidates, file.path(cfg$paths$output_dir, "main_text_candidate_triads.csv"))

cat(
  "Assigned bounded role hypotheses for ", nrow(roles),
  " triads; ", nrow(main_candidates), " meet the prespecified main-text screen.\n",
  sep = ""
)
