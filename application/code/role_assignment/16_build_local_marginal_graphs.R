#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "13_eraca_config_controlled_data.R"))
source(file.path(script_dir, "helpers_eraca.R"))

cfg <- eraca_config()
ensure_eraca_directories(cfg)

uvmr_path <- file.path(cfg$paths$output_dir, "uvmr_primary_edges.csv")
if (!file.exists(uvmr_path)) eraca_stop("UVMR primary edge file not found: ", uvmr_path)
uvmr <- utils::read.csv(uvmr_path, stringsAsFactors = FALSE, check.names = FALSE)
required_uvmr <- c("task_id", "exposure_name", "outcome_name", "edge_status")
missing_uvmr <- setdiff(required_uvmr, names(uvmr))
if (length(missing_uvmr)) eraca_stop("UVMR edge file lacks: ", paste(missing_uvmr, collapse = ", "))

literature <- read_literature_edges(cfg$paths$literature_edges)
triads <- build_triad_registry(cfg)
graphs <- vector("list", nrow(triads))
names(graphs) <- triads$triad_id
ledger_rows <- list()

edge_pairs <- function(nodes) {
  do.call(rbind, lapply(nodes, function(from) {
    data.frame(from = from, to = setdiff(nodes, from), stringsAsFactors = FALSE)
  }))
}

lookup_mr <- function(from, to) {
  row <- uvmr[uvmr$exposure_name == from & uvmr$outcome_name == to, , drop = FALSE]
  if (!nrow(row)) return(list(status = "not_estimable", task_id = NA_character_, p = NA_real_, p_fdr = NA_real_))
  list(
    status = row$edge_status[[1]],
    task_id = row$task_id[[1]],
    p = if ("pval" %in% names(row)) row$pval[[1]] else NA_real_,
    p_fdr = if ("p_fdr" %in% names(row)) row$p_fdr[[1]] else NA_real_
  )
}

lookup_literature <- function(candidate, outcome, from, to) {
  if (!nrow(literature)) return(list(status = NA_character_, tier = NA_character_, reference = NA_character_))
  row <- literature[
    literature$candidate == candidate & literature$outcome == outcome &
      literature$from == from & literature$to == to,
    , drop = FALSE
  ]
  if (!nrow(row)) return(list(status = NA_character_, tier = NA_character_, reference = NA_character_))
  list(
    status = row$direction_status[[1]],
    tier = if ("evidence_tier" %in% names(row)) row$evidence_tier[[1]] else NA_character_,
    reference = if ("reference_id" %in% names(row)) row$reference_id[[1]] else NA_character_
  )
}

for (i in seq_len(nrow(triads))) {
  triad <- triads[i, , drop = FALSE]
  nodes <- c("BMI", triad$candidate, triad$outcome)
  adjacency <- matrix(0L, 3L, 3L, dimnames = list(nodes, nodes))
  status_matrix <- matrix("not_estimable", 3L, 3L, dimnames = list(nodes, nodes))
  diag(status_matrix) <- "self"
  pairs <- edge_pairs(nodes)
  local_rows <- list()

  for (j in seq_len(nrow(pairs))) {
    from <- pairs$from[[j]]
    to <- pairs$to[[j]]
    mr <- lookup_mr(from, to)
    lit <- lookup_literature(triad$candidate, triad$outcome, from, to)
    consensus <- combine_edge_evidence(mr$status, lit$status)
    status_matrix[from, to] <- consensus
    adjacency[from, to] <- as.integer(edge_is_supported(consensus))
    local_rows[[j]] <- data.frame(
      triad_id = triad$triad_id,
      candidate = triad$candidate,
      outcome = triad$outcome,
      from = from,
      to = to,
      uvmr_task_id = mr$task_id,
      mr_status = mr$status,
      mr_p = mr$p,
      mr_p_fdr = mr$p_fdr,
      literature_status = lit$status,
      literature_tier = lit$tier,
      literature_reference = lit$reference,
      consensus_status = consensus,
      edge_present = edge_is_supported(consensus),
      stringsAsFactors = FALSE
    )
  }
  local_ledger <- do.call(rbind, local_rows)
  ledger_rows[[triad$triad_id]] <- local_ledger
  graphs[[triad$triad_id]] <- list(
    triad = triad,
    nodes = nodes,
    adjacency = adjacency,
    status = status_matrix,
    cyclic = has_directed_cycle(adjacency),
    edge_ledger = local_ledger
  )
}

ledger <- do.call(rbind, ledger_rows)
if (nrow(ledger) != 168L * 6L) {
  eraca_stop("Expected 1008 marginal edge-ledger rows; found ", nrow(ledger), ".")
}

saveRDS(graphs, file.path(cfg$paths$graphs, "marginal_graphs_168.rds"), compress = "xz")
write_csv_atomic(ledger, file.path(cfg$paths$graphs, "marginal_edge_ledger.csv"))
write_csv_atomic(
  data.frame(
    triad_id = names(graphs),
    cyclic = vapply(graphs, `[[`, logical(1), "cyclic"),
    n_edges = vapply(graphs, function(x) sum(x$adjacency), integer(1)),
    stringsAsFactors = FALSE
  ),
  file.path(cfg$paths$graphs, "marginal_graph_qc.csv")
)

cat("Constructed 168 local marginal graphs and 1008 directed edge records.\n")
