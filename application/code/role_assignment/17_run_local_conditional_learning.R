#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "13_eraca_config_controlled_data.R"))
source(file.path(script_dir, "helpers_eraca.R"))

require_eraca_packages(c("TwoSampleMR", "MVMR", "ieugwasr", "data.table"))
cfg <- eraca_config()
ensure_eraca_directories(cfg)

args <- commandArgs(trailingOnly = TRUE)
aggregate_only <- "--aggregate-only" %in% args || eraca_env_flag("ERACA_AGGREGATE_ONLY", FALSE)
overwrite <- "--overwrite" %in% args || eraca_env_flag("ERACA_OVERWRITE", FALSE)

marginal_path <- file.path(cfg$paths$graphs, "marginal_graphs_168.rds")
if (!file.exists(marginal_path)) eraca_stop("Marginal graph file not found: ", marginal_path)
marginal_graphs <- readRDS(marginal_path)
triads <- build_triad_registry(cfg)

checkpoint_dir <- file.path(cfg$paths$graphs, "conditional_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

select_triads <- function(dat) {
  id <- Sys.getenv("ERACA_TRIAD_ID", unset = "")
  if (nzchar(id)) return(dat[dat$triad_id == id, , drop = FALSE])
  start <- suppressWarnings(as.integer(Sys.getenv("ERACA_TRIAD_START", unset = "")))
  end <- suppressWarnings(as.integer(Sys.getenv("ERACA_TRIAD_END", unset = "")))
  if (is.finite(start) || is.finite(end)) {
    if (!is.finite(start)) start <- 1L
    if (!is.finite(end)) end <- nrow(dat)
    return(dat[seq.int(max(1L, start), min(nrow(dat), end)), , drop = FALSE])
  }
  dat
}

trait_panel <- function(trait, target_snps) {
  assoc <- load_association_cache(trait, cfg)
  own <- load_instrument(trait, cfg)
  panel <- rbind(assoc, own)
  panel <- panel[panel$SNP %in% target_snps, , drop = FALSE]
  panel <- panel[!duplicated(panel$SNP), , drop = FALSE]
  panel
}

joint_instrument_snps <- function(exposures) {
  rows <- do.call(rbind, lapply(exposures, function(trait) {
    iv <- load_instrument(trait, cfg)
    data.frame(SNP = iv$SNP, pval = iv$pval, trait = trait, stringsAsFactors = FALSE)
  }))
  rows <- rows[order(rows$pval), , drop = FALSE]
  rows <- rows[!duplicated(rows$SNP), , drop = FALSE]
  if (nrow(rows) <= 1L) return(rows$SNP)
  clumped <- try(
    ieugwasr::ld_clump_local(
      dat = data.frame(rsid = rows$SNP, pval = rows$pval, id = "joint"),
      clump_kb = cfg$thresholds$clump_kb,
      clump_r2 = cfg$thresholds$clump_r2,
      clump_p = 1,
      bfile = cfg$paths$bfile,
      plink_bin = cfg$paths$plink
    ),
    silent = TRUE
  )
  if (inherits(clumped, "try-error") || is.null(clumped) || !nrow(clumped)) {
    eraca_stop("Joint IV clumping failed for exposures: ", paste(exposures, collapse = ", "))
  }
  clumped$rsid
}

make_mv_data <- function(exposures, outcome, target_snps = NULL) {
  if (is.null(target_snps)) target_snps <- joint_instrument_snps(exposures)
  exposure_long <- do.call(rbind, lapply(exposures, function(trait) {
    panel <- trait_panel(trait, target_snps)
    format_tsmr_data(panel, "exposure", trait)
  }))
  outcome_panel <- trait_panel(outcome, target_snps)
  outcome_dat <- format_tsmr_data(outcome_panel, "outcome", outcome)
  mv_dat <- TwoSampleMR::mv_harmonise_data(exposure_long, outcome_dat, harmonise_strictness = 2)
  if (is.null(mv_dat$exposure_beta) || nrow(mv_dat$exposure_beta) < 2L) {
    eraca_stop("Fewer than two SNPs remain after multivariable harmonisation.")
  }
  mv_dat
}

extract_conditional_f <- function(mv_dat) {
  formatted <- MVMR::format_mvmr(
    BXGs = mv_dat$exposure_beta,
    BYG = mv_dat$outcome_beta,
    seBXGs = mv_dat$exposure_se,
    seBYG = mv_dat$outcome_se,
    RSID = rownames(mv_dat$exposure_beta)
  )
  strength <- try(MVMR::strength_mvmr(formatted, gencov = 0), silent = TRUE)
  if (inherits(strength, "try-error") || is.null(strength)) {
    return(setNames(rep(NA_real_, ncol(mv_dat$exposure_beta)), colnames(mv_dat$exposure_beta)))
  }
  strength <- as.data.frame(strength)
  numeric_cols <- names(strength)[vapply(strength, is.numeric, logical(1))]
  if (!length(numeric_cols)) {
    return(setNames(rep(NA_real_, ncol(mv_dat$exposure_beta)), colnames(mv_dat$exposure_beta)))
  }

  # MVMR::strength_mvmr() currently returns one row (F-statistic) with one
  # numeric column per exposure (for example, exposure1 and exposure2).
  # Map that column-oriented result back to the actual exposure labels before
  # applying extract_trait_stat().  The previous implementation selected only
  # the first numeric column and named it with the row name, which made every
  # BMI/candidate lookup return NA while the log still showed valid F values.
  n_exposures <- ncol(mv_dat$exposure_beta)
  if (nrow(strength) == 1L && length(numeric_cols) >= n_exposures) {
    values <- as.numeric(strength[1, numeric_cols[seq_len(n_exposures)], drop = TRUE])
    names(values) <- colnames(mv_dat$exposure_beta)
    return(values)
  }

  values <- strength[[numeric_cols[[1]]]]
  names(values) <- rownames(strength) %||% colnames(mv_dat$exposure_beta)[seq_along(values)]
  values
}

match_trait_column <- function(label, columns) {
  idx <- match(label, columns)
  if (is.na(idx)) idx <- match(sanitize_trait_id(label), sanitize_trait_id(columns))
  idx
}

extract_trait_stat <- function(values, label) {
  if (is.null(values) || !length(values)) return(NA_real_)
  idx <- match_trait_column(label, names(values))
  if (is.na(idx)) return(NA_real_)
  as.numeric(values[[idx]])
}

run_direct_edge_mvmr <- function(from, to, adjuster, triad_id) {
  result <- tryCatch({
    mv_dat <- make_mv_data(c(from, adjuster), to)
    fit <- TwoSampleMR::mv_multiple(mv_dat)$result
    direct <- fit[fit$exposure == from, , drop = FALSE]
    if (!nrow(direct)) direct <- fit[1, , drop = FALSE]
    f_stats <- extract_conditional_f(mv_dat)
    conditional_f <- extract_trait_stat(f_stats, from)
    data.frame(
      triad_id = triad_id,
      from = from,
      to = to,
      adjuster = adjuster,
      status = "estimated",
      nsnp = nrow(mv_dat$exposure_beta),
      b = as.numeric(direct$b[[1]]),
      se = as.numeric(direct$se[[1]]),
      pval = as.numeric(direct$pval[[1]]),
      conditional_F = as.numeric(conditional_f),
      direct_edge_supported = is.finite(direct$pval[[1]]) && direct$pval[[1]] < 0.05 &&
        is.finite(conditional_f) && conditional_f >= cfg$thresholds$f_min,
      message = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      triad_id = triad_id, from = from, to = to, adjuster = adjuster,
      status = "not_estimable", nsnp = NA_integer_, b = NA_real_, se = NA_real_,
      pval = NA_real_, conditional_F = NA_real_, direct_edge_supported = NA,
      message = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  result
}

compute_same_set_response <- function(candidate, outcome, triad_id) {
  tryCatch({
    mv_dat <- make_mv_data(c("BMI", candidate), outcome)
    fit <- TwoSampleMR::mv_multiple(mv_dat)$result
    bmi_row <- fit[fit$exposure == "BMI", , drop = FALSE]
    z_row <- fit[fit$exposure == candidate, , drop = FALSE]
    if (!nrow(bmi_row) || !nrow(z_row)) eraca_stop("BMI or candidate coefficient absent from MVMR result.")

    exposure_names <- colnames(mv_dat$exposure_beta)
    x_index <- match_trait_column("BMI", exposure_names)
    z_index <- match_trait_column(candidate, exposure_names)
    if (is.na(x_index) || is.na(z_index)) eraca_stop("MVMR exposure matrix lacks named BMI/candidate columns.")
    x <- mv_dat$exposure_beta[, x_index]
    z <- mv_dat$exposure_beta[, z_index]
    y <- as.numeric(mv_dat$outcome_beta)
    beta_adjusted <- as.numeric(bmi_row$b[[1]])
    alpha_z <- as.numeric(z_row$b[[1]])
    identity <- compute_weighted_response_identity(
      x, z, y, as.numeric(mv_dat$outcome_se), beta_adjusted, alpha_z
    )
    f_stats <- extract_conditional_f(mv_dat)

    data.frame(
      triad_id = triad_id,
      candidate = candidate,
      outcome = outcome,
      status = "estimated",
      nsnp_same_set = length(x),
      beta_unadjusted_same_set = identity$beta_unadjusted,
      beta_adjusted_same_set = beta_adjusted,
      alpha_z = alpha_z,
      alpha_z_se = as.numeric(z_row$se[[1]]),
      alpha_z_pval = as.numeric(z_row$pval[[1]]),
      delta_xz = identity$delta_xz,
      observed_response_same_set = identity$observed_response,
      predicted_response_alpha_delta = identity$predicted_response,
      identity_error = identity$identity_error,
      identity_verified = abs(identity$identity_error) <= cfg$thresholds$response_identity_tolerance,
      conditional_F_BMI = extract_trait_stat(f_stats, "BMI"),
      conditional_F_candidate = extract_trait_stat(f_stats, candidate),
      message = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      triad_id = triad_id, candidate = candidate, outcome = outcome,
      status = "not_estimable", nsnp_same_set = NA_integer_,
      beta_unadjusted_same_set = NA_real_, beta_adjusted_same_set = NA_real_,
      alpha_z = NA_real_, alpha_z_se = NA_real_, alpha_z_pval = NA_real_,
      delta_xz = NA_real_, observed_response_same_set = NA_real_,
      predicted_response_alpha_delta = NA_real_, identity_error = NA_real_,
      identity_verified = FALSE, conditional_F_BMI = NA_real_,
      conditional_F_candidate = NA_real_, message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

make_mrsl_matrices <- function(nodes) {
  target_snps <- joint_instrument_snps(nodes)
  mv_dat <- make_mv_data(nodes[1:2], nodes[3], target_snps = target_snps)
  beta_matrix <- cbind(mv_dat$exposure_beta, mv_dat$outcome_beta)
  se_matrix <- cbind(mv_dat$exposure_se, mv_dat$outcome_se)
  colnames(beta_matrix) <- colnames(se_matrix) <- nodes
  rownames(beta_matrix) <- rownames(se_matrix) <- rownames(mv_dat$exposure_beta)
  iv_indicator <- sapply(nodes, function(trait) {
    as.integer(rownames(beta_matrix) %in% load_instrument(trait, cfg)$SNP)
  })
  colnames(iv_indicator) <- nodes
  rownames(iv_indicator) <- rownames(beta_matrix)
  list(beta = beta_matrix, se = se_matrix, iv = iv_indicator)
}

# MRSL 0.1.0 returns a fixed-length vector of 10 positions from its C++
# topological-sort routine, although the input graph here contains three
# nodes.  Use only the first n_nodes positions and validate them before
# applying the permutation.  The STEP2 iteration remains the package's own
# algorithm; this wrapper only protects the permutation step.
run_mrsl_fixed_topological_order <- function(M_amatt, data_sum_beta, data_sum_se,
                                              beta, cutoff, adj_methods,
                                              use_eggers_step2, vary_mvmr_adj) {
  n_nodes <- ncol(M_amatt)
  nnames <- colnames(M_amatt)
  if (is.null(nnames)) nnames <- as.character(seq_len(n_nodes))

  mrsl_ns <- asNamespace("MRSL")
  dfs_fun <- get("DFS", envir = mrsl_ns)
  step2_fun <- get("STEP2", envir = mrsl_ns)
  order0 <- as.integer(dfs_fun(amattt = M_amatt, n_nodes = n_nodes))
  order0 <- order0[seq_len(n_nodes)]

  if (length(order0) != n_nodes || any(!is.finite(order0)) ||
      any(order0 < 0L | order0 >= n_nodes) || length(unique(order0)) != n_nodes) {
    stop("MRSL returned an invalid topological order for the supplied graph.")
  }

  permutation <- n_nodes - order0
  if (any(permutation < 1L | permutation > n_nodes)) {
    stop("MRSL topological permutation is outside the adjacency-matrix bounds.")
  }

  current <- M_amatt[permutation, permutation, drop = FALSE]
  sorted_names <- nnames[permutation]
  just <- TRUE
  iteration <- 1L
  while (just) {
    iteration <- iteration + 1L
    step <- step2_fun(
      current, data_sum_beta, data_sum_se, beta,
      adj_methods, use_eggers_step2, vary_mvmr_adj
    )
    updated <- step$amatt
    just <- !all(updated == current)
    current <- updated
  }

  dimnames(current) <- list(sorted_names, sorted_names)
  list(amatt = current, iteration = iteration, nnames = sorted_names)
}

run_mrsl_local <- function(graph) {
  if (!cfg$mrsl$enabled) {
    return(list(status = "disabled", primary = graph$adjacency, strategies = list(), message = "MRSL disabled by configuration"))
  }
  if (!requireNamespace("MRSL", quietly = TRUE)) {
    return(list(status = "package_unavailable", primary = graph$adjacency, strategies = list(), message = "MRSL package is not installed"))
  }
  if (graph$cyclic) {
    return(list(status = "cyclic_marginal_graph", primary = graph$adjacency, strategies = list(), message = "MRSL topological pruning not forced for a cyclic marginal graph"))
  }
  if (!sum(graph$adjacency)) {
    return(list(status = "empty_marginal_graph", primary = graph$adjacency, strategies = list(), message = "No marginal edge to prune"))
  }

  tryCatch({
    matrices <- make_mrsl_matrices(graph$nodes)
    strategies <- list()
    primary_method <- as.integer(cfg$mrsl$adj_methods)
    if (length(primary_method) != 1L || is.na(primary_method) || primary_method != 1L) {
      stop("The installed MRSL package does not provide Open_path_adj; only adj_methods=1 is supported.")
    }
    # MRSL 0.1.0 calls an absent Open_path_adj() helper for methods 2 and 3.
    # The prespecified ERACA primary strategy is method 1, so retain that
    # analysis and record the package limitation rather than failing the triad.
    for (method in primary_method) {
      strategies[[as.character(method)]] <- run_mrsl_fixed_topological_order(
        M_amatt = graph$adjacency,
        data_sum_beta = matrices$beta,
        data_sum_se = matrices$se,
        beta = matrices$iv,
        cutoff = cfg$mrsl$cutoff,
        adj_methods = method,
        use_eggers_step2 = cfg$mrsl$use_eggers_step2,
        vary_mvmr_adj = cfg$mrsl$vary_mvmr_adj
      )
    }
    primary_key <- as.character(primary_method)
    primary <- strategies[[primary_key]]$amatt
    primary <- primary[graph$nodes, graph$nodes, drop = FALSE]
    list(
      status = "estimated",
      primary = primary,
      strategies = strategies,
      message = "Primary MRSL adjustment method 1 estimated; methods 2 and 3 were not run because MRSL 0.1.0 does not supply Open_path_adj."
    )
  }, error = function(e) {
    list(status = "error", primary = graph$adjacency, strategies = list(), message = conditionMessage(e))
  })
}

run_one_triad <- function(triad) {
  graph <- marginal_graphs[[triad$triad_id]]
  nodes <- graph$nodes
  present_edges <- which(graph$adjacency != 0, arr.ind = TRUE)
  direct_rows <- list()
  if (nrow(present_edges)) {
    for (i in seq_len(nrow(present_edges))) {
      from <- rownames(graph$adjacency)[present_edges[i, 1]]
      to <- colnames(graph$adjacency)[present_edges[i, 2]]
      adjuster <- setdiff(nodes, c(from, to))[[1]]
      direct_rows[[i]] <- run_direct_edge_mvmr(from, to, adjuster, triad$triad_id)
    }
  }
  direct <- if (length(direct_rows)) do.call(rbind, direct_rows) else data.frame()
  mrsl <- run_mrsl_local(graph)
  response <- compute_same_set_response(triad$candidate, triad$outcome, triad$triad_id)

  consensus_rows <- list()
  pairs <- which(matrix(1L, 3L, 3L) - diag(3L) == 1L, arr.ind = TRUE)
  for (i in seq_len(nrow(pairs))) {
    from <- nodes[pairs[i, 1]]
    to <- nodes[pairs[i, 2]]
    marginal_present <- graph$adjacency[from, to] == 1L
    conditional_present <- if (mrsl$status == "estimated") {
      mrsl$primary[from, to] == 1L
    } else {
      NA
    }
    mvmr_row <- if (nrow(direct)) direct[direct$from == from & direct$to == to, , drop = FALSE] else data.frame()
    consensus_rows[[i]] <- data.frame(
      triad_id = triad$triad_id,
      candidate = triad$candidate,
      outcome = triad$outcome,
      from = from,
      to = to,
      marginal_present = marginal_present,
      mrsl_status = mrsl$status,
      mrsl_conditional_present = conditional_present,
      mvmr_status = if (nrow(mvmr_row)) mvmr_row$status[[1]] else "not_run",
      mvmr_direct_supported = if (nrow(mvmr_row)) mvmr_row$direct_edge_supported[[1]] else NA,
      graph_resolution = if (mrsl$status == "estimated") "conditional" else "marginal_unresolved",
      final_edge_present = if (mrsl$status == "estimated") conditional_present else marginal_present,
      message = mrsl$message,
      stringsAsFactors = FALSE
    )
  }

  list(
    triad = triad,
    direct_edges = direct,
    mrsl = mrsl,
    consensus = do.call(rbind, consensus_rows),
    response = response
  )
}

if (!aggregate_only) {
  selected <- select_triads(triads)
  if (!nrow(selected)) eraca_stop("No triad matched the requested triad filter.")
  errors <- list()
  for (i in seq_len(nrow(selected))) {
    triad <- selected[i, , drop = FALSE]
    checkpoint <- file.path(checkpoint_dir, paste0(triad$triad_id, ".rds"))
    cat(sprintf("[%03d/%03d] BMI - %s - %s\n", i, nrow(selected), triad$candidate, triad$outcome))
    if (file.exists(checkpoint) && !overwrite) {
      cat("  checkpoint reused\n")
      next
    }
    tryCatch({
      result <- run_one_triad(triad)
      saveRDS(result, checkpoint, compress = "xz")
    }, error = function(e) {
      errors[[triad$triad_id]] <<- data.frame(
        triad_id = triad$triad_id, candidate = triad$candidate, outcome = triad$outcome,
        message = conditionMessage(e), stringsAsFactors = FALSE
      )
      cat("  ERROR:", conditionMessage(e), "\n")
    })
  }
  if (length(errors)) {
    path <- file.path(
      cfg$paths$logs,
      paste0("conditional_errors_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid(), ".csv")
    )
    write_csv_atomic(do.call(rbind, errors), path)
    eraca_warn("Some local conditional analyses failed. Review: ", path)
  }

  # A filtered worker writes only its triad checkpoints.  Do not overwrite
  # the aggregate with a partial set when several workers run concurrently.
  if (nrow(selected) < nrow(triads)) {
    cat(
      "Filtered conditional worker completed ", nrow(selected),
      " triads; aggregation was skipped. Run --step=17-aggregate after all workers finish.\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }
}

checkpoint_files <- list.files(checkpoint_dir, pattern = "^T[0-9]{3}_.*\\.rds$", full.names = TRUE)
if (!length(checkpoint_files)) eraca_stop("No conditional-learning checkpoints are available.")
results <- lapply(checkpoint_files, readRDS)

direct_edges <- data.table::rbindlist(lapply(results, `[[`, "direct_edges"), fill = TRUE)
consensus <- data.table::rbindlist(lapply(results, `[[`, "consensus"), fill = TRUE)
responses <- data.table::rbindlist(lapply(results, `[[`, "response"), fill = TRUE)
mrsl_graphs <- setNames(lapply(results, `[[`, "mrsl"), vapply(results, function(x) x$triad$triad_id, character(1)))

write_csv_atomic(as.data.frame(direct_edges), file.path(cfg$paths$graphs, "conditional_edges_mvmr.csv"))
write_csv_atomic(as.data.frame(consensus), file.path(cfg$paths$graphs, "conditional_graph_consensus.csv"))
write_csv_atomic(as.data.frame(responses), file.path(cfg$paths$graphs, "same_set_response_decomposition.csv"))
saveRDS(mrsl_graphs, file.path(cfg$paths$graphs, "mrsl_graphs_168.rds"), compress = "xz")

cat("Aggregated ", length(results), " of 168 local conditional-learning checkpoints.\n", sep = "")
