#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "13_eraca_config_controlled_data.R"))
source(file.path(script_dir, "helpers_eraca.R"))

require_eraca_packages(c("TwoSampleMR", "data.table"))
cfg <- eraca_config()
ensure_eraca_directories(cfg)
validate_eraca_config(cfg, check_files = FALSE)

args <- commandArgs(trailingOnly = TRUE)
aggregate_only <- "--aggregate-only" %in% args || eraca_env_flag("ERACA_AGGREGATE_ONLY", FALSE)
overwrite <- "--overwrite" %in% args || eraca_env_flag("ERACA_OVERWRITE", FALSE)

all_tasks <- build_direction_tasks(cfg)
write_csv_atomic(all_tasks, file.path(cfg$paths$output_dir, "direction_task_registry.csv"))
tasks <- select_task_subset(all_tasks)
if (!nrow(tasks) && !aggregate_only) eraca_stop("No direction task matched the requested task filter.")

checkpoint_path <- function(task_id) {
  file.path(cfg$paths$uvmr_checkpoints, paste0(task_id, ".rds"))
}

extract_heterogeneity_p <- function(heterogeneity) {
  if (is.null(heterogeneity) || !nrow(heterogeneity)) return(NA_real_)
  ivw <- heterogeneity[grepl("Inverse variance weighted", heterogeneity$method), , drop = FALSE]
  if (!nrow(ivw)) ivw <- heterogeneity[1, , drop = FALSE]
  as.numeric(ivw$Q_pval[[1]])
}

run_one_direction <- function(task) {
  exposure_name <- task$exposure[[1]]
  outcome_name <- task$outcome[[1]]
  iv <- load_instrument(exposure_name, cfg)
  outcome_assoc <- load_association_cache(outcome_name, cfg)
  outcome_assoc <- outcome_assoc[outcome_assoc$SNP %in% iv$SNP, , drop = FALSE]
  if (!nrow(outcome_assoc)) eraca_stop("No overlapping SNP association for task ", task$task_id)

  exposure_dat <- format_tsmr_data(iv, "exposure", exposure_name)
  outcome_dat <- format_tsmr_data(outcome_assoc, "outcome", outcome_name)
  dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2)
  dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
  if (!nrow(dat)) eraca_stop("No SNP retained after harmonisation for task ", task$task_id)

  nsnp <- nrow(dat)
  methods <- method_list_for_nsnp(nsnp)
  results <- TwoSampleMR::mr(dat, method_list = methods)
  if (!nrow(results)) eraca_stop("No MR estimate returned for task ", task$task_id)
  results$task_id <- task$task_id
  results$family <- task$family
  results$exposure_name <- exposure_name
  results$outcome_name <- outcome_name
  results$outcome_type <- task$outcome_type
  results$OR <- if (task$outcome_type == "binary") exp(results$b) else NA_real_
  results$OR_lci95 <- if (task$outcome_type == "binary") exp(results$b - 1.96 * results$se) else NA_real_
  results$OR_uci95 <- if (task$outcome_type == "binary") exp(results$b + 1.96 * results$se) else NA_real_

  heterogeneity <- if (nsnp >= 2L) {
    tryCatch(TwoSampleMR::mr_heterogeneity(dat, method_list = methods), error = function(e) data.frame())
  } else data.frame()
  pleiotropy <- if (nsnp >= 3L) {
    tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) data.frame())
  } else data.frame()
  leave_one_out <- if (nsnp >= 3L) {
    tryCatch(TwoSampleMR::mr_leaveoneout(dat), error = function(e) data.frame())
  } else data.frame()
  steiger <- tryCatch(TwoSampleMR::directionality_test(dat), error = function(e) data.frame())

  f_values <- (dat$beta.exposure / dat$se.exposure)^2
  diagnostics <- data.frame(
    task_id = task$task_id,
    family = task$family,
    exposure = exposure_name,
    outcome = outcome_name,
    nsnp = nsnp,
    mean_F = mean(f_values, na.rm = TRUE),
    median_F = stats::median(f_values, na.rm = TRUE),
    min_F = min(f_values, na.rm = TRUE),
    heterogeneity_p = extract_heterogeneity_p(heterogeneity),
    egger_intercept = if (nrow(pleiotropy)) as.numeric(pleiotropy$egger_intercept[[1]]) else NA_real_,
    egger_intercept_p = if (nrow(pleiotropy)) as.numeric(pleiotropy$pval[[1]]) else NA_real_,
    steiger_direction = if (nrow(steiger)) as.character(steiger$correct_causal_direction[[1]]) else NA_character_,
    steiger_p = if (nrow(steiger)) as.numeric(steiger$steiger_pval[[1]]) else NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    task = task,
    methods = results,
    diagnostics = diagnostics,
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy,
    leave_one_out = leave_one_out,
    steiger = steiger,
    harmonised = dat
  )
}

if (!aggregate_only) {
  error_rows <- list()
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    path <- checkpoint_path(task$task_id)
    cat(sprintf("[%03d/%03d] %s -> %s\n", i, nrow(tasks), task$exposure, task$outcome))
    if (file.exists(path) && !overwrite) {
      cat("  checkpoint reused\n")
      next
    }
    tryCatch({
      result <- run_one_direction(task)
      saveRDS(result, path, compress = "xz")
    }, error = function(e) {
      error_rows[[task$task_id]] <<- data.frame(
        task_id = task$task_id, exposure = task$exposure, outcome = task$outcome,
        message = conditionMessage(e), stringsAsFactors = FALSE
      )
      cat("  ERROR:", conditionMessage(e), "\n")
    })
  }
  if (length(error_rows)) {
    error_file <- file.path(
      cfg$paths$logs,
      paste0("uvmr_errors_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid(), ".csv")
    )
    write_csv_atomic(do.call(rbind, error_rows), error_file)
    eraca_warn("Some UVMR tasks failed. Review: ", error_file)
  }

  # A filtered worker must not overwrite the aggregate with a partial FDR
  # universe.  Aggregate explicitly after all scheduled workers finish.
  if (nrow(tasks) < nrow(all_tasks)) {
    cat(
      "Filtered UVMR worker completed ", nrow(tasks),
      " tasks; aggregation was skipped. Run --step=15-aggregate after all workers finish.\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }
}

checkpoint_files <- list.files(cfg$paths$uvmr_checkpoints, pattern = "^D[0-9]{3}_.*\\.rds$", full.names = TRUE)
if (!length(checkpoint_files)) eraca_stop("No UVMR checkpoint files are available for aggregation.")

checkpoints <- lapply(checkpoint_files, readRDS)
all_methods <- data.table::rbindlist(lapply(checkpoints, `[[`, "methods"), fill = TRUE)
diagnostics <- data.table::rbindlist(lapply(checkpoints, `[[`, "diagnostics"), fill = TRUE)

primary_rows <- lapply(checkpoints, function(x) {
  nsnp <- x$diagnostics$nsnp[[1]]
  primary_name <- primary_method_for_nsnp(nsnp)
  hit <- x$methods[x$methods$method == primary_name, , drop = FALSE]
  if (!nrow(hit)) {
    hit <- x$methods[grepl(primary_name, x$methods$method, fixed = TRUE), , drop = FALSE]
  }
  if (!nrow(hit)) hit <- x$methods[1, , drop = FALSE]
  hit[1, , drop = FALSE]
})
primary <- data.table::rbindlist(primary_rows, fill = TRUE)
primary <- merge(primary, diagnostics, by = c("task_id"), all.x = TRUE, suffixes = c("", ".diag"))
primary$p_fdr <- stats::p.adjust(primary$pval, method = "BH")

secondary_concordance <- function(task_id, primary_b) {
  rows <- all_methods[all_methods$task_id == task_id, , drop = FALSE]
  rows <- rows[!grepl("Wald ratio|Inverse variance weighted", rows$method), , drop = FALSE]
  if (!nrow(rows)) return(FALSE)
  any(is.finite(rows$b) & sign(rows$b) == sign(primary_b) & rows$pval < 0.05)
}

primary$secondary_concordant <- mapply(
  secondary_concordance, primary$task_id, primary$b,
  USE.NAMES = FALSE
)
primary$direction_conflict <- vapply(seq_len(nrow(primary)), function(i) {
  rows <- all_methods[all_methods$task_id == primary$task_id[[i]], , drop = FALSE]
  rows <- rows[rows$pval < 0.05 & is.finite(rows$b), , drop = FALSE]
  nrow(rows) > 1L && length(unique(sign(rows$b))) > 1L
}, logical(1))
primary$unstable <- primary$direction_conflict |
  (is.finite(primary$egger_intercept_p) & primary$egger_intercept_p < 0.05)
primary$edge_status <- mapply(
  classify_edge_evidence,
  pval = primary$pval,
  p_fdr = primary$p_fdr,
  b = primary$b,
  nsnp = primary$nsnp,
  mean_f = primary$mean_F,
  secondary_concordant = primary$secondary_concordant,
  unstable = primary$unstable,
  MoreArgs = list(
    f_min = cfg$thresholds$f_min,
    nominal_alpha = cfg$thresholds$nominal_alpha,
    fdr_alpha = cfg$thresholds$fdr_alpha
  ),
  USE.NAMES = FALSE
)
primary$fdr_scope_n <- nrow(primary)
primary$all_388_available <- nrow(primary) == 388L

write_csv_atomic(as.data.frame(all_methods), file.path(cfg$paths$output_dir, "uvmr_all_methods.csv"))
write_csv_atomic(as.data.frame(diagnostics), file.path(cfg$paths$output_dir, "uvmr_diagnostics.csv"))
write_csv_atomic(as.data.frame(primary), file.path(cfg$paths$output_dir, "uvmr_primary_edges.csv"))

cat(
  "UVMR aggregation completed with ", nrow(primary),
  " of 388 direction tasks. FDR scope n=", nrow(primary), ".\n",
  sep = ""
)
