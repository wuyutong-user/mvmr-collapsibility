#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "13_eraca_config_controlled_data.R"))
source(file.path(script_dir, "helpers_eraca.R"))

cfg <- eraca_config()
ensure_eraca_directories(cfg)

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args || eraca_env_flag("ERACA_DRY_RUN", FALSE)
overwrite <- "--overwrite" %in% args || eraca_env_flag("ERACA_OVERWRITE", FALSE)

validate_eraca_config(cfg, check_files = TRUE)
manifest <- resolve_outcome_manifest(cfg, allow_template = TRUE)
traits <- eraca_trait_registry(cfg)

source_path <- function(trait) {
  if (trait == "BMI") return(find_bmi_gwas(cfg))
  if (trait %in% cfg$candidates$name) return(find_candidate_gwas(trait, cfg))
  manifest$zip_path[match(trait, manifest$name)]
}

input_manifest <- traits
input_manifest$source_path <- vapply(traits$name, source_path, character(1))
input_manifest$file_exists <- file.exists(input_manifest$source_path)
input_manifest$instrument_path <- vapply(traits$name, instrument_path, character(1), cfg = cfg)
input_manifest$association_cache_path <- vapply(
  traits$name, association_cache_path, character(1), cfg = cfg
)
write_csv_atomic(input_manifest, file.path(cfg$paths$output_dir, "input_manifest.csv"))

if (dry_run) {
  cat("Dry run completed: 27 trait inputs and output paths validated.\n")
  quit(save = "no", status = 0)
}

instrument_qc <- list()
instrument_errors <- list()
instrument_unavailable <- list()

for (i in seq_len(nrow(traits))) {
  trait <- traits$name[[i]]
  out_path <- instrument_path(trait, cfg)
  cat(sprintf("[%02d/27] Preparing instruments for %s\n", i, trait))

  result <- tryCatch({
    if (file.exists(out_path) && !overwrite) {
      iv <- readRDS(out_path)
      status <- "checkpoint_reused"
    } else {
      significant <- read_trait_associations(
        trait, cfg, p_threshold = cfg$thresholds$instrument_p,
        outcome_manifest = manifest
      )
      iv <- clump_trait_instruments(significant, cfg)
      if (!nrow(iv)) eraca_stop("No genome-wide-significant instrument for ", trait)
      saveRDS(iv, out_path, compress = "xz")
      status <- "prepared"
    }
    data.frame(
      trait = trait,
      status = status,
      nsnp = nrow(iv),
      mean_F = if ("F_stat" %in% names(iv)) mean(iv$F_stat, na.rm = TRUE) else NA_real_,
      median_F = if ("F_stat" %in% names(iv)) stats::median(iv$F_stat, na.rm = TRUE) else NA_real_,
      min_F = if ("F_stat" %in% names(iv)) min(iv$F_stat, na.rm = TRUE) else NA_real_,
      source_path = source_path(trait),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    if (trait %in% cfg$outcomes$name) {
      # An outcome-specific instrument set is only required for reverse UVMR.
      # If outcome SNPs cannot be linked to the local LD reference, preserve
      # forward analyses and mark reverse directions as not estimable.
      empty_iv <- empty_instrument_table(trait, "binary")
      saveRDS(empty_iv, out_path, compress = "xz")
      instrument_unavailable[[trait]] <<- data.frame(
        stage = "instrument", trait = trait,
        status = "outcome_instrument_unavailable",
        message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      data.frame(
        trait = trait,
        status = "outcome_instrument_unavailable",
        nsnp = 0L,
        mean_F = NA_real_, median_F = NA_real_, min_F = NA_real_,
        source_path = source_path(trait),
        stringsAsFactors = FALSE
      )
    } else {
      instrument_errors[[trait]] <<- data.frame(
        stage = "instrument", trait = trait,
        message = conditionMessage(e), stringsAsFactors = FALSE
      )
      NULL
    }
  })
  if (!is.null(result)) instrument_qc[[trait]] <- result
}

if (length(instrument_qc)) {
  write_csv_atomic(do.call(rbind, instrument_qc), file.path(cfg$paths$output_dir, "instrument_qc.csv"))
}
if (length(instrument_errors)) {
  write_csv_atomic(do.call(rbind, instrument_errors), file.path(cfg$paths$logs, "instrument_errors.csv"))
  eraca_stop(
    "Instrument preparation failed for ", length(instrument_errors),
    " trait(s). Review logs/instrument_errors.csv and rerun."
  )
}
if (length(instrument_unavailable)) {
  write_csv_atomic(
    do.call(rbind, instrument_unavailable),
    file.path(cfg$paths$logs, "instrument_unavailable.csv")
  )
  cat(
    "Outcome instrument sets unavailable for ", length(instrument_unavailable),
    " outcome(s); forward analyses will proceed and reverse directions will be marked not estimable.\n",
    sep = ""
  )
}

tasks <- build_direction_tasks(cfg)
cache_qc <- list()
cache_errors <- list()

for (i in seq_len(nrow(traits))) {
  trait <- traits$name[[i]]
  out_path <- association_cache_path(trait, cfg)
  cat(sprintf("[%02d/27] Building association cache for %s\n", i, trait))

  result <- tryCatch({
    if (file.exists(out_path) && !overwrite) {
      assoc <- readRDS(out_path)
      status <- "checkpoint_reused"
    } else {
      connected_exposures <- unique(tasks$exposure[tasks$outcome == trait])
      target_snps <- unique(load_instrument(trait, cfg)$SNP)
      for (exposure in connected_exposures) {
        target_snps <- union(target_snps, load_instrument(exposure, cfg)$SNP)
      }
      assoc <- read_trait_associations(
        trait, cfg, target_snps = target_snps, outcome_manifest = manifest
      )
      own_iv <- load_instrument(trait, cfg)
      # Instrument tables contain F_stat, whereas association tables do not.
      # Bind by column name and fill absent fields rather than using base rbind.
      assoc <- data.table::rbindlist(
        list(assoc, own_iv), fill = TRUE, use.names = TRUE
      )
      assoc <- as.data.frame(assoc, stringsAsFactors = FALSE)
      assoc <- assoc[order(match(assoc$SNP, target_snps)), , drop = FALSE]
      assoc <- assoc[!duplicated(assoc$SNP), , drop = FALSE]
      saveRDS(assoc, out_path, compress = "xz")
      status <- "prepared"
    }
    data.frame(
      trait = trait, status = status, cached_snps = nrow(assoc),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cache_errors[[trait]] <<- data.frame(
      stage = "association_cache", trait = trait, message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    NULL
  })
  if (!is.null(result)) cache_qc[[trait]] <- result
}

if (length(cache_qc)) {
  write_csv_atomic(do.call(rbind, cache_qc), file.path(cfg$paths$output_dir, "association_cache_qc.csv"))
}
if (length(cache_errors)) {
  write_csv_atomic(do.call(rbind, cache_errors), file.path(cfg$paths$logs, "association_cache_errors.csv"))
  eraca_stop(
    "Association-cache preparation failed for ", length(cache_errors),
    " trait(s). Review logs/association_cache_errors.csv and rerun."
  )
}

cat("Instrument and association-cache preparation completed for all 27 traits.\n")
