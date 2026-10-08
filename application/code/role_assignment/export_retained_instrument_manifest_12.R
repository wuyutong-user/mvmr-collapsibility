#!/usr/bin/env Rscript

# Export the exact SNP lineage used by the strict same-set MVMR response for
# the 12 resolved Independent-cause triads. This script is intentionally an
# export-only step: it does not alter the role assignments or rerun the full
# application pipeline.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "13_eraca_config_controlled_data.R"))
source(file.path(script_dir, "helpers_eraca.R"))

require_eraca_packages(c("TwoSampleMR", "MVMR", "ieugwasr", "data.table"))
cfg <- eraca_config()
ensure_eraca_directories(cfg)

triad_ids <- c(
  "T043_BMI_DBP_CAD", "T044_BMI_DBP_HF", "T047_BMI_DBP_AAA",
  "T048_BMI_DBP_TAA", "T050_BMI_DBP_PE", "T051_BMI_DBP_SAH",
  "T052_BMI_DBP_ICH", "T053_BMI_DBP_DVT", "T054_BMI_DBP_AVS",
  "T055_BMI_DBP_AF", "T103_BMI_TG_AAA", "T110_BMI_TG_AVS"
)

triads <- build_triad_registry(cfg)
selected <- triads[match(triad_ids, triads$triad_id), , drop = FALSE]
if (any(is.na(selected$triad_id))) {
  eraca_stop("One or more requested triads were not found in the registry.")
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
    data.frame(SNP = iv$SNP, pval = iv$pval, trait = trait,
      stringsAsFactors = FALSE)
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
    ), silent = TRUE
  )
  if (inherits(clumped, "try-error") || is.null(clumped) || !nrow(clumped)) {
    eraca_stop("Joint IV clumping failed for exposures: ",
      paste(exposures, collapse = ", "))
  }
  clumped$rsid
}

make_mv_data_with_sources <- function(exposures, outcome) {
  target_snps <- joint_instrument_snps(exposures)
  exposure_panels <- lapply(exposures, function(trait) trait_panel(trait, target_snps))
  exposure_long <- do.call(rbind, lapply(seq_along(exposures), function(i) {
    format_tsmr_data(exposure_panels[[i]], "exposure", exposures[[i]])
  }))
  outcome_panel <- trait_panel(outcome, target_snps)
  outcome_dat <- format_tsmr_data(outcome_panel, "outcome", outcome)
  mv_dat <- TwoSampleMR::mv_harmonise_data(
    exposure_long, outcome_dat, harmonise_strictness = 2
  )
  if (is.null(mv_dat$exposure_beta) || nrow(mv_dat$exposure_beta) < 2L) {
    eraca_stop("Fewer than two SNPs remain after multivariable harmonisation.")
  }
  list(
    target_snps = target_snps,
    exposure_long = exposure_long,
    outcome_dat = outcome_dat,
    mv_dat = mv_dat,
    instrument_sets = lapply(exposures, function(trait) load_instrument(trait, cfg))
  )
}

matrix_value <- function(x, snp, exposure) {
  if (is.null(x)) return(NA_real_)
  if (is.null(dim(x))) return(NA_real_)
  if (!(snp %in% rownames(x)) || !(exposure %in% colnames(x))) return(NA_real_)
  as.numeric(x[snp, exposure])
}

vector_value <- function(x, snp) {
  if (is.null(x)) return(NA_real_)
  if (!is.null(names(x)) && snp %in% names(x)) return(as.numeric(x[[snp]]))
  NA_real_
}

source_numeric <- function(dat, snp, column) {
  if (is.null(dat) || !nrow(dat) || !(column %in% names(dat))) return(NA_real_)
  hit <- dat[dat$SNP == snp, column, drop = TRUE]
  if (!length(hit)) return(NA_real_)
  as.numeric(hit[[1]])
}

p_from_beta_se <- function(beta, se) {
  if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NA_real_)
  2 * stats::pnorm(-abs(beta / se))
}

source_value <- function(dat, snp, column) {
  if (is.null(dat) || !nrow(dat) || !(column %in% names(dat))) return(NA_character_)
  hit <- dat[dat$SNP == snp, column, drop = TRUE]
  if (!length(hit)) return(NA_character_)
  as.character(hit[[1]])
}

rows <- list()
errors <- list()
for (i in seq_len(nrow(selected))) {
  triad <- selected[i, , drop = FALSE]
  tid <- triad$triad_id[[1]]
  candidate <- triad$candidate[[1]]
  outcome <- triad$outcome[[1]]
  cat(sprintf("[%02d/%02d] Exporting %s\n", i, nrow(selected), tid))
  result <- tryCatch({
    x <- make_mv_data_with_sources(c("BMI", candidate), outcome)
    mv <- x$mv_dat
    final_snps <- rownames(mv$exposure_beta)
    bmi_iv <- x$instrument_sets[[1]]
    z_iv <- x$instrument_sets[[2]]
    bmi_ids <- bmi_iv$SNP
    z_ids <- z_iv$SNP
    exp_names <- colnames(mv$exposure_beta)
    if (!("BMI" %in% exp_names) || !(candidate %in% exp_names)) {
      eraca_stop("Harmonised exposure matrix lacks BMI or candidate column.")
    }
    do.call(rbind, lapply(x$target_snps, function(snp) {
      bmi_source <- x$exposure_long[x$exposure_long$exposure == "BMI", , drop = FALSE]
      z_source <- x$exposure_long[x$exposure_long$exposure == candidate, , drop = FALSE]
      # Use the harmonised matrices for retained rows; retain source-level
      # association values for post-clump rows removed at outcome harmonisation.
      beta_bmi <- matrix_value(mv$exposure_beta, snp, "BMI")
      if (!is.finite(beta_bmi)) beta_bmi <- source_numeric(bmi_source, snp, "beta.exposure")
      se_bmi <- matrix_value(mv$exposure_se, snp, "BMI")
      if (!is.finite(se_bmi)) se_bmi <- source_numeric(bmi_source, snp, "se.exposure")
      beta_z <- matrix_value(mv$exposure_beta, snp, candidate)
      if (!is.finite(beta_z)) beta_z <- source_numeric(z_source, snp, "beta.exposure")
      se_z <- matrix_value(mv$exposure_se, snp, candidate)
      if (!is.finite(se_z)) se_z <- source_numeric(z_source, snp, "se.exposure")
      beta_y <- vector_value(mv$outcome_beta, snp)
      if (!is.finite(beta_y)) beta_y <- source_numeric(x$outcome_dat, snp, "beta.outcome")
      se_y <- vector_value(mv$outcome_se, snp)
      if (!is.finite(se_y)) se_y <- source_numeric(x$outcome_dat, snp, "se.outcome")
      bmi_flag <- snp %in% bmi_ids
      z_flag <- snp %in% z_ids
      source_origin <- if (bmi_flag && z_flag) {
        "BMI_and_candidate_instrument_sets"
      } else if (bmi_flag) {
        "BMI_instrument_set_only"
      } else if (z_flag) {
        "candidate_instrument_set_only"
      } else {
        "post_union_clump_only"
      }
      data.frame(
        triad_id = tid, candidate = candidate, outcome = outcome, SNP = snp,
        final_retained = snp %in% final_snps,
        BMI_selection_flag = bmi_flag,
        candidate_selection_flag = z_flag,
        SNP_source_selection_origin = source_origin,
        beta_BMI = beta_bmi, se_BMI = se_bmi,
        p_BMI = p_from_beta_se(beta_bmi, se_bmi),
        beta_Z = beta_z, se_Z = se_z,
        p_Z = p_from_beta_se(beta_z, se_z),
        beta_Y = beta_y, se_Y = se_y,
        p_Y = p_from_beta_se(beta_y, se_y),
        effect_allele_BMI = source_value(x$exposure_long[x$exposure_long$exposure == "BMI", , drop = FALSE], snp, "effect_allele.exposure"),
        other_allele_BMI = source_value(x$exposure_long[x$exposure_long$exposure == "BMI", , drop = FALSE], snp, "other_allele.exposure"),
        effect_allele_Z = source_value(x$exposure_long[x$exposure_long$exposure == candidate, , drop = FALSE], snp, "effect_allele.exposure"),
        other_allele_Z = source_value(x$exposure_long[x$exposure_long$exposure == candidate, , drop = FALSE], snp, "other_allele.exposure"),
        effect_allele_Y = source_value(x$outcome_dat, snp, "effect_allele.outcome"),
        other_allele_Y = source_value(x$outcome_dat, snp, "other_allele.outcome"),
        post_union_clump_retained = TRUE,
        outcome_harmonised_retained = snp %in% final_snps,
        final_same_set_row = snp %in% final_snps,
        stringsAsFactors = FALSE
      )
    }))
  }, error = function(e) {
    errors[[tid]] <<- data.frame(
      triad_id = tid, candidate = candidate, outcome = outcome,
      message = conditionMessage(e), stringsAsFactors = FALSE
    )
    NULL
  })
  if (!is.null(result)) rows[[tid]] <- result
}

if (length(rows)) {
  write_csv_atomic(
    do.call(rbind, rows),
    file.path(cfg$paths$output_dir, "retained_instrument_manifest_12.csv")
  )
}
if (length(errors)) {
  write_csv_atomic(
    do.call(rbind, errors),
    file.path(cfg$paths$logs, "retained_instrument_manifest_12_errors.csv")
  )
  eraca_stop("Manifest export failed for ", length(errors), " triad(s).")
}

cat("Retained-instrument manifest export completed for ", length(rows),
  " triads.\n", sep = "")
