#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

eraca_stop <- function(...) stop(paste0(...), call. = FALSE)
eraca_warn <- function(...) warning(paste0(...), call. = FALSE, immediate. = TRUE)

sanitize_trait_id <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

ensure_eraca_directories <- function(cfg) {
  dirs <- unlist(cfg$paths[c(
    "output_dir", "instruments", "association_cache", "uvmr_checkpoints",
    "graphs", "logs"
  )], use.names = FALSE)
  for (path in dirs) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) eraca_stop("Cannot create output directory: ", path)
  }
  invisible(dirs)
}

require_eraca_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    eraca_stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Install them on the server before running this step."
    )
  }
  invisible(TRUE)
}

eraca_trait_registry <- function(cfg) {
  bmi <- cfg$exposure
  bmi$role_family <- "exposure"
  candidates <- cfg$candidates
  candidates$role_family <- "candidate"
  outcomes <- cfg$outcomes[, c("name", "source_kind", "analysis_type")]
  outcomes$raw_id <- NA_character_
  outcomes$overlapping_measure <- FALSE
  outcomes$role_family <- "outcome"

  keep <- c(
    "name", "raw_id", "source_kind", "analysis_type",
    "overlapping_measure", "role_family"
  )
  if (!"overlapping_measure" %in% names(bmi)) bmi$overlapping_measure <- FALSE
  do.call(rbind, lapply(list(bmi, candidates, outcomes), function(x) x[, keep]))
}

validate_eraca_config <- function(cfg, check_files = TRUE) {
  if (nrow(cfg$exposure) != 1L || cfg$exposure$name[[1]] != "BMI") {
    eraca_stop("Configuration must contain exactly one BMI exposure row.")
  }
  if (nrow(cfg$candidates) != 12L) eraca_stop("Expected 12 candidate traits.")
  if (nrow(cfg$outcomes) != 14L) eraca_stop("Expected 14 outcomes.")
  if (anyDuplicated(cfg$candidates$name)) eraca_stop("Candidate trait names are duplicated.")
  if (anyDuplicated(cfg$outcomes$name)) eraca_stop("Outcome names are duplicated.")

  numeric_thresholds <- unlist(cfg$thresholds, use.names = TRUE)
  if (any(!is.finite(numeric_thresholds))) eraca_stop("All thresholds must be finite.")
  if (cfg$thresholds$instrument_p <= 0 || cfg$thresholds$instrument_p >= 1) {
    eraca_stop("instrument_p must be between 0 and 1.")
  }

  if (check_files) {
    required <- c(cfg$paths$variants, cfg$paths$plink)
    missing <- required[!file.exists(required)]
    if (length(missing)) eraca_stop("Missing required files: ", paste(missing, collapse = "; "))
    if (!file.exists(paste0(cfg$paths$bfile, ".bed"))) {
      eraca_stop("LD reference .bed file not found for prefix: ", cfg$paths$bfile)
    }
    if (!dir.exists(cfg$paths$candidate_dir)) {
      eraca_stop("Candidate GWAS directory not found: ", cfg$paths$candidate_dir)
    }
    if (!dir.exists(cfg$paths$outcome_dir)) {
      eraca_stop("Outcome ZIP directory not found: ", cfg$paths$outcome_dir)
    }
  }
  invisible(TRUE)
}

build_direction_tasks <- function(cfg) {
  trait_type <- setNames(
    eraca_trait_registry(cfg)$analysis_type,
    eraca_trait_registry(cfg)$name
  )
  rows <- list()
  add_pair <- function(a, b, family) {
    list(
      data.frame(family = family, exposure = a, outcome = b, stringsAsFactors = FALSE),
      data.frame(family = family, exposure = b, outcome = a, stringsAsFactors = FALSE)
    )
  }

  for (z in cfg$candidates$name) rows <- c(rows, add_pair("BMI", z, "bmi_candidate"))
  for (z in cfg$candidates$name) {
    for (y in cfg$outcomes$name) rows <- c(rows, add_pair(z, y, "candidate_outcome"))
  }
  for (y in cfg$outcomes$name) rows <- c(rows, add_pair("BMI", y, "bmi_outcome"))

  tasks <- do.call(rbind, rows)
  tasks$exposure_type <- unname(trait_type[tasks$exposure])
  tasks$outcome_type <- unname(trait_type[tasks$outcome])
  tasks$task_id <- sprintf(
    "D%03d_%s_to_%s",
    seq_len(nrow(tasks)), sanitize_trait_id(tasks$exposure), sanitize_trait_id(tasks$outcome)
  )
  tasks <- tasks[, c(
    "task_id", "family", "exposure", "outcome", "exposure_type", "outcome_type"
  )]
  rownames(tasks) <- NULL
  tasks
}

build_triad_registry <- function(cfg) {
  triads <- expand.grid(
    candidate = cfg$candidates$name,
    outcome = cfg$outcomes$name,
    stringsAsFactors = FALSE
  )
  candidate_order <- setNames(seq_len(nrow(cfg$candidates)), cfg$candidates$name)
  outcome_order <- setNames(seq_len(nrow(cfg$outcomes)), cfg$outcomes$name)
  triads <- triads[order(candidate_order[triads$candidate], outcome_order[triads$outcome]), ]
  triads$triad_id <- sprintf(
    "T%03d_BMI_%s_%s", seq_len(nrow(triads)),
    sanitize_trait_id(triads$candidate), sanitize_trait_id(triads$outcome)
  )
  triads$overlapping_measure <- cfg$candidates$overlapping_measure[
    match(triads$candidate, cfg$candidates$name)
  ]
  triads <- triads[, c("triad_id", "candidate", "outcome", "overlapping_measure")]
  rownames(triads) <- NULL
  triads
}

find_bmi_gwas <- function(cfg) {
  paths <- unique(cfg$paths$bmi_raw_candidates[nzchar(cfg$paths$bmi_raw_candidates)])
  hits <- paths[file.exists(paths)]
  if (!length(hits)) {
    eraca_stop("BMI GWAS not found. Checked: ", paste(paths, collapse = "; "))
  }
  hits[[1]]
}

find_candidate_gwas <- function(trait, cfg) {
  row <- cfg$candidates[cfg$candidates$name == trait, , drop = FALSE]
  if (nrow(row) != 1L) eraca_stop("Unknown candidate trait: ", trait)
  files <- list.files(cfg$paths$candidate_dir, full.names = TRUE, pattern = "\\.bgz$")
  hits <- files[grepl(row$raw_id[[1]], basename(files), fixed = TRUE)]
  if (length(hits) != 1L) {
    eraca_stop(
      "Expected one GWAS file for ", trait, " (", row$raw_id[[1]], "); found ",
      length(hits), "."
    )
  }
  hits[[1]]
}

write_outcome_manifest_template <- function(cfg, zip_files, path) {
  template <- cfg$outcomes[, c("name", "full_name")]
  template$zip_file <- ""
  template$available_zip_files <- paste(basename(zip_files), collapse = " | ")
  utils::write.csv(template, path, row.names = FALSE, na = "")
  invisible(path)
}

resolve_outcome_manifest <- function(cfg, allow_template = TRUE) {
  manifest_path <- cfg$paths$outcome_manifest
  if (file.exists(manifest_path)) {
    manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
    needed <- c("name", "zip_file")
    if (!all(needed %in% names(manifest))) {
      eraca_stop("Outcome manifest requires columns: name, zip_file")
    }
    manifest <- merge(cfg$outcomes[, c("name", "full_name")], manifest, by = "name", all.x = TRUE)
    blank_zip <- is.na(manifest$zip_file) | !nzchar(trimws(manifest$zip_file))
    if (any(blank_zip)) {
      eraca_stop("Outcome manifest has blank zip_file entries: ", manifest_path)
    }
    manifest$zip_path <- ifelse(
      grepl("^/", manifest$zip_file), manifest$zip_file,
      file.path(cfg$paths$outcome_dir, manifest$zip_file)
    )
    missing <- manifest$zip_path[!file.exists(manifest$zip_path)]
    if (length(missing)) eraca_stop("Outcome ZIP files missing: ", paste(missing, collapse = "; "))
    return(manifest)
  }

  zip_files <- list.files(cfg$paths$outcome_dir, pattern = "\\.zip$", full.names = TRUE)
  inferred <- character(nrow(cfg$outcomes))
  for (i in seq_len(nrow(cfg$outcomes))) {
    hits <- zip_files[grepl(cfg$outcomes$zip_regex[[i]], basename(zip_files), ignore.case = TRUE)]
    if (length(hits) == 1L) inferred[[i]] <- basename(hits)
  }
  if (all(nzchar(inferred)) && !anyDuplicated(inferred)) {
    manifest <- cfg$outcomes[, c("name", "full_name")]
    manifest$zip_file <- inferred
    manifest$zip_path <- file.path(cfg$paths$outcome_dir, inferred)
    return(manifest)
  }

  if (allow_template) {
    dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
    write_outcome_manifest_template(cfg, zip_files, manifest_path)
  }
  eraca_stop(
    "Could not uniquely infer all 14 outcome ZIP files. A template was written to ",
    manifest_path, ". Fill zip_file for each outcome and rerun."
  )
}

load_variants_table <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  objects <- ls(env)
  if (!length(objects)) eraca_stop("No object found in variants RData: ", path)
  variants <- env[[objects[[1]]]]
  as.data.frame(variants, stringsAsFactors = FALSE)
}

attach_variant_metadata <- function(dat, variants, source_label) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  if (nrow(dat) != nrow(variants)) {
    eraca_stop(
      "Row mismatch between ", source_label, " and variants table: ",
      nrow(dat), " versus ", nrow(variants), "."
    )
  }
  # The UKB GWAS files already contain rsid.  The earlier preparation scripts
  # therefore appended variants1[,-1]; retain that behavior while avoiding
  # duplicate columns when the metadata table uses the same field names.
  add_names <- setdiff(names(variants), names(dat))
  if (length(add_names)) dat <- cbind(dat, variants[, add_names, drop = FALSE])
  dat
}

standardize_ukb_associations <- function(dat, trait, analysis_type) {
  required <- c("rsid", "beta", "se", "alt", "ref", "pval")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    eraca_stop("UKB GWAS for ", trait, " lacks columns: ", paste(missing, collapse = ", "))
  }
  get_col <- function(name, default = NA) if (name %in% names(dat)) dat[[name]] else rep(default, nrow(dat))
  out <- data.frame(
    SNP = as.character(dat$rsid),
    beta = suppressWarnings(as.numeric(dat$beta)),
    se = suppressWarnings(as.numeric(dat$se)),
    effect_allele = toupper(as.character(dat$alt)),
    other_allele = toupper(as.character(dat$ref)),
    eaf = suppressWarnings(as.numeric(get_col("minor_AF"))),
    pval = suppressWarnings(as.numeric(dat$pval)),
    samplesize = suppressWarnings(as.numeric(get_col("n_complete_samples"))),
    chr = as.character(get_col("chr")),
    pos = suppressWarnings(as.numeric(get_col("pos"))),
    trait = trait,
    analysis_type = analysis_type,
    stringsAsFactors = FALSE
  )
  out <- out[
    nzchar(out$SNP) & is.finite(out$beta) & is.finite(out$se) & out$se > 0 &
      is.finite(out$pval), , drop = FALSE
  ]
  out[!duplicated(out$SNP), , drop = FALSE]
}

standardize_outcome_associations <- function(dat, trait) {
  required <- c("ID", "OR", "LOG(OR)_SE", "P", "A1", "ALT", "REF")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    eraca_stop("Outcome GWAS for ", trait, " lacks columns: ", paste(missing, collapse = ", "))
  }
  get_col <- function(name, default = NA) if (name %in% names(dat)) dat[[name]] else rep(default, nrow(dat))
  odds_ratio <- suppressWarnings(as.numeric(dat$OR))
  effect <- toupper(as.character(dat$A1))
  alt <- toupper(as.character(dat$ALT))
  ref <- toupper(as.character(dat$REF))
  other <- ifelse(effect == alt, ref, alt)
  out <- data.frame(
    SNP = as.character(dat$ID),
    beta = log(odds_ratio),
    se = suppressWarnings(as.numeric(dat$`LOG(OR)_SE`)),
    effect_allele = effect,
    other_allele = other,
    eaf = suppressWarnings(as.numeric(get_col("A1_FREQ"))),
    pval = suppressWarnings(as.numeric(dat$P)),
    samplesize = suppressWarnings(as.numeric(get_col("OBS_CT"))),
    chr = as.character(get_col("#CHROM")),
    pos = suppressWarnings(as.numeric(get_col("POS"))),
    trait = trait,
    analysis_type = "binary",
    stringsAsFactors = FALSE
  )
  out <- out[
    nzchar(out$SNP) & is.finite(odds_ratio) & odds_ratio > 0 &
      is.finite(out$beta) & is.finite(out$se) & out$se > 0 & is.finite(out$pval),
    , drop = FALSE
  ]
  out[!duplicated(out$SNP), , drop = FALSE]
}

read_ukb_trait <- function(trait, cfg, target_snps = NULL, p_threshold = NULL) {
  require_eraca_packages("data.table")
  path <- if (identical(trait, "BMI")) find_bmi_gwas(cfg) else find_candidate_gwas(trait, cfg)
  dat <- data.table::fread(path, showProgress = FALSE)
  variants <- load_variants_table(cfg$paths$variants)
  dat <- attach_variant_metadata(dat, variants, basename(path))
  if (!is.null(target_snps)) dat <- dat[dat$rsid %in% target_snps, , drop = FALSE]
  if (!is.null(p_threshold)) {
    p_numeric <- suppressWarnings(as.numeric(dat$pval))
    dat <- dat[is.finite(p_numeric) & p_numeric < p_threshold, , drop = FALSE]
  }
  registry <- eraca_trait_registry(cfg)
  type <- registry$analysis_type[match(trait, registry$name)]
  standardize_ukb_associations(dat, trait, type)
}

read_outcome_zip <- function(zip_path, trait, target_snps = NULL, p_threshold = NULL) {
  require_eraca_packages("data.table")
  contents <- utils::unzip(zip_path, list = TRUE)
  members <- contents$Name[grepl("\\.csv$", contents$Name, ignore.case = TRUE)]
  if (!length(members)) eraca_stop("No CSV members found in: ", zip_path)
  pieces <- vector("list", length(members))
  temp_root <- tempfile(pattern = paste0("eraca_", sanitize_trait_id(trait), "_"))
  dir.create(temp_root, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)

  for (i in seq_along(members)) {
    utils::unzip(zip_path, files = members[[i]], exdir = temp_root)
    extracted <- file.path(temp_root, members[[i]])
    if (!file.exists(extracted)) {
      hits <- list.files(temp_root, basename(members[[i]]), recursive = TRUE, full.names = TRUE)
      if (!length(hits)) eraca_stop("Failed to extract ZIP member: ", members[[i]])
      extracted <- hits[[1]]
    }
    part <- as.data.frame(data.table::fread(extracted, showProgress = FALSE), stringsAsFactors = FALSE)
    if (!is.null(target_snps) && "ID" %in% names(part)) {
      part <- part[part$ID %in% target_snps, , drop = FALSE]
    }
    if (!is.null(p_threshold) && "P" %in% names(part)) {
      p_numeric <- suppressWarnings(as.numeric(part$P))
      part <- part[is.finite(p_numeric) & p_numeric < p_threshold, , drop = FALSE]
    }
    pieces[[i]] <- part
    unlink(extracted, force = TRUE)
  }
  combined <- data.table::rbindlist(pieces, fill = TRUE)
  standardize_outcome_associations(as.data.frame(combined), trait)
}

read_trait_associations <- function(trait, cfg, target_snps = NULL, p_threshold = NULL,
                                    outcome_manifest = NULL) {
  if (trait %in% c("BMI", cfg$candidates$name)) {
    return(read_ukb_trait(trait, cfg, target_snps, p_threshold))
  }
  manifest <- outcome_manifest %||% resolve_outcome_manifest(cfg)
  row <- manifest[manifest$name == trait, , drop = FALSE]
  if (nrow(row) != 1L) eraca_stop("Outcome is absent from manifest: ", trait)
  read_outcome_zip(row$zip_path[[1]], trait, target_snps, p_threshold)
}

eraca_ld_cache <- new.env(parent = emptyenv())

load_ld_reference_table <- function(cfg) {
  key <- gsub("[^A-Za-z0-9]", "_", cfg$paths$bfile)
  if (exists(key, envir = eraca_ld_cache, inherits = FALSE)) {
    return(get(key, envir = eraca_ld_cache, inherits = FALSE))
  }
  bim_path <- paste0(cfg$paths$bfile, ".bim")
  if (!file.exists(bim_path)) return(NULL)
  require_eraca_packages("data.table")
  bim <- data.table::fread(
    bim_path, header = FALSE, showProgress = FALSE, data.table = FALSE
  )
  if (ncol(bim) < 4L) eraca_stop("PLINK BIM file has fewer than four columns: ", bim_path)
  reference <- data.frame(
    ref_chr = as.character(bim[[1]]),
    ref_id = as.character(bim[[2]]),
    ref_pos = suppressWarnings(as.numeric(bim[[4]])),
    stringsAsFactors = FALSE
  )
  assign(key, reference, envir = eraca_ld_cache)
  reference
}

normalise_chr_for_ld <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("^CHR", "", x)
}

clump_trait_instruments <- function(dat, cfg) {
  require_eraca_packages(c("ieugwasr", "data.table"))
  dat <- dat[is.finite(dat$pval) & dat$pval < cfg$thresholds$instrument_p, , drop = FALSE]
  if (!nrow(dat)) return(dat)
  dat <- dat[order(dat$pval), , drop = FALSE]
  dat <- dat[!duplicated(dat$SNP), , drop = FALSE]
  if (nrow(dat) > 1L) {
    reference <- load_ld_reference_table(cfg)
    dat$clump_id <- dat$SNP
    if (!is.null(reference)) {
      direct_match <- dat$SNP %in% reference$ref_id
      dat_key <- paste(
        normalise_chr_for_ld(dat$chr),
        suppressWarnings(as.character(as.integer(dat$pos))),
        sep = ":"
      )
      ref_key <- paste(
        normalise_chr_for_ld(reference$ref_chr),
        suppressWarnings(as.character(as.integer(reference$ref_pos))),
        sep = ":"
      )
      coordinate_id <- reference$ref_id[match(dat_key, ref_key)]
      use_coordinate <- !direct_match & !is.na(coordinate_id) & nzchar(coordinate_id)
      dat$clump_id[use_coordinate] <- coordinate_id[use_coordinate]
      matched <- dat$clump_id %in% reference$ref_id
      if (!any(matched)) {
        eraca_stop(
          "No significant SNP for ", unique(dat$trait),
          " matched the PLINK reference by rsID or chromosome-position."
        )
      }
      if (any(!matched)) {
        eraca_warn(
          sum(!matched), " significant SNP(s) for ", unique(dat$trait),
          " were absent from the PLINK reference and were excluded before clumping."
        )
        dat <- dat[matched, , drop = FALSE]
      }
      dat <- dat[!duplicated(dat$clump_id), , drop = FALSE]
    }
    clump_input <- data.frame(rsid = dat$clump_id, pval = dat$pval, id = dat$trait)
    clumped <- try(
      ieugwasr::ld_clump_local(
        dat = clump_input,
        clump_kb = cfg$thresholds$clump_kb,
        clump_r2 = cfg$thresholds$clump_r2,
        clump_p = 1,
        bfile = cfg$paths$bfile,
        plink_bin = cfg$paths$plink
      ),
      silent = TRUE
    )
    if (inherits(clumped, "try-error") || is.null(clumped) || !nrow(clumped)) {
      eraca_stop("Local LD clumping failed for trait: ", unique(dat$trait))
    }
    dat <- dat[dat$clump_id %in% clumped$rsid, , drop = FALSE]
  }
  dat$clump_id <- NULL
  dat$F_stat <- (dat$beta / dat$se)^2
  dat
}

instrument_path <- function(trait, cfg) {
  file.path(cfg$paths$instruments, paste0(sanitize_trait_id(trait), ".rds"))
}

association_cache_path <- function(trait, cfg) {
  file.path(cfg$paths$association_cache, paste0(sanitize_trait_id(trait), ".rds"))
}

empty_instrument_table <- function(trait, analysis_type = "binary") {
  data.frame(
    SNP = character(), beta = numeric(), se = numeric(),
    effect_allele = character(), other_allele = character(), eaf = numeric(),
    pval = numeric(), samplesize = numeric(), chr = character(), pos = numeric(),
    trait = character(), analysis_type = character(), F_stat = numeric(),
    stringsAsFactors = FALSE
  )
}

load_instrument <- function(trait, cfg) {
  path <- instrument_path(trait, cfg)
  if (!file.exists(path)) eraca_stop("Instrument file missing for ", trait, ": ", path)
  readRDS(path)
}

load_association_cache <- function(trait, cfg) {
  path <- association_cache_path(trait, cfg)
  if (!file.exists(path)) eraca_stop("Association cache missing for ", trait, ": ", path)
  readRDS(path)
}

format_tsmr_data <- function(dat, type = c("exposure", "outcome"), label = NULL) {
  require_eraca_packages("TwoSampleMR")
  type <- match.arg(type)
  if (!nrow(dat)) eraca_stop("Cannot format an empty association table.")
  label <- label %||% unique(dat$trait)[[1]]
  formatted <- TwoSampleMR::format_data(
    dat,
    type = type,
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize",
    chr_col = "chr", pos_col = "pos"
  )
  formatted[[type]] <- label
  formatted[[paste0("id.", type)]] <- sanitize_trait_id(label)
  formatted
}

primary_method_for_nsnp <- function(nsnp) {
  if (!is.finite(nsnp) || nsnp < 1L) return(NA_character_)
  if (nsnp == 1L) return("Wald ratio")
  "Inverse variance weighted"
}

method_list_for_nsnp <- function(nsnp) {
  if (nsnp < 1L) return(character())
  if (nsnp == 1L) return("mr_wald_ratio")
  if (nsnp == 2L) return("mr_ivw")
  c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
}

compute_weighted_response_identity <- function(x, z, y, se_y, beta_adjusted, alpha_z) {
  values <- list(x = x, z = z, y = y, se_y = se_y)
  lengths <- vapply(values, length, integer(1))
  if (length(unique(lengths)) != 1L || lengths[[1]] < 1L) {
    eraca_stop("x, z, y and se_y must have the same positive length.")
  }
  keep <- is.finite(x) & is.finite(z) & is.finite(y) & is.finite(se_y) & se_y > 0
  if (!any(keep)) eraca_stop("No finite SNP association remains for the response identity.")
  x <- x[keep]
  z <- z[keep]
  y <- y[keep]
  w <- 1 / se_y[keep]^2
  denominator <- sum(w * x^2)
  if (!is.finite(denominator) || denominator <= 0) {
    eraca_stop("Weighted BMI denominator must be positive.")
  }
  beta_unadjusted <- sum(w * x * y) / denominator
  delta_xz <- sum(w * x * z) / denominator
  observed <- beta_unadjusted - beta_adjusted
  predicted <- alpha_z * delta_xz
  list(
    beta_unadjusted = beta_unadjusted,
    delta_xz = delta_xz,
    observed_response = observed,
    predicted_response = predicted,
    identity_error = observed - predicted
  )
}

classify_edge_evidence <- function(pval, p_fdr, b, nsnp, mean_f,
                                   secondary_concordant = FALSE,
                                   unstable = FALSE,
                                   f_min = 10,
                                   nominal_alpha = 0.05,
                                   fdr_alpha = 0.05) {
  if (!is.finite(pval) || !is.finite(b) || nsnp < 1L) return("not_estimable")
  if (isTRUE(unstable)) return("unstable")
  if (!is.finite(mean_f) || mean_f < f_min) return("weak_instrument")
  if (is.finite(p_fdr) && p_fdr < fdr_alpha && nsnp >= 3L && isTRUE(secondary_concordant)) {
    return("robust_supported")
  }
  if (pval < nominal_alpha && nsnp >= 2L) return("sensitivity_qualified")
  if (pval < nominal_alpha) return("nominal")
  "not_supported"
}

edge_is_supported <- function(status) {
  status %in% c(
    "robust_supported", "sensitivity_qualified", "nominal",
    "literature_supported", "conditional_supported"
  )
}

combine_edge_evidence <- function(mr_status, literature_status = NA_character_) {
  lit <- literature_status %||% NA_character_
  if (length(lit) == 0L || is.na(lit) || !nzchar(lit)) return(mr_status)
  lit_supported <- lit %in% c("supported", "robust_supported", "literature_supported")
  lit_against <- lit %in% c("contradicted", "against")
  mr_supported <- edge_is_supported(mr_status)
  if ((lit_supported && mr_status %in% c("unstable", "not_supported")) ||
      (lit_against && mr_supported)) return("unstable")
  if (lit_supported && mr_supported) return("robust_supported")
  if (lit_supported) return("literature_supported")
  mr_status
}

assign_role_from_edges <- function(edge_states, overlapping_measure = FALSE) {
  required_names <- c("X_to_Z", "Z_to_X", "Z_to_Y", "Y_to_Z")
  missing <- setdiff(required_names, names(edge_states))
  if (length(missing)) eraca_stop("Missing role edge states: ", paste(missing, collapse = ", "))
  if (isTRUE(overlapping_measure)) {
    return(list(role = "overlapping_measure", confidence = "unclassified",
                reason = "BMI and candidate are overlapping adiposity measures"))
  }
  present <- vapply(edge_states[required_names], edge_is_supported, logical(1))
  names(present) <- required_names
  if ((present[["X_to_Z"]] && present[["Z_to_X"]]) ||
      (present[["Z_to_Y"]] && present[["Y_to_Z"]])) {
    return(list(role = "feedback_or_cycle", confidence = "unclassified",
                reason = "At least one trait pair has supported edges in both directions"))
  }
  if (present[["X_to_Z"]] && present[["Z_to_Y"]]) {
    return(list(role = "mediator_out_of_library", confidence = "unclassified",
                reason = "The graph contains BMI -> candidate -> outcome"))
  }
  if (present[["Z_to_X"]] && present[["Z_to_Y"]]) {
    return(list(role = "confounder", confidence = "provisional",
                reason = "Candidate has supported arrows to BMI and outcome"))
  }
  if (present[["X_to_Z"]] && present[["Y_to_Z"]]) {
    return(list(role = "collider", confidence = "provisional",
                reason = "BMI and outcome have supported arrows to candidate"))
  }
  if (present[["Z_to_Y"]] && !present[["X_to_Z"]] && !present[["Z_to_X"]] &&
      !present[["Y_to_Z"]]) {
    return(list(role = "independent_cause", confidence = "provisional",
                reason = "Candidate-to-outcome edge is supported without a supported BMI-candidate edge"))
  }
  if (present[["Z_to_X"]] && !present[["Z_to_Y"]] && !present[["Y_to_Z"]]) {
    return(list(role = "upstream_surrogate_of_exposure", confidence = "provisional",
                reason = "Candidate-to-BMI edge is supported without a supported candidate-to-outcome edge"))
  }
  if (present[["Y_to_Z"]] && !present[["X_to_Z"]] && !present[["Z_to_X"]] &&
      !present[["Z_to_Y"]]) {
    return(list(role = "downstream_surrogate_of_outcome", confidence = "provisional",
                reason = "Outcome-to-candidate edge is supported without a supported BMI-candidate edge"))
  }
  if (present[["X_to_Z"]] && !present[["Z_to_Y"]] && !present[["Y_to_Z"]]) {
    return(list(role = "downstream_surrogate_of_exposure", confidence = "provisional",
                reason = "BMI-to-candidate edge is supported without a supported candidate-to-outcome edge"))
  }
  if (!any(present)) {
    return(list(role = "insufficient_evidence", confidence = "unclassified",
                reason = "No required direction edge has sufficient support"))
  }
  list(role = "ambiguous", confidence = "unclassified",
       reason = "Supported edges do not uniquely match a prespecified role template")
}

new_adjacency <- function(nodes = c("BMI", "candidate", "outcome")) {
  matrix(0L, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
}

has_directed_cycle <- function(adjacency) {
  n <- nrow(adjacency)
  state <- integer(n)
  visit <- function(v) {
    if (state[[v]] == 1L) return(TRUE)
    if (state[[v]] == 2L) return(FALSE)
    state[[v]] <<- 1L
    children <- which(adjacency[v, ] != 0)
    for (child in children) if (visit(child)) return(TRUE)
    state[[v]] <<- 2L
    FALSE
  }
  for (v in seq_len(n)) if (state[[v]] == 0L && visit(v)) return(TRUE)
  FALSE
}

select_task_subset <- function(tasks) {
  task_id <- Sys.getenv("ERACA_TASK_ID", unset = "")
  if (nzchar(task_id)) return(tasks[tasks$task_id == task_id, , drop = FALSE])
  start <- suppressWarnings(as.integer(Sys.getenv("ERACA_TASK_START", unset = "")))
  end <- suppressWarnings(as.integer(Sys.getenv("ERACA_TASK_END", unset = "")))
  if (is.finite(start) || is.finite(end)) {
    if (!is.finite(start)) start <- 1L
    if (!is.finite(end)) end <- nrow(tasks)
    start <- max(1L, start)
    end <- min(nrow(tasks), end)
    return(tasks[seq.int(start, end), , drop = FALSE])
  }
  tasks
}

read_literature_edges <- function(path) {
  if (!file.exists(path)) return(data.frame())
  dat <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("candidate", "outcome", "from", "to", "direction_status")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    eraca_stop("Literature-edge file lacks columns: ", paste(missing, collapse = ", "))
  }
  dat
}

write_csv_atomic <- function(dat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- paste0(path, ".tmp_", Sys.getpid())
  utils::write.csv(dat, temp, row.names = FALSE, na = "")
  if (!file.rename(temp, path)) {
    unlink(temp, force = TRUE)
    eraca_stop("Failed to write output atomically: ", path)
  }
  invisible(path)
}
