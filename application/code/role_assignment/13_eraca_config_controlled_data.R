#!/usr/bin/env Rscript

eraca_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

eraca_env_flag <- function(name, default = FALSE) {
  value <- tolower(eraca_env(name, if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

eraca_config_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (!is.null(candidate) && grepl("13_eraca_config_controlled_data\\.R$", candidate)) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  normalizePath("13_eraca_config_controlled_data.R", mustWork = FALSE)
}

eraca_candidates <- function() {
  data.frame(
    name = c(
      "WC", "HC", "SBP", "DBP", "TC", "HDL-C", "LDL-C", "TG",
      "Glucose", "HbA1c", "Smoking", "Alcohol"
    ),
    raw_id = c(
      "48_irnt", "49_irnt", "4080_irnt", "4079_irnt", "30690_irnt",
      "30760_irnt", "30780_irnt", "30870_irnt", "30740_irnt",
      "30750_irnt", "20116_2", "20117_2"
    ),
    source_kind = rep("ukb_bgz", 12),
    analysis_type = c(rep("continuous", 10), "binary_liability", "binary_liability"),
    overlapping_measure = c(TRUE, TRUE, rep(FALSE, 10)),
    stringsAsFactors = FALSE
  )
}

eraca_outcomes <- function() {
  data.frame(
    name = c(
      "CAD", "HF", "IS", "TIA", "AAA", "TAA", "PVD", "PE", "SAH",
      "ICH", "DVT", "AVS", "AF", "HTN"
    ),
    full_name = c(
      "Coronary artery disease", "Heart failure", "Ischemic stroke",
      "Transient ischemic attack", "Abdominal aortic aneurysm",
      "Thoracic aortic aneurysm", "Peripheral vascular disease",
      "Pulmonary embolism", "Subarachnoid hemorrhage",
      "Intracerebral hemorrhage", "Deep vein thrombosis",
      "Aortic valve stenosis", "Atrial fibrillation", "Arterial hypertension"
    ),
    zip_regex = c(
      "(^|[^A-Za-z])CAD([^A-Za-z]|$)|coronary.*artery",
      "(^|[^A-Za-z])HF([^A-Za-z]|$)|heart.*failure",
      "(^|[^A-Za-z])IS([^A-Za-z]|$)|ischemic.*stroke",
      "(^|[^A-Za-z])TIA([^A-Za-z]|$)|transient.*ischemic",
      "(^|[^A-Za-z])AAA([^A-Za-z]|$)|abdominal.*aortic",
      "(^|[^A-Za-z])TAA([^A-Za-z]|$)|thoracic.*aortic",
      "(^|[^A-Za-z])PVD([^A-Za-z]|$)|peripheral.*vascular",
      "(^|[^A-Za-z])PE([^A-Za-z]|$)|pulmonary.*embolism",
      "(^|[^A-Za-z])SAH([^A-Za-z]|$)|subarachnoid.*hemorrhage",
      "(^|[^A-Za-z])ICH([^A-Za-z]|$)|intracerebral.*hemorrhage",
      "(^|[^A-Za-z])DVT([^A-Za-z]|$)|deep.*vein.*thrombosis",
      "(^|[^A-Za-z])AVS([^A-Za-z]|$)|aortic.*valve.*stenosis",
      "(^|[^A-Za-z])AF([^A-Za-z]|$)|atrial.*fibrillation",
      "(^|[^A-Za-z])HTN([^A-Za-z]|$)|hypertension"
    ),
    source_kind = rep("outcome_zip", 14),
    analysis_type = rep("binary", 14),
    stringsAsFactors = FALSE
  )
}

eraca_config <- function() {
  config_path <- eraca_config_file()
  code_dir <- dirname(config_path)
  package_root_default <- normalizePath(file.path(code_dir, "..", ".."), mustWork = FALSE)
  package_root <- eraca_env("ERACA_PACKAGE_ROOT", package_root_default)
  controlled_root <- eraca_env("ERACA_CONTROLLED_ROOT", file.path(package_root, "controlled_data"))
  output_dir <- eraca_env("ERACA_OUTPUT_DIR", file.path(controlled_root, "ERACA_RESULT"))

  list(
    version = "0.1.0",
    exposure = data.frame(
      name = "BMI", raw_id = "21001_irnt", source_kind = "ukb_bgz",
      analysis_type = "continuous", stringsAsFactors = FALSE
    ),
    candidates = eraca_candidates(),
    outcomes = eraca_outcomes(),
    paths = list(
      package_root = package_root,
      code_dir = code_dir,
      controlled_root = controlled_root,
      candidate_dir = eraca_env("ERACA_CANDIDATE_DIR", file.path(controlled_root, "cov_data")),
      bmi_raw_candidates = c(
        eraca_env("ERACA_BMI_GWAS", ""),
        file.path(controlled_root, "exp_data", "21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz"),
        file.path(controlled_root, "cov_data", "21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")
      ),
      variants = eraca_env("ERACA_VARIANTS_RDATA", file.path(controlled_root, "variants1.Rdata")),
      outcome_dir = eraca_env("ERACA_OUTCOME_DIR", file.path(controlled_root, "outcome_gwas")),
      outcome_manifest = eraca_env(
        "ERACA_OUTCOME_MANIFEST",
        file.path(controlled_root, "eraca_outcome_manifest.csv")
      ),
      plink = eraca_env("ERACA_PLINK", unname(Sys.which("plink"))),
      bfile = eraca_env(
        "ERACA_LD_BFILE",
        file.path(controlled_root, "ld_reference", "reference")
      ),
      output_dir = output_dir,
      instruments = file.path(output_dir, "instruments"),
      association_cache = file.path(output_dir, "association_cache"),
      uvmr_checkpoints = file.path(output_dir, "uvmr_checkpoints"),
      graphs = file.path(output_dir, "graphs"),
      logs = file.path(output_dir, "logs"),
      literature_edges = eraca_env(
        "ERACA_LITERATURE_EDGES",
        file.path(code_dir, "literature_edge_evidence.csv")
      ),
      application_response = eraca_env(
        "ERACA_APPLICATION_RESPONSE",
        file.path(package_root, "data", "application_derived", "Fig6_source_data_MR.csv")
      )
    ),
    thresholds = list(
      instrument_p = as.numeric(eraca_env("ERACA_INSTRUMENT_P", "5e-8")),
      clump_kb = as.numeric(eraca_env("ERACA_CLUMP_KB", "10000")),
      clump_r2 = as.numeric(eraca_env("ERACA_CLUMP_R2", "0.001")),
      f_min = as.numeric(eraca_env("ERACA_F_MIN", "10")),
      nominal_alpha = as.numeric(eraca_env("ERACA_NOMINAL_ALPHA", "0.05")),
      fdr_alpha = as.numeric(eraca_env("ERACA_FDR_ALPHA", "0.05")),
      response_identity_tolerance = as.numeric(
        eraca_env("ERACA_IDENTITY_TOLERANCE", "1e-8")
      )
    ),
    mrsl = list(
      enabled = eraca_env_flag("ERACA_ENABLE_MRSL", TRUE),
      cutoff = as.numeric(eraca_env("ERACA_MRSL_CUTOFF", "0.05")),
      adj_methods = as.integer(eraca_env("ERACA_MRSL_ADJ_METHOD", "1")),
      use_eggers_step2 = as.integer(eraca_env("ERACA_MRSL_EGGER", "0")),
      vary_mvmr_adj = as.integer(eraca_env("ERACA_MRSL_VARY_ADJ", "0"))
    )
  )
}
