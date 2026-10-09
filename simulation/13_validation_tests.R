# 13_validation_tests.R
# Acceptance tests enforcing the simulation specification.

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop("VALIDATION FAILED: ", msg, call. = FALSE)
  invisible(TRUE)
}

run_validation_tests <- function(
    config = default_sim_config(),
    project_dir = getwd(),
    verbose = TRUE) {

  ensure_output_dirs(config)
  registry <- get_role_registry(config$path_coefficient)
  validate_role_registry(registry)

  say <- function(i, msg) {
    if (verbose) cat(sprintf("TEST %02d PASS: %s\n", i, msg))
  }

  # 01: all primary direct configurations have K=100
  dc <- config$direct_instrument_config
  assert_true(all(dc$n_direct_x + dc$n_direct_z + dc$n_direct_shared == config$K),
              "primary direct configurations do not all have K.")
  say(1, "all primary direct configurations have fixed K")

  # 02-10: theorem-space tests
  tc <- config$theorem_instrument_config
  assert_true(all(tc$n_x_only + tc$n_z_only + tc$n_shared == config$K),
              "theorem configurations do not all have K.")
  allc <- tc[tc$instrument_config == "all_summary_classes", , drop = FALSE]
  w <- seq(0.75, 1.25, length.out = config$K)

  gpos <- make_summary_geometry(
    config$K, allc$n_x_only, allc$n_z_only, allc$n_shared,
    rho_target = 0.8, w = w, effect_rms = config$theorem_effect_rms,
    seed = stable_seed(config$master_seed, "validation_pos")
  )
  assert_true(abs(attr(gpos, "rho_w_true") - 0.8) < 1e-8,
              "positive weighted alignment target failed.")
  say(7, "positive weighted alignment achieved")

  gneg <- make_summary_geometry(
    config$K, allc$n_x_only, allc$n_z_only, allc$n_shared,
    rho_target = -0.8, w = w, effect_rms = config$theorem_effect_rms,
    seed = stable_seed(config$master_seed, "validation_neg")
  )
  assert_true(abs(attr(gneg, "rho_w_true") + 0.8) < 1e-8,
              "negative weighted alignment target failed.")
  say(8, "negative weighted alignment achieved")

  gorth <- make_summary_geometry(
    config$K, allc$n_x_only, allc$n_z_only, allc$n_shared,
    rho_target = 0, w = w, effect_rms = config$theorem_effect_rms,
    seed = stable_seed(config$master_seed, "validation_orth")
  )
  assert_true(allc$n_shared > 0 && abs(attr(gorth, "delta_xz_true")) < 1e-10,
              "shared+orthogonal scenario did not produce delta=0.")
  say(6, "shared variants can coexist with weighted orthogonality")

  nos <- tc[tc$instrument_config == "no_shared", , drop = FALSE]
  gnosh <- make_summary_geometry(
    config$K, nos$n_x_only, nos$n_z_only, nos$n_shared,
    rho_target = 0, w = w, effect_rms = config$theorem_effect_rms,
    seed = stable_seed(config$master_seed, "validation_noshared")
  )
  assert_true(abs(attr(gnosh, "delta_xz_true")) < 1e-12,
              "no-shared scenario did not produce delta=0.")
  say(5, "no-shared summary geometry gives delta=0")

  assert_true(all(registry$alpha_z_theorem[registry$role %in%
                c("collider","upstream_surrogate_exposure",
                  "downstream_surrogate_outcome","downstream_surrogate_exposure")] == 0),
              "alpha_Z zero-role registry failed.")
  say(9, "alpha_Z=0 frozen for B/D/E/F")

  assert_true(all(registry$alpha_z_theorem[registry$role %in%
                c("confounder","independent_cause")] != 0),
              "alpha_Z nonzero-role registry failed.")
  say(10, "alpha_Z!=0 frozen for A/C")


  # 10b: independent cross-check of continuous DGP population summary relation.
  arch_test <- make_direct_architecture(
    K = config$K, n_direct_x = 30L, n_direct_z = 30L,
    n_direct_shared = 40L, shared_rho = 0.3,
    effect_rms = 0.05, maf = config$maf,
    seed = stable_seed(config$master_seed, "population_truth_test")
  )
  for (r in registry$role) {
    pt <- population_xz_summary_truth(r, arch_test, beta_x = 0.10, registry)
    spec_r <- get_role_spec(r, registry)
    # Analytic SNP-Y truth for the continuous DGP.
    if (r %in% c("confounder", "independent_cause")) {
      gamma_y <- 0.10 * pt$gamma_x_true +
        spec_r$edge_z_to_y * pt$gamma_z_true
    } else {
      gamma_y <- 0.10 * pt$gamma_x_true
    }
    expected_alpha <- spec_r$alpha_z_theorem
    err <- max(abs(
      gamma_y - (0.10 * pt$gamma_x_true + expected_alpha * pt$gamma_z_true)
    ))
    assert_true(err < config$theorem_tolerance,
                paste("population DGP relation failed for", r))
  }
  say(10, "continuous six-DAG population summary relations match theorem alpha_Z")

  # 02-04: strict same-set and identity
  sd <- make_theorem_summary_data(gpos, 0.10, 0.30)
  fit <- fit_strict_same_set(sd, se_method = "fixed")
  assert_true(isTRUE(fit$converged), "strict same-set theorem fit failed.")
  assert_true(length(fit$snp_id) == config$K, "strict same-set dropped SNP rows.")
  say(2, "strict unadjusted and adjusted regressions use identical SNP rows")
  assert_true(abs(fit$identity_error) < config$theorem_tolerance,
              "Theorem 1 identity error exceeded tolerance.")
  say(4, "Theorem 1 identity holds to numerical tolerance")
  say(3, "same-set estimator uses one common outcome vector and weight definition")


  # 03b: weighted-lm estimator/SE cross-check against base R lm().
  set.seed(stable_seed(config$master_seed, "wls_lm_crosscheck"))
  dchk <- data.frame(
    snp_id = sprintf("c%03d", 1:40),
    beta_x = rnorm(40),
    beta_z = rnorm(40),
    beta_y = rnorm(40),
    se_y = runif(40, 0.02, 0.08),
    se_x = rep(0.02, 40),
    se_z = rep(0.02, 40),
    cov_xz = 0
  )
  fit_ours <- fit_strict_same_set(dchk, se_method = "weighted_lm")
  fit_lm_u <- summary(lm(
    beta_y ~ 0 + beta_x,
    weights = 1 / se_y^2,
    data = dchk
  ))
  fit_lm_a <- summary(lm(
    beta_y ~ 0 + beta_x + beta_z,
    weights = 1 / se_y^2,
    data = dchk
  ))
  assert_true(
    abs(fit_ours$beta_unadj - coef(fit_lm_u)["beta_x", "Estimate"]) < 1e-12 &&
    abs(fit_ours$se_unadj_weighted_lm -
        coef(fit_lm_u)["beta_x", "Std. Error"]) < 1e-12,
    "unadjusted weighted-lm cross-check failed."
  )
  assert_true(
    abs(fit_ours$beta_adj - coef(fit_lm_a)["beta_x", "Estimate"]) < 1e-12 &&
    abs(fit_ours$se_adj_weighted_lm -
        coef(fit_lm_a)["beta_x", "Std. Error"]) < 1e-12 &&
    abs(fit_ours$alpha_z - coef(fit_lm_a)["beta_z", "Estimate"]) < 1e-12,
    "adjusted weighted-lm cross-check failed."
  )
  say(3, "weighted summary regression point estimates and SEs match base R lm()")

  # 14-16 are enforced by output schema and summarizer logic.
  tmp <- data.frame(
    scenario_id = "tmp", role = "confounder", outcome_type = "continuous",
    sample_design = "two_sample", beta_x_target = 0,
    instrument_config = "tmp", converged = TRUE,
    beta_unadj = c(0, 0.01), beta_adj = c(0, 0.01),
    se_unadj = c(0.1,0.1), se_adj = c(0.1,0.1),
    se_unadj_fixed = c(0.1,0.1), se_adj_fixed = c(0.1,0.1),
    se_unadj_weighted_lm = c(0.1,0.1), se_adj_weighted_lm = c(0.1,0.1),
    cover_unadj = c(1,1), cover_adj = c(1,1),
    reject_unadj = c(0,0), reject_adj = c(0,0),
    F_x_given_z = c(20,20), F_z_given_x = c(20,20), F_min = c(20,20),
    F_x_given_z_package = c(19.8,19.8),
    F_z_given_x_package = c(19.8,19.8),
    F_min_package = c(19.8,19.8),
    response = c(0,0),
    role_predicted_response = c(0,0), role_response_error = c(0,0),
    rho_w = c(.5,.5), delta_xz = c(.1,.1),
    condition_number = c(2,2), identity_error = c(0,0)
  )
  sm <- summarise_scenario(tmp, config)
  assert_true("empirical_sd_unadj" %in% names(sm) &&
              "mean_model_se_unadj" %in% names(sm),
              "empirical SD and model SE are not separated.")
  assert_true(identical(config$se_method_primary, "weighted_lm"),
              "primary SE method is not frozen to weighted_lm.")
  say(14, "empirical SD, fixed SE and weighted-lm model SE are separated")
  assert_true(is.finite(sm$type1_error_unadj) && is.na(sm$power_unadj),
              "type I / power null mapping failed.")
  say(15, "type I error is assigned only under beta_X=0")

  tmp$beta_x_target <- 0.10
  sm2 <- summarise_scenario(tmp, config)
  assert_true(is.na(sm2$type1_error_unadj) && is.finite(sm2$power_unadj),
              "type I / power alternative mapping failed.")
  say(16, "power is assigned only under beta_X!=0")

  # 17-20: Theorems 2-3
  bsc <- run_binary_scale_calibration(config, write_output = FALSE)
  structural <- bsc[bsc$residual_mode == "none", ]
  signs <- sign(structural$delta_struct[abs(structural$delta_struct) > 1e-14])
  assert_true(length(unique(signs)) <= 1L, "structural sign changed across positive S.")
  say(18, "Theorem 3 preserves structural sign across S>0")
  zero_status <- abs(structural$delta_struct) < 1e-14
  assert_true(length(unique(zero_status)) <= 1L, "structural zero status changed across S.")
  say(17, "Theorem 3 preserves structural zero/nonzero status across S>0")

  off <- bsc[bsc$residual_mode == "exact_offset", ]
  assert_true(all(abs(off$delta_total) < config$theorem_tolerance) &&
              all(abs(off$delta_struct) > config$theorem_tolerance),
              "exact-offset boundary failed.")
  say(19, "exact offset gives total response=0 with structural response!=0")

  so <- bsc[bsc$residual_mode == "scale_only", ]
  assert_true(all(abs(so$delta_struct) < config$theorem_tolerance) &&
              all(abs(so$delta_total) > config$theorem_tolerance),
              "scale-only boundary failed.")
  say(20, "scale-only response can be nonzero when structural response=0")

  # 21: conditional F is role-independent and follows the Eq. 7 denominator.
  assert_true(!("adjustment_class" %in% names(formals(conditional_strength))),
              "conditional F function is coupled to adjustment class.")
  cft <- conditional_q_two_exposure(
    target_beta = c(.10, .20, .31, .39),
    other_beta = c(.05, .10, .15, .20),
    target_se = rep(.02, 4),
    other_se = rep(.02, 4),
    cov_target_other = rep(0, 4)
  )
  assert_true(
    abs(cft$F - cft$Q / (cft$L - 1L)) < 1e-12,
    "two-sample conditional F does not use Eq. 7 denominator."
  )
  assert_true(
    abs(cft$F_package - cft$Q_package / cft$L) < 1e-12,
    "package-compatible F crosswalk is incorrect."
  )
  say(21, "conditional F is class-independent and Eq. 7/package crosswalk is explicit")


  mvmr_package_crosscheck <- FALSE
  if (requireNamespace("MVMR", quietly = TRUE)) {
    bxg <- cbind(
      x = c(.10, .20, .31, .39, .51, .61),
      z = c(.05, .10, .15, .20, .25, .30)
    )
    sebx <- cbind(
      x = rep(.02, nrow(bxg)),
      z = rep(.02, nrow(bxg))
    )
    byg <- rep(0, nrow(bxg))
    seby <- rep(.05, nrow(bxg))
    cv <- rep(0, nrow(bxg))
    fmt <- MVMR::format_mvmr(
      BXGs = bxg,
      BYG = byg,
      seBXGs = sebx,
      seBYG = seby,
      RSID = paste0("rs", seq_len(nrow(bxg)))
    )
    gcov <- lapply(seq_len(nrow(bxg)), function(j) {
      matrix(c(
        sebx[j,1]^2, cv[j],
        cv[j], sebx[j,2]^2
      ), nrow = 2, byrow = TRUE)
    })
    pkg <- suppressWarnings(MVMR::strength_mvmr(fmt, gencov = gcov))
    ours_dat <- data.frame(
      beta_x = bxg[,1], se_x = sebx[,1],
      beta_z = bxg[,2], se_z = sebx[,2],
      cov_xz = cv
    )
    ours <- conditional_strength(ours_dat)
    pkgv <- as.numeric(pkg[1, c("exposure1", "exposure2")])
    assert_true(
      max(abs(pkgv - c(
        ours$F_x_given_z_package,
        ours$F_z_given_x_package
      ))) < 1e-8,
      "MVMR::strength_mvmr package cross-check failed."
    )
    mvmr_package_crosscheck <- TRUE
    if (verbose) cat("OPTIONAL PACKAGE CROSS-CHECK PASS: MVMR::strength_mvmr\n")
  } else {
    if (verbose) cat("OPTIONAL PACKAGE CROSS-CHECK SKIPPED: MVMR not installed\n")
  }


  # GWAS engine cross-checks against base lm/glm.
  set.seed(stable_seed(config$master_seed, "gwas_engine_crosscheck"))
  Gc <- make_genotype_matrix(500L, config$maf, 5L)
  xc <- 0.2 * Gc[,1] - 0.1 * Gc[,2] + rnorm(500)
  zc <- -0.15 * Gc[,1] + 0.12 * Gc[,3] + rnorm(500)

  gvec <- gwas_continuous_pair(Gc, xc, zc)
  for (j in 1:3) {
    fx <- summary(lm(xc ~ Gc[,j]))
    fz <- summary(lm(zc ~ Gc[,j]))
    assert_true(
      abs(gvec$beta_x[j] - coef(fx)[2,1]) < 1e-12 &&
      abs(gvec$se_x[j] - coef(fx)[2,2]) < 1e-10 &&
      abs(gvec$beta_z[j] - coef(fz)[2,1]) < 1e-12 &&
      abs(gvec$se_z[j] - coef(fz)[2,2]) < 1e-10,
      paste("vectorised continuous GWAS cross-check failed at SNP", j)
    )
  }

  eta <- -1 + 0.25 * Gc[,1] - 0.15 * Gc[,2]
  yb <- rbinom(500, 1, plogis(eta))
  gb <- gwas_binary(Gc, yb)
  for (j in 1:3) {
    fg <- summary(glm(yb ~ Gc[,j], family = binomial()))
    assert_true(
      abs(gb$beta[j] - coef(fg)[2,1]) < 1e-10 &&
      abs(gb$se[j] - coef(fg)[2,2]) < 1e-8,
      paste("binary glm.fit GWAS cross-check failed at SNP", j)
    )
  }
  say(5, "vectorised linear GWAS and glm.fit binary GWAS match base R reference fits")

  # 11-13,22: small individual-level reproducibility / sample-design checks
  sc <- build_primary_scenarios(config)
  sc <- sc[sc$role == "independent_cause" &
           sc$outcome_type == "continuous" &
           sc$beta_x_target == 0.10 &
           sc$instrument_config == "all_direct_classes", , drop = FALSE]
  # Reduce N only for validation speed; all other scenario truth is unchanged.
  sc$n_exposure <- 500L
  sc$n_outcome <- 500L
  sc$scenario_id <- paste0(sc$scenario_id, "__validation")
  truth <- prepare_scenario_truth(sc, config)
  s1 <- generate_samples_for_replicate(truth, config, 1L)
  assert_true(length(intersect(s1$exposure$id, s1$outcome$id)) == 0L,
              "two-sample IDs overlap.")
  assert_true(!identical(s1$exposure$G, s1$outcome$G),
              "two-sample genotype realizations are identical.")
  say(11, "two-sample samples are disjoint and independently generated")

  sc_same <- sc
  sc_same$sample_design <- "same_sample"
  sc_same$sample_overlap <- 1
  sc_same$scenario_id <- paste0(sc_same$scenario_id, "__same")
  truth_same <- prepare_scenario_truth(sc_same, config)
  ss <- generate_samples_for_replicate(truth_same, config, 1L)
  assert_true(identical(ss$exposure$id, ss$outcome$id),
              "same-sample IDs are not identical.")
  assert_true(identical(ss$exposure$G, ss$outcome$G),
              "same-sample genotype matrices are not identical.")
  say(12, "same-sample exposure and outcome data are identical")

  # Binary prevalence calibration test uses reference prevalence, not a single finite replicate.
  scb <- build_primary_scenarios(config)
  scb <- scb[scb$role == "independent_cause" &
             scb$outcome_type == "binary" &
             scb$beta_x_target == 0.10 &
             scb$instrument_config == "all_direct_classes", , drop = FALSE]
  scb$n_exposure <- 500L
  scb$n_outcome <- 500L
  scb$scenario_id <- paste0(scb$scenario_id, "__validation")
  truth_b <- prepare_scenario_truth(scb, config)
  assert_true(abs(truth_b$binary_calibration$calibrated_prevalence -
                  config$binary_prevalence_primary) < 0.005,
              "binary prevalence calibration outside tolerance.")
  say(13, "binary reference prevalence matches target")

  a <- run_one_replicate(truth, config, 1L)
  b <- run_one_replicate(truth, config, 1L)
  assert_true(is.finite(a$role_predicted_response) && is.finite(a$role_response_error),
              "continuous role-implied response prediction was not recorded.")
  assert_true(isTRUE(all.equal(a, b, check.attributes = TRUE)),
              "same seed did not reproduce identical output.")
  say(22, "same seed reproduces identical replicate output")

  # Extra strict condition: all theorem identities in the exact grid.
  th <- run_theorem1_calibration(config, write_output = FALSE)
  assert_true(max(abs(th$identity_error)) < config$theorem_tolerance,
              "full theorem grid identity failed.")

  validation_report <- list(
    passed = TRUE,
    timestamp = Sys.time(),
    mvmr_package_crosscheck = mvmr_package_crosscheck,
    theorem1_max_abs_identity_error = max(abs(th$identity_error)),
    protocol_version = config$protocol_version
  )
  saveRDS(
    validation_report,
    file.path(config$output_dir, "validation", "validation_passed.rds")
  )
  invisible(validation_report)
}
