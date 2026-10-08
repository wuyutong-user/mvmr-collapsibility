> Historical simulation audit. See the top-level README and public-release verification for current checks.

# Cross-validation audit — MVMR simulation v1.1

## Scope

This audit checks the implementation against the frozen simulation protocol, the six Fig.1 structures, and the current manuscript theory. It separates checks that can be completed without an R runtime from checks that must run on the target server.

## 1. Theory-to-code mapping

### Theorem 1

**Requirement:** unadjusted and Z-included no-intercept weighted summary regressions must use the same instrument rows and weighting scheme; the exact response identity is

`beta_unadj - beta_adj = alpha_Z * delta_XZ`.

**Implementation:**
- `06_ivw_estimators.R::fit_strict_same_set()` fits both regressions to identical SNP IDs and a common outcome-based weight vector.
- `10_run_six_dag_primary.R` requires all K=100 prespecified SNPs to have finite GWAS estimates; a failed binary SNP fit fails the replicate rather than silently changing K.
- every finite-sample replicate stores the algebraic response, `alpha_Z * delta_XZ`, and the identity error.
- `08_theorem_calibration.R` performs exact population-summary calibration.

**Independent numerical cross-check:** maximum absolute identity error = **4.9960e-16** in an independent Python implementation, versus protocol tolerance 1e-10.

### Role-implied response versus algebraic identity

The code now explicitly distinguishes two concepts that must not be conflated:

1. **algebraic identity:** fitted `alpha_Z_hat * delta_hat`, which must equal the fitted coefficient response by linear algebra;
2. **role-implied prediction:** for continuous outcomes, the Fig.1 DGP supplies the population `alpha_Z` and population SNP-X/SNP-Z associations, so the code separately computes the theorem-implied structural response from the prespecified role.

The latter is used to assess whether the six-DAG simulation reproduces the role prediction; the former is an implementation identity check.

### Six Fig.1 roles

`01_role_registry.R` freezes exactly six roles. The continuous population DGP was independently checked for the relation

`gamma_Y = beta_X * gamma_X + alpha_Z * gamma_Z`,

with `alpha_Z=0.30` for confounder and independent cause and `alpha_Z=0` for collider and the three surrogate roles.

**Independent cross-check:** maximum algebraic residual across all six roles = **0**.

### Theorems 2 and 3

`09_binary_scale_calibration.R` directly constructs

`ell = S^-1(beta_X* x + alpha_Z* z) + r_scale`

and evaluates scale factors 0.75, 1.00 and 1.50 under no residual, reinforcing residual, partial offset, exact offset and scale-only cases.

**Independent numerical cross-check:** maximum decomposition error = **4.9960e-16**. Exact-offset cases had total response below 5e-16 while the structural response was non-zero; scale-only cases had structural response zero and total response 0.05.

## 2. Instrument definitions and geometry

The code separates:
- `direct_class`: direct genetic source in the individual-level DGP;
- `summary_class`: population SNP-X/SNP-Z association support used by the theorem.

This prevents a direct-to-Z SNP from being incorrectly labelled Z-only when a pathway such as Z→X makes it summary-level shared.

Exact theorem-space geometry is constructed with weighted Gram-Schmidt operations rather than relying on random sign cancellation.

**Independent numerical cross-check:** across the frozen 30/30/40, 0/60/40, 60/0/40 and 50/50/0 configurations and target weighted correlations +0.80, 0, -0.80 and +0.95 where feasible, maximum absolute correlation-target error = **3.3307e-16**. No-shared and shared-but-orthogonal designs produced zero weighted projection to floating-point precision.

## 3. Sample design

Primary finite-sample simulations are truly non-overlapping two-sample analyses:
- one exposure GWAS sample generates SNP-X and SNP-Z estimates;
- a separately generated outcome GWAS sample generates SNP-Y estimates.

The 100% same-sample sensitivity reuses the exact exposure population as the outcome population. Common-random-number keys intentionally pair same-sample, two-sample, beta, strength, alignment, N and prevalence contrasts so that controlled comparisons are not confounded by unrelated random architecture changes.

The genetic architecture seed is fixed by instrument configuration and does not change when only RMS strength or alignment is changed.

## 4. Binary DGP

The common latent factor U explicitly enters the logistic outcome linear predictor in every relevant role. This corrects the legacy binary code, which generated U for X/Z but omitted it from `plogis()`.

The scenario-level intercept is calibrated to the target prevalence using a large reference population and is then frozen across Monte Carlo replicates. The realistic binary GWAS uses SNP-wise logistic regression and no rare-outcome approximation.

The realistic binary analysis is labelled as target-relative total-logOR behaviour; it does not claim that structural and scale components can be separately recovered from marginal GWAS statistics.

## 5. Precision and inference

The implementation stores separately:
- empirical Monte Carlo variance;
- empirical Monte Carlo SD;
- replicate-level model SE;
- fixed-effect SE;
- weighted-lm SE;
- variance ratio `Var(adj)/Var(unadj)`;
- relative efficiency;
- RMSE;
- 95% coverage;
- Monte Carlo SE of mean estimates/bias.

This corrects the legacy code, where replicate SD was named `sebeta`.

Empirical type I error is reported only when beta_X=0; empirical power is reported only when beta_X is non-zero. Wilson 95% Monte Carlo intervals are stored for rejection proportions.

## 6. Weighted-regression implementation

The primary finite-sample SE is the weighted-lm SE corresponding to `summary(lm(..., weights=1/seY^2))`; fixed-effect SE is retained separately. The primary p-value uses a normal reference.

**Independent Python/statsmodels cross-check:** maximum difference in WLS coefficients or weighted-lm SE = **2.2204e-16**.

The server validation suite additionally compares the custom R weighted-regression implementation directly with base R `lm`.

## 7. Conditional F

The primary conditional-strength statistic implements Sanderson, Spiller & Bowden (2021), Eq.7:

`F_TS,k = Q_xk / [L-(K-1)]`.

For two exposures, the denominator is L-1. The nuisance coefficient is obtained by no-intercept IVW, and the SNP-X/SNP-Z covariance is retained because X and Z are estimated in the same exposure GWAS sample.

For crosswalk/reproducibility, the code also returns the current public `MVMR::strength_mvmr` convention: no-intercept OLS nuisance coefficient and Q/L. If MVMR is installed, server preflight compares the package output numerically with the implementation.

Conditional F has no route to the adjustment-class registry and is used only for reliability/estimability.

## 8. 100-core server execution audit

Version 1.1 is revised for the user's 100-core server:
- one persistent PSOCK cluster is created and reused across scenarios;
- replicate tasks use load-balanced scheduling;
- worker BLAS/OpenMP thread counts are forced to one, preventing nested 100×oversubscription;
- replicate seeds are deterministic and independent of scheduling order;
- output is checkpointed in chunks (default 100);
- reruns with `SIM_RESUME=1` validate and reuse complete checkpoints;
- final scenario RDS files are written atomically.

The 100-core implementation assumes all requested cores are on one node. A multi-node scheduler allocation would require an explicit host list.

## 9. Static source audit completed here

At the time of this historical audit, R/Rscript was unavailable
in the audit environment. The following checks were performed there. Later local
R checks are recorded separately under `verification/`; this historical record
does not describe the current release environment:

- **21/21 R files:** balanced strings, parentheses, brackets and braces passed a separate static scanner;
- **server shell script:** `bash -n` passed;
- exact weighted-geometry check: passed;
- Theorem 1 identity check: passed;
- six-DAG continuous population relation: passed;
- Theorems 2–3 decomposition/boundary checks: passed;
- weighted WLS coefficient/SE cross-check against statsmodels: passed.

The numerical results are also saved in `validation_python_crosscheck.json`.

## 10. Mandatory R-side publication gate

A literal guarantee that R code will execute on the user's specific server cannot be made without executing that server's R interpreter and installed package versions. To close that gap, production is programmatically gated by `97_preflight.R`.

Before the full run it must:
1. verify the frozen source MD5 manifest;
2. parse every R file with the actual server R parser;
3. run all formal validation tests;
4. compare custom WLS against base R `lm`;
5. compare vectorised linear GWAS against base R `lm`;
6. compare binary `glm.fit` against base R `glm`;
7. compare the MVMR-package conditional-F crosswalk when `REQUIRE_MVMR=1`;
8. run exact Theorem 1 and Theorems 2–3 calibrations;
9. run continuous and binary finite-sample smoke tests;
10. record the exact frozen strength calibration used.

`99_run_all.R` refuses to start production if preflight was not rerun after the frozen strength pilot.

## 11. Manuscript consistency consequence

The revised code follows the planned revised Methods, not the legacy simulation prose. After final results are generated, the main-text simulation section and Supplementary simulation methods must be updated to describe:
- strict same-set primary analysis;
- non-overlapping two-sample primary simulation;
- fixed K=100 instrument configurations;
- direct-source versus summary-association classes;
- prevalence-calibrated binary DGP with U in the logit;
- precision/variance/coverage outputs;
- Eq.7 conditional F plus package crosswalk;
- exact Theorem 2–3 common-scale calibration;
- sample-overlap, geometry, N, prevalence, beta=0.05 and strength sensitivities.

The current Supplementary text still describes the previous same-sample design, old instrument selection and zero logistic intercept, so it should not be retained unchanged after this revised simulation is adopted.

## Overall audit judgement

**Design/theory consistency:** PASS after v1.1 revisions.  
**Independent mathematical/numerical cross-validation:** PASS.  
**Static source/shell audit:** PASS.  
**Server R-runtime validation:** PENDING until `97_preflight.R` returns `PREFLIGHT PASS` on the target server.  
**Production authorization:** only after the post-pilot preflight passes.
