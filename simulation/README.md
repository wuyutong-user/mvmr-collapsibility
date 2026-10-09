# MVMR Simulation Protocol v1.1

This project implements the revised simulation study for the six-role summary-data MVMR collapsibility/adjustment framework. It is designed for a one-node server with up to 100 CPU cores and uses a persistent PSOCK cluster, deterministic checkpoints and reproducible replicate-specific seeds.

## What is separated by design

The code intentionally separates three scientific tasks.

1. **Exact Theorem 1 calibration (`08_theorem_calibration.R`)**
   - summary-association level;
   - strict identical SNP set, outcome vector and weighting matrix;
   - exact check of `beta_unadj - beta_adj = alpha_Z * delta_XZ`;
   - exact positive, negative and orthogonal weighted association geometries.

2. **Six-DAG finite-sample simulation (`10_run_six_dag_primary.R`)**
   - six Fig.1 roles;
   - primary non-overlapping two-sample design;
   - strict same-set IVW/MVMR-IVW;
   - response, bias/target-relative error, empirical variance/SD, model SE, variance ratio, RMSE, coverage, empirical type I error, empirical power, conditional F and conditioning diagnostics;
   - practical different-instrument-set estimates retained separately within the same replicate.

3. **Exact binary structural+scale calibration (`09_binary_scale_calibration.R`)**
   - directly constructs `ell = S^-1(beta_X* x + alpha_Z* z) + r_scale`;
   - checks the structural-plus-scale decomposition;
   - checks positive-scale sign and zero/non-zero invariance;
   - includes exact-offset and scale-only boundary cases.

## Critical terminology

`direct_class` and `summary_class` are not interchangeable.

- `direct_class`: direct genetic source in the individual-level DGP (`direct_X`, `direct_Z`, `direct_shared`).
- `summary_class`: population summary-association class used by the theorem (`X_only`, `Z_only`, `shared`).

For example, a direct-to-Z SNP can become summary-level shared when the DAG contains `Z -> X`. The code records these concepts separately.

## Six roles

- `confounder`
- `collider`
- `independent_cause`
- `upstream_surrogate_exposure`
- `downstream_surrogate_outcome`
- `downstream_surrogate_exposure`

The role is always prespecified. No simulation function assigns a role from coefficient movement.

## Primary sample design

- Exposure GWAS sample: N=10,000; estimates SNP-X and SNP-Z.
- Outcome GWAS sample: independent N=10,000; estimates SNP-Y.
- Primary overlap: 0%.
- Same-sample sensitivity: 100% overlap.
- X and Z share the exposure GWAS sample, so the SNP-X/SNP-Z sampling covariance is retained for conditional-F calculations.

## Primary instrument configurations

Total K remains 100 in every primary configuration.

| configuration | direct-X | direct-Z | direct-shared |
|---|---:|---:|---:|
| all_direct_classes | 30 | 30 | 40 |
| no_direct_x | 0 | 60 | 40 |
| no_direct_z | 60 | 0 | 40 |
| no_direct_shared | 50 | 50 | 0 |

The exact theorem module has the corresponding population summary-association configurations.

## Weighted regression and SE convention

The point estimator is the no-intercept weighted summary regression required by the theorem. For finite-sample inference, the primary model-based SE is the weighted-lm SE returned by the same algebra as `summary(lm(..., weights=1/seY^2))`. Fixed-effect SEs are retained as additional diagnostics. The primary p-value uses the normal reference, matching the current `TwoSampleMR::mv_multiple` convention.

Empirical Monte Carlo SD is a different quantity and is stored separately.

## Conditional F

The primary statistic follows Sanderson, Spiller & Bowden (2021), Eq.7:

`F_TS,k = Q_xk / [L-(K-1)]`.

With two exposures this is `Q/(L-1)`. The nuisance coefficient is estimated by no-intercept IVW, and SNP-X/SNP-Z sampling covariance is included.

For auditability, the code also returns a package-compatible crosswalk reproducing the current public `MVMR::strength_mvmr` convention (OLS nuisance coefficient and `Q/L`). If the MVMR package is installed, server preflight compares the crosswalk numerically against the package.

Conditional F is a reliability diagnostic only; it never changes an adjustment class.

## Binary outcome

The individual-level binary DGP includes the common U term in the logistic outcome model. The primary event prevalence is calibrated to approximately 10% at scenario level and the calibrated intercept is then frozen across replicates. SNP-Y associations are obtained by SNP-wise logistic regression without a rare-outcome approximation.

The realistic binary simulation reports target-relative total-logOR behavior. It does not claim to recover separate structural and scale components. Those components are tested only in the controlled summary-level Theorems 2–3 module.

## Strength pilot and freeze

The production grid is intentionally gated by a strength pilot. The pilot uses the same latent genetic architecture and common random numbers across candidate RMS values, so changing RMS changes strength rather than the underlying random architecture. It selects the RMS closest to median minimum conditional-F targets 25, 10 and 5 and saves the choices to:

`simulation_output/validation/strength_calibration_frozen.rds`

Production loads this frozen file automatically. Do not retune the RMS after viewing the final simulation results.

## 100-core HPC implementation

The previous scripts recreated a cluster inside each parameter setting. Version 1.1 instead creates one persistent PSOCK cluster and reuses it across scenarios. Each worker is forced to one BLAS/OpenMP thread to prevent nested oversubscription.

Replicates are checkpointed in chunks (default 100). If a job is interrupted, rerunning with `SIM_RESUME=1` reuses valid checkpoints; replicate-specific deterministic seeds ensure that load balancing does not change results.

**Assumption:** the command below assumes the 100 cores belong to one node. If a scheduler gives 100 cores across multiple nodes, the PSOCK host list needs to be adapted to the scheduler allocation.

## Recommended R environment

Use R 4.3 or later if possible. Production code itself uses base R/statistics/parallel. The `MVMR` package is recommended for the independent preflight cross-check of conditional F.

If needed:

```bash
Rscript -e 'if (!requireNamespace("MVMR", quietly=TRUE)) install.packages("MVMR", repos=c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))'
```

## Complete publication run on a 100-core one-node server

The easiest command is:

```bash
unzip mvmr_simulation_v1_1.zip
cd mvmr_simulation_v1_1
bash RUN_SERVER_100_CORES.sh
```

The wrapper performs the required phases in order:

1. parse + validation + theorem + smoke preflight;
2. conditional-strength pilot;
3. preflight again using the frozen pilot calibration;
4. full primary six-DAG grid;
5. sample-overlap/alignment/N/prevalence/beta=0.05 sensitivities;
6. strong/borderline/weak conditional-strength sensitivity.

Defaults can be overridden, for example:

```bash
SIM_CORES=100 SIM_CHUNK_SIZE=100 PILOT_NREP=100 REQUIRE_MVMR=1 bash RUN_SERVER_100_CORES.sh
```

## Manual commands

```bash
cd mvmr_simulation_v1_1
mkdir -p simulation_output/logs

export SIM_CORES=100
export SIM_CHUNK_SIZE=100
export SIM_RESUME=1
export PILOT_NREP=100
export REQUIRE_MVMR=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript 97_preflight.R 2>&1 | tee simulation_output/logs/01_preflight_before_pilot.log
Rscript 98_pilot_strength.R 2>&1 | tee simulation_output/logs/02_strength_pilot.log
Rscript 97_preflight.R 2>&1 | tee simulation_output/logs/03_preflight_after_pilot.log
RUN_FULL_SIM=1 Rscript 99_run_all.R 2>&1 | tee simulation_output/logs/04_primary_full.log
Rscript 95_run_sensitivities_all.R 2>&1 | tee simulation_output/logs/05_sensitivities.log
Rscript 94_run_strength_sensitivity.R 2>&1 | tee simulation_output/logs/06_strength_sensitivity.log
```

For a detached primary run after phases 1–3 have passed:

```bash
nohup env RUN_FULL_SIM=1 SIM_CORES=100 SIM_CHUNK_SIZE=100 SIM_RESUME=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  Rscript 99_run_all.R \
  > simulation_output/logs/04_primary_full.log 2>&1 &
echo $!
```

## Production gates

`97_preflight.R` does all of the following before production is permitted:

- parses every R file;
- runs the formal validation suite;
- independently cross-checks weighted-regression estimates/SEs against base R `lm`;
- cross-checks vectorized continuous GWAS against base R `lm`;
- cross-checks binary `glm.fit` associations against base R `glm`;
- optionally/when required cross-checks conditional F against `MVMR::strength_mvmr`;
- runs the exact Theorem 1 and Theorems 2–3 modules;
- runs small continuous and binary finite-sample smoke tests;
- records the exact frozen strength calibration in the audit.

`99_run_all.R` refuses to start if the preflight was not rerun after the frozen strength pilot.

## Outputs

- primary raw: `simulation_output/raw/primary/`
- sensitivity raw: `simulation_output/raw/sensitivity/`
- primary summary: `simulation_output/summary/primary_scenario_summary.csv`
- sensitivity summaries: `simulation_output/summary/sensitivity_*_scenario_summary.csv`
- exact theorem outputs: `simulation_output/theorem/`
- frozen strength pilot and audits: `simulation_output/validation/`
- logs: `simulation_output/logs/`

## Acceptance rule

Any code modification must preserve the formal validation tests in `13_validation_tests.R`. A change that breaks same-set identity, exact theorem geometry, sample-design separation, binary scale decomposition, GWAS/weighted-regression cross-checks, conditional-F crosswalk, reproducibility or type-I/power mapping is not an acceptable implementation of the frozen protocol.
