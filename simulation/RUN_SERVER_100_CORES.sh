#!/usr/bin/env bash
set -euo pipefail

# Complete one-node 100-core publication run.
# Override with e.g. SIM_CORES=96 bash RUN_SERVER_100_CORES.sh

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

export SIM_CORES="${SIM_CORES:-100}"
export SIM_CHUNK_SIZE="${SIM_CHUNK_SIZE:-100}"
export SIM_RESUME="${SIM_RESUME:-1}"
export PILOT_NREP="${PILOT_NREP:-100}"
export REQUIRE_MVMR="${REQUIRE_MVMR:-1}"

# Prevent nested BLAS/OpenMP oversubscription inside 100 R workers.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p simulation_output/logs

cat "=== Phase 1: preflight before strength pilot ==="
Rscript 97_preflight.R 2>&1 | tee simulation_output/logs/01_preflight_before_pilot.log

cat "=== Phase 2: prespecified conditional-strength pilot ==="
Rscript 98_pilot_strength.R 2>&1 | tee simulation_output/logs/02_strength_pilot.log

cat "=== Phase 3: preflight again using frozen pilot calibration ==="
Rscript 97_preflight.R 2>&1 | tee simulation_output/logs/03_preflight_after_pilot.log

cat "=== Phase 4: primary six-DAG production grid ==="
RUN_FULL_SIM=1 Rscript 99_run_all.R 2>&1 | tee simulation_output/logs/04_primary_full.log

cat "=== Phase 5: prespecified sample/alignment/N/prevalence/effect-size sensitivities ==="
Rscript 95_run_sensitivities_all.R 2>&1 | tee simulation_output/logs/05_sensitivities.log

cat "=== Phase 6: strong/borderline/weak conditional-strength sensitivity ==="
Rscript 94_run_strength_sensitivity.R 2>&1 | tee simulation_output/logs/06_strength_sensitivity.log

cat "=== ALL REQUESTED SIMULATIONS COMPLETE ==="
cat "Primary summary: simulation_output/summary/primary_scenario_summary.csv"
cat "Sensitivity summaries: simulation_output/summary/sensitivity_*_scenario_summary.csv"
cat "Preflight audit: simulation_output/validation/preflight_audit.txt"
