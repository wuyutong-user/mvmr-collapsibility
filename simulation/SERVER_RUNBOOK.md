# Server runbook — 100 cores, one node

## Recommended complete run

```bash
unzip mvmr_simulation_v1_1.zip
cd mvmr_simulation_v1_1
bash RUN_SERVER_100_CORES.sh
```

The wrapper uses 100 cores by default, checkpoints every 100 replicates, resumes completed checkpoints, requires the MVMR package cross-check by default, and runs preflight → strength pilot → post-pilot preflight → primary grid → standard sensitivities → strength sensitivities.

## If MVMR is not installed

```bash
Rscript -e 'if (!requireNamespace("MVMR", quietly=TRUE)) install.packages("MVMR", repos=c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))'
```

Then run the wrapper.

## Full manual sequence

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

Because `set -o pipefail` is used in the wrapper, a failed R phase stops the complete run rather than being hidden by `tee`.

## Resume after interruption

Simply rerun:

```bash
bash RUN_SERVER_100_CORES.sh
```

The strength pilot is deterministic, and production/sensitivity scenario chunks are reused when valid. For primary production alone:

```bash
RUN_FULL_SIM=1 SIM_CORES=100 SIM_CHUNK_SIZE=100 SIM_RESUME=1 Rscript 99_run_all.R
```

## Detached primary job

After both preflights and the strength pilot have passed:

```bash
nohup env RUN_FULL_SIM=1 SIM_CORES=100 SIM_CHUNK_SIZE=100 SIM_RESUME=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  Rscript 99_run_all.R \
  > simulation_output/logs/04_primary_full.log 2>&1 &
echo $!
```

## Files to inspect before accepting a production run

```bash
cat simulation_output/validation/preflight_audit.txt
cat simulation_output/validation/strength_rms_recommendation.csv
head simulation_output/summary/primary_scenario_summary.csv
```

The post-pilot preflight must report `PREFLIGHT PASS` and the production runner must finish with `PRIMARY PRODUCTION RUN COMPLETE.`
