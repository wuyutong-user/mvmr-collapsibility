> Controlled-data pipeline. The paths below are examples to replace on an authorised system.
> The bundled literature-edge CSV is a starter template. This pipeline does not by itself
> reproduce the final display roles or the corrected Collider adjustment decisions;
> use the reporting scripts described in the top-level README for the submitted register.

# ERACA server execution guide

## 1. What this folder runs

The scripts implement 168 local `BMI–candidate trait–cardiovascular outcome` analyses:

1. prepare instruments and association caches for 27 traits;
2. run 388 bidirectional UVMR direction tasks;
3. build 168 marginal three-node graphs;
4. run local MVMR and, when available, official MRSL pruning;
5. assign bounded role hypotheses;
6. evaluate the same-instrument-set coefficient identity and merge the existing application response.

Controlled GWAS files remain under the server paths. Do not upload them to a public repository.

## 2. Upload layout

Upload the whole `reproducibility_package` folder to the server. The examples below assume:

```bash
export ERACA_PACKAGE_ROOT=/path/to/Supplementary_Code_1/application
cd "$ERACA_PACKAGE_ROOT"
```

The default controlled paths reproduce the earlier scripts:

```bash
export ERACA_CONTROLLED_ROOT=/path/to/controlled_data
export ERACA_CANDIDATE_DIR=/path/to/controlled_data/cov_data
export ERACA_OUTCOME_DIR=/path/to/outcome_gwas
export ERACA_PLINK=/path/to/plink
export ERACA_LD_BFILE=/path/to/ld_reference/reference
export ERACA_OUTPUT_DIR=/path/to/controlled_data/ERACA_RESULT
```

Local LD clumping first matches SNP IDs directly to the PLINK `.bim` file. If
the outcome GWAS uses rsIDs but the reference panel uses another variant-ID
scheme, the code falls back to chromosome–position matching and retains the
original GWAS SNP ID in the saved instrument table.

If an outcome still cannot be matched to the LD reference, its outcome
instrument file is recorded as unavailable. Forward analyses with that outcome
as the dependent variable continue; reverse UVMR directions using that outcome
as the exposure are retained as `not_estimable` rather than being assigned
unclumped SNPs.

If the BMI file is in another location, set:

```bash
export ERACA_BMI_GWAS=/absolute/path/to/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
```

## 3. R packages

Use R 4.3 or later. Check the existing server environment first:

```bash
Rscript -e 'pkgs <- c("data.table","TwoSampleMR","ieugwasr","MVMR","igraph"); print(setNames(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE), pkgs))'
```

Install missing standard dependencies in a user library:

```bash
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); install.packages(c("data.table","MVMR","igraph","remotes"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages(c("TwoSampleMR","ieugwasr"), repos=c("https://mrcieu.r-universe.dev","https://cloud.r-project.org"))'
```

Install the official MRSL implementation:

```bash
Rscript -e 'remotes::install_github("hhoulei/MRSL", upgrade="never")'
```

MRSL requires a working compiler plus BLAS/LAPACK. If MRSL is unavailable, local MVMR still runs; the graph is labelled `marginal_unresolved` and role confidence is downgraded. The pipeline never silently substitutes a fabricated conditional graph.

## 4. Outcome ZIP manifest

The code first tries to match ZIP filenames using outcome names and abbreviations. If the files are generically named, such as `AAA outcome1.zip`, run:

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=14 --dry-run
```

When automatic matching fails, the script writes:

```text
/path/to/controlled_data/eraca_outcome_manifest.csv
```

Fill the `zip_file` column for all 14 rows using filenames from `/path/to/outcome_gwas`, then repeat the dry run. The required columns are:

```csv
name,full_name,zip_file
AAA,Abdominal aortic aneurysm,actual_file_name.zip
```

## 5. Literature edge prior

`code/role_assignment/literature_edge_evidence.csv` is a structured input, not an
automatically generated literature review. It currently contains a small set of
starter rows to demonstrate the format. Before confirmatory role assignment,
replace or expand these rows using the reviewed external evidence for the
candidate–outcome pairs. Missing literature rows are treated as missing prior
evidence; they do not create a causal edge.

## 6. Run tests first

```bash
Rscript code/role_assignment/tests/test_eraca_synthetic_triads.R
```

The final line must be:

```text
All ERACA synthetic-triad tests passed.
```

These fixtures test task counts and graph-to-role rules. They are software tests, not additional scientific simulations.

## 7. Validate inputs

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=14 --dry-run
```

Expected message:

```text
Dry run completed: 27 trait inputs and output paths validated.
```

## 8. Prepare instruments and cached SNP associations

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=14
```

This reads each large GWAS source only during preparation and writes reusable RDS files under:

```text
$ERACA_OUTPUT_DIR/instruments/
$ERACA_OUTPUT_DIR/association_cache/
```

Existing checkpoints are reused. Add `--overwrite` only when the source data, phenotype mapping or clumping settings have changed.

## 9. Smoke test BMI–LDL-C–CAD

Run the six required marginal directions:

```bash
for id in \
  D013_BMI_to_LDL_C \
  D014_LDL_C_to_BMI \
  D193_LDL_C_to_CAD \
  D194_CAD_to_LDL_C \
  D361_BMI_to_CAD \
  D362_CAD_to_BMI
do
  ERACA_TASK_ID="$id" Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15
done
```

Aggregate the available direction results and build the local graph:

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15-aggregate
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=16
ERACA_TRIAD_ID=T085_BMI_LDL_C_CAD Rscript code/role_assignment/19_run_eraca_pipeline.R --step=17
```

For this partial smoke test, role assignment requires:

```bash
ERACA_ALLOW_PARTIAL=true Rscript code/role_assignment/19_run_eraca_pipeline.R --step=18
```

Inspect:

```text
$ERACA_OUTPUT_DIR/uvmr_primary_edges.csv
$ERACA_OUTPUT_DIR/graphs/conditional_graph_consensus.csv
$ERACA_OUTPUT_DIR/graphs/same_set_response_decomposition.csv
$ERACA_OUTPUT_DIR/role_assignments_168.csv
```

## 10. Full serial run

After the smoke test succeeds:

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=all --skip-tests
```

The run is checkpointed. Repeating the command resumes completed work.

## 11. Batch UVMR on a scheduler

Split the 388 directions into ranges. Example with four jobs:

```bash
ERACA_TASK_START=1   ERACA_TASK_END=97  Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15
ERACA_TASK_START=98  ERACA_TASK_END=194 Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15
ERACA_TASK_START=195 ERACA_TASK_END=291 Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15
ERACA_TASK_START=292 ERACA_TASK_END=388 Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15
```

After all jobs finish:

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=15-aggregate
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=16
```

Do not run the aggregation step concurrently with jobs that are still writing checkpoints.

## 12. Batch the 168 local graphs

Example with four jobs:

```bash
ERACA_TRIAD_START=1   ERACA_TRIAD_END=42  Rscript code/role_assignment/19_run_eraca_pipeline.R --step=17
ERACA_TRIAD_START=43  ERACA_TRIAD_END=84  Rscript code/role_assignment/19_run_eraca_pipeline.R --step=17
ERACA_TRIAD_START=85  ERACA_TRIAD_END=126 Rscript code/role_assignment/19_run_eraca_pipeline.R --step=17
ERACA_TRIAD_START=127 ERACA_TRIAD_END=168 Rscript code/role_assignment/19_run_eraca_pipeline.R --step=17
```

Then aggregate and assign roles:

```bash
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=17-aggregate
Rscript code/role_assignment/19_run_eraca_pipeline.R --step=18
```

## 13. Main outputs

| File | Meaning |
|---|---|
| `direction_task_registry.csv` | 388 prespecified directions |
| `uvmr_primary_edges.csv` | FDR-adjusted primary edge evidence |
| `graphs/marginal_edge_ledger.csv` | Literature and UVMR evidence kept separately |
| `graphs/conditional_edges_mvmr.csv` | Adjust-all local MVMR direct-edge diagnostics |
| `graphs/conditional_graph_consensus.csv` | MRSL-pruned or explicitly unresolved graph edges |
| `graphs/same_set_response_decomposition.csv` | Fixed-instrument-set theorem identity |
| `role_assignments_168.csv` | Bounded role hypotheses and confidence |
| `role_response_concordance_168.csv` | Role-specific response evaluation |
| `main_text_candidate_triads.csv` | Prespecified screen for main-text examples |

## 14. Interpretation safeguards

- `not_supported` is not proof of no edge.
- `marginal_unresolved` is not a conditional causal graph.
- `feedback_or_cycle`, `mediator_out_of_library`, `overlapping_measure`, `ambiguous` and `insufficient_evidence` are valid final outputs.
- Use `observed_response_same_set` for the algebraic summary-regression identity. For the 14 binary outcomes, the structural and log-OR scale components are not separately identified; therefore `response_concordance` remains `descriptive_binary_total_response_scale_sensitive`. The legacy Figure 6 response is retained separately because its unadjusted and adjusted estimates may use different SNP sets.
- Do not use the observed response to revise the frozen role label.
