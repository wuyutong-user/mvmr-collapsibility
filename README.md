# MVMR collapsibility

Code and aggregate results accompanying *A collapsibility-based framework for
candidate-trait adjustment in multivariable Mendelian randomization using
summary data*.

## Download and run

Run the commands below from the repository root or the extracted code folder.
The [v1.0.0 manuscript archive](https://github.com/wuyutong-user/mvmr-collapsibility/releases/tag/v1.0.0)
is retained as the version cited in the manuscript. This maintenance candidate
improves code organisation and path handling; it retains the archived data.

## Aggregate application results

[Supplementary_Data_1.csv](./Supplementary_Data_1.csv) contains 168 application
triads and eight fields: triad ID, candidate trait, outcome, working structural
class, mediator-like flag, adjustment label, observed logOR response and minimum
conditional F. It contains aggregate results, rather than participant-level
or SNP-level GWAS data. The CSV is also included in the code archive.

To reconstruct and check this export using Python 3:

```sh
python application/remap_collider_reporting.py
python application/export_reporting.py
```

The scripts reproduce the archived reporting revision and verify the result
against `application/expected/Supplementary_Data_1.csv`. The adjustment counts
are Necessary 5, Unnecessary 129, Overadjustment 22 and Unclassified 12.
They do not re-estimate MR effects or validate retained-instrument assumptions.

## Code entry points

| Task | Entry point | Inputs |
|---|---|---|
| Reproduce the aggregate results CSV | `application/export_reporting.py` | Bundled reporting register |
| Reproduce the archived Collider reporting revision | `application/remap_collider_reporting.py` | Bundled pre-revision register |
| Check the simulation | `simulation/97_preflight.R` | Bundled configuration and calibration records |
| Run a new controlled-data application | `application/code/role_assignment/19_run_eraca_pipeline.R` | Authorised GWAS inputs and LD reference |
| Locate figure sources | `plotting/README.md` | Figure-specific archived inputs or regenerated replicates |

The original preparation scripts `application/code/09_*` through `12_*` are
retained as historical sources. Use the documented role-assignment pipeline
for a new application run. Historical reporting postprocessors are under
`application/provenance/`.

## Simulation

`simulation/` contains the simulation v1.1 sources, archived scenario summaries,
exact calibrations and the strength-pilot freeze. The simulation R sources are
unchanged in this maintenance copy.

Run the preflight in a working copy; it creates validation outputs:

```sh
cd simulation
Rscript 97_preflight.R
```

See `simulation/README.md` for configuration and production-run instructions.
Historical `nominal_adjustment_class` fields remain as simulation metadata;
`COLLIDER_INTERPRETATION_20261006.md` documents the reporting interpretation.

## Controlled-data application

See `application/code/role_assignment/README_SERVER.md` and
`environment/paths.env.example` for setup. A full run requires authorised
GWAS inputs, the relevant evidence register, PLINK and an appropriate LD
reference. The bundled literature-edge CSV is a starter template, rather than
the complete evidence register used for the final reporting view.

Direct-MVMR P/F support is recorded separately from MRSL graph pruning and
the final edge states. Participant-level data, controlled SNP-association
caches, credentials and third-party package libraries are not distributed.
No full application estimation or production simulation was rerun for this
maintenance copy.

## Figures and environments

`plotting/` contains a slim set of sources and inputs for numerical main
Figures 2, 3 and 5. See `plotting/README.md` for commands and dependencies.
The Figure 3 input combines the required fields of 24 archived binary
simulation outputs without re-estimation. The slim set excludes PowerPoint
schematics, obsolete figure scripts and preview exporters.

The top-level ZIP is still the earlier full maintenance bundle, including
supplementary plotting sources. Its reconciliation is deferred. The archived
Figure 4 schematic in that ZIP needs its workflow wording synchronised with
the manuscript before a final release. Figures 1 and 4 are schematic drawings;
their PowerPoint sources are outside the slim numerical plotting directory.

The historical application session records R 4.3.2 on Ubuntu 22.04.3 LTS;
loaded TwoSampleMR, MVMR and MRSL versions were not recorded. Local check
records are separate from that historical environment. The simulation uses
base R, stats and parallel; figure scripts list their additional dependencies.

## Provenance and verification

`CODE_PROVENANCE.md` describes the source and maintenance records. Existing source manifests and verification records are historical
records. `MAINTENANCE_VERIFICATION_20261008.json` records the new checks;
`release_manifest_sha256.csv` records the files in the earlier ZIP bundle.

Within the source bundle, `source_manifest.csv` compares original and public
source bytes, excluding the three manifests themselves.
`archive_manifest_sha256.csv` covers bundle files except itself and the release
manifest. `release_manifest_sha256.csv` covers bundle files except itself. The prior manifests are preserved under
`application/provenance/maintenance_baseline_20261008/`.

The companion HTML belongs to the v1.0.0 archive. Its embedded ZIP and
`restore_code_from_html.py` restore that version, rather than this maintenance
copy. Restoration checks file hashes and does not execute analysis code.

The top-level `Supplementary_Code_1.zip` retains the earlier full maintenance
bundle and has not been rebuilt during the plotting-directory reduction. The
source manifests describe that ZIP's contents; the current slim plotting
folder has its own `plotting/files_sha256.csv`.
`repository_files_sha256.csv` covers the current repository candidate,
including the retained ZIP, except the repository manifest itself.
