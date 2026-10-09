# MVMR collapsibility

Code and aggregate results accompanying *A collapsibility-based framework for
candidate-trait adjustment in multivariable Mendelian randomization using
summary data*.

## Download and run

Run the commands below from the repository root or the extracted code folder.
The [v1.0.0 manuscript archive](https://github.com/wuyutong-user/mvmr-collapsibility/releases/tag/v1.0.0)
is retained as the version cited in the manuscript. This maintenance candidate
improves code organisation and path handling; it retains the archived data.
The repository source folders and `Supplementary_Code_1.zip` contain the same
source bundle. Extract the ZIP and enter its `Supplementary_Code_1` folder
before running the commands below.

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
| Reproduce numerical main figures | `plotting/README.md` | Bundled summaries and frozen plotting inputs |

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

The repository and ZIP both contain the same slim numerical plotting sources.
PowerPoint schematics, obsolete plotting scripts, supplementary plotting
scripts and preview exporters are outside this bundle. Figures 1 and 4 are
schematic drawings; their editable sources are outside this numerical code
bundle. See the manuscript's figures for those diagrams.

The historical application session records R 4.3.2 on Ubuntu 22.04.3 LTS;
loaded TwoSampleMR, MVMR and MRSL versions were not recorded. Local check
records are separate from that historical environment. The simulation uses
base R, stats and parallel; figure scripts list their additional dependencies.

## Provenance and verification

`CODE_PROVENANCE.md` describes the source and maintenance records. Earlier
verification records document the original release and prior maintenance checks. `MAINTENANCE_VERIFICATION_20261008.json` records source checks and the
9 October 2026 ZIP synchronisation; `release_manifest_sha256.csv` records the
files in the current source bundle.

The source bundle excludes the top-level ZIP and `repository_files_sha256.csv`
to avoid including an archive or its repository checksum list inside itself.
`source_manifest.csv` records current source bytes and retained original hashes
where available; new derived inputs have no single original-file hash.
It excludes the three bundle manifests themselves.
`archive_manifest_sha256.csv` covers bundle files except itself and the release
manifest. `release_manifest_sha256.csv` covers bundle files except itself.
Historical manifests are preserved under
`application/provenance/maintenance_baseline_20261008/`.

The companion HTML belongs to the v1.0.0 archive. Its embedded ZIP and
`restore_code_from_html.py` restore that version, rather than this maintenance
copy. Restoration checks file hashes and does not execute analysis code.

`Supplementary_Code_1.zip` is the synchronised download copy of the source
bundle, under a single `Supplementary_Code_1/` folder. Its files match the
repository source files byte for byte. `plotting/files_sha256.csv` covers the
slim plotting directory except that manifest itself.
`repository_files_sha256.csv` covers the repository files, including the ZIP,
except the repository manifest itself.
