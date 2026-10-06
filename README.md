# Supplementary Code 1

Public release prepared 6 October 2026 for the IJE manuscript.

The complete source code and reproducibility records are in
[Supplementary_Code_1.zip](./Supplementary_Code_1.zip). Download and extract the
archive, then follow the instructions below from the extracted
`Supplementary_Code_1` directory.

[Supplementary_Data_1.csv](./Supplementary_Data_1.csv) contains the 168-row
aggregate application results. The controlled inputs required for a full
application rerun are described below.

## Current reporting version

The application contains 168 BMI–candidate-trait–outcome triads. The corrected
adjustment counts are Necessary 5, Unnecessary 129, Overadjustment 22 and
Unclassified 12. Nine Collider decisions were changed to Unclassified. Working
roles, estimates, coefficient responses, conditional F statistics and
Independent-cause decisions were retained.

Run these commands from the package directory with Python 3:

```sh
python application/remap_collider_reporting.py
python application/export_reporting.py
```

The second command writes `application/check_output/Supplementary_Data_1.csv`
and compares it with the supplied expected CSV. These commands reproduce the
reporting export from archived aggregate results, rather than re-estimating MR
effects or validating retained-instrument assumptions.

## Simulation

`simulation/` contains the simulation v1.1 sources, archived scenario summaries,
exact calibrations and the strength-pilot freeze. All original simulation R
files are unchanged. The original MD5 manifest is retained as
`SOURCE_MANIFEST_MD5_original.tsv`; `SOURCE_MANIFEST_MD5.tsv` records the public
files, including the dated clarification of the historical audit document.

Run the preflight in a working copy, because it creates validation outputs:

```sh
cd simulation
Rscript 97_preflight.R
```

Use `RUN_SERVER_100_CORES.sh` only for a new production run with suitable compute
resources. Historical `nominal_adjustment_class` fields are preserved as
simulation metadata; the corrected Collider interpretation is described in
`COLLIDER_INTERPRETATION_20261006.md`.

## Application inputs and analysis pipeline

`application/code/role_assignment/19_run_eraca_pipeline.R` is the controlled-data
pipeline. See its `README_SERVER.md` and `environment/paths.env.example` for
input configuration. Replace example paths with authorised input locations.
The literature-edge CSV is a starter template, not the complete evidential
register used for the final reporting view. Historical reporting postprocessors
are in `application/provenance/` and are not current analysis entry points.

Participant-level data, controlled SNP-association caches, the LD reference,
credentials and third-party package libraries are not included. A complete
application rerun needs authorised GWAS inputs, the relevant evidence register,
PLINK and the appropriate LD reference. No full GWAS, application estimation or
production simulation rerun was performed for this public release.

## Figures

See `plotting/README.md` for current and historical figure sources. Defaults use
the package layout and can be overridden with environment variables. Plots
that need raw Monte Carlo replicates require those outputs to be regenerated;
only the archived summary and plotting data listed in this package are supplied.

## Environments and verification

The archived application session reports R 4.3.2 on Ubuntu 22.04.3 LTS but does
not record the loaded versions of TwoSampleMR, MVMR and MRSL. Those historical
versions remain unavailable. Local check records and the Figure 5 rendering
session are separate records, not a reconstruction of that environment.
The simulation uses base R, stats and parallel. See the relevant scripts for
optional cross-check and plotting dependencies.

`PUBLIC_RELEASE_VERIFICATION.json` describes public-release checks.
`verification_results.json` and `verification/` retain the preceding revision's
checks. Source-copy hashes are historical records; current bytes are recorded
in `release_manifest_sha256.csv`. Path/documentation edits and retained files
are listed in `public_source_provenance.csv`.

## AI assistance

ChatGPT/Codex assisted manuscript revisions, the reporting revision for the nine
Collider labels, edits and exports of Figures 4 and 5, and preparation of this
public package. The historical simulation audit also records ChatGPT assistance.
These statements describe known use; they do not establish the origin of every
historical script. AI assistance does not replace scientific verification or
author responsibility for the code, analysis and submitted findings.

## HTML copy

The companion HTML provides readable source and an exact embedded ZIP.
`restore_code_from_html.py` can restore all package files from the HTML and
checks each SHA256 hash. Restoration does not execute analysis code.

Top-level `source_manifest.csv` compares original and public source bytes,
excluding the three manifest files themselves. `archive_manifest_sha256.csv`
covers the public files apart from that manifest and the release manifest.
`release_manifest_sha256.csv` covers every package file except itself.
Original manifests are retained under `application/provenance/original_*.csv`.
