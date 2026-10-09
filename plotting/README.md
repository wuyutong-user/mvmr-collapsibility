# Main-figure plotting sources

This directory contains the R sources and inputs for numerical main Figures 2,
3 and 5. PowerPoint schematics, obsolete figure scripts, preview exporters and
supplementary plotting scripts are outside this slim directory.

## Run from the repository root

```sh
Rscript plotting/R/Figure2.R
Rscript plotting/R/Figure3.R
```

For Figure 5, change to the directory containing its adjacent input:

```sh
cd plotting/figure5
Rscript Figure5.R
```

Figure 2 reads three archived CSV summaries. Figure 3 reads two archived CSV
summaries and `plotting_data/Figure3_binary_replicates.rds`, which combines the
eight required columns of the 24 archived binary primary simulation outputs
(1000 replicates each). Responses are copied without re-estimation and checked
against the frozen summary means. Figure 5 reads the archived 168-triad CSV;
the decision counts are Necessary 5, Unnecessary 129, Overadjustment 22 and
Unclassified 12.

## Dependencies

Figure 2 requires readr, dplyr, tidyr, ggplot2, patchwork, scales and jsonlite.
Its panel-alignment check uses the bundled helper and Python 3. Set
`NATURE_FIGURE_PYTHON` to the Python executable if it is not on PATH.

Figure 3 requires readr, dplyr, ggplot2, patchwork and scales. Its SVG export
also uses svglite. Figure 5 requires ggplot2, dplyr, ggrepel, patchwork, svglite
and jsonlite. Available Cairo fonts affect glyph rendering. ragg and svglite
provide optional additional Figure 2 exports.

`MVMR_PROJECT_ROOT`, `MVMR_SIM_OUTPUT_ROOT`, `MVMR_FIG2_DATA_ROOT`,
`MVMR_FIG2_OUTPUT_ROOT` and `MVMR_FIG3_OUTPUT_ROOT` can override applicable
defaults. Output folders are created when rendering.

## Scope and verification

This preparation changes file organisation, input-loading paths and stale
figure-number messages. Numerical plotting calculations, colours and display
settings are retained. No simulation or application estimation was rerun.
R syntax, source preservation and input-grid/summary agreement were checked;
full image rendering was not repeated in this preparation. Required plotting
packages must be installed before running the commands above.

The top-level `Supplementary_Code_1.zip` contains this same slim plotting
directory. The ZIP and repository sources are synchronised.
`files_sha256.csv` covers this directory except itself;
`repository_files_sha256.csv` covers the repository files except itself.
