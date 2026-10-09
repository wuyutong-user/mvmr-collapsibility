# Figure 2: continuous-outcome simulations, panels a-c.
# Reads archived plotting summaries and Monte Carlo outputs from simulation v1.1.
# Configure MVMR_SIM_OUTPUT_ROOT or MVMR_FIG2_DATA_ROOT for an external run.

# ----------------------------
# 0. Packages
# ----------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required packages before rendering Figure 2: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ----------------------------
# 1. Paths
# ----------------------------
# Package-relative input and output locations can be overridden with
# the documented environment variables.

PUBLIC_SCRIPT_DIR <- local({
  frames <- sys.frames()
  files <- Filter(Negate(is.null), lapply(frames, function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- if (length(files)) tail(files, 1)[[1]] else if (length(arg)) {
    sub("^--file=", "", arg[[1]])
  } else stop("Source this script by filename or run it with Rscript.")
  dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
})
source(file.path(PUBLIC_SCRIPT_DIR, "..", "helpers", "public_paths.R"))
PUBLIC_PACKAGE_ROOT <- normalizePath(file.path(PUBLIC_SCRIPT_DIR, "..", ".."),
                                     winslash = "/", mustWork = TRUE)
PUBLIC_PLOT_ROOT <- dirname(PUBLIC_SCRIPT_DIR)
SIM_OUTPUT_ROOT <- ije_simulation_output(PUBLIC_PACKAGE_ROOT)
PLOT_DATA_DEFAULT <- file.path(PUBLIC_PLOT_ROOT, "plotting_data")

FIG2_DATA_ROOT_ENV <- Sys.getenv("MVMR_FIG2_DATA_ROOT", unset = "")
LEGACY_FIG4_DATA_ROOT_ENV <- Sys.getenv("MVMR_FIG4_DATA_ROOT", unset = "")

required_files <- c(
  "FigS5_continuous_response_calibration_data.csv",
  "FigS6_continuous_target_relative_distribution_data.csv",
  "FigS8_continuous_precision_reliability_data.csv"
)

candidate_data_roots <- unique(Filter(
  nzchar,
  c(
    FIG2_DATA_ROOT_ENV,
    LEGACY_FIG4_DATA_ROOT_ENV,
    PLOT_DATA_DEFAULT,
    getwd()
  )
))

root_has_all_files <- function(root) {
  all(file.exists(file.path(root, required_files)))
}

matching_roots <- candidate_data_roots[
  vapply(candidate_data_roots, root_has_all_files, logical(1))
]

if (length(matching_roots) == 0) {
  stop(
    paste0(
      "Could not locate all Figure 2 input files.\n",
      "Set MVMR_FIG2_DATA_ROOT to the directory containing:\n  ",
      paste(required_files, collapse = "\n  "),
      "\nCandidate roots checked:\n  ",
      paste(candidate_data_roots, collapse = "\n  ")
    ),
    call. = FALSE
  )
}

FIG2_DATA_ROOT <- matching_roots[[1]]

FIG2_OUTPUT_ROOT <- Sys.getenv(
  "MVMR_FIG2_OUTPUT_ROOT",
  unset = file.path(
    dirname(SIM_OUTPUT_ROOT),
    "main_figures_v2",
    "Figure2"
  )
)

output_stem <- "Figure2_continuous_outcomes_v2"

SOURCE_DATA_DIR <- file.path(FIG2_OUTPUT_ROOT, "source_data", output_stem)
dir.create(FIG2_OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

message("Figure 2 data root:   ", FIG2_DATA_ROOT)
message("Figure 2 output root: ", FIG2_OUTPUT_ROOT)

# ----------------------------
# 2. Read data
# ----------------------------

response_dat <- read_csv(
  file.path(
    FIG2_DATA_ROOT,
    "FigS5_continuous_response_calibration_data.csv"
  ),
  show_col_types = FALSE
)

distribution_dat <- read_csv(
  file.path(
    FIG2_DATA_ROOT,
    "FigS6_continuous_target_relative_distribution_data.csv"
  ),
  show_col_types = FALSE
)

precision_dat <- read_csv(
  file.path(
    FIG2_DATA_ROOT,
    "FigS8_continuous_precision_reliability_data.csv"
  ),
  show_col_types = FALSE
)

# ----------------------------
# 3. Frozen ordering / labels
# ----------------------------

ROLE_LEVELS <- c(
  "confounder",
  "collider",
  "independent_cause",
  "upstream_surrogate_exposure",
  "downstream_surrogate_outcome",
  "downstream_surrogate_exposure"
)

ROLE_LABELS <- c(
  confounder = "Confounder",
  collider = "Collider",
  independent_cause = "Independent cause",
  upstream_surrogate_exposure = "Upstream surrogate\nof exposure",
  downstream_surrogate_outcome = "Downstream surrogate\nof outcome",
  downstream_surrogate_exposure = "Downstream surrogate\nof exposure"
)

CONFIG_LEVELS <- c(
  "all_direct_classes",
  "no_direct_shared",
  "no_direct_x",
  "no_direct_z"
)

CONFIG_LABELS <- c(
  all_direct_classes = "All direct\nclasses",
  no_direct_shared = "No direct-\nshared",
  no_direct_x = "No direct-\nX",
  no_direct_z = "No direct-\nZ"
)

# Supporting panels use compact display labels to protect final-size readability.
# The full formal configuration names remain in panel a and in the caption.
CONFIG_LABELS_COMPACT <- c(
  "All direct\nclasses" = "All classes",
  "No direct-\nshared" = "No shared",
  "No direct-\nX" = "No X",
  "No direct-\nZ" = "No Z"
)

CONFIG_LEGEND_LABELS <- c(
  "All direct classes" = "All classes",
  "No direct-shared" = "No shared",
  "No direct-X" = "No X",
  "No direct-Z" = "No Z"
)

# ------------------------------------------------------------
# Manuscript palette v7
# ------------------------------------------------------------
# Prespecified figure palette:
#   #CC247C  #E95351  #F7A24F  #FBEB66
#   #4EA660  #79CAFB  #5292F7  #AA77E9
# The figure keeps neutral grids/borders for readability, and uses this
# palette for all substantive data encodings.
COL_MAGENTA <- "#CC247C"
COL_CORAL   <- "#E95351"
COL_ORANGE  <- "#F7A24F"
COL_YELLOW  <- "#FBEB66"
COL_GREEN   <- "#4EA660"
COL_SKY     <- "#79CAFB"
COL_BLUE    <- "#5292F7"
COL_PURPLE  <- "#AA77E9"
COL_BLACK   <- "#222222"
COL_DARK    <- "#4D4D4D"
COL_GRID    <- "#ECEFF1"
COL_STRIP   <- "#F4F6F8"
COL_BORDER  <- "#C9CED3"

ESTIMATOR_COLOURS <- c(
  "Unadjusted IVW" = COL_CORAL,
  "Adjusted MVMR-IVW" = COL_BLUE
)

QUANTITY_COLOURS <- c(
  "Observed mean" = COL_BLUE,
  "Population prediction" = COL_ORANGE
)

QUANTITY_FILLS <- c(
  "Observed mean" = COL_BLUE,
  "Population prediction" = "white"
)

# Panel c uses redundant colour + shape encoding.
# Filled shapes 21–24 retain distinct silhouettes while allowing
# a dark outline, improving readability at final journal size.
CONFIG_SHAPES <- c(
  "All direct classes" = 21,
  "No direct-shared" = 24,
  "No direct-X" = 22,
  "No direct-Z" = 23
)

# Configuration colours drawn directly from the v6 palette.
CONFIG_FILLS <- c(
  "All direct classes" = COL_MAGENTA,
  "No direct-shared" = COL_GREEN,
  "No direct-X" = COL_BLUE,
  "No direct-Z" = COL_PURPLE
)

# ----------------------------
# 4. Data integrity checks
# ----------------------------

assert_columns <- function(dat, cols, object_name) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

assert_columns(
  response_dat,
  c(
    "role", "outcome_type", "beta_x_target", "instrument_config",
    "mean_response", "mean_role_predicted_response",
    "mean_role_response_error", "rmse_role_response_error", "n_converged"
  ),
  "response_dat"
)

assert_columns(
  distribution_dat,
  c(
    "role", "instrument_config", "beta_x_target",
    "estimator", "target_relative_error"
  ),
  "distribution_dat"
)

assert_columns(
  precision_dat,
  c(
    "role", "outcome_type", "beta_x_target",
    "instrument_config", "empirical_sd_ratio"
  ),
  "precision_dat"
)

# Figure 2 uses continuous outcomes and beta_X = 0.10.
panel_a_raw <- response_dat |>
  filter(
    outcome_type == "continuous",
    beta_x_target == 0.10
  )

panel_b_raw <- distribution_dat |>
  filter(
    role == "confounder",
    beta_x_target == 0.10
  )

panel_c_raw <- precision_dat |>
  filter(
    outcome_type == "continuous",
    beta_x_target == 0.10
  )

if (nrow(panel_a_raw) != 24L) {
  stop("Panel a expected 24 rows (6 roles x 4 configurations).", call. = FALSE)
}

if (nrow(panel_b_raw) != 16000L) {
  stop(
    "Panel b expected 16,000 rows (4 configurations x 2 estimators x 2,000 replicates).",
    call. = FALSE
  )
}

if (nrow(panel_c_raw) != 24L) {
  stop("Panel c expected 24 rows (6 roles x 4 configurations).", call. = FALSE)
}

if (!setequal(unique(panel_a_raw$role), ROLE_LEVELS)) {
  stop("Unexpected role set in Panel a data.", call. = FALSE)
}

if (!setequal(unique(panel_a_raw$instrument_config), CONFIG_LEVELS)) {
  stop("Unexpected instrument-configuration set in Panel a data.", call. = FALSE)
}

# ----------------------------
# 5. Figure-specific theme
# ----------------------------
# This is intentionally independent of all previous manuscript plotting themes.

theme_fig4 <- function(base_size = 8.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(
        colour = COL_GRID,
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = COL_BORDER,
        fill = NA,
        linewidth = 0.32
      ),
      strip.background = element_rect(
        fill = COL_STRIP,
        colour = COL_BORDER,
        linewidth = 0.32
      ),
      strip.text = element_text(
        colour = COL_BLACK,
        face = "plain",
        size = base_size - 0.1,
        margin = margin(2.4, 2.4, 2.4, 2.4)
      ),
      axis.title = element_text(
        colour = COL_BLACK,
        size = base_size
      ),
      axis.text = element_text(
        colour = COL_BLACK,
        size = base_size - 0.7
      ),
      axis.ticks = element_line(
        colour = COL_DARK,
        linewidth = 0.30
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(
        colour = COL_BLACK,
        size = base_size - 0.5,
        face = "plain"
      ),
      legend.text = element_text(
        colour = COL_BLACK,
        size = base_size - 0.7
      ),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(3.6, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ----------------------------
# 5a. Monte Carlo interval helper
# ----------------------------

add_response_mc_ci <- function(dat, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)

  dat |>
    mutate(
      response_error_second_moment =
        pmax(
          rmse_role_response_error^2 -
            mean_role_response_error^2,
          0
        ),
      response_sd =
        sqrt(
          (n_converged / pmax(n_converged - 1, 1)) *
            response_error_second_moment
        ),
      mcse_mean_response =
        response_sd / sqrt(n_converged),
      response_ci_lower =
        mean_response - z * mcse_mean_response,
      response_ci_upper =
        mean_response + z * mcse_mean_response
    )
}

# ----------------------------
# 6. Panel a
# Continuous coefficient response across six causal structures
# ----------------------------

panel_a_wide <- panel_a_raw |>
  add_response_mc_ci() |>
  mutate(
    role = factor(role, levels = ROLE_LEVELS),
    role_label = factor(
      ROLE_LABELS[as.character(role)],
      levels = unname(ROLE_LABELS[ROLE_LEVELS])
    ),
    instrument_config = factor(
      instrument_config,
      levels = CONFIG_LEVELS
    ),
    config_index = as.numeric(instrument_config),
    config_label = factor(
      CONFIG_LABELS[as.character(instrument_config)],
      levels = unname(CONFIG_LABELS[CONFIG_LEVELS])
    ),
    observed_x = config_index - 0.08,
    predicted_x = config_index + 0.08
  )

# Long-form source table is retained for source-data export and legend semantics.
panel_a <- bind_rows(
  panel_a_wide |>
    transmute(
      role,
      role_label,
      instrument_config,
      config_label,
      x = observed_x,
      quantity = "Observed mean",
      response = mean_response,
      response_ci_lower,
      response_ci_upper,
      mcse_mean_response
    ),
  panel_a_wide |>
    transmute(
      role,
      role_label,
      instrument_config,
      config_label,
      x = predicted_x,
      quantity = "Population prediction",
      response = mean_role_predicted_response,
      response_ci_lower = NA_real_,
      response_ci_upper = NA_real_,
      mcse_mean_response = NA_real_
    )
) |>
  mutate(
    quantity = factor(
      quantity,
      levels = c("Observed mean", "Population prediction")
    )
  )

p_a <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    colour = "#6E6E6E"
  ) +
  geom_errorbar(
    data = panel_a |>
      filter(quantity == "Observed mean"),
    aes(
      x = x,
      ymin = response_ci_lower,
      ymax = response_ci_upper
    ),
    width = 0.060,
    linewidth = 0.60,
    colour = COL_BLUE
  ) +
  geom_point(
    data = panel_a,
    aes(
      x = x,
      y = response,
      shape = quantity,
      fill = quantity,
      colour = quantity
    ),
    size = 2.15,
    stroke = 0.75
  ) +
  facet_wrap(
    ~ role_label,
    ncol = 3,
    scales = "fixed"
  ) +
  scale_shape_manual(
    values = c(
      "Observed mean" = 21,
      "Population prediction" = 21
    )
  ) +
  scale_fill_manual(values = QUANTITY_FILLS) +
  scale_colour_manual(values = QUANTITY_COLOURS) +
  scale_x_continuous(
    breaks = 1:4,
    labels = unname(CONFIG_LABELS[CONFIG_LEVELS]),
    limits = c(0.65, 4.35),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.05, 0.10, 0.15, 0.20),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  coord_cartesian(
    ylim = c(-0.012, 0.220),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = expression(
      "Coefficient response, " * Delta[X]
    ),
    shape = NULL,
    fill = NULL,
    colour = NULL
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        fill = unname(QUANTITY_FILLS),
        colour = unname(QUANTITY_COLOURS),
        linewidth = 0
      )
    ),
    fill = "none",
    colour = "none"
  ) +
  theme_fig4() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 7.0,
      lineheight = 0.95,
      margin = margin(t = 2)
    )
  )

# ----------------------------
# 7. Panel b
# Confounder target-relative error before and after adjustment
# ----------------------------

panel_b <- panel_b_raw |>
  mutate(
    instrument_config = factor(
      instrument_config,
      levels = CONFIG_LEVELS
    ),
    config_label = factor(
      CONFIG_LABELS[as.character(instrument_config)],
      levels = unname(CONFIG_LABELS[CONFIG_LEVELS])
    ),
    estimator = recode(
      estimator,
      "Unadjusted" = "Unadjusted IVW",
      "Adjusted" = "Adjusted MVMR-IVW"
    ),
    estimator = factor(
      estimator,
      levels = c("Unadjusted IVW", "Adjusted MVMR-IVW")
    )
  )

panel_b_means <- panel_b |>
  group_by(instrument_config, config_label, estimator) |>
  summarise(
    mean_estimation_error = mean(target_relative_error, na.rm = TRUE),
    .groups = "drop"
  )

p_b <- ggplot(
  panel_b,
  aes(
    x = config_label,
    y = target_relative_error,
    fill = estimator
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    colour = "#6E6E6E"
  ) +
  geom_boxplot(
    width = 0.66,
    colour = COL_DARK,
    position = position_dodge2(
      width = 0.76,
      preserve = "single"
    ),
    linewidth = 0.38,
    outlier.size = 0.38,
    outlier.alpha = 0.45,
    fatten = 1.25
  ) +
  # geom_point(
  #   data = panel_b_means,
  #   aes(
  #     x = config_label,
  #     y = mean_estimation_error,
  #     group = estimator
  #   ),
  #   inherit.aes = FALSE,
  #   position = position_dodge(width = 0.76),
  #   shape = 23,
  #   size = 1.85,
  #   stroke = 0.55,
  #   fill = "white",
  #   colour = COL_BLACK
  # ) +
  scale_fill_manual(
    values = ESTIMATOR_COLOURS,
    labels = c(
      "Unadjusted IVW" = "Unadjusted",
      "Adjusted MVMR-IVW" = "Adjusted"
    )
  ) +
  scale_x_discrete(
    labels = CONFIG_LABELS_COMPACT
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.05),
    expand = expansion(mult = c(0.04, 0.05))
  ) +
  coord_cartesian(
    ylim = c(-0.12, 0.30)
  ) +
  labs(
    x = NULL,
    y = "Estimation error",
    fill = NULL
  ) +
  theme_fig4() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.margin = margin(0, 0, 1.5, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 7.0,
      lineheight = 0.95,
      margin = margin(t = 2)
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title = NULL
    )
  )

# ----------------------------
# 8. Panel c
# Empirical SD ratio after adjustment
# ----------------------------

panel_c <- panel_c_raw |>
  transmute(
    role = factor(role, levels = ROLE_LEVELS),
    role_label = factor(
      ROLE_LABELS[as.character(role)],
      levels = rev(unname(ROLE_LABELS[ROLE_LEVELS]))
    ),
    config_label = factor(
      config_label,
      levels = c(
        "All direct classes",
        "No direct-shared",
        "No direct-X",
        "No direct-Z"
      )
    ),
    empirical_sd_ratio
  )

p_c <- ggplot(
  panel_c,
  aes(
    x = empirical_sd_ratio,
    y = role_label,
    shape = config_label,
    fill = config_label
  )
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.35,
    linetype = "dashed",
    colour = "#6E6E6E"
  ) +
  geom_point(
    colour = COL_DARK,
    size = 2.90,
    stroke = 0.62,
    position = position_dodge(width = 0.48)
  ) +
  scale_shape_manual(
    values = CONFIG_SHAPES,
    labels = CONFIG_LEGEND_LABELS
  ) +
  scale_fill_manual(
    values = CONFIG_FILLS,
    labels = CONFIG_LEGEND_LABELS,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(1.00, 1.05, 1.10, 1.15, 1.20),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  coord_cartesian(
    xlim = c(0.95, 1.20),
    clip = "off"
  ) +
  labs(
    x = expression(
      "Empirical SD ratio, " *
      SD[adjusted] / SD[unadjusted]
    ),
    y = NULL,
    shape = NULL
  ) +
  theme_fig4() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.margin = margin(0, 0, 1.5, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    axis.text.y = element_text(
      size = 7.1,
      lineheight = 0.95
    )
  ) +
  guides(
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      order = 1,
      title = NULL,
      override.aes = list(
        fill = unname(CONFIG_FILLS),
        colour = rep(COL_DARK, 4),
        size = rep(3.0, 4),
        stroke = rep(0.62, 4)
      )
    )
  )

# ----------------------------
# 10. Assemble Figure 2
# ----------------------------
# Panel a spans the full width as the hero panel.
# Panels b and c share the second row. Their legends remain local to each
# panel so estimator and instrument-configuration encodings stay distinct.

figure_design <- "
AA
BC
"

# Keep panel-specific legends separate.
# This is more robust across patchwork versions and also preserves
# the distinct semantics of estimator and genetic-configuration legends.

p_b <- p_b +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p_c <- p_c +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

fig2_core <- wrap_plots(
  p_a,
  p_b,
  p_c,
  design = figure_design
) +
  plot_layout(
    widths = c(1.00, 1.20),
    heights = c(1.64, 1.26)
  )

fig2 <- fig2_core +
  plot_annotation(
    tag_levels = "a"
  )

# Render-time patchwork alignment gate.
panel_alignment_r <- Sys.getenv(
  "NATURE_FIGURE_PANEL_ALIGNMENT_R",
  unset = file.path(PUBLIC_PACKAGE_ROOT, "plotting", "helpers", "panel_alignment.R")
)
panel_alignment_audit <- Sys.getenv(
  "NATURE_FIGURE_PANEL_ALIGNMENT_AUDIT",
  unset = file.path(PUBLIC_PACKAGE_ROOT, "plotting", "helpers", "audit_panel_alignment.py")
)
panel_alignment_python <- Sys.getenv(
  "NATURE_FIGURE_PYTHON",
  unset = if (nzchar(Sys.which("python3"))) Sys.which("python3") else Sys.which("python")
)

if (!file.exists(panel_alignment_r)) {
  stop(
    paste0("Panel-alignment helper not found: ", panel_alignment_r),
    call. = FALSE
  )
}

source(panel_alignment_r)
write_patchwork_panel_layout(
  plot = fig2_core,
  manifest_path = file.path(
    FIG2_OUTPUT_ROOT,
    paste0(output_stem, ".alignment-layout.json")
  ),
  width_in = 7.20,
  height_in = 6.45,
  panel_ids = NULL,
  exemptions = list()
)

# ----------------------------
# 11. Save source data used by Figure 2
# ----------------------------

write_csv(
  panel_a |>
    arrange(role, instrument_config, quantity),
  file.path(
    SOURCE_DATA_DIR,
    "Figure2a_continuous_response_source.csv"
  )
)

write_csv(
  panel_b |>
    arrange(instrument_config, estimator),
  file.path(
    SOURCE_DATA_DIR,
    "Figure2b_confounder_error_distribution_source.csv"
  )
)

write_csv(
  panel_b_means |>
    arrange(instrument_config, estimator),
  file.path(
    SOURCE_DATA_DIR,
    "Figure2b_confounder_mean_error_source.csv"
  )
)

write_csv(
  panel_c |>
    arrange(role, config_label),
  file.path(
    SOURCE_DATA_DIR,
    "Figure2c_empirical_sd_ratio_source.csv"
  )
)

# ----------------------------
# 12. Export
# ----------------------------

FIG_WIDTH_IN  <- 7.20
FIG_HEIGHT_IN <- 6.45

pdf_file <- file.path(
  FIG2_OUTPUT_ROOT,
  paste0(output_stem, ".pdf")
)

png_file <- file.path(
  FIG2_OUTPUT_ROOT,
  paste0(output_stem, ".png")
)

if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    pdf_file,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    family = "sans",
    bg = "white"
  )
} else {
  grDevices::pdf(
    pdf_file,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    useDingbats = FALSE,
    bg = "white"
  )
}
print(fig2)
grDevices::dev.off()

# Prefer ragg for high-resolution raster output when installed.
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = png_file,
    plot = fig2,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    bg = "white"
  )

  tiff_file <- file.path(
    FIG2_OUTPUT_ROOT,
    paste0(output_stem, ".tiff")
  )

  ggsave(
    filename = tiff_file,
    plot = fig2,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    units = "in",
    dpi = 600,
    device = ragg::agg_tiff,
    compression = "lzw",
    bg = "white"
  )
} else {
  ggsave(
    filename = png_file,
    plot = fig2,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# SVG is useful for final vector editing but is optional.
if (requireNamespace("svglite", quietly = TRUE)) {
  svg_file <- file.path(
      FIG2_OUTPUT_ROOT,
      paste0(output_stem, ".svg")
  )
  svglite::svglite(
    svg_file,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    bg = "white"
  )
  print(fig2)
  grDevices::dev.off()
}

message("Figure 2 rendering complete.")
message("PDF: ", pdf_file)
message("PNG: ", png_file)
message("Source data: ", SOURCE_DATA_DIR)

# Optional interactive display when run in RStudio.
print(fig2)
