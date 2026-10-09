# Figure 3 v10: binary response distributions and the three-facet
# structural--residual decomposition. Panel-a fills match Figure 2.
#
# Panel a shows replicate distributions. Panel b plots deterministic
# calibration points against the positive common scale S; joined points
# show the trend across the three prespecified S values, not uncertainty.
#
# Panel-a inputs are the archived binary simulation responses, combined
# without re-estimation. Plot calculations and display settings are retained.

required_packages <- c("readr", "dplyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

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
project_root <- ije_simulation_root(PUBLIC_PACKAGE_ROOT)
raw_file <- file.path(PUBLIC_PLOT_ROOT, "plotting_data", "Figure3_binary_replicates.rds")
summary_file <- file.path(
  PUBLIC_PLOT_ROOT, "plotting_data",
  "FigS9_binary_total_logOR_response_data.csv"
)
scale_file <- file.path(
  PUBLIC_PLOT_ROOT, "plotting_data",
  "FigS4_binary_scale_calibration_data.csv"
)
output_dir <- Sys.getenv(
  "MVMR_FIG3_OUTPUT_ROOT",
  unset = file.path(project_root, "main_figures_v2", "Figure3")
)
output_stem <- "Figure3_binary_outcomes_v10_classic_decomposition"
source_dir <- file.path(output_dir, "source_data", output_stem)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

role_levels <- c(
  "confounder", "collider", "independent_cause",
  "upstream_surrogate_exposure", "downstream_surrogate_outcome",
  "downstream_surrogate_exposure"
)
role_labels <- c(
  confounder = "Confounder",
  collider = "Collider",
  independent_cause = "Independent cause",
  upstream_surrogate_exposure = "Upstream surrogate\nof exposure",
  downstream_surrogate_outcome = "Downstream surrogate\nof outcome",
  downstream_surrogate_exposure = "Downstream surrogate\nof exposure"
)
config_levels <- c(
  "all_direct_classes", "no_direct_shared", "no_direct_x", "no_direct_z"
)
config_labels <- c(
  all_direct_classes = "All direct\nclasses",
  no_direct_shared = "No direct-\nshared",
  no_direct_x = "No direct-X",
  no_direct_z = "No direct-Z"
)
mode_levels <- c(
  "No scale residual", "Reinforce", "Partial offset",
  "Exact offset", "Scale only"
)
mode_labels <- c(
  "No scale residual" = "Scale = 0 (no residual)",
  "Reinforce" = "Same sign (reinforced)",
  "Partial offset" = "Partial cancellation",
  "Exact offset" = "Complete cancellation",
  "Scale only" = "Structural = 0"
)

# Panel-a fills are the exact four direct-source colours used in Figure 2.
# Panel b uses a separate three-colour component key.
config_fills <- c(
  "All direct\nclasses" = "#CC247C",
  "No direct-\nshared" = "#4EA660",
  "No direct-X" = "#5292F7",
  "No direct-Z" = "#AA77E9"
)
col_box_edge <- "#38434C"
col_struct <- "#659ABA"
col_scale <- "#E08232"
col_total <- "#B0101A"
col_text <- "#222222"
col_grid <- "#E9EEF2"
col_border <- "#C9CED3"
col_strip <- "#F4F6F8"              # match Figure 2 facet headers

if (!file.exists(raw_file)) {
  stop("Figure 3 binary replicate input was not found: ", raw_file)
}
panel_a_raw <- readRDS(raw_file)
needed <- c(
  "role", "outcome_type", "sample_design", "beta_x_target",
  "instrument_config", "replicate", "response", "converged"
)
if (!is.data.frame(panel_a_raw) || !all(needed %in% names(panel_a_raw))) {
  stop("Unexpected Figure 3 binary replicate structure")
}
if (anyNA(panel_a_raw[needed]) || any(!panel_a_raw$converged) ||
    any(!is.finite(panel_a_raw$response))) {
  stop("Incomplete binary replicates or non-finite responses")
}
scenario_counts <- table(interaction(
  panel_a_raw$role, panel_a_raw$instrument_config, drop = TRUE
))
if (length(scenario_counts) != 24L || any(scenario_counts != 1000L) ||
    anyDuplicated(panel_a_raw[c("role", "instrument_config", "replicate")])) {
  stop("Expected 24 binary scenarios with 1000 unique replicates each")
}
panel_a_raw <- panel_a_raw[, needed, drop = FALSE]
if (
  nrow(panel_a_raw) != 24000L ||
  !setequal(unique(panel_a_raw$role), role_levels) ||
  !setequal(unique(panel_a_raw$instrument_config), config_levels) ||
  !all(panel_a_raw$outcome_type == "binary") ||
  !all(panel_a_raw$sample_design == "two_sample") ||
  !all(panel_a_raw$beta_x_target == 0.10)
) {
  stop("Binary panel-a scenario grid failed validation")
}

panel_a <- panel_a_raw |>
  mutate(
    role = factor(role, levels = role_levels),
    role_label = factor(
      role_labels[as.character(role)],
      levels = unname(role_labels[role_levels])
    ),
    instrument_config = factor(instrument_config, levels = config_levels),
    config_label = factor(
      config_labels[as.character(instrument_config)],
      levels = unname(config_labels[config_levels])
    )
  )

# Check the distribution-derived means against the frozen summary dataset.
summary_dat <- readr::read_csv(summary_file, show_col_types = FALSE) |>
  filter(
    outcome_type == "binary",
    sample_design == "two_sample",
    beta_x_target == 0.10
  ) |>
  select(role, instrument_config, mean_response)
summary_check <- panel_a |>
  group_by(role, instrument_config) |>
  summarise(
    n = n(),
    mean_from_raw = mean(response),
    .groups = "drop"
  ) |>
  left_join(
    summary_dat,
    by = c("role", "instrument_config")
  )
if (
  nrow(summary_check) != 24L ||
  any(summary_check$n != 1000L) ||
  anyNA(summary_check$mean_response) ||
  max(abs(summary_check$mean_from_raw -
          summary_check$mean_response)) > 1e-10
) {
  stop("Raw and frozen-summary binary responses disagree")
}

scale_dat <- readr::read_csv(scale_file, show_col_types = FALSE)
needed_scale <- c(
  "simulation_module", "S", "residual_mode", "delta_struct",
  "delta_scale_observed", "delta_total", "decomposition_error"
)
if (!all(needed_scale %in% names(scale_dat))) {
  stop("Scale-calibration input is missing required fields")
}
panel_b <- scale_dat |>
  filter(simulation_module == "binary_common_scale_exact") |>
  mutate(
    S_facet = factor(
      sprintf("S = %.2f", S),
      levels = c("S = 0.75", "S = 1.00", "S = 1.50")
    ),
    mode_label = factor(
      mode_labels[residual_mode],
      levels = rev(unname(mode_labels[mode_levels]))
    ),
    delta_struct_plot = if_else(abs(delta_struct) < 1e-12, 0, delta_struct),
    delta_scale_plot = if_else(
      abs(delta_scale_observed) < 1e-12, 0, delta_scale_observed
    ),
    delta_total_plot = if_else(abs(delta_total) < 1e-12, 0, delta_total)
  )
if (
  nrow(panel_b) != 15L ||
  anyNA(panel_b$mode_label) ||
  anyNA(panel_b$S_facet) ||
  max(abs(panel_b$delta_struct_plot + panel_b$delta_scale_plot -
          panel_b$delta_total_plot)) > 1e-12 ||
  max(abs(panel_b$decomposition_error)) > 1e-12
) {
  stop("Binary common-scale decomposition failed validation")
}

theme_figure3 <- function(base_size = 8.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = col_grid, linewidth = 0.28),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = col_border, fill = NA, linewidth = 0.35
      ),
      strip.background = element_rect(
        fill = col_strip, colour = col_border, linewidth = 0.35
      ),
      strip.text = element_text(
        colour = col_text, size = 8.1,
        margin = margin(2.2, 2.2, 2.2, 2.2)
      ),
      axis.text = element_text(colour = col_text, size = 7.6),
      axis.title = element_text(colour = col_text, size = 8.6),
      axis.ticks = element_line(colour = "#5B6770", linewidth = 0.3),
      legend.title = element_blank(),
      legend.text = element_text(colour = col_text, size = 7.7),
      legend.key = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

# a: six DAG-role distributions, with one shared scale and no near-zero zoom.
# The four box fills identify the direct genetic-source configurations.
p_a <- ggplot(panel_a, aes(x = config_label, y = response, fill = config_label)) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    linewidth = 0.4, colour = "#69737A"
  ) +
  geom_boxplot(
    width = 0.60,
    colour = col_box_edge,
    linewidth = 0.38,
    outlier.shape = 16,
    outlier.size = 0.28,
    outlier.alpha = 0.25
  ) +
  stat_summary(
    fun = mean, geom = "point",
    shape = 16, size = 1.15, colour = col_text
  ) +
  facet_wrap(~ role_label, ncol = 3, scales = "fixed") +
  scale_fill_manual(values = config_fills, guide = "none") +
  scale_y_continuous(
    breaks = c(-0.10, 0, 0.10, 0.20, 0.30),
    labels = scales::label_number(accuracy = 0.01),
    limits = c(-0.15, 0.37),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Simulated binary-outcome response",
    x = NULL,
    y = "Total marginal logOR response"
  ) +
  theme_figure3() +
  theme(
    panel.spacing = grid::unit(2.3, "mm"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      size = 6.7, lineheight = 0.92,
      margin = margin(t = 2)
    ),
    axis.text.y = element_text(size = 7.1),
    axis.title.y = element_text(margin = margin(r = -8), vjust = -5),
    plot.title = element_text(
      face = "bold", size = 8.8, colour = col_text,
      margin = margin(b = 3)
    ),
    plot.margin = margin(4, 3, 4, 3)
  )

# b: from zero to the structural component, then to the total response.
# Arrow directions retain the sign of both additive components.
struct_segments <- panel_b |>
  filter(abs(delta_struct_plot) > 1e-12)
scale_segments <- panel_b |>
  filter(abs(delta_scale_plot) > 1e-12)
component_colours <- c(
  "Structural component" = "#619CC3",
  "Residual-scale component" = "#E88B3A",
  "Total response" = "#BE1725"
)
p_b <- ggplot(panel_b, aes(y = mode_label)) +
  geom_vline(
    xintercept = 0, linetype = "dashed",
    linewidth = 0.36, colour = "#69737A"
  ) +
  geom_segment(
    data = struct_segments,
    aes(
      x = 0, xend = delta_struct_plot,
      y = mode_label, yend = mode_label,
      colour = "Structural component"
    ),
    linewidth = 0.76, lineend = "round",
    arrow = grid::arrow(length = grid::unit(1.18, "mm"), type = "closed")
  ) +
  geom_segment(
    data = scale_segments,
    aes(
      x = delta_struct_plot, xend = delta_total_plot,
      y = mode_label, yend = mode_label,
      colour = "Residual-scale component"
    ),
    linewidth = 0.76, lineend = "round",
    arrow = grid::arrow(length = grid::unit(1.18, "mm"), type = "closed")
  ) +
  geom_point(
    data = panel_b,
    aes(x = delta_total_plot, y = mode_label, colour = "Total response"),
    shape = 18, size = 2.45
  ) +
  facet_grid(. ~ S_facet) +
  scale_colour_manual(
    values = component_colours,
    breaks = names(component_colours),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-0.30, -0.15, 0, 0.15, 0.30, 0.45),
    labels = scales::label_number(accuracy = 0.01),
    limits = c(-0.34, 0.53),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(limits = rev(unname(mode_labels[mode_levels]))) +
  labs(
    title = "Structural–residual decomposition under positive common scaling",
    x = "logOR coefficient response", y = NULL
  ) +
  guides(
    colour = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(
        linewidth = c(0.8, 0.8, 0),
        shape = c(NA, NA, 18),
        size = c(NA, NA, 2.5)
      )
    )
  ) +
  theme_figure3() +
  theme(
    panel.grid.major.y = element_line(
      colour = "#E8EDF0", linewidth = 0.22, linetype = "dashed"
    ),
    panel.grid.major.x = element_line(colour = "#E8EDF0", linewidth = 0.22),
    axis.text.x = element_text(size = 7.0),
    axis.text.y = element_text(size = 7.2),
    panel.spacing.x = grid::unit(2.5, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(8, "mm"),
    legend.key.height = grid::unit(4, "mm"),
    legend.spacing.x = grid::unit(1.7, "mm"),
    legend.margin = margin(3, 6, 3, 6),
    legend.background = element_rect(
      fill = "white", colour = "#D7DFE4", linewidth = 0.25
    ),
    legend.text = element_text(size = 7.1, colour = col_text),
    plot.title = element_text(
      face = "bold", size = 8.8, colour = col_text,
      margin = margin(b = 4)
    ),
    plot.margin = margin(3, 3, 3, 3)
  )

fig3 <- (p_a / p_b) +
  plot_layout(heights = c(1.47, 1.00)) +
  plot_annotation(tag_levels = "a")

readr::write_csv(
  panel_a |>
    arrange(role, instrument_config, replicate),
  file.path(source_dir, "Figure3a_binary_replicate_response_source.csv")
)
readr::write_csv(
  summary_check |>
    arrange(role, instrument_config),
  file.path(source_dir, "Figure3a_binary_mean_check.csv")
)
readr::write_csv(
  panel_b |>
    arrange(S, residual_mode),
  file.path(source_dir, "Figure3b_binary_scale_calibration_source.csv")
)

width_in <- 7.20
height_in <- 7.10
pdf_file <- file.path(output_dir, paste0(output_stem, ".pdf"))
png_file <- file.path(output_dir, paste0(output_stem, ".png"))
tiff_file <- file.path(output_dir, paste0(output_stem, ".tiff"))
svg_file <- file.path(output_dir, paste0(output_stem, ".svg"))

if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    pdf_file, width = width_in, height = height_in,
    family = "sans", bg = "white"
  )
} else {
  grDevices::pdf(
    pdf_file, width = width_in, height = height_in,
    useDingbats = FALSE, bg = "white"
  )
}
print(fig3)
grDevices::dev.off()

if (requireNamespace("ragg", quietly = TRUE)) {
  ggplot2::ggsave(
    png_file, fig3, width = width_in, height = height_in,
    units = "in", dpi = 600, device = ragg::agg_png, bg = "white"
  )
  ggplot2::ggsave(
    tiff_file, fig3, width = width_in, height = height_in,
    units = "in", dpi = 600, device = ragg::agg_tiff,
    compression = "lzw", bg = "white"
  )
} else {
  ggplot2::ggsave(
    png_file, fig3, width = width_in, height = height_in,
    units = "in", dpi = 600, bg = "white"
  )
  ggplot2::ggsave(
    tiff_file, fig3, width = width_in, height = height_in,
    units = "in", dpi = 600, compression = "lzw", bg = "white"
  )
}
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(
    svg_file, width = width_in, height = height_in, bg = "white"
  )
  print(fig3)
  grDevices::dev.off()
}

message("Figure 3 v10 three-facet structural-residual decomposition complete")
message("PDF: ", pdf_file)
message("PNG: ", png_file)
message("Source data: ", source_dir)
