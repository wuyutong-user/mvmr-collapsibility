plot_library <- Sys.getenv("IJE_PLOTTING_LIBRARY", "")
if (nzchar(plot_library)) .libPaths(c(plot_library, .libPaths()))
library(patchwork)
# Figure 5: UK Biobank application.
# Panels show working structural roles, target-specific adjustment decisions,
# and observed total marginal logOR response with conditional instrument strength.
# Reads the adjacent frozen 168-triad input; no effects are re-estimated.
# Use a UTF-8 locale where available for labels and portable path handling.
try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)

required_pkgs <- c("ggplot2", "dplyr", "ggrepel", "svglite", "jsonlite")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

library(ggplot2)
library(dplyr)
library(ggrepel)

plot_data_dir <- "."
out_dir <- "."
data_file <- "Figure5_plotting_input_168.csv"
if (!file.exists(data_file)) {
  stop("Figure 5 plotting data were not found: ", data_file)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -------------------------------------------------------------------------
# 1. Read frozen application data and define display contracts
# -------------------------------------------------------------------------

dat <- read.csv(
  data_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA"),
  fileEncoding = "UTF-8-BOM"
)
# Older R builds retain the UTF-8 BOM in the first CSV header.
names(dat)[1] <- sub("^\\ufeff", "", names(dat)[1])

required_cols <- c(
  "triad_id", "candidate", "outcome", "raw_role", "reporting_role",
  "display_role_key", "display_role", "display_role_reason",
  "mediator_like_display_marker", "adjustment_decision",
  "display_adjustment_decision", "display_decision_reason",
  "observed_response_same_set", "F_min", "identity_error"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (nrow(dat) != 168) {
  stop("Expected 168 triads, found ", nrow(dat), ".")
}
if (anyDuplicated(dat[c("candidate", "outcome")]) > 0) {
  stop("Duplicate candidate-outcome combinations detected.")
}
if (length(unique(dat$candidate)) != 12 || length(unique(dat$outcome)) != 14) {
  stop("The application grid is not a complete 12 x 14 candidate-outcome grid.")
}

role_order <- c(
  "Confounder",
  "Collider",
  "Independent cause",
  "Upstream surrogate of exposure",
  "Downstream surrogate of outcome",
  "Downstream surrogate of exposure"
)

role_labels <- c(
  confounder = "Confounder",
  collider = "Collider",
  independent_cause = "Independent cause",
  upstream_surrogate_of_exposure = "Upstream surrogate of exposure",
  downstream_surrogate_of_outcome = "Downstream surrogate of outcome",
  downstream_surrogate_of_exposure = "Downstream surrogate of exposure"
)

role_display <- c(
  "Confounder" = "Confounder",
  "Collider" = "Collider",
  "Independent cause" = "Independent cause",
  "Upstream surrogate of exposure" = "Upstream surrogate\nof exposure",
  "Downstream surrogate of outcome" = "Downstream surrogate\nof outcome",
  "Downstream surrogate of exposure" = "Downstream surrogate\nof exposure"
)

# Distinct colours encode structural roles and target-specific decisions.
role_colors <- c(
  "Confounder" = "#7b95c6",
  "Collider" = "#49c2d9",
  "Independent cause" = "#67a583",
  "Upstream surrogate of exposure" = "#a2c986",
  "Downstream surrogate of outcome" = "#fded95",
  "Downstream surrogate of exposure" = "#ffc1a6"
)

decision_order <- c("Necessary", "Unnecessary", "Overadjustment", "Unclassified")
decision_colors <- c(
  "Necessary" = "#67a583",
  "Unnecessary" = "#d0e2c0",
  "Overadjustment" = "#f59c7c",
  "Unclassified" = "#a8c3d3"
)

candidate_order <- c(
  "WC", "HC", "SBP", "DBP", "HDL-C", "LDL-C", "TC", "TG",
  "Glucose", "HbA1c", "Smoking", "Alcohol"
)
candidate_display <- c(
  "WC" = "Waist\ncircumference",
  "HC" = "Hip\ncircumference",
  "SBP" = "Systolic blood\npressure",
  "DBP" = "Diastolic blood\npressure",
  "HDL-C" = "HDL cholesterol",
  "LDL-C" = "LDL cholesterol",
  "TC" = "Total cholesterol",
  "TG" = "Triglycerides",
  "Glucose" = "Glucose",
  "HbA1c" = "HbA1c",
  "Smoking" = "Smoking",
  "Alcohol" = "Alcohol"
)

outcome_order <- c(
  "AAA", "AF", "AVS", "CAD", "DVT", "HF", "HTN",
  "ICH", "IS", "PE", "PVD", "SAH", "TAA", "TIA"
)

dat <- dat %>%
  mutate(
    structural_role = case_when(
      display_role_key == "confounder" ~ "Confounder",
      display_role_key == "collider" ~ "Collider",
      display_role_key == "independent_cause" ~ "Independent cause",
      display_role_key == "upstream_surrogate_of_exposure" ~
        "Upstream surrogate of exposure",
      display_role_key == "downstream_surrogate_of_outcome" ~
        "Downstream surrogate of outcome",
      display_role_key == "downstream_surrogate_of_exposure" ~
        "Downstream surrogate of exposure",
      TRUE ~ NA_character_
    ),
    structural_role = factor(structural_role, levels = role_order),
    adjustment_decision = factor(
      adjustment_decision,
      levels = decision_order
    ),
    F_min = as.numeric(F_min),
    observed_response_same_set = as.numeric(observed_response_same_set),
    identity_error = as.numeric(identity_error),
    F_status = ifelse(
      !is.na(F_min) & F_min < 10,
      "F_min < 10",
      "F_min >= 10"
    ),
    candidate_plot = factor(
      candidate,
      levels = rev(candidate_order),
      labels = rev(unname(candidate_display[candidate_order]))
    ),
    outcome_plot = factor(outcome, levels = outcome_order),
    role_plot = factor(structural_role, levels = rev(role_order)),
    mediator_like_display_marker = as.logical(mediator_like_display_marker),
    example_number = case_when(
      candidate == "Glucose" & outcome == "AAA" ~ "1",
      candidate == "HDL-C" & outcome == "HF" ~ "2",
      TRUE ~ NA_character_
    )
  )

role_counts <- table(factor(dat$structural_role, levels = role_order))
decision_counts <- table(
  factor(dat$adjustment_decision, levels = decision_order)
)
F_counts <- table(
  factor(dat$F_status, levels = c("F_min >= 10", "F_min < 10"))
)
role_order_c <- role_order

expected_roles <- as.integer(c(5, 9, 28, 23, 0, 103))
expected_decisions <- as.integer(c(5, 129, 22, 12))
cat("Structural role counts observed: ", paste(as.integer(role_counts), collapse = ", "), "\n", sep = "")
cat("Adjustment decision counts observed: ", paste(as.integer(decision_counts), collapse = ", "), "\n", sep = "")
if (!identical(as.integer(role_counts), expected_roles)) {
  stop("Structural-role counts do not match the frozen application table.")
}
if (!identical(as.integer(decision_counts), expected_decisions)) {
  stop("Adjustment-decision counts do not match the frozen application table.")
}
if (sum(dat$mediator_like_display_marker) != 25) {
  stop("Expected 25 mediator-like display markers.")
}
if (sum(F_counts) != 168) {
  stop("Conditional-F status does not cover all 168 triads.")
}
if (max(abs(dat$identity_error), na.rm = TRUE) > 1e-10) {
  stop("Same-set identity error exceeds the plotting-data tolerance.")
}

role_structure <- c(
  "Confounder" = "Z -> X; Z -> Y",
  "Collider" = "X -> Z <- Y",
  "Independent cause" = "Z -> Y",
  "Upstream surrogate of exposure" = "Z -> X -> Y",
  "Downstream surrogate of outcome" = "X -> Y -> Z",
  "Downstream surrogate of exposure" = "X -> Z"
)
role_legend_labels <- c(
  "Confounder\n(n = 5)",
  "Collider\n(n = 9)",
  "Independent cause\n(n = 28)",
  "Upstream\nsurrogate of\nexposure (n = 23)",
  "Downstream\nsurrogate of\noutcome (n = 0)",
  "Downstream\nsurrogate of\nexposure (n = 103)"
)
role_legend_title <- "Working DAG roles\nM = mediator-like pathway"
decision_legend_labels <- paste0(
  decision_order,
  "\n(n = ",
  as.integer(decision_counts[decision_order]),
  ")"
)

# Keep a copy of the exact plotting data with display fields for traceability.
write.csv(
  dat,
  file.path(out_dir, "Figure5_collider_revision_20261006_source_data_168.csv"),
  row.names = FALSE,
  na = ""
)

cat("\nFigure 4 simplified display data audit\n")
cat("------------------------\n")
cat("Triads: ", nrow(dat), "\n", sep = "")
cat("Structural roles: ", paste(role_counts, collapse = ", "), "\n", sep = "")
cat("Adjustment decisions: ", paste(decision_counts, collapse = ", "), "\n", sep = "")
cat("Conditional F: ", paste(F_counts, collapse = ", "), "\n", sep = "")
cat("Maximum identity error: ", format(max(abs(dat$identity_error), na.rm = TRUE), scientific = TRUE), "\n", sep = "")

# -------------------------------------------------------------------------
# 2. Shared theme and panel helper settings
# -------------------------------------------------------------------------

font_family <- "sans"

theme_application <- theme_classic(base_size = 7.0, base_family = "sans") +
  theme(
    text = element_text(family = "sans", colour = "#20252a"),
    plot.title.position = "panel",
    plot.title = element_text(
      family = font_family, size = 9.2, face = "bold", hjust = 0,
      margin = margin(b = 0)
    ),
    plot.subtitle = element_text(
      family = font_family, size = 6.7, colour = "#69727a", hjust = 0,
      margin = margin(t = 0, b = 2)
    ),
    plot.caption = element_text(
      family = font_family, size = 5.6, colour = "#69727a", hjust = 0,
      margin = margin(t = 2)
    ),
    axis.text = element_text(family = font_family, size = 6.0, colour = "#20252a"),
    axis.title = element_text(family = font_family, size = 6.9, colour = "#20252a"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#51585e"),
    axis.line = element_line(linewidth = 0.3, colour = "#51585e"),
    legend.title = element_text(family = font_family, size = 6.3, face = "bold"),
    legend.text = element_text(family = font_family, size = 5.8),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.key.size = grid::unit(4.2, "mm"),
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.spacing.y = grid::unit(0.4, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )

example_dat <- dat %>%
  filter(!is.na(example_number)) %>%
  mutate(
    example_label = ifelse(
      example_number == "1",
      "\u2460 BMI\u2013Glucose\u2013AAA",
      "\u2461 BMI\u2013HDL cholesterol\u2013HF"
    )
  )

example_label_family <- "Arial Unicode MS"

# -------------------------------------------------------------------------
# 3. Panel a - working structural role map
# -------------------------------------------------------------------------

p_a <- ggplot(dat, aes(x = outcome_plot, y = candidate_plot)) +
  geom_tile(
    aes(fill = structural_role),
    colour = "#8c9398",
    linewidth = 0.22
  ) +
  geom_text(
    data = subset(dat, mediator_like_display_marker),
    aes(label = "M"),
    inherit.aes = TRUE,
    size = 1.65,
    fontface = "bold",
    colour = "#1b2025",
    family = font_family,
    hjust = 1.55,
    vjust = 1.65
  ) +
  scale_fill_manual(
    values = role_colors,
    breaks = role_order,
    drop = FALSE,
    labels = role_legend_labels,
    name = role_legend_title,
    guide = guide_legend(
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      order = 1,
      override.aes = list(colour = "#8c9398")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed(clip = "off") +
  labs(
    title = "a  Simplified working DAG-role map",
    subtitle = "Bidirectional MR and local conditional evidence",
    x = "Outcome",
    y = NULL
  ) +
  theme_application +
  theme(
    axis.text.x = element_text(family = font_family, angle = 45, vjust = 1, hjust = 1, size = 5.7),
    axis.text.y = element_text(family = font_family, size = 5.7, lineheight = 0.88),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(family = font_family, size = 5.3, face = "bold", lineheight = 0.9),
    legend.text = element_text(family = font_family, size = 5.0, lineheight = 0.86),
    legend.key.width = grid::unit(3.2, "mm"),
    legend.key.height = grid::unit(2.6, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm"),
    legend.spacing.y = grid::unit(0.2, "mm")
  )

# -------------------------------------------------------------------------
# 4. Panel b - target-specific adjustment decision map
# -------------------------------------------------------------------------

p_b <- ggplot(dat, aes(x = outcome_plot, y = candidate_plot)) +
  geom_tile(
    aes(fill = adjustment_decision),
    colour = "#8c9398",
    linewidth = 0.22
  ) +
  scale_fill_manual(
    values = decision_colors,
    breaks = decision_order,
    labels = decision_legend_labels,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(
      nrow = 1,
      ncol = 4,
      byrow = TRUE,
      order = 1,
      override.aes = list(colour = "#8c9398")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed(clip = "off") +
  labs(
    title = "b  Target-specific adjustment decision map",
    subtitle = "Four-class display for the BMI total-effect target",
    x = "Outcome",
    y = NULL
  ) +
  theme_application +
  theme(
    axis.text.x = element_text(family = font_family, angle = 45, vjust = 1, hjust = 1, size = 5.7),
    axis.text.y = element_text(family = font_family, size = 5.7, lineheight = 0.88),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.key.width = grid::unit(3.2, "mm"),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.spacing.x = grid::unit(0.4, "mm"),
    legend.text = element_text(family = font_family, size = 5.0, lineheight = 0.88)
  )

# -------------------------------------------------------------------------
# 5. Panel c - response and conditional instrument strength
# -------------------------------------------------------------------------

p_c <- ggplot(dat, aes(x = observed_response_same_set, y = role_plot)) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.42,
    colour = "#555d63"
  ) +
  geom_boxplot(
    aes(fill = structural_role),
    width = 0.36,
    alpha = 0.32,
    outlier.shape = NA,
    linewidth = 0.30,
    colour = "#626a70"
  ) +
  geom_point(
    aes(colour = structural_role, shape = F_status),
    position = position_jitter(width = 0, height = 0.105, seed = 606),
    size = 1.9,
    stroke = 0.58,
    alpha = 0.88
  ) +
  ggrepel::geom_text_repel(
    data = example_dat,
    aes(
      x = observed_response_same_set,
      y = role_plot,
      label = example_label
    ),
    inherit.aes = FALSE,
    seed = 606,
    size = 2.0,
    colour = "#20252a",
    family = example_label_family,
    direction = "x",
    nudge_x = 0.035,
    hjust = 0,
    box.padding = 0.22,
    point.padding = 0.16,
    min.segment.length = 0,
    segment.colour = "#626a70",
    segment.size = 0.28,
    max.overlaps = Inf
  ) +
  scale_fill_manual(values = role_colors, breaks = role_order, drop = FALSE, guide = "none") +
  scale_colour_manual(
    values = role_colors,
    breaks = role_order,
    labels = role_legend_labels,
    drop = FALSE,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("F_min >= 10" = 16, "F_min < 10" = 1),
    breaks = c("F_min >= 10", "F_min < 10"),
    labels = c("F_min >= 10", "F_min < 10"),
    name = "Conditional instrument strength"
  ) +
  scale_y_discrete(
    limits = rev(role_order_c),
    labels = role_display[rev(role_order_c)],
    drop = FALSE
  ) +
  labs(
    title = "c  Observed total marginal logOR response",
    subtitle = "Same-set BMI-only IVW versus candidate-adjusted MVMR-IVW",
    x = "Observed total marginal logOR response",
    y = NULL
  ) +
  theme_application +
  theme(
    axis.text.y = element_text(family = font_family, size = 6.0, lineheight = 0.88),
    axis.text.x = element_text(family = font_family, angle = 45, vjust = 1, hjust = 1, size = 6.0),
    panel.grid.major.x = element_line(colour = "#e3e6e8", linewidth = 0.28),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.97, 0.96),
    legend.justification = c("right", "top"),
    legend.background = element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = 0.90),
      colour = "#d2d7da",
      linewidth = 0.25
    ),
    legend.margin = margin(2, 2, 2, 2),
    legend.title = element_text(family = font_family, size = 6.1, face = "bold", lineheight = 0.92),
    legend.text = element_text(family = font_family, size = 5.8),
    legend.key.width = grid::unit(5.0, "mm"),
    legend.key.height = grid::unit(3.4, "mm"),
    plot.title = element_text(
      family = font_family, size = 9.2, face = "bold", hjust = 0,
      margin = margin(l = -15, b = 0)
    ),
    plot.margin = margin(4, 28, 4, 4)
  )

# -------------------------------------------------------------------------
# 6. Assemble the final 5:5 + full-width layout
# -------------------------------------------------------------------------

fig_width_mm = 183
fig_height_mm = 165

fig_grob <- patchwork::patchworkGrob(
  (p_a | p_b) / p_c + patchwork::plot_layout(heights = c(1.1, 0.9))
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.draw(fig_grob)
}

base_name <- file.path(out_dir, "Figure5_collider_revision_20261006")

# Cairo PDF is used here so the circled-number case labels can use the
# installed Arial Unicode MS fallback while preserving vector text.
grDevices::cairo_pdf(
  paste0(base_name, ".pdf"),
  width = fig_width_mm / 25.4,
  height = fig_height_mm / 25.4,
  family = "sans",
  onefile = TRUE
)
draw_figure()
grDevices::dev.off()

grDevices::png(
  paste0(base_name, ".png"),
  width = fig_width_mm / 25.4,
  height = fig_height_mm / 25.4,
  units = "in",
  res = 600,
  type = "cairo",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

grDevices::tiff(
  paste0(base_name, ".tiff"),
  width = fig_width_mm / 25.4,
  height = fig_height_mm / 25.4,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

old_wd <- getwd()
setwd(out_dir)
svglite::svglite(
  "Figure5_collider_revision_20261006.svg",
  width = fig_width_mm / 25.4,
  height = fig_height_mm / 25.4
)
draw_figure()
grDevices::dev.off()
setwd(old_wd)

# Layout manifest records the intended physical composition. The top row is
# explicitly 5:5; panel c spans the full figure width.
layout_manifest <- list(
  figure = "Figure 5 Collider revision — 6 October 2026",
  width_mm = fig_width_mm,
  height_mm = fig_height_mm,
  panels = list(
    a = list(row = 1, column = 1, width_fraction = 0.5, role = "structural role map"),
    b = list(row = 1, column = 2, width_fraction = 0.5, role = "adjustment decision map"),
    c = list(row = 2, column = "1:2", width_fraction = 1.0, role = "response and reliability; conditional-F legend inside panel")
  ),
  counts = list(
    necessary = as.integer(decision_counts["Necessary"]),
    unnecessary = as.integer(decision_counts["Unnecessary"]),
    overadjustment = as.integer(decision_counts["Overadjustment"]),
    unclassified = as.integer(decision_counts["Unclassified"])
  ),
  conditional_F = list(
    F_min_ge_10 = as.integer(F_counts["F_min >= 10"]),
    F_min_lt_10 = as.integer(F_counts["F_min < 10"])
  )
)
jsonlite::write_json(
  layout_manifest,
  paste0(base_name, ".layout.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  paste0(base_name, ".sessionInfo.txt")
)

cat("\nSaved Figure 5 outputs:\n")
cat("  ", normalizePath(paste0(base_name, ".pdf"), winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(paste0(base_name, ".png"), winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(paste0(base_name, ".tiff"), winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(paste0(base_name, ".svg"), winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  Plotting data: ", normalizePath(data_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
