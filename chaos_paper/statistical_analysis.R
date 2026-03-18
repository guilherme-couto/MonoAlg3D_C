#!/usr/bin/env Rscript

# ==========================================================
# --- ENVIRONMENT CONFIGURATION (WSL/LINUX) ---
# ==========================================================
local_lib <- Sys.getenv("R_LIBS_USER")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(local_lib)

required_packages <- c("mgcv", "ggplot2", "dplyr", "brglm2", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ==========================================================
# --- GLOBAL VARIABLES & STYLING ---
# ==========================================================
# Change these values to adapt the script to different analyses
REF_FIBROSIS   <- "Stochastic" # The baseline for Firth's Regression
MIN_EVENTS_GAM <- 5            # Minimum arrhythmic events to fit a GAM curve

# Order for plotting
FIB_ORDER <- c("Stochastic", "Compact", "Diffuse", "Interstitial", "Patchy")

# Color palette map
FIB_COLORS <- c(
  "Stochastic"   = "#444444",
  "Compact"      = "#0000a2",
  "Diffuse"      = "#50ad9f",
  "Interstitial" = "#e9c716",
  "Patchy"       = "#bc272d"
)

# Helper function to capitalize first letters
capitalize <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}

# ==========================================================
# --- DATA LOADING & PREPARATION ---
# ==========================================================
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) stop("Error: Provide output directory.")
input_dir <- args[1]
file_path <- file.path(input_dir, "data_for_gam_r.csv")

if (!file.exists(file_path)) stop(paste("Error: File not found at", file_path))

cat("--- STARTING GENERIC ANALYSIS ---\n")
data <- read.csv(file_path)

# Clean, Capitalize, and Factorize Fibrosis Types
data$fibrosis_type <- capitalize(trimws(as.character(data$fibrosis_type)))
data$fibrosis_type <- factor(data$fibrosis_type, levels = FIB_ORDER)

# Set Reference Level
if (REF_FIBROSIS %in% levels(data$fibrosis_type)) {
  data$fibrosis_type <- relevel(data$fibrosis_type, ref = REF_FIBROSIS)
} else {
  stop(paste("Error: Reference level", REF_FIBROSIS, "not found in data."))
}

# ==========================================================
# --- STEP 1: AUTOMATIC DIAGNOSIS ---
# ==========================================================
event_counts <- data %>%
  group_by(fibrosis_type) %>%
  summarise(Total_Successes = sum(successes), .groups = 'drop')

active_types <- event_counts %>% filter(Total_Successes >= MIN_EVENTS_GAM) %>% pull(fibrosis_type)
rare_types   <- event_counts %>% filter(Total_Successes < MIN_EVENTS_GAM) %>% pull(fibrosis_type)

cat("\n[DIAGNOSIS]\n")
cat("Active Types (GAM Curve):", paste(active_types, collapse=", "), "\n")
cat("Rare Types (Global Risk Only):", paste(rare_types, collapse=", "), "\n")

# ==========================================================
# --- STEP 2: GLOBAL RISK ANALYSIS (FIRTH/BRGLM) ---
# ==========================================================
cat("\n=== RUNNING FIRTH'S REGRESSION (BRGLM2) ===\n")
total_success  <- sum(data$successes)
total_failures <- sum(data$failures)

robust_model <- NULL
odds_ratios  <- NULL

if (total_success > 0 && total_failures > 0) {
  robust_model <- glm(cbind(successes, failures) ~ fibrosis_type,
                      family = binomial(link = "logit"),
                      data = data, method = "brglmFit")

  odds_ratios <- exp(cbind(OR = coef(robust_model), confint(robust_model)))
} else {
  cat("WARNING: Monotone dataset (100% or 0% arrhythmia). Firth skipped.\n")
}

# ==========================================================
# --- STEP 3: SHAPE ANALYSIS & PEAK DETECTION (GAM) ---
# ==========================================================
cat("\n=== RUNNING GAM (CURVES) ===\n")
model_gam <- NULL
gam_plot_data <- data.frame()
peaks <- data.frame()

if (length(active_types) > 0) {
  data_gam <- data %>% filter(fibrosis_type %in% active_types) %>% droplevels()

  tryCatch({
    model_gam <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5),
             family = quasibinomial(link = "logit"),
             data = data_gam)

    # Generate continuous predictions for plotting and peak detection
    grid <- expand.grid(
      fibrosis_type = active_types,
      density = seq(min(data_gam$density), max(data_gam$density), length.out = 200)
    )
    p <- predict(model_gam, newdata = grid, type = "link", se.fit = TRUE)
    grid$prob <- model_gam$family$linkinv(p$fit)
    grid$ci_lower <- model_gam$family$linkinv(p$fit - 1.96 * p$se.fit)
    grid$ci_upper <- model_gam$family$linkinv(p$fit + 1.96 * p$se.fit)

    gam_plot_data <- grid

    # Peak Detection (Req 3: Max probability density)
    peaks <- gam_plot_data %>%
      filter(prob > 0) %>%
      group_by(fibrosis_type) %>%
      slice_max(order_by = prob, n = 1, with_ties = FALSE) %>%
      select(Fibrosis_Pattern = fibrosis_type, Peak_Density = density, Max_Probability = prob) %>%
      arrange(Fibrosis_Pattern)

  }, error = function(e) {
    cat("Error fitting GAM:", e$message, "\n")
  })
} else {
  cat("GAM SKIPPED (No types with enough events).\n")
}

# Add flat lines for rare types (for plotting consistency)
if (length(rare_types) > 0) {
  rare_data <- expand.grid(
    fibrosis_type = factor(rare_types, levels = FIB_ORDER),
    density = seq(min(data$density), max(data$density), length.out = 200),
    prob = 0, ci_lower = 0, ci_upper = 0
  )
  gam_plot_data <- bind_rows(gam_plot_data, rare_data)
}
gam_plot_data$fibrosis_type <- factor(gam_plot_data$fibrosis_type, levels = FIB_ORDER)

# ==========================================================
# --- REPORT GENERATION ---
# ==========================================================
cat("\nExporting Statistical Report...\n")
summary_file <- file.path(input_dir, "statistical_analysis_report.txt")
sink(summary_file)

cat("==========================================================\n")
cat("       AUTOMATED STATISTICAL FRAMEWORK (GAM + FIRTH)      \n")
cat("==========================================================\n\n")

cat("1. DIAGNOSTIC SUMMARY (Absolute Arrhythmia Incidence)\n")
print(as.data.frame(event_counts))
cat("\nAnalysis Strategy:\n")
cat(" - Global Risk Comparison: Firth's Bias Reduction\n")
cat(" - Density Curve Shape: GAM Splines (Active Types Only)\n\n")

cat("----------------------------------------------------------\n")
cat("2. GLOBAL RISK HIERARCHY (FIRTH RESULTS)\n")
if (!is.null(odds_ratios)) {
  cat(paste("Reference Baseline:", REF_FIBROSIS, "(OR = 1.0)\n\n"))
  print(odds_ratios)
  cat("\n(Significance p-values):\n")
  print(coef(summary(robust_model))[,4])
} else {
  cat("Model not fitted (Monotone data).\n")
}

cat("\n----------------------------------------------------------\n")
cat("3. DENSITY DEPENDENCE & PEAK VULNERABILITY (GAM RESULTS)\n")
if (nrow(peaks) > 0) {
  cat("Estimated Peak Vulnerability (Density at max probability):\n")
  print(as.data.frame(peaks))
  cat("\n[GAM Summary]\n")
  print(summary(model_gam))
} else {
  cat("No density curves modeled (Insufficient events).\n")
}
sink()

# ==========================================================
# --- PLOTTING ---
# ==========================================================
cat("\nGenerating Plots...\n")

# PLOT 1: GAM CURVES
if (nrow(gam_plot_data) > 0) {
  max_y_val <- max(gam_plot_data$ci_upper, na.rm = TRUE)
  y_limit_top <- min(max(max_y_val * 1.1, 0.1), 1.0)

  p_gam <- ggplot(gam_plot_data, aes(x = density, y = prob, color = fibrosis_type, fill = fibrosis_type)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
    scale_color_manual(values = FIB_COLORS) +
    scale_fill_manual(values = FIB_COLORS) +
    labs(title = "Statistical Model of Arrhythmia Risk", x = "Fibrosis Density", y = "Probability of Sustained Reentry") +
    theme_minimal(base_family = "sans", base_size = 16) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme(
      axis.text = element_text(size = 14, color = "black", face = "bold"),
      axis.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.position = c(0.15, 0.8),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.box.background = element_rect(color = "black")
    ) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(0, y_limit_top))

  ggsave(file.path(input_dir, "statistical_model_fit.png"), plot = p_gam, width = 9, height = 6, dpi = 300)
}

# PLOT 2: FOREST PLOT (Odds Ratios)
if (!is.null(odds_ratios)) {
  or_df <- as.data.frame(odds_ratios)
  or_df$fibrosis_type <- rownames(or_df)
  or_df <- or_df %>% filter(fibrosis_type != "(Intercept)")
  or_df$fibrosis_type <- gsub("fibrosis_type", "", or_df$fibrosis_type)

  ref_row <- data.frame(OR = 1.0, `2.5 %` = 1.0, `97.5 %` = 1.0, fibrosis_type = REF_FIBROSIS, check.names = FALSE)
  or_df <- rbind(or_df, ref_row)

  # Y-axis order (Top to Bottom mapping)
  or_df$fibrosis_type <- factor(or_df$fibrosis_type, levels = rev(FIB_ORDER))
  colnames(or_df)[c(2,3)] <- c("lower_ci", "upper_ci")

  p_forest <- ggplot(or_df, aes(x = OR, y = fibrosis_type, color = fibrosis_type)) +
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 1) +
    geom_pointrange(aes(xmin = lower_ci, xmax = upper_ci), size = 1.2, linewidth = 1.2) +
    scale_color_manual(values = FIB_COLORS) +
    scale_x_log10(breaks = c(0.01, 0.1, 0.5, 1.0, 2.0), labels = scales::label_number(accuracy = 0.01)) +
    labs(title = "Relative Arrhythmia Risk", subtitle = "Odds Ratios with 95% Confidence Intervals",
         x = "Odds Ratio (Log Scale)\n<-- Protective Effect | Increased Risk -->", y = "Fibrosis Pattern") +
    theme_minimal(base_family = "sans", base_size = 14) +
    theme(
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 16, color = "black", face="bold"),
      axis.title.x = element_text(size = 16),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(input_dir, "odds_ratio_forest_plot.png"), plot = p_forest, width = 8, height = 5, dpi = 300)
}

# PLOT 3: OVERALL INCIDENCE BAR PLOT (Absolute Reentries)
event_counts$fibrosis_type <- factor(event_counts$fibrosis_type, levels = FIB_ORDER)

p_bar <- ggplot(event_counts, aes(x = fibrosis_type, y = Total_Successes, fill = fibrosis_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = FIB_COLORS) +
  labs(title = "Overall Arrhythmia Incidence",
       x = "Fibrosis Pattern",
       y = "Total Number of Sustained Reentries") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )
ggsave(file.path(input_dir, "overall_incidence_barplot.png"), plot = p_bar, width = 8, height = 6, dpi = 300)

cat("GAM + Firth Analysis Completed Successfully.\n")
