#!/usr/bin/env Rscript

# --- CONFIGURAÇÃO DE BIBLIOTECA LOCAL (WSL/LINUX) ---
local_lib <- Sys.getenv("R_LIBS_USER")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(local_lib)

# --- PACKAGES ---
required_packages <- c("mgcv", "ggplot2", "dplyr", "brglm2", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# --- INPUT CONFIGURATION ---
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) stop("Error: Provide output directory.")
input_dir <- args[1]
file_path <- file.path(input_dir, "data_for_gam_r.csv")

if (!file.exists(file_path)) stop("Error: File data_for_gam_r.csv not found.")

cat("--- STARTING GENERIC ANALYSIS ---\n")
data <- read.csv(file_path)
data$fibrosis_type <- as.factor(data$fibrosis_type)

# Define Reference
if ("diffuse" %in% levels(data$fibrosis_type)) {
  data$fibrosis_type <- relevel(data$fibrosis_type, ref = "diffuse")
}

# --- STEP 1: AUTOMATIC DIAGNOSIS ---
event_counts <- data %>%
  group_by(fibrosis_type) %>%
  summarise(Total_Successes = sum(successes))

print(as.data.frame(event_counts))

# Cutoff criterion for GAM curve (less than 5 points makes curve plotting impossible)
MIN_EVENTS_FOR_GAM <- 5
active_types <- event_counts %>% filter(Total_Successes >= MIN_EVENTS_FOR_GAM) %>% pull(fibrosis_type)
rare_types <- event_counts %>% filter(Total_Successes < MIN_EVENTS_FOR_GAM) %>% pull(fibrosis_type)

cat("\n[DIAGNOSIS]\n")
cat("Active Types (GAM Curve):", paste(active_types, collapse=", "), "\n")
cat("Rare Types (Global Risk Only):", paste(rare_types, collapse=", "), "\n")

# --- STEP 2: GLOBAL RISK ANALYSIS (FIRTH/BRGLM) ---
# Always run for all types. This is the basis for statistical comparison.
cat("\n=== RUNNING FIRTH'S REGRESSION (BRGLM2) ===\n")

# Check if there is at least 1 success and 1 failure in the entire dataset
total_success <- sum(data$successes)
total_failures <- sum(data$failures)

if (total_success > 0 && total_failures > 0) {
  # Model: Arrhythmia ~ Type (Ignoring density to avoid singularity in rare types)
  robust_model <- glm(cbind(successes, failures) ~ fibrosis_type,
                      family = binomial(link = "logit"),
                      data = data,
                      method = "brglmFit") # Firth's method

  odds_ratios <- exp(cbind(OR = coef(robust_model), confint(robust_model)))
  firth_summary <- summary(robust_model)
  print(odds_ratios)
} else {
  cat("WARNING: Monotone dataset (100% arrhythmia or 0% arrhythmia). Firth skipped.\n")
  robust_model <- NULL
  odds_ratios <- NULL
}

# --- STEP 3: SHAPE ANALYSIS (GAM) ---
# Run ONLY if there are enough active types.

model_gam <- NULL
if (length(active_types) > 0) {
  cat("\n=== RUNNING GAM (CURVES) ===\n")

  # Filter only those with data for curve
  data_gam <- data %>% filter(fibrosis_type %in% active_types)
  data_gam$fibrosis_type <- droplevels(data_gam$fibrosis_type)

  # Try running GAM. If error (e.g., too few points), catch with tryCatch
  tryCatch({
    model_gam <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5),
                 family = quasibinomial(link = "logit"),
                 data = data_gam)
    print(summary(model_gam))
  }, error = function(e) {
    cat("Error fitting GAM (probably insufficient data for splines):", e$message, "\n")
  })
} else {
  cat("\n=== GAM SKIPPED (No types with enough events for curve) ===\n")
}

# --- REPORT GENERATION ---
summary_file <- file.path(input_dir, "statistical_analysis_report.txt")
sink(summary_file)

cat("==========================================================\n")
cat("       AUTOMATED STATISTICAL FRAMEWORK (GAM + FIRTH)      \n")
cat("==========================================================\n\n")

cat("1. DIAGNOSTIC SUMMARY\n")
print(as.data.frame(event_counts))
cat("\nAnalysis Strategy:\n")
cat(" - Global Risk Comparison: Firth's Bias Reduction (All Types)\n")
cat(" - Density Curve Shape: GAM Splines (Only Types with >", MIN_EVENTS_FOR_GAM, "events)\n\n")

cat("----------------------------------------------------------\n")
cat("2. GLOBAL RISK HIERARCHY (FIRTH RESULTS)\n")
if (!is.null(odds_ratios)) {
  print(odds_ratios)
  cat("\n(Significance):\n")
  print(coef(summary(robust_model))[,4])
} else {
  cat("Model not fitted (Monotone data).\n")
}

cat("\n----------------------------------------------------------\n")
cat("3. DENSITY DEPENDENCE (GAM RESULTS)\n")
if (!is.null(model_gam)) {
  print(summary(model_gam))
} else {
  cat("No density curves modeled (Insufficient events).\n")
}
sink()

# ==========================================================
# --- PLOTTING ---
# ==========================================================
cat("\nGenerating Generic Plots...\n")
plot_data <- NULL

# 1. GAM data (if exists)
if (!is.null(model_gam)) {
  new_data_gam <- expand.grid(
    fibrosis_type = unique(data_gam$fibrosis_type),
    density = seq(min(data_gam$density), max(data_gam$density), length.out = 200)
  )
  preds <- predict(model_gam, newdata = new_data_gam, type = "link", se.fit = TRUE)
  new_data_gam$prob <- model_gam$family$linkinv(preds$fit)
  new_data_gam$ci_lower <- model_gam$family$linkinv(preds$fit - 1.96 * preds$se.fit)
  new_data_gam$ci_upper <- model_gam$family$linkinv(preds$fit + 1.96 * preds$se.fit)
  plot_data <- new_data_gam
}

# 2. "Flat" data for rare types
if (length(rare_types) > 0) {
  cat("Adding flat lines for rare types:", paste(rare_types, collapse=", "), "\n")
  rare_data <- expand.grid(
    fibrosis_type = factor(rare_types, levels=levels(data$fibrosis_type)),
    density = seq(min(data$density), max(data$density), length.out = 200),
    prob = 0, ci_lower = 0, ci_upper = 0
  )

  if (is.null(plot_data)) {
    plot_data <- rare_data
  } else {
    plot_data <- bind_rows(plot_data, rare_data)
  }
}

# Only plot if there is something to plot
if (!is.null(plot_data)) {

  max_y_val <- max(plot_data$ci_upper, na.rm = TRUE)
  y_limit_top <- min(max(max_y_val * 1.1, 0.1), 1.0)

  p <- ggplot(plot_data, aes(x = density, y = prob, color = fibrosis_type, fill = fibrosis_type)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("compact"="#0000a2", "diffuse"="#50ad9f", "interstitial"="#e9c716", "patchy"="#bc272d")) +
    scale_fill_manual(values = c("compact"="#0000a2", "diffuse"="#50ad9f", "interstitial"="#e9c716", "patchy"="#bc272d")) +

    labs(title = "Statistical Model of Arrhythmia Risk",
         x = "Fibrosis Density", y = "Probability of Sustained Reentry") +

    labs(title = "Statistical Model of Arrhythmia Risk",
         x = "Fibrosis Density", y = "Probability of Sustained Reentry") +

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

  ggsave(file.path(input_dir, "statistical_model_fit.png"), plot = p, width = 9, height = 6, dpi = 300)
}

# --- FOREST PLOT (Odds Ratios) ---
if (!is.null(odds_ratios)) {
  cat("Generating Forest Plot...\n")

  or_df <- as.data.frame(odds_ratios)
  or_df$fibrosis_type <- rownames(or_df)

  or_df <- or_df %>% filter(fibrosis_type != "(Intercept)")
  or_df$fibrosis_type <- gsub("fibrosis_type", "", or_df$fibrosis_type)

  # Add reference (Diffuse = 1.0)
  ref_row <- data.frame(OR = 1.0, `2.5 %` = 1.0, `97.5 %` = 1.0, fibrosis_type = "diffuse", check.names = FALSE)
  or_df <- rbind(or_df, ref_row)

  or_df$fibrosis_type <- factor(or_df$fibrosis_type, levels = rev(c("compact", "diffuse", "interstitial", "patchy")))

  colnames(or_df)[2] <- "lower_ci"
  colnames(or_df)[3] <- "upper_ci"

  p_forest <- ggplot(or_df, aes(x = OR, y = fibrosis_type, color = fibrosis_type)) +

    geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 1) +

    geom_pointrange(aes(xmin = lower_ci, xmax = upper_ci), size = 1.2, linewidth = 1.2) +

    scale_color_manual(values = c("compact"="#0000a2", "diffuse"="#50ad9f", "interstitial"="#e9c716", "patchy"="#bc272d")) +

    scale_x_log10(breaks = c(0.01, 0.1, 0.5, 1.0, 2.0), labels = scales::label_number(accuracy = 0.01)) +

    labs(title = "Relative Arrhythmia Risk",
         subtitle = "Odds Ratios with 95% Confidence Intervals",
         x = "Odds Ratio (Log Scale)\n<-- Protective Effect | Increased Risk -->",
         y = "") +

    theme_minimal(base_family = "sans", base_size = 14) +
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 16, color = "black", face="bold"),
          axis.title.x = element_text(size = 16),
          legend.position = "none",
          panel.grid.minor = element_blank())

  ggsave(file.path(input_dir, "odds_ratio_forest_plot.png"), plot = p_forest, width = 8, height = 5, dpi = 300)
}

cat("GAM + Firth Analysis Completed Successfully.\n")
