#!/usr/bin/env Rscript

# ==========================================================
# --- ENVIRONMENT CONFIGURATION ---
# ==========================================================
local_lib <- Sys.getenv("R_LIBS_USER")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(local_lib)

required_packages <- c("mgcv", "ggplot2", "dplyr", "brglm2", "scales", "RColorBrewer")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ==========================================================
# --- GLOBAL CONFIGURATION (EDIT HERE) ---
# ==========================================================
# Default scenarios (used if NO arguments are passed via terminal)
SCENARIOS <- list(
  "Global Fibrosis" = "./2D/full/post_processing/data_for_gam_r.csv",
  "Focal Lesion"    = "./2D/ellipse/post_processing/data_for_gam_r.csv"
  # Add more here in the future if running directly from RStudio:
  # "Global 3D"     = "./3D/full/post_processing/data_for_gam_r.csv"
)

# Colors for scenarios (Will auto-generate if a new scenario is added and not listed here)
SCENARIO_COLORS <- c(
  "Global Fibrosis" = "#1f77b4",
  "Focal Lesion"    = "#ff7f0e",
  "Global 3D"       = "#2ca02c"
)

# Define the baseline for statistical comparisons
REF_FIBROSIS <- "Diffuse"
REF_SCENARIO <- "Global Fibrosis" # Must match one of the scenario names

# Fibrosis visual properties
FIB_ORDER <- c("Stochastic", "Compact", "Diffuse", "Interstitial", "Patchy")
FIB_COLORS <- c(
  "Stochastic"   = "#444444",
  "Compact"      = "#0000a2",
  "Diffuse"      = "#50ad9f",
  "Interstitial" = "#e9c716",
  "Patchy"       = "#bc272d"
)

# Output settings
OUTPUT_DIR <- "scenario_comparison_outputs"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(paste("Created output directory:", OUTPUT_DIR, "\n"))
}

# Helper function to capitalize first letters
capitalize <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}

# ==========================================================
# --- COMMAND LINE ARGUMENT PARSING (GENERIC N-SCENARIOS) ---
# ==========================================================
# Accepts arguments in the format: "Scenario Name|path/to/file.csv"
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  SCENARIOS <- list()
  for (arg in args) {
    parts <- strsplit(arg, "\\|")[[1]]
    if (length(parts) == 2) {
      SCENARIOS[[trimws(parts[1])]] <- trimws(parts[2])
    } else {
      stop(paste("Invalid argument format:", arg, "\nExpected format: 'Scenario Name|path/to/file.csv'"))
    }
  }
}

# Ensure Reference Scenario exists in the current run
if (!REF_SCENARIO %in% names(SCENARIOS)) {
  cat(paste("WARNING: Reference Scenario '", REF_SCENARIO, "' not found in inputs.\n"))
  REF_SCENARIO <- names(SCENARIOS)[1]
  cat(paste("Defaulting reference scenario to the first provided:", REF_SCENARIO, "\n"))
}

# ==========================================================
# --- DATA LOADING & COMBINATION ---
# ==========================================================
cat("Loading and combining datasets...\n")
data_list <- list()

for (scen_name in names(SCENARIOS)) {
  file_path <- SCENARIOS[[scen_name]]
  if (file.exists(file_path)) {
    tmp_data <- read.csv(file_path)
    tmp_data$Scenario <- scen_name

    # String Cleanup & Capitalization
    tmp_data$fibrosis_type <- capitalize(trimws(as.character(tmp_data$fibrosis_type)))
    data_list[[scen_name]] <- tmp_data
    cat(paste(" -> Loaded:", scen_name, "\n"))
  } else {
    stop(paste("Error: File not found for", scen_name, "at", file_path))
  }
}

data_combined <- do.call(rbind, data_list)

# Factorize Levels
scenario_levels <- names(SCENARIOS)
data_combined$Scenario <- factor(data_combined$Scenario, levels = scenario_levels)
data_combined$fibrosis_type <- factor(data_combined$fibrosis_type, levels = FIB_ORDER)

if (REF_FIBROSIS %in% levels(data_combined$fibrosis_type)) {
  data_combined$fibrosis_type <- relevel(data_combined$fibrosis_type, ref = REF_FIBROSIS)
} else {
  stop(paste("Error: Reference fibrosis '", REF_FIBROSIS, "' not found in combined data."))
}

data_combined$Scenario <- relevel(data_combined$Scenario, ref = REF_SCENARIO)

# Create interaction variable (Type | Scenario) for Absolute Hierarchy
data_combined$type_scenario <- as.factor(paste(data_combined$fibrosis_type, data_combined$Scenario, sep = "|"))
ref_level_absolute <- paste(REF_FIBROSIS, REF_SCENARIO, sep = "|")
data_combined$type_scenario <- relevel(data_combined$type_scenario, ref = ref_level_absolute)

# ==========================================================
# --- STATISTICAL MODELING ---
# ==========================================================
cat(paste("\nRunning Combined Firth's Model (Absolute Reference:", ref_level_absolute, ")...\n"))
robust_model <- glm(cbind(successes, failures) ~ type_scenario,
                    family = binomial(link = "logit"),
                    data = data_combined, method = "brglmFit")

ors <- exp(cbind(OR = coef(robust_model), confint(robust_model)))
or_df <- as.data.frame(ors)
or_df$type_scenario <- rownames(or_df)
or_df <- or_df %>% filter(type_scenario != "(Intercept)")
or_df$type_scenario <- gsub("type_scenario", "", or_df$type_scenario)

ref_row <- data.frame(OR = 1.0, `2.5 %` = 1.0, `97.5 %` = 1.0, type_scenario = ref_level_absolute, check.names = FALSE)
or_df <- rbind(or_df, ref_row)

split_names <- do.call(rbind, strsplit(as.character(or_df$type_scenario), "\\|"))
or_df$fibrosis_type <- factor(split_names[, 1], levels = rev(FIB_ORDER))
or_df$Scenario <- factor(split_names[, 2], levels = scenario_levels)
colnames(or_df)[c(2,3)] <- c("lower_ci", "upper_ci")

cat("Running Interaction Model...\n")
interaction_model <- glm(cbind(successes, failures) ~ fibrosis_type * Scenario,
                         family = binomial(link = "logit"),
                         data = data_combined, method = "brglmFit")

# ==========================================================
# --- GAM MODELING (PER SCENARIO) ---
# ==========================================================
cat("Running GAMs for Curve Comparison...\n")

fit_gam_and_predict <- function(df, scenario_name) {
  event_counts <- df %>% group_by(fibrosis_type) %>% summarise(Total = sum(successes), .groups="drop")
  active_types <- event_counts %>% filter(Total >= 5) %>% pull(fibrosis_type)
  rare_types <- event_counts %>% filter(Total < 5) %>% pull(fibrosis_type)

  preds_all <- data.frame()
  model <- NULL

  if (length(active_types) > 0) {
    df_active <- df %>% filter(fibrosis_type %in% active_types) %>% droplevels()
    model <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5),
                 family = quasibinomial(link = "logit"), data = df_active)

    grid <- expand.grid(
      fibrosis_type = active_types,
      density = seq(min(df$density), max(df$density), length.out = 200)
    )
    p <- predict(model, newdata = grid, type = "link", se.fit = TRUE)
    grid$prob <- model$family$linkinv(p$fit)
    grid$lower_ci <- model$family$linkinv(p$fit - 1.96 * p$se.fit)
    grid$upper_ci <- model$family$linkinv(p$fit + 1.96 * p$se.fit)
    grid$Scenario <- scenario_name
    preds_all <- rbind(preds_all, grid)
  }

  if (length(rare_types) > 0) {
    grid_rare <- expand.grid(
      fibrosis_type = rare_types,
      density = seq(min(df$density), max(df$density), length.out = 200),
      prob = 0, lower_ci = 0, upper_ci = 0, Scenario = scenario_name
    )
    preds_all <- rbind(preds_all, grid_rare)
  }
  return(list(predictions = preds_all, model = model))
}

gam_results <- list()
gam_plot_data <- data.frame()

for (scen in scenario_levels) {
  df_scen <- data_combined %>% filter(Scenario == scen)
  res <- fit_gam_and_predict(df_scen, scen)
  gam_results[[scen]] <- res$model
  gam_plot_data <- rbind(gam_plot_data, res$predictions)
}

gam_plot_data$Scenario <- factor(gam_plot_data$Scenario, levels = scenario_levels)
gam_plot_data$fibrosis_type <- factor(gam_plot_data$fibrosis_type, levels = FIB_ORDER)

# Find Peak Vulnerability Densities
peaks <- gam_plot_data %>%
  filter(prob > 0) %>%
  group_by(Scenario, fibrosis_type) %>%
  slice_max(order_by = prob, n = 1, with_ties = FALSE) %>%
  select(Scenario, Fibrosis_Pattern = fibrosis_type, Peak_Density = density, Max_Probability = prob) %>%
  arrange(Scenario, Fibrosis_Pattern)


# ==========================================================
# --- TEXT REPORT EXPORT ---
# ==========================================================
cat("Exporting Statistical Report...\n")
sink(file.path(OUTPUT_DIR, "scenario_comparison_report.txt"))
cat("==========================================================\n")
cat("          SCENARIO COMPARISON STATISTICAL REPORT          \n")
cat("==========================================================\n\n")

cat("--- MODEL 1: ABSOLUTE RISK HIERARCHY (FIRTH) ---\n")
cat(paste("Absolute Baseline:", ref_level_absolute, "(OR = 1.0)\n\n"))
print(or_df[, c("fibrosis_type", "Scenario", "OR", "lower_ci", "upper_ci")])

cat("\n----------------------------------------------------------\n")
cat("--- MODEL 2: SCENARIO-TYPE INTERACTION ANALYSIS (FIRTH) ---\n")
print(summary(interaction_model))

cat("\n----------------------------------------------------------\n")
cat("--- MODEL 3: GAM DENSITY DEPENDENCE & WINDOW SHIFT ---\n")
cat("Peak Vulnerability Density (Right-Shift Analysis):\n")
print(as.data.frame(peaks))

for (scen in scenario_levels) {
  cat(paste("\n[GAM Summary:", scen, "]\n"))
  if(!is.null(gam_results[[scen]])) print(summary(gam_results[[scen]])) else cat("No GAM fitted.\n")
}
sink()


# ==========================================================
# --- PLOTTING ---
# ==========================================================
# Dynamic missing color handling for scenarios
missing_scen <- setdiff(scenario_levels, names(SCENARIO_COLORS))
if (length(missing_scen) > 0) {
  extra_colors <- brewer.pal(max(3, length(missing_scen)), "Set2")
  names(extra_colors) <- missing_scen
  SCENARIO_COLORS <- c(SCENARIO_COLORS, extra_colors)
}

# 1. FOREST PLOT
p_forest <- ggplot(or_df, aes(x = OR, y = fibrosis_type, color = Scenario)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray20", linewidth = 1) +
  geom_pointrange(aes(xmin = lower_ci, xmax = upper_ci), position = position_dodge(width = 0.5), size = 1.2, linewidth = 1.2) +
  scale_color_manual(values = SCENARIO_COLORS) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5, 1.0, 2.0), labels = scales::label_number(accuracy = 0.001)) +
  labs(title = "Comparative Arrhythmia Risk Across Scenarios",
       x = "Odds Ratio (Log Scale)\n<-- Protective Effect | Increased Risk -->", y = "Fibrosis Pattern") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(axis.text.x = element_text(size = 14, color = "black", face="bold"),
        axis.text.y = element_text(size = 16, color = "black", face="bold"),
        axis.title.x = element_text(size = 16), legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 14), panel.grid.minor = element_blank())
ggsave(file.path(OUTPUT_DIR, "forest_plot_comparison.png"), plot = p_forest, width = 10, height = 6, dpi = 300)

# 2. INCIDENCE BAR PLOT
incidence_df <- data_combined %>% group_by(fibrosis_type, Scenario) %>% summarise(total_success = sum(successes), .groups = 'drop')
incidence_df$Scenario <- factor(incidence_df$Scenario, levels = scenario_levels)
incidence_df$fibrosis_type <- factor(incidence_df$fibrosis_type, levels = FIB_ORDER)

p_bar <- ggplot(incidence_df, aes(x = fibrosis_type, y = total_success, fill = Scenario)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color="black") +
  scale_fill_manual(values = SCENARIO_COLORS) +
  labs(title = "Overall Arrhythmia Incidence Across Scenarios", x = "Fibrosis Pattern", y = "Total Number of Sustained Reentries") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(axis.text.x = element_text(size = 14, color = "black", face="bold"),
        axis.text.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16), legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 14), panel.grid.major.x = element_blank())
ggsave(file.path(OUTPUT_DIR, "overall_incidence_comparison.png"), plot = p_bar, width = 9, height = 6, dpi = 300)

# 3. SIDE-BY-SIDE GAM CURVES
cat("Generating Side-by-Side GAM Curves...\n")
max_y_val <- min(max(gam_plot_data$upper_ci, na.rm = TRUE) * 1.10, 1.0)
num_cols <- min(3, length(scenario_levels)) # Auto-adjust layout for N scenarios

p_gam_comp <- ggplot(gam_plot_data, aes(x = density, y = prob, color = fibrosis_type, fill = fibrosis_type)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  facet_wrap(~ Scenario, ncol = num_cols) +
  scale_color_manual(values = FIB_COLORS) +
  scale_fill_manual(values = FIB_COLORS) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(0, max_y_val)) +
  labs(title = "Statistical Models of Arrhythmia Risk by Scenario",
       subtitle = "Comparison of vulnerable windows and peak probabilities",
       x = "Fibrosis Density", y = "Probability of Sustained Reentry") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(
    strip.text = element_text(size = 14, face = "bold", color="black"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.spacing = unit(2, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
ggsave(file.path(OUTPUT_DIR, "gam_curves_comparison.png"), plot = p_gam_comp, width = 6 * num_cols, height = 6, dpi = 300)

# 4. DUMBBELL PLOT (PEAK DENSITY SHIFT)
cat("Generating Peak Density Shift Plot (Dumbbell)...\n")
if (length(scenario_levels) >= 2) {
  # We only plot active types (those that exist in the peaks dataframe)
  peaks$Fibrosis_Pattern <- factor(peaks$Fibrosis_Pattern, levels = rev(FIB_ORDER))

  p_dumbbell <- ggplot(peaks, aes(x = Peak_Density, y = Fibrosis_Pattern)) +
    # The line connecting the scenarios
    geom_line(aes(group = Fibrosis_Pattern), color = "gray50", linewidth = 1.5) +
    # The points for each scenario
    geom_point(aes(color = Scenario), size = 6) +
    scale_color_manual(values = SCENARIO_COLORS) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    labs(title = "Shift in Peak Vulnerability Density",
         x = "Fibrosis Density at Peak Probability",
         y = "Fibrosis Pattern") +
    theme_minimal(base_family = "sans", base_size = 16) +
    theme(axis.text.x = element_text(size = 14, color = "black", face="bold"),
          axis.text.y = element_text(size = 16, color = "black", face="bold"),
          axis.title.x = element_text(size = 16),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 14),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed", color = "gray80"))

  ggsave(file.path(OUTPUT_DIR, "peak_density_shift.png"), plot = p_dumbbell, width = 9, height = 5, dpi = 300)
}

cat("Success! All plots and text reports saved.\n")
