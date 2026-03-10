#!/usr/bin/env Rscript

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
if (length(args) < 2) {
  stop("Usage: Rscript compare_scenarios.R <path_csv_full> <path_csv_ellipse>")
}
file_full <- args[1]
file_ellipse <- args[2]

cat("Loading data...\n")
data_full <- read.csv(file_full)
data_ellipse <- read.csv(file_ellipse)

# 1. STRING CLEANUP (Avoid hidden spaces)
data_full$fibrosis_type <- trimws(as.character(data_full$fibrosis_type))
data_ellipse$fibrosis_type <- trimws(as.character(data_ellipse$fibrosis_type))

# Add scenario labels
data_full$Scenario <- "Global Fibrosis"
data_ellipse$Scenario <- "Focal Lesion"

# Combine datasets
data_combined <- rbind(data_full, data_ellipse)
data_combined$fibrosis_type <- as.factor(data_combined$fibrosis_type)
data_combined$Scenario <- as.factor(data_combined$Scenario)

# Setup Levels and Order
scenario_levels <- c("Global Fibrosis", "Focal Lesion")
visual_order_x <- c("compact", "diffuse", "interstitial", "patchy")
data_combined$Scenario <- factor(data_combined$Scenario, levels = scenario_levels)

if ("diffuse" %in% levels(data_combined$fibrosis_type)) {
  data_combined$fibrosis_type <- relevel(data_combined$fibrosis_type, ref = "diffuse")
} else {
  stop("Error: 'diffuse' level not found in data.")
}

data_combined$type_scenario <- as.factor(paste(data_combined$fibrosis_type, data_combined$Scenario, sep = "|"))
ref_level <- "diffuse|Global Fibrosis"
data_combined$type_scenario <- relevel(data_combined$type_scenario, ref = ref_level)


# ==========================================================
# --- STATISTICAL MODELING ---
# ==========================================================

cat("\nRunning Combined Firth's Model (Reference: Diffuse Global)...\n")
robust_model <- glm(cbind(successes, failures) ~ type_scenario,
                    family = binomial(link = "logit"),
                    data = data_combined, method = "brglmFit")

ors <- exp(cbind(OR = coef(robust_model), confint(robust_model)))
or_df <- as.data.frame(ors)
or_df$type_scenario <- rownames(or_df)
or_df <- or_df %>% filter(type_scenario != "(Intercept)")
or_df$type_scenario <- gsub("type_scenario", "", or_df$type_scenario)
ref_row <- data.frame(OR = 1.0, `2.5 %` = 1.0, `97.5 %` = 1.0, type_scenario = ref_level, check.names = FALSE)
or_df <- rbind(or_df, ref_row)

split_names <- do.call(rbind, strsplit(as.character(or_df$type_scenario), "\\|"))
or_df$fibrosis_type <- factor(split_names[, 1], levels = rev(visual_order_x))
or_df$Scenario <- factor(split_names[, 2], levels = scenario_levels)
colnames(or_df)[c(2,3)] <- c("lower_ci", "upper_ci")

cat("Running Interaction Model...\n")
interaction_model <- glm(cbind(successes, failures) ~ fibrosis_type * Scenario,
                         family = binomial(link = "logit"),
                         data = data_combined, method = "brglmFit")

# --- GAM MODELING (PER SCENARIO) ---
cat("Running GAMs for Curve Comparison...\n")

fit_gam_and_predict <- function(df, scenario_name) {
  # Find active types
  event_counts <- df %>% group_by(fibrosis_type) %>% summarise(Total = sum(successes))
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

df_global <- data_combined %>% filter(Scenario == "Global Fibrosis")
gam_full <- fit_gam_and_predict(df_global, "Global Fibrosis")

df_local <- data_combined %>% filter(Scenario == "Focal Lesion")
gam_ellipse <- fit_gam_and_predict(df_local, "Focal Lesion")

gam_plot_data <- rbind(gam_full$predictions, gam_ellipse$predictions)
gam_plot_data$Scenario <- factor(gam_plot_data$Scenario, levels = scenario_levels)
gam_plot_data$fibrosis_type <- factor(gam_plot_data$fibrosis_type, levels = visual_order_x)

# Find Peak Vulnerability Densities
peaks <- gam_plot_data %>%
  filter(prob > 0) %>% # Ignore flat zero lines for rare types
  group_by(Scenario, fibrosis_type) %>%
  slice_max(order_by = prob, n = 1, with_ties = FALSE) %>%
  select(Scenario, fibrosis_type, Peak_Density = density, Max_Probability = prob) %>%
  arrange(Scenario, fibrosis_type)


# ==========================================================
# --- TEXT REPORT EXPORT ---
# ==========================================================
cat("Exporting Statistical Report...\n")
sink("scenario_comparison_report.txt")
cat("==========================================================\n")
cat("          SCENARIO COMPARISON STATISTICAL REPORT          \n")
cat("==========================================================\n\n")

cat("--- MODEL 1: ABSOLUTE RISK HIERARCHY (FIRTH) ---\n")
print(or_df[, c("fibrosis_type", "Scenario", "OR", "lower_ci", "upper_ci")])

cat("\n----------------------------------------------------------\n")
cat("--- MODEL 2: SCENARIO-TYPE INTERACTION ANALYSIS (FIRTH) ---\n")
print(summary(interaction_model))

cat("\n----------------------------------------------------------\n")
cat("--- MODEL 3: GAM DENSITY DEPENDENCE & WINDOW SHIFT ---\n")
cat("Peak Vulnerability Density (Right-Shift Analysis):\n")
print(as.data.frame(peaks))

cat("\n[GAM Summary: Global Fibrosis]\n")
if(!is.null(gam_full$model)) print(summary(gam_full$model)) else cat("No GAM fitted.\n")

cat("\n[GAM Summary: Focal Lesion]\n")
if(!is.null(gam_ellipse$model)) print(summary(gam_ellipse$model)) else cat("No GAM fitted.\n")
sink()


# ==========================================================
# --- PLOTTING ---
# ==========================================================
scenario_colors <- c("Global Fibrosis" = "#1f77b4", "Focal Lesion" = "#ff7f0e")
fibrosis_colors <- c("compact"="#0000a2", "diffuse"="#50ad9f", "interstitial"="#e9c716", "patchy"="#bc272d")

# 1. FOREST PLOT
p_forest <- ggplot(or_df, aes(x = OR, y = fibrosis_type, color = Scenario)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray20", linewidth = 1) +
  geom_pointrange(aes(xmin = lower_ci, xmax = upper_ci), position = position_dodge(width = 0.5), size = 1.2, linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5, 1.0, 2.0), labels = scales::label_number(accuracy = 0.001)) +
  labs(title = "Protective Effect of Border Zone and Healthy Tissue (2D Domain)",
       x = "Odds Ratio (Log Scale)\n<-- Protective Effect | Increased Risk -->", y = "Fibrosis Pattern") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(axis.text.x = element_text(size = 14, color = "black", face="bold"),
        axis.text.y = element_text(size = 16, color = "black", face="bold"),
        axis.title.x = element_text(size = 16), legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 14), panel.grid.minor = element_blank())
ggsave("forest_plot_comparison.png", plot = p_forest, width = 10, height = 6, dpi = 300)

# 2. INCIDENCE BAR PLOT
incidence_df <- data_combined %>% group_by(fibrosis_type, Scenario) %>% summarise(total_success = sum(successes), .groups = 'drop')
incidence_df$Scenario <- factor(incidence_df$Scenario, levels = scenario_levels)
incidence_df$fibrosis_type <- factor(incidence_df$fibrosis_type, levels = visual_order_x)

p_bar <- ggplot(incidence_df, aes(x = fibrosis_type, y = total_success, fill = Scenario)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color="black") +
  scale_fill_manual(values = scenario_colors) +
  labs(title = "Overall Arrhythmia Incidence by Scenario", x = "Fibrosis Pattern", y = "Total Number of Sustained Reentries") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(axis.text.x = element_text(size = 14, color = "black", face="bold"),
        axis.text.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16), legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 14), panel.grid.major.x = element_blank())
ggsave("overall_incidence_comparison.png", plot = p_bar, width = 9, height = 6, dpi = 300)

# 3. SIDE-BY-SIDE GAM CURVES (VULNERABLE WINDOW SHIFT)
cat("Generating Side-by-Side GAM Curves...\n")
max_y_val <- min(max(gam_plot_data$upper_ci, na.rm = TRUE) * 1.10, 1.0) # Synchronized dynamic Y limit

p_gam_comp <- ggplot(gam_plot_data, aes(x = density, y = prob, color = fibrosis_type, fill = fibrosis_type)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +

  # FACET WRAP is the magic here: creates side-by-side plots with identically scaled axes
  facet_wrap(~ Scenario, ncol = 2) +

  scale_color_manual(values = fibrosis_colors) +
  scale_fill_manual(values = fibrosis_colors) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(0, max_y_val)) +

  labs(title = "Statistical Models of Arrhythmia Risk by Scenario",
       subtitle = "Side-by-side comparison and shift in vulnerable window",
       x = "Fibrosis Density", y = "Probability of Sustained Reentry") +

  theme_minimal(base_family = "sans", base_size = 16) +
  theme(
    strip.text = element_text(size = 14, face = "bold", color="black"), # Facet titles
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.spacing = unit(2, "lines"), # Adds space between the two charts
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("gam_curves_comparison.png", plot = p_gam_comp, width = 12, height = 6, dpi = 300)

cat("Success! All plots and text reports saved.\n")
