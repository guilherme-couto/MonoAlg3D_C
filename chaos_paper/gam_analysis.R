#!/usr/bin/env Rscript

# Check and install necessary packages
if (!require("mgcv", quietly = TRUE)) install.packages("mgcv", repos = "http://cran.us.r-project.org")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "http://cran.us.r-project.org")

library(mgcv)
library(ggplot2)
library(dplyr)

# --- Configuration ---
# You need to pass the directory where the CSV is located or set it manually
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Please provide the path to the 'post_processing' directory created by the Python script.")
}
input_dir <- args[1]

input_dir <- "./2D/full/post_processing"
file_path <- file.path(input_dir, "data_for_gam_r_window.csv")

if (!file.exists(file_path)) {
  stop(paste("File not found:", file_path))
}

cat("Loading Filtered Data (Window):", file_path, "\n")
data <- read.csv(file_path)
data$fibrosis_type <- as.factor(data$fibrosis_type)

# Relevel to Diffuse (Reference)
data$fibrosis_type <- relevel(data$fibrosis_type, ref = "diffuse")

# --- GAM MODEL (Vulnerable Window) ---
cat("Fitting GAM on Vulnerable Window (0.1 - 0.55)...\n")
model <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5),
             family = quasibinomial(link = "logit"),
             data = data)
# Note: Reduzi k=10 para k=5, pois temos menos pontos no eixo X agora.

# --- REPORT ---
summary_file <- file.path(input_dir, "gam_window_report.txt")
sink(summary_file)
cat("=== GAM ANALYSIS: VULNERABLE WINDOW (0.1 - 0.55) ===\n")
cat("NOTE: Zeros outside this range were excluded to improve model fit.\n\n")
print(summary(model))
sink()

# --- PLOT ---
new_data <- expand.grid(
  fibrosis_type = unique(data$fibrosis_type),
  density = seq(min(data$density), max(data$density), length.out = 200)
)
preds <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
new_data$prob <- model$family$linkinv(preds$fit)
new_data$ci_lower <- model$family$linkinv(preds$fit - 1.96 * preds$se.fit)
new_data$ci_upper <- model$family$linkinv(preds$fit + 1.96 * preds$se.fit)

p <- ggplot(new_data, aes(x = density, y = prob, color = fibrosis_type, fill = fibrosis_type)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("compact"="#1f77b4", "diffuse"="#ff7f0e", "interstitial"="#2ca02c", "patchy"="#d62728")) +
  scale_fill_manual(values = c("compact"="#1f77b4", "diffuse"="#ff7f0e", "interstitial"="#2ca02c", "patchy"="#d62728")) +
  labs(title = "Statistical Model (Vulnerable Window Only)",
       x = "Fibrosis Density", y = "Probability") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 0.3))

ggsave(file.path(input_dir, "gam_window_fit.png"), plot = p, width = 8, height = 6, dpi = 300)
cat("Done.\n")
