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
REF_FIBROSIS   <- "Uniform"
REAL_REF       <- "Diffuse" # Reference when Uniform is removed from GAM
MIN_EVENTS_GAM <- 10
GAM_PADDING    <- 0.05 # Padding exclusive for GAM limits

FIB_ORDER <- c("Uniform", "Compact", "Diffuse", "Interstitial", "Patchy")
ANGLE_COL <- "fiber_angle_deg"

FIB_COLORS <- c(
  "Uniform"   = "#444444",
  "Compact"      = "#0000a2",
  "Diffuse"      = "#50ad9f",
  "Interstitial" = "#e9c716",
  "Patchy"       = "#bc272d"
)

ANGLE_COLORS <- c(
  "0"  = "#4477aa",
  "30" = "#66ccee",
  "60" = "#228833",
  "90" = "#cc3311"
)

capitalize <- function(x) { paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2))) }

# ==========================================================
# --- DATA LOADING & ARGS ---
# ==========================================================
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop("Error: Provide output directory, dim, and geom.")
input_dir <- args[1]
dim_arg   <- args[2]
geom_arg  <- args[3]

# Auto Scenario Naming
scenario_name <- paste("Scenario", paste0("(", dim_arg, "-", geom_arg, ")"))
if(dim_arg == "2D" && geom_arg == "full") scenario_name <- "Scenario I"
if(dim_arg == "2D" && geom_arg == "ellipse") scenario_name <- "Scenario II"
if(dim_arg == "3D" && geom_arg == "full") scenario_name <- "Scenario III"
if(dim_arg == "3D" && geom_arg == "ellipse") scenario_name <- "Scenario IV"

file_path <- file.path(input_dir, "data_for_gam_r.csv")
file_path_angles <- file.path(input_dir, "data_angles_for_gam_r.csv")

if (!file.exists(file_path)) stop(paste("Error: File not found at", file_path))

cat("--- STARTING STATISTICAL ANALYSIS FOR", scenario_name, "---\n")

data <- read.csv(file_path)
data$fibrosis_type <- capitalize(trimws(as.character(data$fibrosis_type)))
data$fibrosis_type <- factor(data$fibrosis_type, levels = FIB_ORDER)

if (REF_FIBROSIS %in% levels(data$fibrosis_type)) {
  data$fibrosis_type <- relevel(data$fibrosis_type, ref = REF_FIBROSIS)
} else { stop("Reference level not found.") }

has_angles <- file.exists(file_path_angles)
if (has_angles) {
  data_angles <- read.csv(file_path_angles)
  if (ANGLE_COL %in% colnames(data_angles)) { colnames(data_angles)[colnames(data_angles) == ANGLE_COL] <- "angle" }
  data_angles$fibrosis_type <- factor(capitalize(trimws(as.character(data_angles$fibrosis_type))), levels = FIB_ORDER)
  data_angles$angle <- as.factor(data_angles$angle)
}

# ==========================================================
# --- STEP 1: AUTOMATIC DIAGNOSIS & GAM WINDOW DEFINITION ---
# ==========================================================
event_counts <- data %>% group_by(fibrosis_type) %>% summarise(Total_Successes = sum(successes), .groups = 'drop')
active_types <- event_counts %>% filter(Total_Successes >= MIN_EVENTS_GAM) %>% pull(fibrosis_type)
rare_types   <- event_counts %>% filter(Total_Successes < MIN_EVENTS_GAM) %>% pull(fibrosis_type)

# Calculate dynamic GAM window based strictly on active types
gam_min <- 0; gam_max <- 1
if(length(active_types) > 0) {
    active_data <- data %>% filter(fibrosis_type %in% active_types & successes > 0)
    abs_min <- min(data$density)
    abs_max <- max(data$density)

    if(nrow(active_data) > 0) {
        core_min <- min(active_data$density)
        core_max <- max(active_data$density)
        gam_min <- max(abs_min, core_min - GAM_PADDING)
        gam_max <- min(abs_max, core_max + GAM_PADDING)
        gam_min <- round(gam_min, 2); gam_max <- round(gam_max, 2)
    }
}
# Filter data for GAM models to use only the calculated window
data_gam_filtered <- data %>% filter(density >= gam_min & density <= gam_max)

cat("\n[DIAGNOSIS]\n")
cat("Active Types (GAM Curve):", paste(active_types, collapse=", "), "\n")
cat("Rare Types (Global Risk Only):", paste(rare_types, collapse=", "), "\n")
cat("GAM Window:", paste(gam_min, gam_max), "\n")

# ==========================================================
# --- STEP 2: GLOBAL RISK ANALYSIS ---
# ==========================================================
cat("=== RUNNING FIRTH'S REGRESSION ===\n")
robust_model <- NULL; odds_ratios <- NULL
if (sum(data$successes) > 0) {
  robust_model <- glm(cbind(successes, failures) ~ fibrosis_type, family = binomial(link = "logit"), data = data, method = "brglmFit")
  odds_ratios <- exp(cbind(OR = coef(robust_model), confint(robust_model)))
}

# ==========================================================
# --- STEP 3: SHAPE ANALYSIS & PEAK DETECTION (GAM) ---
# ==========================================================
cat("=== RUNNING GAM (OVERALL CURVES) ===\n")
gam_plot_data_all <- data.frame(); gam_plot_data_real <- data.frame()
peaks_all <- data.frame(); peaks_real <- data.frame()
model_all <- NULL; model_real <- NULL

# 1. Fit for all active types (including Uniform)
if (length(active_types) > 0) {
  df_all <- data_gam_filtered %>% filter(fibrosis_type %in% active_types) %>% droplevels()
  tryCatch({
    model_all <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5), family = quasibinomial, data = df_all)
    grid_all <- expand.grid(fibrosis_type = active_types, density = seq(min(df_all$density), max(df_all$density), length.out = 200))
    p_all <- predict(model_all, newdata = grid_all, type = "link", se.fit = TRUE)
    grid_all$prob <- model_all$family$linkinv(p_all$fit); grid_all$ci_lower <- model_all$family$linkinv(p_all$fit - 1.96 * p_all$se.fit); grid_all$ci_upper <- model_all$family$linkinv(p_all$fit + 1.96 * p_all$se.fit)
    gam_plot_data_all <- grid_all
    peaks_all <- grid_all %>% filter(prob > 0) %>% group_by(fibrosis_type) %>% slice_max(prob, n=1, with_ties=F)
  }, error = function(e) {})
}

# 2. Fit for realistic types (Excluding Uniform)
real_types <- active_types[active_types != REF_FIBROSIS]
if (length(real_types) > 0) {
  df_real <- data_gam_filtered %>% filter(fibrosis_type %in% real_types) %>% droplevels()
  if(REAL_REF %in% levels(df_real$fibrosis_type)) df_real$fibrosis_type <- relevel(df_real$fibrosis_type, ref = REAL_REF)
  tryCatch({
    model_real <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5), family = quasibinomial, data = df_real)
    grid_real <- expand.grid(fibrosis_type = real_types, density = seq(min(df_real$density), max(df_real$density), length.out = 200))
    p_real <- predict(model_real, newdata = grid_real, type = "link", se.fit = TRUE)
    grid_real$prob <- model_real$family$linkinv(p_real$fit); grid_real$ci_lower <- model_real$family$linkinv(p_real$fit - 1.96 * p_real$se.fit); grid_real$ci_upper <- model_real$family$linkinv(p_real$fit + 1.96 * p_real$se.fit)
    gam_plot_data_real <- grid_real
    peaks_real <- grid_real %>% filter(prob > 0) %>% group_by(fibrosis_type) %>% slice_max(prob, n=1, with_ties=F)
  }, error = function(e) {})
}

# Append rare types
if (length(rare_types) > 0) {
  rare_data <- expand.grid(fibrosis_type = factor(rare_types, levels = FIB_ORDER), density = seq(gam_min, gam_max, length.out = 200), prob = 0, ci_lower = 0, ci_upper = 0)
  if(nrow(gam_plot_data_all) > 0) gam_plot_data_all <- bind_rows(gam_plot_data_all, rare_data)
  rare_data_real <- rare_data %>% filter(fibrosis_type != REF_FIBROSIS)
  if(nrow(rare_data_real) > 0 && nrow(gam_plot_data_real) > 0) gam_plot_data_real <- bind_rows(gam_plot_data_real, rare_data_real)
}

if(nrow(gam_plot_data_all) > 0) gam_plot_data_all$fibrosis_type <- factor(gam_plot_data_all$fibrosis_type, levels = FIB_ORDER)
if(nrow(gam_plot_data_real) > 0) gam_plot_data_real$fibrosis_type <- factor(gam_plot_data_real$fibrosis_type, levels = FIB_ORDER)

# ==========================================================
# --- STEP 4: ANGLE-SPECIFIC GAM MODELING ---
# ==========================================================
gam_angle_plot_data_all <- data.frame(); gam_angle_plot_data_real <- data.frame()
if (has_angles) {
  cat("=== RUNNING GAM (ANGLE-SPECIFIC CURVES) ===\n")
  data_angles_filtered <- data_angles %>% filter(density >= gam_min & density <= gam_max)
  angles_present <- sort(unique(as.numeric(as.character(data_angles_filtered$angle))))

  for (ang in angles_present) {
    df_ang <- data_angles_filtered %>% filter(angle == as.character(ang))
    ev_counts <- df_ang %>% group_by(fibrosis_type) %>% summarise(Total = sum(successes), .groups='drop')
    act_types <- ev_counts %>% filter(Total >= MIN_EVENTS_GAM) %>% pull(fibrosis_type)
    rar_types <- ev_counts %>% filter(Total < MIN_EVENTS_GAM) %>% pull(fibrosis_type)

    fit_angle_gam <- function(types, ref=NULL) {
      if(length(types)==0) return(data.frame())
      df_a <- df_ang %>% filter(fibrosis_type %in% types) %>% droplevels()
      if(!is.null(ref) && ref %in% levels(df_a$fibrosis_type)) df_a$fibrosis_type <- relevel(df_a$fibrosis_type, ref=ref)
      tryCatch({
        m <- gam(cbind(successes, failures) ~ fibrosis_type + s(density, by=fibrosis_type, k=5), family = quasibinomial, data = df_a)
        grid <- expand.grid(fibrosis_type = types, density = seq(min(df_a$density), max(df_a$density), length.out = 200))
        p <- predict(m, newdata = grid, type = "link", se.fit = TRUE)
        grid$prob <- m$family$linkinv(p$fit); grid$ci_lower <- m$family$linkinv(p$fit - 1.96*p$se.fit); grid$ci_upper <- m$family$linkinv(p$fit + 1.96*p$se.fit)
        grid$angle <- as.character(ang)
        return(grid)
      }, error=function(e) return(data.frame()))
    }

    gam_angle_plot_data_all <- bind_rows(gam_angle_plot_data_all, fit_angle_gam(act_types))
    if (length(rar_types) > 0) gam_angle_plot_data_all <- bind_rows(gam_angle_plot_data_all, expand.grid(fibrosis_type=rar_types, density=seq(gam_min, gam_max, length.out=200), prob=0, ci_lower=0, ci_upper=0, angle=as.character(ang)))

    act_real <- act_types[act_types != REF_FIBROSIS]; rar_real <- rar_types[rar_types != REF_FIBROSIS]
    gam_angle_plot_data_real <- bind_rows(gam_angle_plot_data_real, fit_angle_gam(act_real, ref=REAL_REF))
    if (length(rar_real) > 0) gam_angle_plot_data_real <- bind_rows(gam_angle_plot_data_real, expand.grid(fibrosis_type=rar_real, density=seq(gam_min, gam_max, length.out=200), prob=0, ci_lower=0, ci_upper=0, angle=as.character(ang)))
  }
  if(nrow(gam_angle_plot_data_all)>0) { gam_angle_plot_data_all$fibrosis_type <- factor(gam_angle_plot_data_all$fibrosis_type, levels = FIB_ORDER); gam_angle_plot_data_all$angle <- factor(gam_angle_plot_data_all$angle, levels = as.character(angles_present)) }
  if(nrow(gam_angle_plot_data_real)>0) { gam_angle_plot_data_real$fibrosis_type <- factor(gam_angle_plot_data_real$fibrosis_type, levels = FIB_ORDER); gam_angle_plot_data_real$angle <- factor(gam_angle_plot_data_real$angle, levels = as.character(angles_present)) }
}

# ==========================================================
# --- REPORT GENERATION ---
# ==========================================================
summary_file <- file.path(input_dir, "statistical_analysis_report.txt")
sink(summary_file)

cat("==========================================================\n")
cat(paste("       STATISTICAL FRAMEWORK -", toupper(scenario_name), "    \n"))
cat("==========================================================\n\n")

cat("1. DIAGNOSTIC SUMMARY & GAM WINDOW\n")
print(as.data.frame(event_counts))
cat("\nGAM Filtering Logic:\n")
cat(paste(" - Active Types (fitted):", paste(active_types, collapse=", "), "\n"))
cat(paste(" - Final GAM Analysis Window (with 0.05 padding): [", gam_min, ",", gam_max, "]\n"))

cat("\n----------------------------------------------------------\n")
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
cat("3. GAM RESULTS (MODEL: ALL PATTERNS)\n")
if(!is.null(model_all)) { print(summary(model_all)) }

cat("\n----------------------------------------------------------\n")
cat("4. GAM RESULTS (MODEL: REALISTIC PATTERNS ONLY)\n")
if (nrow(peaks_real) > 0) {
  cat("Estimated Peak Vulnerability:\n")
  print(as.data.frame(peaks_real))
  cat("\n[GAM Summary]\n")
  if(!is.null(model_real)) print(summary(model_real))
} else {
  cat("No realistic patterns to analyze.\n")
}
sink()

# ==========================================================
# --- PLOTTING ---
# ==========================================================
cat("=== PLOTTING ===\n")
plot_gam_curves <- function(df, title, filename) {
  if(nrow(df) == 0) return()
  y_top <- min(max(df$ci_upper, na.rm = TRUE) * 1.1, 1.0)
  p <- ggplot(df, aes(x=density, y=prob, color=fibrosis_type, fill=fibrosis_type)) +
    geom_line(linewidth=1.2) + geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.2, color=NA) +
    scale_color_manual(values=FIB_COLORS) + scale_fill_manual(values=FIB_COLORS) +
    labs(title=title, x="Fibrosis Density", y="Probability of Sustained Reentry") +
    theme_minimal(base_family="sans", base_size=16) +
    scale_x_continuous(labels=scales::label_number(accuracy=0.01)) + scale_y_continuous(labels=scales::label_number(accuracy=0.01)) +
    coord_cartesian(ylim = c(0, y_top)) +
    theme(axis.text=element_text(size=14, color="black", face="bold"), axis.title=element_text(size=16),
          legend.title=element_blank(), legend.text=element_text(size=14),
          legend.position=c(0.15, 0.8), legend.background=element_rect(fill="white", color="black")) +
    guides(fill="none")
  ggsave(file.path(input_dir, filename), plot=p, width=9, height=6, dpi=300)
}

# PLOT 1: GAM CURVES
plot_gam_curves(gam_plot_data_all, paste(scenario_name, "- Statistical Model"), "statistical_model_fit_all.png")
plot_gam_curves(gam_plot_data_real, paste(scenario_name, "- Statistical Model of Arrhythmia Risk"), "statistical_model_fit_realistic.png")

# PLOT 2: FOREST PLOT (Odds Ratios)
if (!is.null(odds_ratios)) {
  or_df <- as.data.frame(odds_ratios)
  or_df$fibrosis_type <- gsub("fibrosis_type", "", rownames(or_df))
  or_df <- or_df %>% filter(fibrosis_type != "(Intercept)" & fibrosis_type != REF_FIBROSIS)
  or_df$fibrosis_type <- factor(or_df$fibrosis_type, levels = rev(FIB_ORDER[FIB_ORDER != REF_FIBROSIS]))
  colnames(or_df)[c(2,3)] <- c("lower_ci", "upper_ci")
  p_forest <- ggplot(or_df, aes(x = OR, y = fibrosis_type, color = fibrosis_type)) +
    geom_vline(aes(xintercept = 1.0, linetype = "Baseline"), color = FIB_COLORS[REF_FIBROSIS], linewidth = 1.2) +
    geom_pointrange(aes(xmin = lower_ci, xmax = upper_ci), size = 1.2, linewidth = 1.2) +
    scale_color_manual(values = FIB_COLORS, guide="none") +
    scale_linetype_manual(name="", values=c("Baseline"="dashed"), labels=paste("Baseline:", REF_FIBROSIS)) +
    scale_x_log10(breaks = c(0.01, 0.1, 0.5, 1.0, 2.0), labels = scales::label_number(accuracy = 0.01)) +
    labs(title = paste(scenario_name, "- Relative Arrhythmia Risk"), subtitle = paste("Odds Ratios with 95% Confidence Intervals\nBaseline:", REF_FIBROSIS, "(OR = 1.0)"),
         x = "Odds Ratio (Log Scale)\n<-- Protective Effect | Increased Risk -->", y = "") +
    theme_minimal(base_family = "sans", base_size = 14) +
    theme(axis.text.x = element_text(size=14, color="black", face="bold"), axis.text.y = element_text(size=16, color="black", face="bold"),
          axis.title = element_text(size=16), legend.position = "none", panel.grid.minor = element_blank())
  ggsave(file.path(input_dir, "odds_ratio_forest_plot.png"), plot = p_forest, width = 8, height = 5, dpi = 300)
}

# PLOT 3: OVERALL INCIDENCE BAR PLOT
ref_val <- event_counts %>% filter(fibrosis_type == REF_FIBROSIS) %>% pull(Total_Successes)
if(length(ref_val) == 0) ref_val <- 0
bar_data <- event_counts %>% filter(fibrosis_type != REF_FIBROSIS)
bar_data$fibrosis_type <- factor(bar_data$fibrosis_type, levels = FIB_ORDER[FIB_ORDER != REF_FIBROSIS])

p_bar <- ggplot(bar_data, aes(x = fibrosis_type, y = Total_Successes, fill = fibrosis_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_hline(aes(yintercept = ref_val, linetype = "Baseline"), color = FIB_COLORS[REF_FIBROSIS], linewidth = 1.2) +
  scale_fill_manual(values = FIB_COLORS, guide = "none") +
  scale_linetype_manual(name = "", values = c("Baseline" = "dashed"), labels = paste0("Baseline: ", REF_FIBROSIS, " (", ref_val, " reentries)")) +
  labs(title = paste(scenario_name, "- Overall Arrhythmia Incidence"), x = "", y = "Total Number of Sustained Reentries") +
  theme_minimal(base_family = "sans", base_size = 16) +
  theme(axis.text.x = element_text(size=14, color="black", face="bold"), axis.text.y = element_text(size=14, color="black", face="bold"),
        axis.title = element_text(size=16), legend.position = "top", legend.justification = "left", legend.text = element_text(size=14, color=FIB_COLORS[REF_FIBROSIS]),
        panel.grid.major.x = element_blank())
ggsave(file.path(input_dir, "overall_incidence_barplot.png"), plot = p_bar, width = 8, height = 6, dpi = 300)

# ANGLE SPECIFIC PLOTS
if (has_angles) {
  plot_gam_angles <- function(df, title, filename) {
    if(nrow(df) == 0) return()
    y_top <- min(max(df$ci_upper, na.rm=T) * 1.1, 1.0)
    p <- ggplot(df, aes(x=density, y=prob, color=fibrosis_type, fill=fibrosis_type)) +
      geom_line(linewidth=1.2) + geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.2, color=NA) +
      facet_wrap(~ angle, ncol=2, labeller=labeller(angle=function(x) paste0("Fiber Orientation: ", x, "°"))) +
      scale_color_manual(values=FIB_COLORS) + scale_fill_manual(values=FIB_COLORS) +
      coord_cartesian(ylim=c(0, y_top)) +
      labs(title=title, x="Fibrosis Density", y="Probability of Sustained Reentry") +
      theme_minimal(base_family="sans", base_size=16) +
      scale_x_continuous(labels=scales::label_number(accuracy=0.01)) + scale_y_continuous(labels=scales::label_number(accuracy=0.01)) +
      theme(strip.text=element_text(size=14, face="bold", color="black"), axis.text=element_text(size=14, color="black", face="bold"),
            axis.title=element_text(size=16), legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=14),
            panel.spacing=unit(2, "lines"), panel.border=element_rect(color="black", fill=NA, linewidth=0.5)) + guides(fill="none")
    ggsave(file.path(input_dir, filename), plot=p, width=12, height=8, dpi=300)
  }

  plot_gam_angles(gam_angle_plot_data_all, paste(scenario_name, "- Arrhythmia Risk by Fiber Orientation (All)"), "statistical_model_fit_by_angle_all.png")
  plot_gam_angles(gam_angle_plot_data_real, paste(scenario_name, "- Arrhythmia Risk by Fiber Orientation (Realistic)"), "statistical_model_fit_by_angle_realistic.png")

  event_counts_angle <- data_angles %>% group_by(fibrosis_type, angle) %>% summarise(Total_Successes = sum(successes), .groups='drop')
  p_bar_angle <- ggplot(event_counts_angle, aes(x=fibrosis_type, y=Total_Successes, fill=angle)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7, color="black") +
    scale_fill_manual(values=ANGLE_COLORS, name="Fiber Angle (°)") +
    labs(title=paste(scenario_name, "- Arrhythmia Incidence by Pattern and Fiber Orientation"), x="", y="Total Number of Sustained Reentries") +
    theme_minimal(base_family="sans", base_size=16) +
    theme(axis.text.x=element_text(size=14, color="black", face="bold"), axis.text.y=element_text(size=14, color="black", face="bold"),
          axis.title=element_text(size=16), legend.position="top", legend.title=element_text(size=14, face="bold"),
          legend.text=element_text(size=14), panel.grid.major.x=element_blank())
  ggsave(file.path(input_dir, "overall_incidence_barplot_by_angle.png"), plot=p_bar_angle, width=10, height=6, dpi=300)
}

cat("Statiscal Analysis Completed Successfully.\n")
