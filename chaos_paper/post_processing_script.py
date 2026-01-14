import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Binomial

# --- Configuration & Style ---
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    'figure.max_open_warning': 0,
    'axes.spines.top': False,
    'axes.spines.right': False
})

FIBROSIS_PALETTE = {
    'compact': '#1f77b4',      # Blue
    'diffuse': '#ff7f0e',      # Orange
    'interstitial': '#2ca02c', # Green
    'patchy': '#d62728'        # Red
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Statistical inference and analysis of cardiac simulations.")
    parser.add_argument('--root_dir', type=str, default='.', help="Root directory containing simulation folders.")
    parser.add_argument('--dim', type=str, required=True, choices=['2D', '3D'], help="Simulation dimension.")
    parser.add_argument('--geom', type=str, required=True, choices=['full', 'ellipse'], help="Geometry scenario.")
    parser.add_argument('--threshold', type=float, default=1000.0, help="Sustained arrhythmia time threshold (ms).")
    return parser.parse_args()

def load_data(file_path):
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        sys.exit(1)
    print(f"Loading data from: {file_path}")
    return pd.read_csv(file_path)

def calculate_basic_stats(df, threshold):
    """
    Calculates Probability, SE, and CI for every configuration.
    """
    df['is_sustained'] = (df['final_time_ms'] >= threshold).astype(int)

    stats = df.groupby(['fibrosis_type', 'fiber_angle_deg', 'density'])['is_sustained'].agg(['mean', 'count', 'std']).reset_index()
    stats.rename(columns={'mean': 'prob', 'count': 'n_samples'}, inplace=True)

    # Standard Error and CI
    stats['se'] = np.sqrt((stats['prob'] * (1 - stats['prob'])) / stats['n_samples'])
    stats['ci_95'] = 1.96 * stats['se']
    stats['ci_lower'] = (stats['prob'] - stats['ci_95']).clip(lower=0.0)
    stats['ci_upper'] = (stats['prob'] + stats['ci_95']).clip(upper=1.0)

    return stats

def run_logistic_regression(df_raw):
    """
    Performs a Logistic Regression to determine the Odds Ratio of each fibrosis type.
    Formula: is_sustained ~ density + C(fibrosis_type) + C(fiber_angle_deg)
    """
    print("Running Logistic Regression Model...")

    # We set one fibrosis type as reference (baseline). Let's pick 'compact' or the first one alphabetically.
    # The model handles this, but let's be explicit if needed.

    formula = 'is_sustained ~ density + C(fibrosis_type) + C(fiber_angle_deg)'

    try:
        model = smf.logit(formula=formula, data=df_raw).fit(disp=0)
        return model
    except Exception as e:
        print(f"Warning: Logistic regression failed (possibly due to separation or lack of variance). {e}")
        return None

def calculate_auc(stats_df):
    """
    Calculates Area Under the Curve (AUC) for Probability vs Density
    using the trapezoidal rule. This represents 'Total Arrhythmic Burden'.
    """
    auc_results = []

    for f_type in stats_df['fibrosis_type'].unique():
        subset = stats_df[stats_df['fibrosis_type'] == f_type]
        # Calculate AUC for each angle, then average, or average probs then AUC?
        # Let's do AUC per angle to get variance.

        for angle in subset['fiber_angle_deg'].unique():
            angle_data = subset[subset['fiber_angle_deg'] == angle].sort_values('density')
            auc = np.trapezoid(angle_data['prob'], angle_data['density'])
            auc_results.append({
                'fibrosis_type': f_type,
                'fiber_angle_deg': angle,
                'auc': auc
            })

    return pd.DataFrame(auc_results)

def interpolate_smooth(x, y, num_points=200):
    """PCHIP Interpolation for smooth curves without shading."""
    if len(x) < 2: return x, y
    x_new = np.linspace(x.min(), x.max(), num_points)
    f_prob = PchipInterpolator(x, y)
    y_smooth = np.clip(f_prob(x_new), 0, 1)
    return x_new, y_smooth

def plot_trend_overview(stats_df, output_dir):
    print("Generating clean trend overview...")
    angles = sorted(stats_df['fiber_angle_deg'].unique())
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    axes = axes.flatten()

    # Dynamic Y limit
    y_max = stats_df['prob'].max()
    y_limit = y_max * 1.15 if y_max > 0 else 0.1

    for i, angle in enumerate(angles):
        ax = axes[i]
        subset = stats_df[stats_df['fiber_angle_deg'] == angle]

        for f_type in FIBROSIS_PALETTE.keys():
            data = subset[subset['fibrosis_type'] == f_type].sort_values('density')
            if data.empty: continue

            x_s, y_s = interpolate_smooth(data['density'].values, data['prob'].values)
            ax.plot(x_s, y_s, color=FIBROSIS_PALETTE[f_type], label=f_type, linewidth=2.5)

        ax.set_title(f"Angle: {angle}°", fontweight='bold')
        ax.set_ylim(0, y_limit)
        if i >= 2: ax.set_xlabel("Fibrosis Density")
        if i % 2 == 0: ax.set_ylabel("Prob. Sustained Reentry")

        if i == 0: ax.legend(frameon=False, loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "trend_overview.png"), dpi=300)
    plt.close()

def plot_aggregated_comparison(stats_df, output_dir):
    print("Generating clean aggregated comparison...")

    agg_stats = stats_df.groupby(['fibrosis_type', 'density'])['prob'].mean().reset_index()

    plt.figure(figsize=(8, 6))
    y_max = agg_stats['prob'].max()
    y_limit = y_max * 1.15 if y_max > 0 else 0.1

    for f_type in FIBROSIS_PALETTE.keys():
        data = agg_stats[agg_stats['fibrosis_type'] == f_type].sort_values('density')
        x_s, y_s = interpolate_smooth(data['density'].values, data['prob'].values)
        plt.plot(x_s, y_s, color=FIBROSIS_PALETTE[f_type], label=f_type, linewidth=3)

    plt.title("Aggregated Arrhythmogenicity (All Angles)")
    plt.xlabel("Fibrosis Density")
    plt.ylabel("Mean Prob. Sustained Reentry")
    plt.ylim(0, y_limit)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "aggregated_comparison.png"), dpi=300)
    plt.close()

def plot_heatmaps(stats_df, output_dir):
    print("Generating heatmaps...")
    global_max = stats_df['prob'].max()
    vmax = np.ceil(global_max * 10) / 10.0
    if vmax == 0: vmax = 0.1

    # Combined plot only (as requested implicitly by focus on comparisons)
    fibrosis_types = sorted(list(FIBROSIS_PALETTE.keys()))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    cbar_ax = fig.add_axes([.91, .3, .02, .4])

    for i, f_type in enumerate(fibrosis_types):
        subset = stats_df[stats_df['fibrosis_type'] == f_type]
        pivot = subset.pivot(index="fiber_angle_deg", columns="density", values="prob")

        sns.heatmap(pivot, ax=axes[i], annot=False, cmap="rocket_r",
                    vmin=0, vmax=vmax, cbar=(i==0), cbar_ax=None if i else cbar_ax)
        axes[i].set_title(f_type.capitalize(), fontweight='bold')
        axes[i].set_ylabel("Angle (°)")
        axes[i].invert_yaxis()

    cbar_ax.set_ylabel('Prob. Sustained Reentry', rotation=270, labelpad=20)
    fig.suptitle("Comparative Vulnerability Heatmaps", fontsize=16)
    plt.subplots_adjust(right=0.89)
    plt.savefig(os.path.join(output_dir, "heatmap_combined.png"), dpi=300)
    plt.close()

def generate_full_report(stats_df, df_raw, logistic_model, auc_df, output_dir):
    print("Writing full statistical report...")

    # 1. Save detailed CSV with SE and CI for EVERY configuration
    csv_path = os.path.join(output_dir, "full_statistical_data.csv")
    stats_df.to_csv(csv_path, index=False, float_format='%.5f')

    # 2. Generate Text Report
    report_path = os.path.join(output_dir, "analysis_report.txt")
    with open(report_path, "w") as f:
        f.write("=== ADVANCED STATISTICAL ANALYSIS REPORT ===\n\n")
        f.write("PART 1: LOGISTIC REGRESSION INFERENCE\n")
        f.write("-------------------------------------\n")
        f.write("Model: Sustained_Reentry ~ Density + Fibrosis_Type + Fiber_Angle\n")
        f.write("Interpretation: Odds Ratio (OR) > 1 indicates increased risk relative to reference.\n")
        f.write("P-value < 0.05 indicates statistical significance.\n\n")

        if logistic_model:
            # Extract Odds Ratios and P-values
            params = logistic_model.params
            conf = logistic_model.conf_int()
            conf['OR'] = params
            conf.columns = ['2.5%', '97.5%', 'OR']
            # Exponentiate to get OR
            conf = np.exp(conf)
            pvalues = logistic_model.pvalues

            results = pd.concat([conf, pvalues], axis=1)
            results.rename(columns={0: 'p-value'}, inplace=True)

            f.write(results.to_string())
            f.write("\n\n")
            f.write("Summary:\n")
            # Simple interpretation of the most significant finding
            sig_vars = results[results['p-value'] < 0.05].index.tolist()
            f.write(f"Significant predictors (p<0.05): {', '.join(sig_vars)}\n")
        else:
            f.write("Logistic regression could not converge.\n")

        f.write("\nPART 2: ARRHYTHMIC BURDEN (AUC ANALYSIS)\n")
        f.write("----------------------------------------\n")
        f.write("Total area under probability curve (measure of overall risk across all densities).\n\n")

        auc_summary = auc_df.groupby('fibrosis_type')['auc'].agg(['mean', 'std']).sort_values('mean', ascending=False)
        f.write(auc_summary.to_string())
        f.write("\n\n")

        f.write("PART 3: UNCERTAINTY OVERVIEW (Standard Errors)\n")
        f.write("----------------------------------------------\n")
        f.write("Summary of Standard Errors (SE) by Fibrosis Type:\n\n")
        se_summary = stats_df.groupby('fibrosis_type')['se'].describe()[['mean', 'min', 'max', 'std']]
        f.write(se_summary.to_string())
        f.write("\n\n")
        f.write("NOTE: For the full table of SE and CI for every single configuration,\n")
        f.write(f"please refer to the file: {os.path.basename(csv_path)}\n")

def main():
    args = parse_arguments()

    base_path = os.path.join(args.root_dir, args.dim, args.geom)
    input_csv = os.path.join(base_path, "analysis", "simulation_results.csv")
    output_dir = os.path.join(base_path, "post_processing_v3")
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load Data
    df_raw = load_data(input_csv)

    # 2. Basic Stats (Prob, SE, CI)
    df_stats = calculate_basic_stats(df_raw, args.threshold)

    # 3. Advanced Stats (Logistic Regression & AUC)
    logistic_model = run_logistic_regression(df_raw)
    auc_df = calculate_auc(df_stats)

    # 4. Visualization (Clean curves, no shading)
    plot_trend_overview(df_stats, output_dir)
    plot_aggregated_comparison(df_stats, output_dir)
    plot_heatmaps(df_stats, output_dir)

    # 5. Reporting
    generate_full_report(df_stats, df_raw, logistic_model, auc_df, output_dir)

    print(f"\nAnalysis completed. Results saved to: {output_dir}")

if __name__ == "__main__":
    main()
