import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from pathlib import Path
from scipy.interpolate import PchipInterpolator

# ==========================================================
# --- GLOBAL SETTINGS & AESTHETICS ---
# ==========================================================
# Match R aesthetics: 14pt bold for ticks, 16pt normal for axis labels/titles
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.labelweight': 'normal',
    'axes.titlesize': 16,
    'axes.titleweight': 'normal',
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'figure.max_open_warning': 0,
    'axes.spines.top': False,
    'axes.spines.right': False
})

# Capitalized names to match the R scripts perfectly
FIBROSIS_PALETTE = {
    'Stochastic': '#444444',
    'Compact': '#0000a2',
    'Diffuse': '#50ad9f',
    'Interstitial': '#e9c716',
    'Patchy': '#bc272d'
}

CATEGORY_ORDER = ['Stochastic', 'Compact', 'Diffuse', 'Interstitial', 'Patchy']

# Define the exact name of the angle column in your CSV
ANGLE_COL = 'fiber_angle_deg'
TARGET_ANGLES = [0, 30, 60, 90]

# ==========================================================
# --- HELPER FUNCTIONS ---
# ==========================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description="Analysis: Window-Focused Metrics & Interpolation")
    parser.add_argument('--root_dir', type=str, default='.', help="Root directory")
    parser.add_argument('--dim', type=str, required=True, choices=['2D', '3D'])
    parser.add_argument('--geom', type=str, required=True, choices=['full', 'ellipse'])
    parser.add_argument('--threshold', type=float, default=1000.0, help="Sustained reentry threshold (ms)")
    parser.add_argument('--force_min', type=float, default=None, help="Force window min density")
    parser.add_argument('--force_max', type=float, default=None, help="Force window max density")
    parser.add_argument('--output_dir', type=str, default='post_processing')
    parser.add_argument('--pad_left', type=float, default=0.0, help="Extra density margin before the active window")
    parser.add_argument('--pad_right', type=float, default=0.0, help="Extra density margin after the active window")
    return parser.parse_args()

def apply_bold_ticks(ax):
    """Helper to apply bold weight to tick labels matching R plots."""
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')

def detect_window(df, force_min=None, force_max=None):
    """Detect boundaries based on the presence of arrhythmia (>0)."""
    density_activity = df.groupby('density')['is_sustained'].sum()
    active_densities = density_activity[density_activity > 0].index.tolist()

    if not active_densities:
        print("WARNING: No arrhythmias detected in the dataset!")
        return 0.05, 0.9 # Default safe fallback

    auto_min = min(active_densities)
    auto_max = max(active_densities)

    final_min = force_min if force_min is not None else auto_min
    final_max = force_max if force_max is not None else auto_max

    return final_min, final_max

def add_pchip_curves(ax, df_viz, min_w, max_w):
    """Helper to plot PCHIP interpolated curves on a given Matplotlib axis."""
    for f_type in CATEGORY_ORDER:
        if f_type not in df_viz['fibrosis_type'].values:
            continue

        subset = df_viz[df_viz['fibrosis_type'] == f_type].sort_values('density')
        if subset.empty:
            continue

        x, y = subset['density'].values, subset['prob'].values

        # PCHIP requires at least 3 points to interpolate smoothly
        if len(x) > 2:
            x_new = np.linspace(x.min(), x.max(), 300)
            y_new = np.clip(PchipInterpolator(x, y)(x_new), 0, 1)
        else:
            x_new, y_new = x, y

        color = FIBROSIS_PALETTE[f_type]
        ax.plot(x_new, y_new, color=color, label=f_type, linewidth=3)
        ax.fill_between(x_new, 0, y_new, color=color, alpha=0.1)

    ax.set_xlim(min_w - 0.02, max_w + 0.02)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    apply_bold_ticks(ax)

# ==========================================================
# --- MAIN PROCESSING & PLOTTING ---
# ==========================================================
def main():
    args = parse_arguments()

    # Pathlib usage for clean cross-platform paths
    base_path = Path(args.root_dir) / args.dim / args.geom
    csv_path = base_path / "analysis" / "simulation_results.csv"
    out_dir = base_path / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not csv_path.exists():
        print(f"Error: Dataset not found at {csv_path}")
        sys.exit(1)

    print(f"Loading data from {csv_path}...")
    df_raw = pd.read_csv(csv_path)

    # Standardize string formatting (Capitalized first letter)
    if df_raw['fibrosis_type'].dtype == object:
        df_raw['fibrosis_type'] = df_raw['fibrosis_type'].str.strip().str.capitalize()

    # Binary outcome classification
    df_raw['is_sustained'] = (df_raw['final_time_ms'] >= args.threshold).astype(int)

    # Detect core window (strictly where arrhythmias happen)
    core_min, core_max = detect_window(df_raw, args.force_min, args.force_max)

    # Expand window with padding
    min_w = max(df_raw['density'].min(), core_min - args.pad_left)
    max_w = min(df_raw['density'].max(), core_max + args.pad_right)

    print(f"Core Active Window: [{core_min:.2f}, {core_max:.2f}]")
    print(f"Padded Analysis Window: [{min_w:.2f}, {max_w:.2f}]")

    # Filter data within the expanded active window
    df_window = df_raw[(df_raw['density'] >= min_w) & (df_raw['density'] <= max_w)].copy()

    # ---------------------------------------------------------
    # DATA AGGREGATION
    # ---------------------------------------------------------
    # 1. Overall Aggregation for Visuals
    df_viz = df_window.groupby(['fibrosis_type', 'density'])['is_sustained'].mean().reset_index()
    df_viz.rename(columns={'is_sustained': 'prob'}, inplace=True)

    # 2. Overall Aggregation for R (Firth/GAM)
    df_r = df_window.groupby(['fibrosis_type', 'density'])['is_sustained'].agg(['count', 'sum']).reset_index()
    df_r.rename(columns={'count': 'trials', 'sum': 'successes'}, inplace=True)
    df_r['failures'] = df_r['trials'] - df_r['successes']

    # 3. Angle-Specific Aggregation for Visuals and R
    if ANGLE_COL in df_window.columns:
        # Aggregation for R
        df_angle_r = df_window.groupby(['fibrosis_type', 'density', ANGLE_COL])['is_sustained'].agg(['count', 'sum']).reset_index()
        df_angle_r.rename(columns={'count': 'trials', 'sum': 'successes'}, inplace=True)
        df_angle_r['failures'] = df_angle_r['trials'] - df_angle_r['successes']

        # Aggregation for Plotting
        df_viz_angle = df_window.groupby(['fibrosis_type', 'density', ANGLE_COL])['is_sustained'].mean().reset_index()
        df_viz_angle.rename(columns={'is_sustained': 'prob'}, inplace=True)
    else:
        print(f"WARNING: Column '{ANGLE_COL}' not found. Angle-specific analysis will be skipped.")
        df_angle_r = None

    # ---------------------------------------------------------
    # PLOTTING
    # ---------------------------------------------------------
    print("Generating plots...")

    # --- PLOT 1: Overall Curves ---
    plt.figure(figsize=(10, 6))
    ax_main = plt.gca()

    max_prob = df_viz['prob'].max()
    y_limit = min(max(max_prob * 1.15, 0.1), 1.0)

    add_pchip_curves(ax_main, df_viz, min_w, max_w)

    plt.title(f"Overall Vulnerability Profile (Window: {core_min:.2f} - {core_max:.2f})")
    plt.xlabel("Fibrosis Density")
    plt.ylabel("Probability of Sustained Reentry")
    plt.ylim(0, y_limit)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "window_curves_overall.png", dpi=300)
    plt.close()

    # --- PLOT 2: Raw Pointplot ---
    plt.figure(figsize=(8, 6))
    ax_pt = sns.pointplot(
        data=df_window, x="fibrosis_type", y="is_sustained", hue="fibrosis_type",
        order=CATEGORY_ORDER, palette=FIBROSIS_PALETTE, legend=False,
        capsize=.1, marker="o", linestyle="none"
    )
    plt.title(f"Risk Intensity (Window: {core_min:.2f} - {core_max:.2f})")
    plt.ylabel("Probability of Sustained Reentry")
    plt.xlabel("Fibrosis Pattern")
    ax_pt.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    apply_bold_ticks(ax_pt)
    plt.tight_layout()
    plt.savefig(out_dir / "window_pointplot.png", dpi=300)
    plt.close()

    # --- PLOT 3: Angle-Specific Curves (2x2 Grid) ---
    if df_angle_r is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
        axes = axes.flatten()

        # Find global max prob across all angles to synchronize Y-axis beautifully
        max_prob_angles = df_viz_angle[df_viz_angle[ANGLE_COL].isin(TARGET_ANGLES)]['prob'].max()
        global_y_limit = min(max(max_prob_angles * 1.15, 0.1), 1.0) if not pd.isna(max_prob_angles) else 1.0

        for i, angle in enumerate(TARGET_ANGLES):
            ax = axes[i]
            df_subset = df_viz_angle[df_viz_angle[ANGLE_COL] == angle]

            if not df_subset.empty:
                add_pchip_curves(ax, df_subset, min_w, max_w)

            ax.set_title(f"Fiber Orientation: {angle}°")
            ax.set_ylim(0, global_y_limit)

            # Subplot label cleanup
            if i >= 2: ax.set_xlabel("Fibrosis Density")
            if i % 2 == 0: ax.set_ylabel("Probability of Reentry")

        # Create a single legend for the entire 2x2 figure at the bottom
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc='lower center', ncol=len(CATEGORY_ORDER), bbox_to_anchor=(0.5, -0.05))

        plt.tight_layout()
        plt.savefig(out_dir / "window_curves_by_angle.png", dpi=300, bbox_inches='tight')
        plt.close()

    # ---------------------------------------------------------
    # EXPORTING CSVS
    # ---------------------------------------------------------
    df_r.to_csv(out_dir / "data_for_gam_r.csv", index=False)
    print(f"Exported overall statistical data: {out_dir / 'data_for_gam_r.csv'}")

    if df_angle_r is not None:
        df_angle_r.to_csv(out_dir / "data_angles_for_gam_r.csv", index=False)
        print(f"Exported angle-specific statistical data: {out_dir / 'data_angles_for_gam_r.csv'}")

    print(f"Post-processing completed successfully!")

if __name__ == "__main__":
    main()
