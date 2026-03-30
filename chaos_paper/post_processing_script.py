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
    'axes.titlesize': 18,
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
    parser.add_argument('--pad_left', type=float, default=0.05, help="Extra density margin before the active window")
    parser.add_argument('--pad_right', type=float, default=0.05, help="Extra density margin after the active window")
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

        # 1. Plot the interpolated smooth curve
        ax.plot(x_new, y_new, color=color, label=f_type, linewidth=3)

        # 2. Plot the shaded area beneath the curve
        ax.fill_between(x_new, 0, y_new, color=color, alpha=0.1)

        # 3. Plot the actual observed data points (discrete markers)
        ax.plot(x, y, marker='o', linestyle='none', color=color, markersize=5, alpha=0.9)

    ax.set_xlim(min_w - 0.02, max_w + 0.02)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    apply_bold_ticks(ax)

def add_raw_lines_curves(ax, df_viz, min_w, max_w):
    """Helper to plot raw data points connected by straight lines."""
    for f_type in CATEGORY_ORDER:
        if f_type not in df_viz['fibrosis_type'].values:
            continue

        subset = df_viz[df_viz['fibrosis_type'] == f_type].sort_values('density')
        if subset.empty:
            continue

        x = subset['density'].values
        y = subset['prob'].values

        color = FIBROSIS_PALETTE[f_type]

        ax.plot(x, y, marker='o', linestyle='-', color=color, label=f_type,
                linewidth=2.5, markersize=6, alpha=0.9)

    ax.set_xlim(min_w - 0.02, max_w + 0.02)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    apply_bold_ticks(ax)

def plot_vulnerability_windows(df_raw, output_dir):
    """Generates a Gantt-style chart showing absolute intervals of arrhythmias with counts."""
    plt.figure(figsize=(10, 5))
    ax = plt.gca()

    # Find the minimum density step in the simulations to define block width
    densities = np.sort(df_raw['density'].unique())
    step = np.min(np.diff(densities)) if len(densities) > 1 else 0.05

    y_ticks = []
    y_labels = []

    # Variables to track the absolute boundaries of the active windows for zooming
    min_window_edge = float('inf')
    max_window_edge = float('-inf')

    # Reverse to plot Stochastic at the top (or matching other plot conventions)
    for i, f_type in enumerate(reversed(CATEGORY_ORDER)):
        y_pos = i * 10
        y_ticks.append(y_pos)
        y_labels.append(f_type)

        color = FIBROSIS_PALETTE.get(f_type, '#333333')

        # Identify all densities where at least one sustained reentry occurred and their counts
        # This creates a pandas Series: density -> reentry count
        active_counts = df_raw[(df_raw['fibrosis_type'] == f_type) & (df_raw['is_sustained'] > 0)].groupby('density')['is_sustained'].sum()

        # Densities sorted
        active_d = active_counts.index.values

        if len(active_d) == 0:
            continue

        # Update global window boundaries for X-axis zoom
        min_window_edge = min(min_window_edge, active_d.min())
        max_window_edge = max(max_window_edge, active_d.max())

        # Logic for grouping adjacent density points and accumulating counts
        intervals = [] # List of tuples: (start_x, width)
        counts = []    # List of total counts per interval

        start_d = active_d[0]
        prev_d = active_d[0]
        current_count = active_counts[start_d] # Start accumulation

        margin = step/5 # Margin to ensure width will overlap density boundaries nicely

        for d in active_d[1:]:
            # If the current density is close to the previous step, it's continuous
            # Using isclose and a loose tolerance for step inconsistencies
            if np.isclose(d - prev_d, step) or (d - prev_d) <= step * 1.05:
                prev_d = d
                current_count += active_counts[d] # Accumulate reentry count
            else:
                # Break in continuity detected, save the interval data
                width = (prev_d - start_d) + margin*2
                intervals.append((start_d - margin, width))
                counts.append(current_count) # Save accumulated count

                # Reset for the new interval
                start_d = d
                prev_d = d
                current_count = active_counts[d]

        # Append the final interval data
        width = (prev_d - start_d) + margin*2
        intervals.append((start_d - margin, width))
        counts.append(current_count)

        # Plot the broken horizontal bars (wide, outlined)
        ax.broken_barh(intervals, (y_pos - 4, 8), facecolors=color, edgecolor='black', linewidth=1.5)

        # Plot text labels inside bars
        # We must add labels for each interval segment individually
        for interval, count in zip(intervals, counts):
            # Calculate exact center coordinates (X is midpoint, Y is y_pos)
            center_x = interval[0] + (interval[1] / 2)

            # Draw bold, white, centered text
            ax.text(center_x, y_pos, str(int(count)),
                    color='white', fontweight='bold', fontsize=14,
                    ha='center', va='center')

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Fibrosis Density")
    ax.set_title("Scenario III - Vulnerability Windows Overview")

    # Set X limits based ONLY on the active window edges +/- 0.05
    if min_window_edge != float('inf'):
        ax.set_xlim(min_window_edge - 0.06, max_window_edge + 0.06)
    else:
        ax.set_xlim(0, 1) # Fallback if no arrhythmias exist

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

    ax.yaxis.grid(False)                   # Remove horizontal grid lines
    ax.spines['left'].set_visible(False)   # Remove the Y-axis spine line
    ax.tick_params(axis='y', length=0)     # Remove the small ticks next to the labels

    # Apply standard bold texts
    apply_bold_ticks(ax)

    plt.tight_layout()
    plt.savefig(output_dir / "vulnerability_windows_gantt.png", dpi=300)
    plt.close()

def plot_combined_curves_and_gantt(df_raw, df_viz, min_w, max_w, output_dir):
    """Generates a combined figure: Top = Raw Curves, Bottom = Gantt Windows."""
    # Create a figure with 2 subplots sharing the X axis. Height ratio 2:1 gives more space to curves.
    fig, (ax_curves, ax_gantt) = plt.subplots(nrows=2, ncols=1, figsize=(12, 10),
                                              sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    # ==========================================
    # --- TOP SUBPLOT: Curves ---
    # ==========================================
    ax_curves.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)
    max_prob = df_viz['prob'].max()
    y_limit = min(max(max_prob * 1.15, 0.1), 1.0)

    # Use the existing function to plot the raw lines
    add_raw_lines_curves(ax_curves, df_viz, min_w, max_w)

    ax_curves.set_title("Scenario III - Empirical Arrhythmia Probability & Vulnerability Windows Overview")
    ax_curves.set_ylabel("Probability of Sustained Reentry")
    ax_curves.set_ylim(0, y_limit)
    ax_curves.legend(loc='upper right')
    # X-axis label is omitted here because it's shared with the bottom plot

    # ==========================================
    # --- BOTTOM SUBPLOT: Gantt ---
    # ==========================================
    densities = np.sort(df_raw['density'].unique())
    step = np.min(np.diff(densities)) if len(densities) > 1 else 0.05

    y_ticks = []
    y_labels = []

    for i, f_type in enumerate(reversed(CATEGORY_ORDER)):
        y_pos = i * 10
        y_ticks.append(y_pos)
        y_labels.append(f_type)

        color = FIBROSIS_PALETTE.get(f_type, '#333333')

        active_counts = df_raw[(df_raw['fibrosis_type'] == f_type) & (df_raw['is_sustained'] > 0)].groupby('density')['is_sustained'].sum()
        active_d = active_counts.index.values

        if len(active_d) == 0:
            continue

        intervals = []
        counts = []
        start_d = active_d[0]
        prev_d = active_d[0]
        current_count = active_counts[start_d]
        margin = step / 5

        for d in active_d[1:]:
            if np.isclose(d - prev_d, step) or (d - prev_d) <= step * 1.05:
                prev_d = d
                current_count += active_counts[d]
            else:
                width = (prev_d - start_d) + margin * 2
                intervals.append((start_d - margin, width))
                counts.append(current_count)

                start_d = d
                prev_d = d
                current_count = active_counts[d]

        width = (prev_d - start_d) + margin * 2
        intervals.append((start_d - margin, width))
        counts.append(current_count)

        ax_gantt.broken_barh(intervals, (y_pos - 4, 8), facecolors=color, edgecolor='black', linewidth=1.5)

        for interval, count in zip(intervals, counts):
            center_x = interval[0] + (interval[1] / 2)
            ax_gantt.text(center_x, y_pos, str(int(count)), color='white', fontweight='bold', fontsize=14, ha='center', va='center')

    ax_gantt.set_yticks(y_ticks)
    ax_gantt.set_yticklabels(y_labels)
    ax_gantt.set_xlabel("Fibrosis Density")

    # Styling Gantt subplot
    ax_gantt.yaxis.grid(False)
    ax_gantt.xaxis.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_gantt.spines['left'].set_visible(False)
    ax_gantt.tick_params(axis='y', length=0)

    # Ensure X-axis limits are exactly bounded by the active window (same as top plot)
    ax_gantt.set_xlim(min_w - 0.02, max_w + 0.02)
    ax_gantt.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

    ax_gantt.text(0.99, 0.99, "* Numbers inside bars indicate\nthe sustained reentry count per interval",
                  transform=ax_gantt.transAxes, fontsize=11, style='italic',
                  color='gray', ha='right', va='center')

    apply_bold_ticks(ax_gantt)

    plt.tight_layout()
    plt.savefig(output_dir / "combined_curves_and_gantt.png", dpi=300)
    plt.close()


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
    min_w = max(df_raw['density'].min(), round(core_min - args.pad_left, 2))
    max_w = min(df_raw['density'].max(), round(core_max + args.pad_right, 2))

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
        df_viz_angle = None

    # ---------------------------------------------------------
    # PLOTTING
    # ---------------------------------------------------------
    print("Generating plots...")

    # --- Gantt Chart for Vulnerability Windows ---
    # This runs before df_window filtering, showing the absolute spectrum of the data
    plot_vulnerability_windows(df_raw, out_dir)

    # --- Combined Curves and Gantt ---
    plot_combined_curves_and_gantt(df_raw, df_viz, min_w, max_w, out_dir)

    # --- Overall Curves ---
    plt.figure(figsize=(10, 6))
    ax_main = plt.gca()
    ax_main.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)

    max_prob = df_viz['prob'].max()
    y_limit = min(max(max_prob * 1.15, 0.1), 1.0)

    # add_pchip_curves(ax_main, df_viz, min_w, max_w)
    add_raw_lines_curves(ax_main, df_viz, min_w, max_w)

    plt.title("Scenario III - Empirical Arrhythmia Probability")
    plt.xlabel("Fibrosis Density")
    plt.ylabel("Probability of Sustained Reentry")
    plt.ylim(0, y_limit)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "window_curves_overall.png", dpi=300)
    plt.close()

    # --- Raw Pointplot ---
    plt.figure(figsize=(8, 6))
    ax_pt = sns.pointplot(
        data=df_window, x="fibrosis_type", y="is_sustained", hue="fibrosis_type",
        order=CATEGORY_ORDER, palette=FIBROSIS_PALETTE, legend=False,
        capsize=.1, marker="o", linestyle="none"
    )
    plt.title("Scenario III - Arrhythmia Risk Intensity")
    plt.ylabel("Probability of Sustained Reentry")
    plt.xlabel("Fibrosis Pattern")
    ax_pt.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    apply_bold_ticks(ax_pt)
    plt.tight_layout()
    plt.savefig(out_dir / "window_pointplot.png", dpi=300)
    plt.close()

    # --- Angle-Specific Curves (2x2 Grid) ---
    if df_angle_r is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
        axes = axes.flatten()

        # Find global max prob across all angles to synchronize Y-axis beautifully
        max_prob_angles = df_viz_angle[df_viz_angle[ANGLE_COL].isin(TARGET_ANGLES)]['prob'].max()
        global_y_limit = min(max(max_prob_angles * 1.15, 0.1), 1.0) if not pd.isna(max_prob_angles) else 1.0

        for i, angle in enumerate(TARGET_ANGLES):
            ax = axes[i]
            ax.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)
            df_subset = df_viz_angle[df_viz_angle[ANGLE_COL] == angle]

            if not df_subset.empty:
                # add_pchip_curves(ax, df_subset, min_w, max_w)
                add_raw_lines_curves(ax, df_subset, min_w, max_w)

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
