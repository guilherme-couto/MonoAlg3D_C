import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    'figure.max_open_warning': 0,
    'axes.spines.top': False,
    'axes.spines.right': False
})

FIBROSIS_PALETTE = {
    'compact': '#1f77b4', 'diffuse': '#ff7f0e',
    'interstitial': '#2ca02c', 'patchy': '#d62728'
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analysis V7: Automated Vulnerable Window Detection")
    parser.add_argument('--root_dir', type=str, default='.', help="Root dir")
    parser.add_argument('--dim', type=str, required=True, choices=['2D', '3D'])
    parser.add_argument('--geom', type=str, required=True, choices=['full', 'ellipse'])
    parser.add_argument('--threshold', type=float, default=1000.0)
    parser.add_argument('--force_min', type=float, default=None, help="Manually set window start")
    parser.add_argument('--force_max', type=float, default=None, help="Manually set window end")

    return parser.parse_args()

def load_data(file_path):
    if not os.path.exists(file_path):
        print(f"Error: Not found {file_path}")
        sys.exit(1)
    return pd.read_csv(file_path)

def detect_window(df, force_min=None, force_max=None):
    """
    Detecta automaticamente onde ocorrem as arritmias.
    Retorna (min_density, max_density).
    """
    # Agrupa por densidade para ver se houve ALGUM sucesso em QUALQUER tipo
    density_activity = df.groupby('density')['is_sustained'].sum()

    # Filtra densidades onde houve pelo menos 1 arritmia
    active_densities = density_activity[density_activity > 0].index.tolist()

    if not active_densities:
        print("WARNING: No arrhythmias detected in the entire dataset!")
        return 0.0, 1.0 # Fallback

    # Define limites automáticos
    auto_min = min(active_densities)
    auto_max = max(active_densities)

    # Aplica override se o usuário pediu
    final_min = force_min if force_min is not None else auto_min
    final_max = force_max if force_max is not None else auto_max

    return final_min, final_max

def process_data(df, threshold, min_w, max_w):
    df['is_sustained'] = (df['final_time_ms'] >= threshold).astype(int)

    # Full Data (Para curvas e AUC global)
    df_viz_full = df.groupby(['fibrosis_type', 'density'])['is_sustained'].mean().reset_index()
    df_viz_full.rename(columns={'is_sustained': 'prob'}, inplace=True)

    # Window Data (Para estatística focada)
    df_window = df[(df['density'] >= min_w) & (df['density'] <= max_w)].copy()

    # Dados crus da janela (para pointplot)
    df_raw_window = df_window

    # Dados agregados da janela (para R)
    df_r = df_window.groupby(['fibrosis_type', 'density'])['is_sustained'].agg(['sum', 'count']).reset_index()
    df_r.rename(columns={'sum': 'successes', 'count': 'trials'}, inplace=True)
    df_r['failures'] = df_r['trials'] - df_r['successes']

    return df_viz_full, df_raw_window, df_r

def plot_density_curves_auc(df_viz_full, output_dir, min_w, max_w):
    print("Generating Density Curves (Full Range) & AUC...")
    plt.figure(figsize=(10, 6))

    auc_list = []
    y_max = df_viz_full['prob'].max()
    ylim = y_max * 1.15 if y_max > 0 else 0.1

    # Sombreado dinâmico baseado na detecção
    plt.axvspan(min_w, max_w, color='gray', alpha=0.1, label=f'Detected Window [{min_w}, {max_w}]')

    for f_type, color in FIBROSIS_PALETTE.items():
        subset = df_viz_full[df_viz_full['fibrosis_type'] == f_type].sort_values('density')
        if subset.empty: continue

        x, y = subset['density'].values, subset['prob'].values

        # PCHIP Interpolation
        if len(x) > 2:
            x_new = np.linspace(x.min(), x.max(), 300)
            y_new = np.clip(PchipInterpolator(x, y)(x_new), 0, 1)
        else:
            x_new, y_new = x, y

        plt.plot(x_new, y_new, color=color, label=f_type, linewidth=3)

        try:
            auc = np.trapezoid(y, x)
        except AttributeError:
            auc = np.trapz(y, x)
        auc_list.append({'fibrosis_type': f_type, 'auc': auc})

    plt.title("Arrhythmia Vulnerability Profile")
    plt.xlabel("Fibrosis Density")
    plt.ylabel("Probability")
    plt.ylim(0, ylim)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "curves_density_full_range.png"), dpi=300)
    plt.close()

    pd.DataFrame(auc_list).sort_values('auc', ascending=False).to_csv(
        os.path.join(output_dir, "auc_results_full.txt"), index=False, sep='\t', float_format='%.4f'
    )

def plot_point_comparison(df_raw_window, output_dir, min_w, max_w):
    print(f"Generating Point Plot (Window {min_w}-{max_w})...")
    plt.figure(figsize=(8, 6))

    sns.pointplot(
        data=df_raw_window,
        x="fibrosis_type",
        y="is_sustained",
        palette=FIBROSIS_PALETTE,
        capsize=.1,
        err_kws={'linewidth': 2},
        linestyle="none",
        marker="o",
        markersize=10
    )

    plt.title(f"Arrhythmia Risk in Detected Window\n(Density {min_w} - {max_w})")
    plt.ylabel("Probability (Mean ± 95% CI)")
    plt.xlabel("Fibrosis Type")
    plt.ylim(0, None)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "comparison_pointplot_window.png"), dpi=300)
    plt.close()

def main():
    args = parse_arguments()
    base = os.path.join(args.root_dir, args.dim, args.geom)
    csv_path = os.path.join(base, "analysis", "simulation_results.csv")
    out_dir = os.path.join(base, "post_processing")
    os.makedirs(out_dir, exist_ok=True)

    # 1. Load
    df_raw = load_data(csv_path)
    df_raw['is_sustained'] = (df_raw['final_time_ms'] >= args.threshold).astype(int) # Ensure flag exists for detection

    # 2. DETECT WINDOW
    min_w, max_w = detect_window(df_raw, args.force_min, args.force_max)

    # Log da janela detectada (importante para reprodutibilidade)
    with open(os.path.join(out_dir, "detected_window_info.txt"), "w") as f:
        f.write(f"Detected Min Density: {min_w}\n")
        f.write(f"Detected Max Density: {max_w}\n")
        if args.force_min or args.force_max:
            f.write("(Manual override was active)\n")

    print(f"\n--- WINDOW DETECTION ---")
    print(f"Analyzing range: [{min_w}, {max_w}]")
    print(f"------------------------\n")

    # 3. Process
    df_viz_full, df_raw_window, df_r = process_data(df_raw, args.threshold, min_w, max_w)

    # 4. Plots
    plot_point_comparison(df_raw_window, out_dir, min_w, max_w)
    plot_density_curves_auc(df_viz_full, out_dir, min_w, max_w)

    # 5. Export R
    df_r.to_csv(os.path.join(out_dir, "data_for_gam_r_window.csv"), index=False)

    print(f"Analysis Complete. Results in: {out_dir}")

if __name__ == "__main__":
    main()
