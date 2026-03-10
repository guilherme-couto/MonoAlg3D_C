import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy.interpolate import PchipInterpolator

# --- SETTINGS ---
sns.set_theme(style="whitegrid", context="talk") # may be "talk" or "notebook" depending on preference
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'figure.max_open_warning': 0,
    'axes.spines.top': False,
    'axes.spines.right': False
})

FIBROSIS_PALETTE = {
    'compact': '#0000a2', 'diffuse': '#50ad9f',
    'interstitial': '#e9c716', 'patchy': '#bc272d'
}

CATEGORY_ORDER = ['compact', 'diffuse', 'interstitial', 'patchy']

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analysis V8: Window-Focused Metrics")
    parser.add_argument('--root_dir', type=str, default='.', help="Root dir")
    parser.add_argument('--dim', type=str, required=True, choices=['2D', '3D'])
    parser.add_argument('--geom', type=str, required=True, choices=['full', 'ellipse'])
    parser.add_argument('--threshold', type=float, default=1000.0)
    parser.add_argument('--force_min', type=float, default=None)
    parser.add_argument('--force_max', type=float, default=None)

    return parser.parse_args()

def load_data(file_path):
    if not os.path.exists(file_path):
        print(f"Error: Not found {file_path}")
        sys.exit(1)
    return pd.read_csv(file_path)

def detect_window(df, force_min=None, force_max=None):
    """
    Detecta limites baseados na presença de arritmia (>0).
    """
    density_activity = df.groupby('density')['is_sustained'].sum()
    active_densities = density_activity[density_activity > 0].index.tolist()

    if not active_densities:
        print("WARNING: No arrhythmias detected!")
        return 0.1, 0.9 # Default safe range

    auto_min = min(active_densities)
    auto_max = max(active_densities)

    # Pequena margem de segurança visual (opcional, para pegar o zero logo antes/depois)
    # Mas para cálculo rigoroso, manteremos estrito ou estenderemos 0.05 se disponível.
    # Vamos manter estrito ao que tem atividade para ser consistente.

    final_min = force_min if force_min is not None else auto_min
    final_max = force_max if force_max is not None else auto_max

    return final_min, final_max

def process_window_data(df, threshold, min_w, max_w):
    df['is_sustained'] = (df['final_time_ms'] >= threshold).astype(int)

    # Filtra APENAS a janela
    df_window = df[(df['density'] >= min_w) & (df['density'] <= max_w)].copy()

    # 1. Dados para Curvas (Médias por ponto na janela)
    df_viz = df_window.groupby(['fibrosis_type', 'density'])['is_sustained'].mean().reset_index()
    df_viz.rename(columns={'is_sustained': 'prob'}, inplace=True)

    # 2. Dados Crús (Para Pointplot)
    df_raw = df_window

    # 3. Dados Agregados (Para R)
    df_r = df_window.groupby(['fibrosis_type', 'density'])['is_sustained'].agg(['sum', 'count']).reset_index()
    df_r.rename(columns={'sum': 'successes', 'count': 'trials'}, inplace=True)
    df_r['failures'] = df_r['trials'] - df_r['successes']

    return df_viz, df_raw, df_r

def calculate_metrics(df_viz, min_w, max_w):
    """
    Calcula AUC e Intensidade Média dentro da janela.
    """
    results = []
    for f_type in FIBROSIS_PALETTE.keys():
        subset = df_viz[df_viz['fibrosis_type'] == f_type].sort_values('density')
        if subset.empty:
            results.append({'fibrosis_type': f_type, 'auc': 0, 'mean_intensity': 0})
            continue

        x, y = subset['density'].values, subset['prob'].values

        # 1. AUC (Total Burden in Window)
        try:
            auc = np.trapezoid(y, x)
        except AttributeError:
            auc = np.trapz(y, x)

        # 2. Mean Intensity (Average probability given we are in the window)
        # Aproximação simples: média aritmética das probabilidades nos pontos amostrados
        mean_intensity = np.mean(y)

        results.append({
            'fibrosis_type': f_type,
            'auc': auc,
            'mean_intensity': mean_intensity
        })

    return pd.DataFrame(results).sort_values('auc', ascending=False)

def plot_auc_bars(metrics_df, output_dir):
    """
    Gráfico de barras para a AUC.
    """
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(
        data=metrics_df,
        x='fibrosis_type',
        y='auc',
        hue='fibrosis_type',
        order=CATEGORY_ORDER,
        palette=FIBROSIS_PALETTE,
        legend=False
    )

    plt.title("Total Arrhythmic Burden (AUC)")
    plt.ylabel("Area Under the Curve", fontsize=16)
    plt.xlabel("Fibrosis Type", fontsize=16)

    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "auc_barplot.png"), dpi=300)
    plt.close()

def plot_focused_analysis(df_viz, df_raw, output_dir, min_w, max_w):

    # --- PLOT 1: Pointplot ---
    plt.figure(figsize=(8, 6))
    ax1 = sns.pointplot(
        data=df_raw,
        x="fibrosis_type",
        y="is_sustained",
        hue="fibrosis_type",
        order=CATEGORY_ORDER,
        palette=FIBROSIS_PALETTE,
        legend=False,
        capsize=.1,
        marker="o",
        linestyle="none"
    )
    plt.title(f"Risk Intensity (Density Window: {min_w:.2f} - {max_w:.2f})")
    plt.ylabel("Probability of Sustained Reentry", fontsize=16)
    plt.xlabel("Fibrosis Type", fontsize=16)
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax1.tick_params(axis='both', labelsize=14)
    ax1.set_xticklabels(ax1.get_xticklabels(), fontweight='bold')
    ax1.set_yticklabels(ax1.get_yticklabels(), fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "window_pointplot.png"), dpi=300)
    plt.close()

    # --- PLOT 2: Curves ---
    plt.figure(figsize=(10, 6))
    max_prob_observed = df_viz['prob'].max()
    y_limit = min(max(max_prob_observed * 1.15, 0.1), 1.0)

    ax2 = plt.gca()

    for f_type in CATEGORY_ORDER:
        if f_type not in df_viz['fibrosis_type'].values: continue
        color = FIBROSIS_PALETTE[f_type]
        subset = df_viz[df_viz['fibrosis_type'] == f_type].sort_values('density')
        if subset.empty: continue

        x, y = subset['density'].values, subset['prob'].values
        if len(x) > 2:
            x_new = np.linspace(x.min(), x.max(), 300)
            y_new = np.clip(PchipInterpolator(x, y)(x_new), 0, 1)
        else:
            x_new, y_new = x, y

        plt.plot(x_new, y_new, color=color, label=f_type, linewidth=3)
        plt.fill_between(x_new, 0, y_new, color=color, alpha=0.1)

    plt.title(f"Vulnerability Profile (Density Window: {min_w:.2f} - {max_w:.2f})")
    plt.xlabel("Fibrosis Density", fontsize=16)
    plt.ylabel("Probability of Sustained Reentry", fontsize=16)
    plt.xlim(min_w - 0.02, max_w + 0.02)
    plt.ylim(0, y_limit)

    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax2.tick_params(axis='both', labelsize=14)
    ax2.set_xticklabels(ax2.get_xticklabels(), fontweight='bold')
    ax2.set_yticklabels(ax2.get_yticklabels(), fontweight='bold')

    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "window_curves.png"), dpi=300)
    plt.close()

def main():
    args = parse_arguments()
    base = os.path.join(args.root_dir, args.dim, args.geom)
    csv_path = os.path.join(base, "analysis", "simulation_results.csv")
    out_dir = os.path.join(base, "post_processing")
    os.makedirs(out_dir, exist_ok=True)

    # 1. Load & Detect
    df_raw = load_data(csv_path)
    df_raw['is_sustained'] = (df_raw['final_time_ms'] >= args.threshold).astype(int)
    min_w, max_w = detect_window(df_raw, args.force_min, args.force_max)

    with open(os.path.join(out_dir, "window_info.txt"), "w") as f:
        f.write(f"Analyzed Window: {min_w} to {max_w}\n")

    print(f"Analyzing Window: [{min_w}, {max_w}]")

    # 2. Process Only Window Data
    df_viz, df_raw_window, df_r = process_window_data(df_raw, args.threshold, min_w, max_w)

    # 3. Calculate Metrics (AUC & Mean Intensity)
    metrics_df = calculate_metrics(df_viz, min_w, max_w)
    metrics_df.to_csv(os.path.join(out_dir, "window_metrics_summary.csv"), index=False, float_format="%.4f")
    plot_auc_bars(metrics_df, out_dir)

    # 4. Plots & Export
    plot_focused_analysis(df_viz, df_raw_window, out_dir, min_w, max_w)
    df_r.to_csv(os.path.join(out_dir, "data_for_gam_r.csv"), index=False)

    print(f"Done. Clean results in: {out_dir}")

if __name__ == "__main__":
    main()
