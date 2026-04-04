import pandas as pd
import glob
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-v0_8-paper")
FIG_WIDTH = 468.0 / 72.27
FIG_HEIGHT = FIG_WIDTH * 0.75
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.titlesize": 14,
    "figure.autolayout": True,
    "savefig.bbox": "tight",
})

TAU_N = 877.75

def run_analysis():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', type=str, default='test', help="Results sub-folder")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_results_path = os.path.join(os.path.dirname(script_dir), "results", args.out_dir)
    
    hold_times = [20, 50, 100, 200, 1550]
    summary_data = []
    all_hits_list = []

    print(f"--- Starting Analysis for {args.out_dir} ---")

    for ht in hold_times:
        search_pattern = os.path.join(base_results_path, f"HT_{ht}_*", "batch_*_rank*.csv")
        files = glob.glob(search_pattern)
        
        if not files:
            print(f" [!] No files for HT {ht}s")
            continue

        cols_to_read = ['code', 't_final', 'x_final', 'z_final', 'e_start', 'e_final'] 
        li = []
        for f in files:
            try:
                li.append(pd.read_csv(f, usecols=cols_to_read))
            except: continue

        if not li: continue
        ht_df = pd.concat(li, axis=0, ignore_index=True)
        
        total = len(ht_df)
        detected_df = ht_df[ht_df['code'] == 0].copy()
        detected_count = len(detected_df)
        survival_rate = detected_count / total if total > 0 else 0
        survival_rate_sim = survival_rate * np.exp(-ht / TAU_N)  # Adjusted for neutron decay
        
        summary_data.append({
            'hold_time': ht,
            'survival_rate': survival_rate,
            'survival_rate_beta': survival_rate_sim,
            'error': np.sqrt(detected_count)/total if total > 0 else 0
        })
        
        detected_df['ht'] = ht
        all_hits_list.append(detected_df)
        print(f" [+] HT {ht}s: Survival = {survival_rate:.4f} ({detected_count}/{total})")

    if not summary_data:
        print(" [Error] No data processed.")
        return

    hits_df = pd.concat(all_hits_list, ignore_index=True) if all_hits_list else None
    summary_df = pd.DataFrame(summary_data)

    # 1: Decay Curve
    plot_decay_curve(summary_df, base_results_path, args.out_dir)

    # 2: Spatial Hit Map (Overall & Comparison)
    plot_spatial_analysis(hits_df, hold_times, base_results_path)

    # 3: Arrival Time Spectrum
    plot_arrival_times(hits_df, hold_times, base_results_path)

    # 4; Energy Conservation Check
    plot_energy_conservation(hits_df, hold_times, base_results_path)

def plot_decay_curve(df, path, label):
    plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    plt.errorbar(df['hold_time'], df['survival_rate_beta'], yerr=df['error'],
                 fmt='o-', markersize=8, capsize=5, elinewidth=1, label=f"Defect: {label}")
    plt.yscale('log')
    plt.xlabel('Holding Time (s)', fontsize=12)
    plt.ylabel('Survival Rate', fontsize=12)
    plt.title('UCN Storage Decay Curve', fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.legend()
    plt.savefig(os.path.join(path, "decay_curve.png"), dpi=300)
    plt.close()

def plot_spatial_analysis(hits_df, hold_times, path):
    fig, axes = plt.subplots(1, 5, figsize=(FIG_WIDTH*5, FIG_HEIGHT), sharey=True)
    
    for i, ht in enumerate(hold_times):
        data = hits_df[hits_df['ht'] == ht]
        ax = axes[i]
        
        if not data.empty:
            hb = ax.hexbin(data['x_final'], data['z_final'], gridsize=40, cmap='inferno', mincnt=1)
            ax.set_title(f"HT: {ht}s (N={len(data)})")
        else:
            ax.set_title(f"HT: {ht}s (No Data)")
        
        ax.set_xlabel('X (m)')
        if i == 0: ax.set_ylabel('Z (m)')

    plt.tight_layout()
    plt.savefig(os.path.join(path, "hit_maps_comparison.png"), dpi=300)
    plt.close()

def plot_arrival_times(hits_df, hold_times, path):
    plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    
    for ht in hold_times:
        data = hits_df[hits_df['ht'] == ht]
        if data.empty: continue

        plt.hist(data['t_final'], bins=50, histtype='step', 
                 linewidth=1.5, label=f'HT {ht}s')

    plt.xlabel('Arrival Time at Dagger (s)')
    plt.ylabel('Neutron Counts')
    plt.title('UCN Arrival Time Spectrum (Deadtime Distribution)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(path, "arrival_times_spectrum.png"), dpi=300)
    plt.close()

def plot_energy_conservation(df, hold_times, path):
    df['e_drift'] = (df['e_final'] - df['e_start']) / df['e_start']
    fig, axes = plt.subplots(1, 5, figsize=(FIG_WIDTH * 5, FIG_HEIGHT), sharey=True)
    
    x_min, x_max = df['e_drift'].min(), df['e_drift'].max()
    for i, ht in enumerate(hold_times):
        data = df[df['ht'] == ht]
        ax = axes[i]
        
        if not data.empty:
            ax.hist(data['e_drift'], bins=50, range=(x_min, x_max), 
                    color='green', alpha=0.7)
            ax.axvline(0, color='red', linestyle='--')
            ax.set_title(f"HT: {ht}s (N={len(data)})")
            ax.set_yscale('log')
        else:
            ax.set_title(f"HT: {ht}s (No Data)")
        ax.set_xlabel('Rel. Energy Drift')
        if i == 0: 
            ax.set_ylabel('Neutron Counts (Log)')

    plt.suptitle('Energy Conservation vs. Hold Time', fontsize=16, y=1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(path, "energy_conservation_comparison.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    run_analysis()