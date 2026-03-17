import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-v0_8-paper")

def analyze_divergence(field_type, file_path):
    df = pd.read_csv(file_path)

    df['error'] = np.sqrt(
        (df['modern_x'] - df['legacy_x'])**2 +
        (df['modern_y'] - df['legacy_y'])**2 +
        (df['modern_z'] - df['legacy_z'])**2
    )

    fig, axes = plt.subplots(4, 1, figsize=(10, 15), sharex=True)

    coords = ['x', 'y', 'z']
    colors = ['tab:blue', 'tab:green', 'tab:orange']

    for i, coord in enumerate(coords):
        axes[i].plot(df['t'], df[f'modern_{coord}'], label=f'Modern {coord.upper()}', alpha=0.8, color=colors[i])
        axes[i].plot(df['t'], df[f'legacy_{coord}'], label=f'Legacy {coord.upper()}', 
                     linestyle='--', color='black', alpha=0.6)
        axes[i].set_ylabel(f'{coord.upper()} Position (m)')
        axes[i].legend(loc='upper right')
        axes[i].grid(True, alpha=0.3)

    axes[3].plot(df['t'], df['error'], color='red', label='Euclidean Error')
    axes[3].set_yscale('log')
    axes[3].set_ylabel('Error (m) - Log Scale')
    axes[3].set_xlabel('Time (s)')
    axes[3].set_title(f'[{field_type.upper()}] Error Divergence Analysis')
    axes[3].grid(True, which="both", ls="-", alpha=0.5)
    axes[3].axhline(y=1e-3, color='gray', linestyle=':', label='1mm Threshold')
    axes[3].legend(loc='upper left')

    try:
        valid_growth = df[(df['error'] > 1e-10) & (df['error'] < 1e-3)]
        if len(valid_growth) > 100:
            t_start = valid_growth['t'].iloc[0]
            t_end = valid_growth['t'].iloc[-1]
            err_start = valid_growth['error'].iloc[0]
            err_end = valid_growth['error'].iloc[-1]
            
            growth_rate = (np.log(err_end) - np.log(err_start)) / (t_end - t_start)
            lyapunov_time = 1 / growth_rate
            print(f"[{field_type}] Estimated Lyapunov Time: {lyapunov_time:.2f} seconds")
        else:
            print(f"[{field_type}] Divergence occurred too fast for Lyapunov estimation.")
    except Exception as e:
        print(f"[{field_type}] Lyapunov estimation failed: {e}")

    plt.tight_layout()

    save_path = f'./results/test/{field_type}_divergence_plot.png'
    plt.savefig(save_path, dpi=300)
    print(f"Figure saved to: {save_path}")
    plt.show()

if __name__ == "__main__":
    analyze_divergence('planar', './tests/files/planar_comparison.csv')
    analyze_divergence('trap', './tests/files/trap_comparison.csv')