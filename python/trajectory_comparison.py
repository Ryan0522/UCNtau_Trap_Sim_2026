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

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    ax1.plot(df['t'], df['modern_x'], label='Modern Engine', alpha=0.8)
    ax1.plot(df['t'], df['legacy_x'], label='Legacy Engine', linestyle='--', color='orange', alpha=0.8)
    ax1.set_ylabel('X Position (m)')
    ax1.set_title('Trajectory Comparison: Modern vs Legacy')
    ax1.legend()
    ax1.grid(True)

    ax2.plot(df['t'], df['error'], color='red', label='Euclidean Distance Error')
    ax2.set_yscale('log')
    ax2.set_ylabel('Error (m) - Log Scale')
    ax2.set_xlabel('Time (s)')
    ax2.set_title('Error Divergence Analysis')
    ax2.grid(True, which="both", ls="-", alpha=0.5)
    ax2.legend()

    try:
        t_start = df[df['error'] > 1e-10]['t'].iloc[0]
        t_end = df[df['error'] > 1e-5]['t'].iloc[0]
        growth_rate = (np.log(1e-5) - np.log(1e-10)) / (t_end - t_start)
        lyapunov_time = 1 / growth_rate
        print(f"Estimated Lyapunov Time: {lyapunov_time:.2f} seconds")
    except:
        print("Error growth is too small or too large to estimate Lyapunov Time.")

    plt.tight_layout()
    if field_type == 'planar':
        plt.savefig('./results/planar_divergence_plot.png')
    else:
        plt.savefig('./results/trap_divergence_plot.png')
    plt.show()

if __name__ == "__main__":
    analyze_divergence('planar', './tests/planar_comparison.csv')
    analyze_divergence('trap', './tests/trap_comparison.csv')