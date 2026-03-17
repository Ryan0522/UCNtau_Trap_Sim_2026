import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def analyze_comparison(modern_csv, legacy_txt):
    df_modern = pd.read_csv(modern_csv)

    legacy_cols = ['count', 't_final', 'e_final', 'code', 'x_final', 'y_final', 'z_final']
    df_legacy = pd.read_csv(legacy_txt, header=None, names=legacy_cols)

    metrics = ['t_final', 'e_final', 'x_final', 'y_final', 'z_final']    

    print(f"{'Metric':<12} | {'KS Stat':<10} | {'P-Value':<10} | {'Status'}")
    print("-" * 65)

    for col in metrics:
        # Kolmogorov-Smirnov
        m_series = pd.to_numeric(df_modern[col], errors='coerce')
        l_series = pd.to_numeric(df_legacy[col], errors='coerce')
        
        m_data = m_series.dropna()
        l_data = l_series.dropna()
        
        if len(m_data) == 0 or len(l_data) == 0:
            print(f"{col:<12} | {'N/A':<10} | {'N/A':<10} | SKIP (Empty Data)")
            continue
        
        ks_stat, p_value = stats.ks_2samp(m_data, l_data)
        
        status = "PASS" if p_value > 0.05 else "FAIL"
        print(f"{col:<12} | {ks_stat:<10.4f} | {p_value:<10.4f} | {status}")

        plt.figure(figsize=(10, 5))
        plt.hist(m_data, bins=50, alpha=0.5, label='Modern', density=True)
        plt.hist(l_data, bins=50, alpha=0.5, label='Legacy', density=True)
        
        title_suffix = " (Flags Only)" if col == 'code' else ""
        plt.title(f"Comparison of {col}{title_suffix}")
        plt.xlabel("Value")
        plt.ylabel("Probability Density")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"./results/test/{col}_dist.png")
        plt.show()

if __name__ == "__main__":
    analyze_comparison('./tests/files/modern_dist_rank0.csv', './tests/files/legacy_dist.txt')