import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def analyze_comparison(modern_csv, legacy_txt):
    df_modern = pd.read_csv(modern_csv)

    legacy_cols = ['count', 't_final', 'e_final', 'z_off', 'x_final', 'y_final', 'z_final']
    df_legacy = pd.read_csv(legacy_txt, header=None, names=legacy_cols)

    metrics = ['t_final', 'e_final', 'z_final']
    
    print(f"{'Metric':<12} | {'KS Stat':<10} | {'P-Value':<10} | {'Status'}")
    print("-" * 50)

    for col in metrics:
        # Kolmogorov-Smirnov
        ks_stat, p_value = stats.ks_2samp(df_modern[col], df_legacy[col])
        
        status = "PASS" if p_value > 0.05 else "FAIL"
        print(f"{col:<12} | {ks_stat:<10.4f} | {p_value:<10.4f} | {status}")

        plt.figure(figsize=(10, 5))
        plt.hist(df_modern[col], bins=50, alpha=0.5, label='Modern', density=True)
        plt.hist(df_legacy[col], bins=50, alpha=0.5, label='Legacy', density=True)
        plt.title(f"Comparison of {col}")
        plt.xlabel("Value")
        plt.ylabel("Probability Density")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.show()

if __name__ == "__main__":
    analyze_comparison('./tests/files/modern_dist_rank0.csv', './tests/files/legacy_dist.txt')