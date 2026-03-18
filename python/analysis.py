import pandas as pd
import glob
import os
import argparse
import matplotlib.pyplot as plt

def run_analysis():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', type=str, default='test', help="The sub-folder in results/")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_results_path = os.path.join(os.path.dirname(script_dir), "results", args.out_dir)
    
    hold_times = [20, 50, 100, 200, 1550]
    summary_data = []

    print(f"Scanning directory: {base_results_path}")

    for ht in hold_times:
        search_pattern = os.path.join(base_results_path, f"HT_{ht}_*", "batch_*_rank*.csv")
        files = glob.glob(search_pattern)
        
        if not files:
            print(f"Warning: No files found for HT {ht}s at {search_pattern}")
            continue

        print(f"Processing HT {ht}s: Found {len(files)} files...")

        cols_to_read = ['code', 't_final', 'e_final'] 
        li = []
        for f in files:
            try:
                df = pd.read_csv(f, usecols=cols_to_read)
                li.append(df)
            except Exception as e:
                print(f"Error reading {f}: {e}")

        ht_df = pd.concat(li, axis=0, ignore_index=True)
        
        total = len(ht_df)
        detected = len(ht_df[ht_df['code'] == 0])
        survival_rate = detected / total if total > 0 else 0
        
        summary_data.append({
            'hold_time': ht,
            'total_neutrons': total,
            'detected_neutrons': detected,
            'survival_rate': survival_rate
        })
        
        output_file = os.path.join(base_results_path, f"integrated_HT_{ht}.csv")
        ht_df.to_csv(output_file, index=False)
        print(f"   Done. Survival Rate: {survival_rate:.4f}. Saved to {output_file}")

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(base_results_path, "final_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"\n[Success] Final summary saved to: {summary_path}")

        # 6. 自動繪圖
        plt.figure(figsize=(8, 6))
        plt.errorbar(summary_df['hold_time'], summary_df['survival_rate'], 
                     fmt='o-', capsize=5, label=f"Defect: {args.out_dir}")
        plt.yscale('log') # 衰減通常看對數座標
        plt.xlabel('Holding Time (s)')
        plt.ylabel('Survival Rate (Detected / Total)')
        plt.title(f'UCN Decay Curve - {args.out_dir}')
        plt.grid(True, which="both", ls="-", alpha=0.5)
        plt.legend()
        
        plot_path = os.path.join(base_results_path, "decay_curve.png")
        plt.savefig(plot_path)
        print(f"[Success] Plot saved to: {plot_path}")
    else:
        print("\n[Error] No data was processed. Check your paths and file names.")

if __name__ == "__main__":
    run_analysis()