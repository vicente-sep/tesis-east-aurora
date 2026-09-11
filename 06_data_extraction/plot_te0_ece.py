import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Paths
base_dir = "/Users/vicentesepulveda/Documents/Tesis_EAST/Extracted_EAST"
csv_dir = os.path.join(base_dir, "CSV_Shots")
ece_dir = os.path.join(base_dir, "ECE")
output_fig = os.path.join(base_dir, "Figures", "Te0_vs_ECE_all_shots.png")

shots = [143064, 143069, 143073, 143074, 143075, 143077, 143079]

fig, axes = plt.subplots(4, 2, figsize=(14, 16))
axes = axes.flatten()

for i, shot in enumerate(shots):
    ax = axes[i]
    csv_file = os.path.join(csv_dir, f"shot_{shot}_data.csv")
    ece_file = os.path.join(ece_dir, f"shot_{shot}_ece_profile.csv")
    
    if os.path.exists(csv_file):
        df_csv = pd.read_csv(csv_file)
        if 'Time (s)' in df_csv.columns and 'Te0 (keV)' in df_csv.columns:
            # Drop nans just in case
            df_plot = df_csv.dropna(subset=['Te0 (keV)'])
            ax.plot(df_plot['Time (s)'], df_plot['Te0 (keV)'], label='Te0 (Thomson)', color='blue', alpha=0.7, linewidth=1.5)
            
    if os.path.exists(ece_file):
        df_ece = pd.read_csv(ece_file)
        te_cols = [c for c in df_ece.columns if c.startswith('Te_ch')]
        if 'Time(s)' in df_ece.columns and te_cols:
            # Proxy for central ECE Te is the max across all channels for each time step
            df_ece['Te_max_ECE'] = df_ece[te_cols].max(axis=1)
            # Remove zeros or negative values which might be bad data
            df_ece['Te_max_ECE'] = df_ece['Te_max_ECE'].replace(0, np.nan)
            df_plot_ece = df_ece.dropna(subset=['Te_max_ECE'])
            ax.plot(df_plot_ece['Time(s)'], df_plot_ece['Te_max_ECE'], label='ECE (Max Te)', color='red', alpha=0.7, linewidth=1.5)
            
    ax.set_title(f"Shot {shot}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Temperature (keV)")
    
    # Plotting range focusing on the flattop mostly, and sensible temperatures
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8) 
    ax.legend()
    ax.grid(True)

# Turn off the last unused subplot
axes[-1].axis('off')
plt.tight_layout()

os.makedirs(os.path.dirname(output_fig), exist_ok=True)
plt.savefig(output_fig, dpi=150)
print(f"Saved figure to {output_fig}")
