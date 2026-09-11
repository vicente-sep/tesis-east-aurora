import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Agg')

base_dir = "/Users/vicentesepulveda/Documents/Tesis_EAST/Extracted_EAST"
ece_dir = os.path.join(base_dir, "ECE")
output_fig = os.path.join(base_dir, "Figures", "ECE_Radial_Profiles_All_Shots.png")

shots = [143064, 143069, 143073, 143074, 143075, 143077, 143079]
ecrh_off_shots = [143064, 143077]
time_slices = [0.5, 2.0, 3.5, 5.0, 6.5, 8.0, 9.5] # 7 slices including early ramp-up and late flat-top/ramp-down

fig, axes = plt.subplots(2, 4, figsize=(22, 10), sharex=True, sharey=True)
axes = axes.flatten()

for i, shot in enumerate(shots):
    ax = axes[i]
    ece_file = os.path.join(ece_dir, f"shot_{shot}_ece_profile.csv")
    if not os.path.exists(ece_file):
        ax.axis('off')
        continue
        
    df = pd.read_csv(ece_file)
    
    # Identify R columns and Te columns
    r_cols = sorted([c for c in df.columns if c.startswith('R_ch')])
    te_cols = sorted([c for c in df.columns if c.startswith('Te_ch')])
    
    # Setup colormap for the time slices
    import matplotlib.cm as cm
    cmap = cm.get_cmap('plasma')
    colors = [cmap(k / max(1, len(time_slices) - 1)) for k in range(len(time_slices))]
    
    # Loop over desired time slices
    for j, t_target in enumerate(time_slices):
        # Find the row closest to t_target
        if df.empty or t_target > df['Time(s)'].max():
            continue
            
        idx = (df['Time(s)'] - t_target).abs().idxmin()
        row = df.loc[idx]
        actual_time = row['Time(s)']
        
        # Extract R and Te arrays
        r_vals = row[r_cols].values.astype(float)
        te_vals = row[te_cols].values.astype(float)
        
        # Clean data: drop NaNs and zeros or extreme outliers
        valid_mask = ~np.isnan(r_vals) & ~np.isnan(te_vals) & (te_vals > 0) & (te_vals < 20)
        
        r_valid = r_vals[valid_mask]
        te_valid = te_vals[valid_mask]
        
        # Sort by R just in case
        sort_idx = np.argsort(r_valid)
        r_valid = r_valid[sort_idx]
        te_valid = te_valid[sort_idx]
        
        ax.plot(r_valid, te_valid, marker='o', color=colors[j], linewidth=2, markersize=5, alpha=0.9, label=f"t = {actual_time:.1f} s")
        
    ecrh_status = "ECRH OFF" if shot in ecrh_off_shots else "ECRH ON"
    title_color = "red" if ecrh_status == "ECRH OFF" else "green"
    
    ax.set_title(f"Shot {shot} ({ecrh_status})", color=title_color, fontweight='bold', fontsize=14)
    ax.set_xlabel("Major Radius R (m)", fontsize=12)
    ax.set_ylabel("Te (keV)", fontsize=12)
    
    # Soft background colors
    if ecrh_status == "ECRH OFF":
        ax.set_facecolor('#fff0f0')
    else:
        ax.set_facecolor('#f0fff0')
        
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Show legend on every plot or just the first?
    # Let's show it on all so it's clear what lines are what, since there are only 4 lines.
    ax.legend(loc='upper right', fontsize=10)
        
# Remove the empty subplot
axes[-1].axis('off')
plt.tight_layout()

os.makedirs(os.path.dirname(output_fig), exist_ok=True)
plt.savefig(output_fig, dpi=150)
print(f"Saved figure to {output_fig}")
