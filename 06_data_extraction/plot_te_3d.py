import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# Set backend to avoid display issues
matplotlib.use('Agg')

base_dir = "/Users/vicentesepulveda/Documents/Tesis_EAST/Extracted_EAST"
ece_dir = os.path.join(base_dir, "ECE")
output_fig = os.path.join(base_dir, "Figures", "Te_3D_ECE_all_shots.png")

shots = [143064, 143069, 143073, 143074, 143075, 143077, 143079]

fig = plt.figure(figsize=(16, 20))

for i, shot in enumerate(shots):
    ece_file = os.path.join(ece_dir, f"shot_{shot}_ece_profile.csv")
    if not os.path.exists(ece_file):
        print(f"Warning: {ece_file} not found")
        continue
        
    df = pd.read_csv(ece_file)
    
    # Filter time between 0 and 10s to focus on the discharge
    df = df[(df['Time(s)'] >= 0) & (df['Time(s)'] <= 10)].copy()
    
    # Subsample data to reduce plotting time and file size for 3D surfaces
    # 10s at 1kHz = 10,000 points. Subsample to ~500 points (every 20th point)
    df = df.iloc[::20, :]
    
    time = df['Time(s)'].values
    
    r_cols = [c for c in df.columns if c.startswith('R_ch')]
    te_cols = [c for c in df.columns if c.startswith('Te_ch')]
    
    # Use mean R for the grid (R is practically constant)
    r_mean = df[r_cols].mean().values
    
    # Clean up Te matrix (interpolate missing data in time, then fill NaNs with 0)
    # This ensures a continuous surface
    te_df = df[te_cols].interpolate(method='linear', axis=0).fillna(0)
    te_matrix = te_df.values
    
    R, T = np.meshgrid(r_mean, time)
    
    ax = fig.add_subplot(4, 2, i+1, projection='3d')
    surf = ax.plot_surface(R, T, te_matrix, cmap='plasma', edgecolor='none')
    
    ax.set_title(f"Shot {shot}")
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Time (s)")
    ax.set_zlabel("Te (keV)")
    ax.set_zlim(0, 8)
    
    # Adjust viewing angle for better perspective
    ax.view_init(elev=30, azim=230)

plt.tight_layout()
os.makedirs(os.path.dirname(output_fig), exist_ok=True)
plt.savefig(output_fig, dpi=150)
print(f"Saved figure to {output_fig}")
