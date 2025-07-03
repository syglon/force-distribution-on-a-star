import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Load data - now expecting 8 columns (without relative energy)
data = np.loadtxt('all_data.dat')
F_x, F_y, F_z, deltaE, deltaL, deltaLx, deltaLy, deltaLz, deltaE_from_potential = data.T

# Updated column names (removed relative energy)
column_names = [
    "Force X (F_x)", "Force Y (F_y)", "Force Z (F_z)",
    "Energy Change (ΔE)",
    "Angular Momentum Change (ΔL)",
    "Angular Momentum X Change (ΔLx)", 
    "Angular Momentum Y Change (ΔLy)",
    "Angular Momentum Z Change (ΔLz)"
    "Energy Change from Potential (ΔE_potential)"
]

# Updated output files (removed relative energy plots)
output_files = [
    "force_x_histogram.png", "force_y_histogram.png", "force_z_histogram.png",
    "deltaE_histogram.png",
    "deltaL_histogram.png",
    "deltaLx_histogram.png", "deltaLy_histogram.png", "deltaLz_histogram.png", "deltaE_potential_histogram.png"
]

# 1. Create basic histograms for all variables
print("Creating basic histograms...")
for col, (data_col, name, fname) in enumerate(zip(data.T, column_names, output_files)):
    filtered_data = data_col[np.isfinite(data_col)]
    lower, upper = np.percentile(filtered_data, [0.5, 99.5])
    filtered_data = filtered_data[(filtered_data >= lower) & (filtered_data <= upper)]
    
    plt.figure(figsize=(10, 6))
    n, bins, patches = plt.hist(filtered_data, bins=100, color='skyblue', 
                               edgecolor='black', alpha=0.7, density=True)
    kde = stats.gaussian_kde(filtered_data)
    x = np.linspace(min(filtered_data), max(filtered_data), 500)
    plt.plot(x, kde(x), 'r-', linewidth=2)
    
    stats_text = (f'Mean = {np.mean(filtered_data):.4e}\n'
                 f'Median = {np.median(filtered_data):.4e}\n'
                 f'Std Dev = {np.std(filtered_data):.4e}\n'
                 f'Skewness = {stats.skew(filtered_data):.2f}\n'
                 f'Kurtosis = {stats.kurtosis(filtered_data):.2f}')
    
    plt.title(f'Distribution of {name} (99% of data)')
    plt.xlabel(name)
    plt.ylabel('Normalized Frequency')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.text(0.97, 0.95, stats_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()

# 2. Angular Momentum Analysis
print("Analyzing angular momentum changes...")
# Signed angular momentum change (projection along initial direction)
L_initial_mag = np.sqrt(deltaLx**2 + deltaLy**2 + deltaLz**2)
valid_indices = L_initial_mag > 1e-10  # Avoid division by zero
deltaL_signed = np.zeros_like(deltaL)
deltaL_signed[valid_indices] = (deltaLx[valid_indices]*deltaLx[valid_indices] + 
                               deltaLy[valid_indices]*deltaLy[valid_indices] + 
                               deltaLz[valid_indices]*deltaLz[valid_indices]) / L_initial_mag[valid_indices]

# Plot signed deltaL distribution
plt.figure(figsize=(10, 6))
filtered = deltaL_signed[np.isfinite(deltaL_signed)]
lower, upper = np.percentile(filtered, [0.5, 99.5])
filtered = filtered[(filtered >= lower) & (filtered <= upper)]

n, bins, patches = plt.hist(filtered, bins=100, color='skyblue', 
                           edgecolor='black', alpha=0.7, density=True)
kde = stats.gaussian_kde(filtered)
x = np.linspace(min(filtered), max(filtered), 500)
plt.plot(x, kde(x), 'r-', linewidth=2)

stats_text = (f'Mean = {np.mean(filtered):.4e}\n'
             f'Median = {np.median(filtered):.4e}\n'
             f'Std Dev = {np.std(filtered):.4e}\n'
             f'Positive changes: {100*np.mean(filtered > 0):.1f}%')

plt.title('Distribution of Signed Angular Momentum Change (99% of data)')
plt.xlabel('Signed ΔL (projection along initial L)')
plt.ylabel('Normalized Frequency')
plt.grid(True, linestyle='--', alpha=0.5)
plt.text(0.97, 0.95, stats_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white', alpha=0.8))
plt.tight_layout()
plt.savefig('deltaL_signed_histogram.png', dpi=150)
plt.close()

# 3. Total Force Magnitude Distribution
print("Analyzing total force magnitude...")
force_mag = np.sqrt(F_x**2 + F_y**2 + F_z**2)

plt.figure(figsize=(10, 6))
filtered_force = force_mag[np.isfinite(force_mag)]
lower, upper = np.percentile(filtered_force, [0.5, 99.5])
filtered_force = filtered_force[(filtered_force >= lower) & (filtered_force <= upper)]

n, bins, patches = plt.hist(filtered_force, bins=100, color='skyblue', 
                           edgecolor='black', alpha=0.7, density=True)
kde = stats.gaussian_kde(filtered_force)
x = np.linspace(min(filtered_force), max(filtered_force), 500)
plt.plot(x, kde(x), 'r-', linewidth=2)

stats_text = (f'Mean = {np.mean(filtered_force):.4e}\n'
             f'Median = {np.median(filtered_force):.4e}\n'
             f'Std Dev = {np.std(filtered_force):.4e}\n'
             f'Skewness = {stats.skew(filtered_force):.2f}\n'
             f'Kurtosis = {stats.kurtosis(filtered_force):.2f}')

plt.title('Distribution of Total Force Magnitude (99% of data)')
plt.xlabel('Force Magnitude (sqrt(Fx²+Fy²+Fz²))')
plt.ylabel('Normalized Frequency')
plt.grid(True, linestyle='--', alpha=0.5)
plt.text(0.97, 0.95, stats_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white', alpha=0.8))
plt.tight_layout()
plt.savefig('total_force_magnitude_histogram.png', dpi=150)
plt.close()
# 4. Energy Change from Potential Analysis
print("Analyzing energy changes from potential...")
plt.figure(figsize=(10, 6))
filtered_deltaE_pot = deltaE_from_potential[np.isfinite(deltaE_from_potential)]
lower, upper = np.percentile(filtered_deltaE_pot, [0.5, 99.5])
filtered_deltaE_pot = filtered_deltaE_pot[(filtered_deltaE_pot >= lower) & (filtered_deltaE_pot <= upper)]

n, bins, patches = plt.hist(filtered_deltaE_pot, bins=100, color='skyblue', 
                          edgecolor='black', alpha=0.7, density=True)
kde = stats.gaussian_kde(filtered_deltaE_pot)
x = np.linspace(min(filtered_deltaE_pot), max(filtered_deltaE_pot), 500)
plt.plot(x, kde(x), 'r-', linewidth=2)

stats_text = (f'Mean = {np.mean(filtered_deltaE_pot):.4e}\n'
             f'Median = {np.median(filtered_deltaE_pot):.4e}\n'
             f'Std Dev = {np.std(filtered_deltaE_pot):.4e}\n'
             f'Skewness = {stats.skew(filtered_deltaE_pot):.2f}\n'
             f'Kurtosis = {stats.kurtosis(filtered_deltaE_pot):.2f}')

plt.title('Distribution of Energy Change from Potential (99% of data)')
plt.xlabel('Energy Change from Potential (ΔE_potential)')
plt.ylabel('Normalized Frequency')
plt.grid(True, linestyle='--', alpha=0.5)
plt.text(0.97, 0.95, stats_text, transform=plt.gca().transAxes,
        fontsize=10, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white', alpha=0.8))
plt.tight_layout()
plt.savefig('deltaE_potential_histogram.png', dpi=150)
plt.close()

print("Analysis complete! Created these plots:")
print("1. Basic histograms for all variables")
print("2. deltaL_signed_histogram.png (direction of angular momentum changes)")
print("3. total_force_magnitude_histogram.png (distribution of total force magnitude)")
print(f"\nSum of deltaE array elements: {np.sum(deltaE):.6e}")
