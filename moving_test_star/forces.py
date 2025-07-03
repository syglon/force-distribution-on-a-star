import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

# Read data from file (only the first column)
data = np.loadtxt('all_data.dat', usecols=(0,))  

# Remove NaN and inf values
data = data[np.isfinite(data)]

# Set bounds for filtering (remove outliers)
lower_bound = np.percentile(data, 0.5)
upper_bound = np.percentile(data, 99.5)

# Create filtered dataset within bounds
filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]

# Define fitting functions
def holtsmark_like(x, a, b, c):
    """Approximate function for Holtsmark-like distribution"""
    return a * np.exp(- (np.abs(x - c) / b)**1.5)

def gaussian(x, a, mu, sigma):
    """Gaussian distribution"""
    return a * np.exp(-0.5 * ((x - mu) / sigma)**2)

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))

# First subplot: Histogram with fits
n, bins, patches = ax1.hist(filtered_data, bins=100,
                            color='skyblue', edgecolor='black',
                            alpha=0.7, density=True, label='Data')
ax1.set_xlabel('Gravitational Force (F_x)')
ax1.set_ylabel('Normalized Frequency')
ax1.set_title('Force Distribution with Fits')
ax1.grid(True, linestyle='--', alpha=0.5)

# Calculate bin centers for fitting
bin_centers = (bins[:-1] + bins[1:]) / 2

# Generate points for smooth curves
x_fit = np.linspace(min(filtered_data), max(filtered_data), 500)

# Try fitting Holtsmark
try:
    params_holtsmark, _ = curve_fit(holtsmark_like, bin_centers, n,
                                   p0=[np.max(n), np.std(filtered_data), np.median(filtered_data)])
    y_holtsmark = holtsmark_like(x_fit, *params_holtsmark)
    ax1.plot(x_fit, y_holtsmark, 'r-', linewidth=2, label='Holtsmark Fit')
    
    # Compute R² for Holtsmark
    y_pred_holtsmark = holtsmark_like(bin_centers, *params_holtsmark)
    ss_res_holtsmark = np.sum((n - y_pred_holtsmark) ** 2)
    ss_tot_holtsmark = np.sum((n - np.mean(n)) ** 2)
    r2_holtsmark = 1 - (ss_res_holtsmark / ss_tot_holtsmark)
except RuntimeError:
    r2_holtsmark = np.nan
    print("Holtsmark fit failed.")

# Try fitting Gaussian
try:
    params_gaussian, _ = curve_fit(gaussian, bin_centers, n,
                                  p0=[np.max(n), np.mean(filtered_data), np.std(filtered_data)])
    y_gaussian = gaussian(x_fit, *params_gaussian)
    ax1.plot(x_fit, y_gaussian, 'g--', linewidth=2, label='Gaussian Fit')
    
    # Compute R² for Gaussian
    y_pred_gaussian = gaussian(bin_centers, *params_gaussian)
    ss_res_gaussian = np.sum((n - y_pred_gaussian) ** 2)
    ss_tot_gaussian = np.sum((n - np.mean(n)) ** 2)
    r2_gaussian = 1 - (ss_res_gaussian / ss_tot_gaussian)
except RuntimeError:
    r2_gaussian = np.nan
    print("Gaussian fit failed.")

ax1.legend()

# Second subplot: Residuals (difference between data and fits)
if not np.isnan(r2_holtsmark) and not np.isnan(r2_gaussian):
    residuals_holtsmark = n - y_pred_holtsmark
    residuals_gaussian = n - y_pred_gaussian
    ax2.plot(bin_centers, residuals_holtsmark, 'r-', linewidth=1, label='Holtsmark Residuals')
    ax2.plot(bin_centers, residuals_gaussian, 'g--', linewidth=1, label='Gaussian Residuals')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_xlabel('Gravitational Force (F_x)')
    ax2.set_ylabel('Residuals')
    ax2.set_title('Residuals (Data - Fit)')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend()

# Compute mean and standard deviation
mean_val = np.mean(filtered_data)
std_val = np.std(filtered_data)

# Add annotations to both plots
annotation = f"Mean = {mean_val:.4f}\nStd Dev = {std_val:.4f}\n"
if not np.isnan(r2_holtsmark):
    annotation += f"Holtsmark $R^2$ = {r2_holtsmark:.4f}\n"
if not np.isnan(r2_gaussian):
    annotation += f"Gaussian $R^2$ = {r2_gaussian:.4f}\n"

for ax in (ax1, ax2):
    ax.text(0.97, 0.95, annotation, transform=ax.transAxes,
            fontsize=12, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()

# Print which fit is better
if not np.isnan(r2_holtsmark) and not np.isnan(r2_gaussian):
    if r2_holtsmark > r2_gaussian:
        print("Holtsmark fit is better (higher R²).")
    else:
        print("Gaussian fit is better (higher R²).")
elif not np.isnan(r2_holtsmark):
    print("Only Holtsmark fit succeeded.")
elif not np.isnan(r2_gaussian):
    print("Only Gaussian fit succeeded.")
else:
    print("Both fits failed.")
