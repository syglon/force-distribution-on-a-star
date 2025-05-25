import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

# Read data from file
data = np.loadtxt('forces.dat')

# Remove NaN and inf values
data = data[np.isfinite(data)]

# Set bounds for filtering (remove outliers)
lower_bound = np.percentile(data, 0.5)
upper_bound = np.percentile(data, 99.5)

# Create filtered dataset within bounds
filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]

# Define fitting function (Holtsmark-like distribution)
def holtsmark_like(x, a, b, c):
    """Approximate function for Holtsmark-like distribution"""
    return a * np.exp(- (np.abs(x - c) / b)**1.5)

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))

# First subplot: Histogram
n, bins, patches = ax1.hist(filtered_data, bins=100,
                            color='skyblue', edgecolor='black',
                            alpha=0.7, density=True)
ax1.set_xlabel('Gravitational Force (F_x)')
ax1.set_ylabel('Normalized Frequency')
ax1.set_title('Filtered Force Distribution (99% of data)')
ax1.grid(True, linestyle='--', alpha=0.5)

# Calculate bin centers for fitting
bin_centers = (bins[:-1] + bins[1:]) / 2

# Second subplot: Curve fit
try:
    # Fit the curve
    params, cov = curve_fit(holtsmark_like, bin_centers, n,
                            p0=[np.max(n), np.std(filtered_data), np.median(filtered_data)])

    # Generate points for smooth curve
    x_fit = np.linspace(min(filtered_data), max(filtered_data), 500)
    y_fit = holtsmark_like(x_fit, *params)

    # Plot the fitted curve
    ax2.plot(x_fit, y_fit, 'r-', linewidth=2)
    ax2.set_xlabel('Gravitational Force (F_x)')
    ax2.set_ylabel('Probability Density')
    ax2.set_title(f'Curve Fit: Holtsmark-like Distribution\na={params[0]:.2f}, b={params[1]:.2f}, c={params[2]:.2f}')
    ax2.grid(True, linestyle='--', alpha=0.5)

    # Compute R^2
    y_pred = holtsmark_like(bin_centers, *params)
    ss_res = np.sum((n - y_pred) ** 2)
    ss_tot = np.sum((n - np.mean(n)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

except RuntimeError:
    # Fallback if fitting fails
    kde = stats.gaussian_kde(filtered_data)
    x_fit = np.linspace(min(filtered_data), max(filtered_data), 500)
    y_fit = kde(x_fit)
    ax2.plot(x_fit, y_fit, 'r-', linewidth=2)
    ax2.set_title('Kernel Density Estimate (Curve Fit Failed)')
    r_squared = np.nan

# Compute mean and standard deviation
mean_val = np.mean(filtered_data)
std_val = np.std(filtered_data)

# Add annotations to both plots
r2_text = f"$R^2$ = {r_squared:.4f}" if not np.isnan(r_squared) else "$R^2$ = N/A"
annotation = f"Mean = {mean_val:.4f}\nStd Dev = {std_val:.4f}\n{r2_text}"

for ax in (ax1, ax2):
    ax.text(0.97, 0.95, annotation, transform=ax.transAxes,
            fontsize=12, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()
