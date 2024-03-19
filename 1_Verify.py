import matplotlib.pyplot as plt
import numpy as np

# Data
overall_quality_factor = [88.1818, 95.7198, 96.9298, 98.3936, 97.6, 99.2, 98.3806, 98.4375, 97.1193]
verify_percentage = [65.38, 88.19, 73.57, 81.11, 86.72, 84.87, 87.04, 84.07, 79.49]
pdb_ids = ['1QP7', '2OF2', '2OFV', '3AC1', '3AD5', '3BYM', '3BYO', '3LCK', '6PDJ']

# Calculate the trend line
z = np.polyfit(overall_quality_factor, verify_percentage, 1)
p = np.poly1d(z)

# Create a colormap based on VERIFY percentage
# Normalization from 0 to 100 for the color bar
colormap = plt.cm.viridis
normalize = plt.Normalize(vmin=0, vmax=100)
colors = [colormap(normalize(value)) for value in verify_percentage]

# Scatter plot
plt.figure(figsize=(10, 6))
scatter = plt.scatter(overall_quality_factor, verify_percentage, color=colors, alpha=0.7)

# Trend line
plt.plot(overall_quality_factor, p(overall_quality_factor), "r--")

# Annotations
for i, pdb_id in enumerate(pdb_ids):
    plt.annotate(pdb_id, (overall_quality_factor[i], verify_percentage[i]), textcoords="offset points", xytext=(0,10), ha='center')

# Plot customization
plt.title('Overall Quality Factor vs VERIFY(%) with PDB IDs')
plt.xlabel('Overall Quality Factor')
plt.ylabel('VERIFY(%)')
plt.grid(True)

# Color bar with correct scale
colorbar = plt.colorbar(scatter, label='VERIFY(%)')
colorbar.set_label('VERIFY(%)')
colorbar.mappable.set_clim(0, 100)

# Save the figure with 300 DPI
plt.savefig('quality_vs_verify_corrected.png', dpi=300)
plt.close()  # Close the figure to prevent it from displaying

# The script assumes that 'quality_vs_verify_corrected.png' will be saved in the current working directory.
