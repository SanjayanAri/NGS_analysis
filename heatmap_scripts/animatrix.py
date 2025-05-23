import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the FASTANI matrix CSV file
# Replace 'fastani_matrix.csv' with the actual filename
data = pd.read_csv('M_ANI.csv', index_col=0)

# Convert non-numeric values to NaN and fill them with 0 (if necessary)
data = data.apply(pd.to_numeric, errors='coerce').fillna(0)

# Create a heatmap with annotations
plt.figure(figsize=(15, 13))  # Adjust the figure size as needed
sns.heatmap(data, annot=True, fmt=".2f", cmap='Reds', linewidths=0.5, linecolor='black', cbar_kws={"shrink": .8})

# Set the title and labels
plt.title('FASTANI Matrix Heatmap', fontsize=16)
plt.xlabel('Genomes', fontsize=14)
plt.ylabel('Genomes', fontsize=14)

# Save the heatmap as a PNG file
plt.savefig('fastani_heatmap.png', dpi=300, bbox_inches='tight')

# Show the heatmap
plt.tight_layout()
plt.show()