import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the presence-absence matrix CSV file
# Replace 'presence_absence_matrix.csv' with your actual file path
data = pd.read_csv('mar_preabs.csv', index_col=0, on_bad_lines='skip')

# Preprocess the data: Convert text IDs to 1 (presence) and '-' to 0 (absence)
binary_data = data.applymap(lambda x: 1 if str(x).strip() != '-' else 0)

# Create a heatmap
plt.figure(figsize=(16, 14))  # Adjust the figure size as needed
sns.heatmap(binary_data, cmap='Reds', cbar_kws={"shrink": .8}, linewidths=0.5, linecolor='black', annot=False)

# Set the title and labels
plt.title('Presence-Absence Heatmap', fontsize=16)
plt.xlabel('Genus', fontsize=14)
plt.ylabel('Clusters', fontsize=14)

# Save the heatmap as a PNG file
plt.savefig('presence_absence_heatmap.png', dpi=300, bbox_inches='tight')

# Show the heatmap
plt.tight_layout()
plt.show()