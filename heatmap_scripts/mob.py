import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet

# Create the data as a DataFrame
data = pd.DataFrame({
    'Samples': ['ISO_NA', 'ISO_SA', 'ISO_LZ', 'ISO_MZ', 'ISO_SV'],
    'Integration, Excision': [12, 34, 1, 0, 0],
    'Repair Replication Recombination': [15, 18, 17, 17, 16],
    'Phage': [11, 11, 12, 14, 12],
    'Stability Transfer Defense': [1, 1, 0, 0, 0],
    'Transfer': [1, 1, 1, 2, 1]
})

# Set 'Samples' as the index
data.set_index('Samples', inplace=True)

# Convert the data to a boolean format for UpSet plotting
boolean_data = data > 0  # Convert values to boolean (True/False)

# Plot the UpSet graph
upset = UpSet(boolean_data, subset_size='count', show_counts=True, orientation='horizontal')
plt.figure(figsize=(10, 6))
upset.plot()
plt.title('UpSet Graph of Mobile Genetic Elements')
plt.savefig('upset_graph.png', dpi=300, bbox_inches='tight')
plt.show()