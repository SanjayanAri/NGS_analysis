import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#load the file in csv
data = pd.read_csv('mar_paper.csv', index_col=0)

# Convert non-numeric values to NaN and fill them with 0
data = data.apply(pd.to_numeric, errors='coerce').fillna(0)


# Create a heatmap
plt.figure(figsize=(14, 10))
sns.heatmap(data, annot=False, cmap='YlGnBu', linewidths=0.5, linecolor='black', cbar_kws={"shrink": .8})

# Set the title and labels
plt.title('Heatmap of virulence gene copies', fontsize=16)
plt.xlabel('Species', fontsize=14)
plt.ylabel('Virulence genes', fontsize=14)

# Save the heatmap as a PNG file
plt.savefig('heatmap.png', dpi=300, bbox_inches='tight')
plt.show()
# Save the heatmap as a PDF file    
plt.savefig('heatmap.pdf', dpi=300, bbox_inches='tight')
# show the heatmap
plt.tight_layout()
plt.show()
