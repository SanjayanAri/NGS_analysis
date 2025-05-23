import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#load the file in csv
data = pd.read_csv('Mobile_OG.csv', index_col=0)
# Create a heatmap
plt.figure(figsize=(16, 12))
sns.heatmap(data, annot=True, cmap='Reds', linewidths=0.5, linecolor='black', cbar_kws={"shrink": .8})

# Set the title and labels
plt.title('Heatmap of mobile genetic elements', fontsize=16)
plt.xlabel('Mobile genetic elements', fontsize=14)
plt.ylabel('Isolates', fontsize=14)

# Save the heatmap as a PNG file
plt.savefig('heatmap.png', dpi=1200, bbox_inches='tight')
plt.show()
# Save the heatmap as a PDF file    
plt.savefig('heatmap.pdf', dpi=1200, bbox_inches='tight')
# show the heatmap
plt.tight_layout()
plt.show()
