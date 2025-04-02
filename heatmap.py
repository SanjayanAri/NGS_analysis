import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# data loading
data = pd.read_csv('VFgenes.csv', index_col=0)

#display check
print(data.head())

#setting figure size
plt.figure(figsize=(10, 8))

#heatmap
sns.heatmap(data, cmap='YlOrRd', annot=True, fmt='d', linewidths=0.5)

#labels
plt.title('virulence genes copies across the strains', fontsize=16)
plt.xlabel('Strains', fontsize=14)
plt.ylabel('Genes', fontsize=14)

plt.show()
