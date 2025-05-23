import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# data loading
data = pd.read_csv('mar_paper.csv', index_col=0)

#display check
print(data.head())

#setting figure size
plt.figure(figsize=(12, 10))

#heatmap
sns.heatmap(data, cmap='Blues', annot=False, fmt='d', linewidths=0.5, linecolor='black', cbar_kws={"shrink": .8})

#labels
plt.title('virulence genes copies across the strains', fontsize=16)
plt.xlabel('Species', fontsize=14)
plt.ylabel('virulence Genes', fontsize=14)

plt.show()
