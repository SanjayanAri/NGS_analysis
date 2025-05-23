import pandas as pd
from bioservices import KEGG
import matplotlib.pyplot as plt
import seaborn as sns

# Load the gene and KO ID data
# Replace 'gene_ko_mapping.txt' with your actual file path
data = pd.read_csv('ISO_SA_ko.csv', sep='\t', header=None, names=['GeneID', 'KOID'])

# Filter out rows where KOID is missing (NaN)
data = data.dropna(subset=['KOID'])

# Initialize KEGG service
kegg = KEGG()

# Fetch pathway information for each KO ID
pathway_counts = {}
for ko_id in data['KOID']:
    try:
        pathways = kegg.get_pathway_by_gene(ko_id)
        if pathways:
            for pathway in pathways:
                pathway_name = pathway['name']
                pathway_counts[pathway_name] = pathway_counts.get(pathway_name, 0) + 1
    except Exception as e:
        print(f"Error fetching pathway for KO ID {ko_id}: {e}")

# Convert pathway counts to a DataFrame for visualization
pathway_df = pd.DataFrame(list(pathway_counts.items()), columns=['Pathway', 'Count'])

# Sort pathways by count
pathway_df = pathway_df.sort_values(by='Count', ascending=False)

# Plot the top pathways
plt.figure(figsize=(12, 8))
sns.barplot(data=pathway_df.head(10), x='Count', y='Pathway', palette='viridis')
plt.title('Top 10 KEGG Pathways', fontsize=16)
plt.xlabel('Gene Count', fontsize=14)
plt.ylabel('Pathway', fontsize=14)
plt.tight_layout()

# Save the plot
plt.savefig('top_kegg_pathways.png', dpi=300)

# Show the plot
plt.show()