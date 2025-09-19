# Step 1: Data Acquisition
import pandas as pd
import numpy as np

# Load the expression data (after downloading and extracting from GEO)
# Example file: "liver_expression.csv"

# Load the expression matrix (e.g., GTEx_Expression_Matrix.csv)
expression_data = pd.read_csv("GTEx_Expression_Matrix.csv", index_col=0)
# Load the metadata file (e.g., GTEx_Metadata.csv)
metadata = pd.read_csv("GTEx_Metadata.csv")
#To isolate liver specific data, load the data into a data frame, and explore the meta data for column like tissue or SMTS. Liver issue might be labeled as liver.
#Then filder the metadata to isolate liver-specific samples and subset the expression ddata by sample ID. Finally, save the isolated liver data for usa
print(metadata.head())
print(metadata.columns)
# Filter rows for liver tissue
liver_samples = metadata[metadata['tissue'] == "Liver"]
# Extract the sample IDs for liver tissue
liver_sample_ids = liver_samples['SampleID'].tolist()
# Subset expression data for liver samples
liver_expression_data = expression_data[liver_sample_ids]
liver_expression_data.to_csv("Liver_Expression_Data.csv")
print(f"Liver expression data dimensions: {liver_expression_data.shape}")
# Extract the sample IDs for liver tissue
liver_sample_ids = liver_samples['SampleID'].tolist()

# Load the expression data (after downloading and extracting from GEO)
data = pd.read_csv("liver_expression.csv", index_col=0)

# Display dataset dimensions
print(f"Dataset dimensions: {data.shape}")

# Apply log transformation to normalize expression levels
data_log = np.log2(data + 1)

# Filter genes with low variance
variance_threshold = 0.5
filtered_data = data_log.loc[data_log.var(axis=1) > variance_threshold]

print(f"Filtered data dimensions: {filtered_data.shape}")

# Step 2: Constructing the Co-Expression Network
from scipy.stats import spearmanr
import networkx as nx

# Calculate the Spearman correlation matrix
correlation_matrix = filtered_data.T.corr(method='spearman')

# Apply a threshold to create an adjacency matrix
threshold = 0.5
adjacency_matrix = (correlation_matrix >= threshold).astype(int)

# Convert the adjacency matrix into a graph
G = nx.from_pandas_adjacency(adjacency_matrix)

print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

# Step 3: Detecting Modules
from networkx.algorithms.community import louvain_communities

# Detect communities using the Louvain algorithm
communities = louvain_communities(G, seed=42)

# Map modules to genes
module_mapping = {node: idx + 1 for idx, community in enumerate(communities) for node in community}
modules = pd.DataFrame({'Gene': list(module_mapping.keys()), 'Module': list(module_mapping.values())})

print(modules.head())

# Step 4: Visualizing the Network
import matplotlib.pyplot as plt

# Assign colors to modules for visualization
module_colors = {module: f"C{module % 10}" for module in modules['Module']}
node_colors = [module_colors[module_mapping[node]] for node in G.nodes()]

# Visualize the network
plt.figure(figsize=(12, 10))
nx.draw(G, with_labels=False, node_color=node_colors, node_size=50, edge_color="gray")
plt.title("Gene Co-Expression Network with Modules")
plt.show()

# Step 5: Biological Interpretation with g:Profiler
from gprofiler import GProfiler

# Initialize the g:Profiler tool
gp = GProfiler(return_dataframe=True)

# Select genes from a specific module (e.g., Module 1)
selected_module = 1
module_genes = modules[modules['Module'] == selected_module]['Gene']

print(f"Genes in selected module: {module_genes.tolist()}")

# Perform enrichment analysis
enrichment_results = gp.profile(organism='hsapiens', query=module_genes.tolist())

# Display top enrichment results
print(enrichment_results[['native', 'p_value', 'term_name', 'source']].head())

# Visualize top terms
significant_results = enrichment_results[enrichment_results['p_value'] < 0.01]

plt.figure(figsize=(10, 6))
plt.barh(significant_results['term_name'][:10], -np.log10(significant_results['p_value'][:10]))
plt.xlabel('-Log10(p-value)')
plt.ylabel('Enriched Terms')
plt.title('Top Enriched Terms')
plt.gca().invert_yaxis()
plt.show()
