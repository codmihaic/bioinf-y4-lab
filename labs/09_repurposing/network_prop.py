#Exploring Network Properties for Drug Repurposing
#Analyze a drug-disease network using the provided dataset to:

#Identify key nodes using centrality metrics.
#Explore shortest paths between specific drugs and diseases.
#Interpret results to gain insights into potential drug repurposing opportunities.
#Step 1: Load and Prepare the Dataset

#Load the dataset from assignment 8 (drug-target, protein-gene, and gene-disease interactions) into Python.
#Construct a graph using NetworkX, ensuring edges are correctly labeled.
import pandas as pd
import networkx as nx

# Load the dataset
data = pd.read_csv("drug_disease_interactions.csv")

# Build a graph
G = nx.Graph()  # Use nx.DiGraph() if directed edges are required
for _, row in data.iterrows():
    G.add_edge(row["Source"], row["Target"], interaction=row["Interaction Type"])
#you must check for missing data and determine whether the graph should be directed or undirected.

#Step 2: Visualize the Network
#Visualize the network to explore its structure.
#Highlight drugs, proteins, genes, and diseases using different colors or shapes.
import matplotlib.pyplot as plt

# Color nodes based on their type
node_colors = []
for node in G.nodes():
    if "Drug" in node:
        node_colors.append("red")
    elif "Disease" in node:
        node_colors.append("blue")
    else:
        node_colors.append("green")

# Visualize the network
plt.figure(figsize=(10, 8))
nx.draw(G, with_labels=True, node_color=node_colors, font_size=8, node_size=500)
plt.show()
# improve the visualization by adding a legend or customizing node sizes.

#Step 3: Calculate Centrality Metrics
#Compute degree and betweenness centrality to identify hubs and bottlenecks.
#Identify and interpret the top 5 nodes for each metric.
# Calculate degree and betweenness centrality
degree_centrality = nx.degree_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G)

# Sort and display top 5 nodes for each metric (students complete this part)
top_degree_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:5]
top_betweenness_nodes = sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)[:5]

#write code to display results clearly (e.g., in a table or plot)

#Step 4: Explore Shortest Paths
#Compute shortest paths between selected drugs and diseases.
#Interpret which drugs are closest to specific diseases and why.
# Shortest paths between specific nodes (complete the node selection)
drug = "Drug_A"
disease = "Disease_1"

# Calculate and display the shortest path
shortest_path = nx.shortest_path(G, source=drug, target=disease)
print(f"Shortest path between {drug} and {disease}: {shortest_path}")
#select specific drug and disease nodes and summarize  results.

#Step 5: Biological Interpretation
#Centrality Results:

#What roles do the top nodes play (e.g., hub proteins, multitarget drugs)?
#Are the top nodes consistent with known biological insights about the dataset?

#Shortest Path Results:

#Which drugs are closest to diseases in terms of path length?
#Do these findings suggest potential candidates for repurposing? Why?

#Discussion Questions

#Network Topology:
#How does the structure (e.g., hubs, clusters) of the network influence drug-disease relationships?

#Biological Relevance:
#Are the top nodes or shortest paths biologically meaningful? Do they align with known roles of drugs, proteins, or genes in disease pathways?

#Limitations:
#What assumptions or limitations in the dataset or network model could affect your conclusions?
