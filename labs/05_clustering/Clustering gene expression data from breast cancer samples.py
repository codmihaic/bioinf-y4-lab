# Hands-on Exercise: Clustering gene expression data from breast cancer samples
# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
#optional visualization libraries
#import numpy as np
#import seaborn as sns


# Load the dataset (direct URL from UCI Repository)
url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data'
columns = ['ID', 'Diagnosis'] + [f'Feature_{i}' for i in range(1, 31)]
df = pd.read_csv(url, header=None, names=columns)

# Drop the ID column and convert 'Diagnosis' to numerical for simplicity
df = df.drop(columns=['ID'])
df['Diagnosis'] = df['Diagnosis'].apply(lambda x: 1 if x == 'M' else 0)

# Display first few rows to understand the dataset
print("Breast Cancer Gene Expression Dataset (first few rows):")
print(df.head())


# Standardize the Data
# Separate features for clustering
X = df.drop(columns=['Diagnosis'])

# Standardize the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Hierarchical Clustering
# Perform hierarchical clustering using 'average' linkage
Z = linkage(X_scaled, method='average')

# Plot the dendrogram
plt.figure(figsize=(12, 8))
dendrogram(Z, labels=df['Diagnosis'].values, leaf_rotation=90, color_threshold=0.7*max(Z[:,2]))
plt.title("Hierarchical Clustering Dendrogram (Average Linkage)")
plt.xlabel("Sample Index")
plt.ylabel("Distance")
plt.show()

# K-means Clustering
# Define the model and fit it to the data
kmeans = KMeans(n_clusters=2, random_state=0)
kmeans_labels = kmeans.fit_predict(X_scaled)

# Add the cluster labels to the original dataframe
df['KMeans_Cluster'] = kmeans_labels

# Visualization
# Reduce data to two dimensions with PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Plot the K-means clusters
plt.figure(figsize=(8, 6))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap='viridis', s=50, alpha=0.7)
plt.title("K-means Clustering (K=2) on PCA-Reduced Data")
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
plt.colorbar(label="Cluster Label")
plt.show()


# DBSCAN for Density-Based Clustering
# Define and fit the DBSCAN model
dbscan = DBSCAN(eps=1.5, min_samples=5)
dbscan_labels = dbscan.fit_predict(X_scaled)

# Add DBSCAN cluster labels to the dataframe
df['DBSCAN_Cluster'] = dbscan_labels
# Visualization
# Plot DBSCAN clusters
plt.figure(figsize=(8, 6))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=dbscan_labels, cmap='plasma', s=50, alpha=0.7)
plt.title("DBSCAN Clustering on PCA-Reduced Data")
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
plt.colorbar(label="Cluster Label")
plt.show()
# Summary of Results
# Show a sample of the dataset with clustering results
print("Sample Clustering Results:")
print(df[['Diagnosis', 'KMeans_Cluster', 'DBSCAN_Cluster']].head(10))
