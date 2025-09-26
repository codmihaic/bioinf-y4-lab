"""
Demo 1 — K-means clustering cu PCA (toy dataset: Breast Cancer WDBC)

"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from pathlib import Path

if __name__ == "__main__":
    # Încărcare dataset
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # Preprocesare
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].apply(lambda x: 1 if x == "M" else 0)

    # Standardizare
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # K-means clustering
    kmeans = KMeans(n_clusters=2, random_state=0)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    # Reducere dimensionalitate pentru vizualizare
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    # Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap="viridis", s=50, alpha=0.7)
    plt.title("K-means Clustering (K=2) pe date PCA")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.colorbar(label="Cluster")
    plt.show()

    # Salvare rezultate (opțional, pentru consistență cu exercițiul)
    output_dir = Path("labs/05_clustering/submissions/demo")
    output_dir.mkdir(parents=True, exist_ok=True)
    df[["Diagnosis", "KMeans_Cluster"]].to_csv(output_dir / "clusters_demo.csv", index=False)
