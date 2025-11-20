from Bio import SeqIO  # pentru citirea FASTA
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

if __name__ == "__main__":
    handle = "codmihaic"

    fasta_path = Path(f"data/work/{handle}/lab04/assign_seq.fasta")
    if not fasta_path.exists():
        raise FileNotFoundError(f"Nu găsesc fișierul: {fasta_path}")

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError("Fișierul FASTA nu conține secvențe.")

    seq_ids = []
    feature_rows = []

    for rec in records:
        seq_ids.append(rec.id)
        seq = str(rec.seq).upper()
        length = len(seq)

        if length == 0:
            gc = freq_A = freq_C = freq_G = freq_T = 0.0
        else:
            count_A = seq.count("A")
            count_C = seq.count("C")
            count_G = seq.count("G")
            count_T = seq.count("T")

            gc = (count_G + count_C) / length
            freq_A = count_A / length
            freq_C = count_C / length
            freq_G = count_G / length
            freq_T = count_T / length

        feature_rows.append(
            [length, gc, freq_A, freq_C, freq_G, freq_T]
        )

    df = pd.DataFrame(
        feature_rows,
        columns=["Length", "GC_Content", "Freq_A", "Freq_C", "Freq_G", "Freq_T"],
    )
    df.insert(0, "Sequence_ID", seq_ids)

    X = df.drop(columns=["Sequence_ID"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)

    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(12, 6))
    dendrogram(Z, labels=df["Sequence_ID"].tolist())
    plt.title("Hierarchical clustering - dendrogram (assign_seq.fasta)")
    plt.xlabel("Secvențe")
    plt.ylabel("Distanță")

    hierarchical_path = output_dir / f"hierarchical_{handle}.png"
    plt.tight_layout()
    plt.savefig(hierarchical_path, dpi=150)
    plt.close()

    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        X_pca[:, 0],
        X_pca[:, 1],
        c=kmeans_labels,
        alpha=0.7
    )
    plt.title("K-means pe assign_seq.fasta - PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.colorbar(scatter, label="Cluster K-means")

    kmeans_path = output_dir / f"kmeans_{handle}.png"
    plt.tight_layout()
    plt.savefig(kmeans_path, dpi=150)
    plt.close()

    sil_kmeans = silhouette_score(X_scaled, kmeans_labels)
    print(f"Silhouette score K-means: {sil_kmeans:.3f}")

    dbscan = DBSCAN(eps=1.5, min_samples=3)
    db_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = db_labels

    plt.figure(figsize=(8, 6))
    scatter_db = plt.scatter(
        X_pca[:, 0],
        X_pca[:, 1],
        c=db_labels,
        alpha=0.7
    )
    plt.title("DBSCAN pe assign_seq.fasta - PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    cbar = plt.colorbar(scatter_db)
    cbar.set_label("Cluster DBSCAN (−1 = noise)")

    dbscan_path = output_dir / f"dbscan_{handle}.png"
    plt.tight_layout()
    plt.savefig(dbscan_path, dpi=150)
    plt.close()

    mask_core = db_labels != -1
    if mask_core.sum() > 1 and len(set(db_labels[mask_core])) > 1:
        sil_dbscan = silhouette_score(X_scaled[mask_core], db_labels[mask_core])
        print(f"Silhouette score DBSCAN (fără noise): {sil_dbscan:.3f}")
    else:
        print("Nu se poate calcula Silhouette pentru DBSCAN (prea puține clustere / puncte).")

    clusters_path = output_dir / f"clusters_{handle}.csv"
    df[["Sequence_ID", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
        clusters_path, index=False
    )
    print(f"Rezultate salvate în: {clusters_path}")
