from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

HANDLE = "codmihaic"

SNP_FILE = Path(f"data/work/codmihaic/lab10/snp_matrix_{HANDLE}.csv")
EXPR_FILE = Path(f"data/work/codmihaic/lab10/expression_matrix_{HANDLE}.csv")
OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_JOINT = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"
OUT_PCA_SNP = OUT_DIR / f"pca_snp_{HANDLE}.png"
OUT_PCA_EXPR = OUT_DIR / f"pca_expr_{HANDLE}.png"
OUT_PCA_JOINT = OUT_DIR / f"pca_joint_{HANDLE}.png"
OUT_PAIRS = OUT_DIR / f"snp_gene_pairs_{HANDLE}.csv"
OUT_CLUSTERS = OUT_DIR / f"clusters_{HANDLE}.csv"
OUT_CTAB = OUT_DIR / f"cluster_vs_subtype_{HANDLE}.csv"

K_MEANS = 3
CORR_THRESHOLD = 0.5

# -----------------------------
# Helpers
# -----------------------------
def zscore_by_feature(df: pd.DataFrame) -> pd.DataFrame:
    mean = df.mean(axis=0)
    std = df.std(axis=0, ddof=1).replace(0, np.nan)
    out = (df - mean) / std
    return out.fillna(0.0)


def ensure_samples_by_features(
    df: pd.DataFrame, kind: str, snp_samples: pd.Index | None = None
) -> pd.DataFrame:
    if kind == "snp":
        return df
    
    if snp_samples is not None:
        col_overlap = len(set(df.columns) & set(snp_samples))
        idx_overlap = len(set(df.index) & set(snp_samples))
        if col_overlap >= idx_overlap:
            return df.T
        else:
            return df 

    def looks_like_samples(index: pd.Index) -> int:
        return sum(str(x).startswith("Sample_") for x in index)
    if looks_like_samples(df.columns) > looks_like_samples(df.index):
        return df.T
    return df


def pca_plot_2d(X: np.ndarray, labels: np.ndarray | None, title: str, out_png: Path) -> None:
    pca = PCA(n_components=2, random_state=0)
    Z = pca.fit_transform(X)

    plt.figure(figsize=(7, 5))
    if labels is None:
        plt.scatter(Z[:, 0], Z[:, 1], s=35)
    else:
        labels = np.asarray(labels)
        for lab in np.unique(labels):
            m = labels == lab
            plt.scatter(Z[m, 0], Z[m, 1], s=35, label=str(lab))
        plt.legend(title="group", fontsize=9)

    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()


def compute_cross_correlation_pairs(
    Xz: pd.DataFrame, 
    Yz: pd.DataFrame,
    thresh: float,
) -> pd.DataFrame:

    n = Xz.shape[0]
    X = Xz.to_numpy(dtype=float)
    Y = Yz.to_numpy(dtype=float)
    corr = (X.T @ Y) / (n - 1) 

    mask = np.abs(corr) > thresh
    snp_idx, gene_idx = np.where(mask)
    rows = []
    for i, j in zip(snp_idx, gene_idx):
        rows.append((Xz.columns[i], Yz.columns[j], float(corr[i, j])))

    out = pd.DataFrame(rows, columns=["snp", "gene", "r"])
    out["abs_r"] = out["r"].abs()
    out = out.sort_values(["abs_r"], ascending=False).drop(columns=["abs_r"]).reset_index(drop=True)
    return out


def main():
    # -------------------------
    # Task 1: Load + harmonize + zscore + concat
    # -------------------------
    snp = pd.read_csv(SNP_FILE, index_col=0)
    expr_raw = pd.read_csv(EXPR_FILE, index_col=0)
    expr = ensure_samples_by_features(expr_raw, kind="expr", snp_samples=snp.index)
    common = snp.index.intersection(expr.index)
    if len(common) == 0:
        raise ValueError("No common sample IDs found between SNP and Expression matrices.")

    snp = snp.loc[common].copy()
    expr = expr.loc[common].copy()
    snp_z = zscore_by_feature(snp)
    expr_z = zscore_by_feature(expr)
    joint = pd.concat([snp_z, expr_z], axis=1)
    joint.to_csv(OUT_JOINT)

    # -------------------------
    # Bonus (and fallback labels): KMeans on joint
    # -------------------------
    km = KMeans(n_clusters=K_MEANS, random_state=0, n_init="auto")
    clusters = km.fit_predict(joint.to_numpy(dtype=float))
    clusters_df = pd.DataFrame(
        {"sample_id": joint.index.astype(str), "cluster": clusters.astype(int)}
    ).set_index("sample_id")
    clusters_df.to_csv(OUT_CLUSTERS)

    # -------------------------
    # Task 2: PCA SNP / Expr / Joint
    # -------------------------
    pca_labels = clusters
    pca_plot_2d(
        snp_z.to_numpy(dtype=float),
        pca_labels,
        title="PCA — SNPs",
        out_png=OUT_PCA_SNP,
    )
    pca_plot_2d(
        expr_z.to_numpy(dtype=float),
        pca_labels,
        title="PCA — Expression",
        out_png=OUT_PCA_EXPR,
    )
    pca_plot_2d(
        joint.to_numpy(dtype=float),
        pca_labels,
        title="PCA — Joint (SNP + Expression)",
        out_png=OUT_PCA_JOINT,
    )

    # -------------------------
    # Task 3: Cross-omics correlation
    # -------------------------
    pairs = compute_cross_correlation_pairs(snp_z, expr_z, thresh=CORR_THRESHOLD)
    pairs.to_csv(OUT_PAIRS, index=False)

    print("DONE!!!")
    print(f"Common samples: {len(common)}")
    print(f"Joint shape: {joint.shape}  (samples x features)")
    print(f"Saved: {OUT_JOINT}")
    print(f"Saved: {OUT_PCA_SNP}")
    print(f"Saved: {OUT_PCA_EXPR}")
    print(f"Saved: {OUT_PCA_JOINT}")
    print(f"Saved: {OUT_PAIRS}   (pairs with |r| > {CORR_THRESHOLD})")
    print(f"Saved: {OUT_CLUSTERS}")


if __name__ == "__main__":
    main()
