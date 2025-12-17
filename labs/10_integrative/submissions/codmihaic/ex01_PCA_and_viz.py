"""
Exercise 10.1 — PCA Single-Omics vs Joint

TODO:
- încărcați SNP și Expression
- normalizați fiecare strat (z-score)
- rulați PCA pe:
    1) strat SNP
    2) strat Expression
    3) strat Joint (concat)
- generați 3 figuri PNG
- comparați vizual distribuția probelor
"""

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

HANDLE = "codmihaic"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# TODO: load, align, normalize, PCA, export figures

def zscore_features(X: pd.DataFrame) -> pd.DataFrame:
    mu = X.mean(axis=0)
    sigma = X.std(axis=0, ddof=0).replace(0, np.nan)
    Z = (X - mu) / sigma
    return Z.fillna(0.0)


def ensure_samples_rows(df: pd.DataFrame) -> pd.DataFrame:
    idx = df.index.astype(str)
    cols = df.columns.astype(str)

    idx_sample_like = np.mean([s.startswith("Sample_") for s in idx]) if len(idx) else 0
    col_sample_like = np.mean([s.startswith("Sample_") for s in cols]) if len(cols) else 0

    if col_sample_like > idx_sample_like:
        return df.T
    return df


def run_pca_and_plot(X: pd.DataFrame, title: str, out_png: Path, n_components: int = 2):
    pca = PCA(n_components=n_components, random_state=42)
    coords = pca.fit_transform(X.values)

    plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1])
    for i, sample in enumerate(X.index):
        plt.annotate(str(sample), (coords[i, 0], coords[i, 1]), fontsize=7, alpha=0.8)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return pca.explained_variance_ratio_


def main():
    if not SNP_CSV.exists():
        raise FileNotFoundError(f"Missing SNP file: {SNP_CSV}")
    if not EXP_CSV.exists():
        raise FileNotFoundError(f"Missing Expression file: {EXP_CSV}")

    snp = pd.read_csv(SNP_CSV, index_col=0)
    expr = pd.read_csv(EXP_CSV, index_col=0)

    snp = ensure_samples_rows(snp)
    expr = ensure_samples_rows(expr)
    snp = snp.apply(pd.to_numeric, errors="coerce")
    snp = snp.fillna(snp.mean(axis=0)).fillna(0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(expr.mean(axis=0)).fillna(0)

    common = snp.index.intersection(expr.index)
    if len(common) < 3:
        raise ValueError(f"Too few common samples after alignment: {len(common)}. Common: {list(common)}")

    snp = snp.loc[common].copy()
    expr = expr.loc[common].copy()
    snp_z = zscore_features(snp)
    expr_z = zscore_features(expr)
    joint = pd.concat([snp_z.add_prefix("SNP_"), expr_z.add_prefix("EXP_")], axis=1)

    ev_snp = run_pca_and_plot(
        snp_z, "PCA — SNP layer", OUT_DIR / f"pca_snp_{HANDLE}.png"
    )
    ev_expr = run_pca_and_plot(
        expr_z, "PCA — Expression layer", OUT_DIR / f"pca_expression_{HANDLE}.png"
    )
    ev_joint = run_pca_and_plot(
        joint, "PCA — Joint (SNP + Expression)", OUT_DIR / f"pca_joint_{HANDLE}.png"
    )

    print("Explained variance ratio (PC1, PC2):")
    print("  SNP:       ", [round(x, 4) for x in ev_snp[:2]])
    print("  Expression:", [round(x, 4) for x in ev_expr[:2]])
    print("  Joint:     ", [round(x, 4) for x in ev_joint[:2]])
    print("\nSaved PNGs in:", OUT_DIR)

    joint_out = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"
    joint.to_csv(joint_out)
    print("Saved joint matrix:", joint_out)


if __name__ == "__main__":
    main()
