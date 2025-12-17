"""
Exercise 10.2 — Identify top SNP–Gene correlations

TODO:
- încărcați matricea integrată multi-omics
- împărțiți rândurile în SNPs vs gene (după indice sau după nume)
- calculați corelații între fiecare SNP și fiecare genă
- filtrați |r| > 0.5
- exportați snp_gene_pairs_<handle>.csv
"""

from pathlib import Path
import pandas as pd
import numpy as np

HANDLE = "codmihaic"
JOINT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/multiomics_concat_{HANDLE}.csv")

OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

# TODO: load joint matrix, compute correlations, export

def main():
    if not JOINT_CSV.exists():
        raise FileNotFoundError(f"Missing joint matrix: {JOINT_CSV}")
    joint = pd.read_csv(JOINT_CSV, index_col=0)

    snp_cols = [c for c in joint.columns if c.startswith("SNP_")]
    gene_cols = [c for c in joint.columns if c.startswith("EXP_")]
    if len(snp_cols) == 0 or len(gene_cols) == 0:
        raise ValueError("Could not identify SNP or gene columns. Check prefixes.")

    snp_mat = joint[snp_cols]
    gene_mat = joint[gene_cols]
    results = []

    for snp in snp_cols:
        for gene in gene_cols:
            r = np.corrcoef(snp_mat[snp], gene_mat[gene])[0, 1]
            if np.isnan(r):
                continue
            if abs(r) > 0.5:
                results.append({
                    "SNP": snp,
                    "Gene": gene.replace("EXP_", ""),
                    "r": r
                })

    res_df = pd.DataFrame(results)
    res_df = res_df.sort_values(by="r", key=np.abs, ascending=False)
    res_df.to_csv(OUT_CSV, index=False)

    print(f"Found {len(res_df)} SNP–gene pairs with |r| > 0.5")
    print("Saved:", OUT_CSV)
    print(res_df.head(10))


if __name__ == "__main__":
    main()