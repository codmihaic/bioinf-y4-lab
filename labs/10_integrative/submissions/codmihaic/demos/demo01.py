"""
Demo 1 — Multi-Omics Alignment (SNPs + Expression)
"""

import pandas as pd
from pathlib import Path

SNP_CSV = Path("data/sample/snp_data.csv")
EXP_CSV = Path("data/sample/expression_data.csv")

df_snp = pd.read_csv(SNP_CSV, index_col=0)
df_exp = pd.read_csv(EXP_CSV, index_col=0)

# Find common samples
common = df_snp.columns.intersection(df_exp.columns)

df_snp = df_snp[common]
df_exp = df_exp[common]

# Z-score normalization
df_snp_norm = (df_snp - df_snp.mean()) / df_snp.std()
df_exp_norm = (df_exp - df_exp.mean()) / df_exp.std()

# Concatenate features
df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=0)

df_joint.to_csv("labs/10_integrative/submissions/codmihaic/multiomics_concat_demo.csv")
print(df_joint.head())
