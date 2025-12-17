import numpy as np
import pandas as pd
from pathlib import Path

HANDLE = "codmihaic"
OUT = Path(f"data/work/{HANDLE}/lab10")
OUT.mkdir(parents=True, exist_ok=True)

np.random.seed(42)
samples = [f"Sample_{i}" for i in range(1, 21)]
n_snps = 600
maf = np.random.uniform(0.05, 0.5, n_snps)

snp_data = []
for p in maf:
    probs = [(1-p)**2, 2*p*(1-p), p**2]
    snp_data.append(np.random.choice([0,1,2], size=len(samples), p=probs))

snp_matrix = pd.DataFrame(
    np.array(snp_data).T,
    index=samples,
    columns=[f"SNP_{i+1}" for i in range(n_snps)]
)

mask = np.random.rand(*snp_matrix.shape) < 0.03
snp_matrix = snp_matrix.mask(mask)
flip_mask = np.random.rand(*snp_matrix.shape) < 0.02
snp_matrix = snp_matrix.where(~flip_mask, 2 - snp_matrix)

out_csv = OUT / f"snp_matrix_{HANDLE}.csv"
snp_matrix.to_csv(out_csv)

print("Saved:", out_csv)
print(snp_matrix.head())
