# GDPR & Data Policy

- Only **public datasets** (e.g. TCGA, GEO) are used in this course.
- Students must not upload clinical or personal data.
- Reports and notebooks must anonymize IDs.
- Projects must clearly state dataset origin and license.

## Data Handling Policy — Bioinformatics & Functional Genomics Lab

This course works with real **bioinformatics datasets** (FASTA/FASTQ, GEO, TCGA, Ensembl, NCBI, etc.).  
To keep the repository clean, fast, and reproducible, we follow these rules:

---

## 1. Local Data (Ignored by Git)
- All files in `data/` are **ignored by Git** (see `.gitignore`).
- Students should save their results (FASTA, GenBank, CSVs, NGS outputs, etc.) under:
```yaml
data/curated/<github_handle>_...
```

- These files are **local only** and not committed.

---

### 2. Public Samples (Versioned in Git)
- The folder `data/sample/` contains **tiny curated files** (≤100 KB) used in:
- CI checks,
- example notebooks/scripts,
- reproducible demos.
- These are the only datasets tracked in the repo.

### 3. Sources of Truth
- All larger data is pulled from **public portals** (Entrez, GEO, TCGA, Ensembl).
- Students must reference the **record IDs/URLs** in their deliverables (e.g., `GSE12345`, `ENSG00000141510`).
- This ensures results are **reproducible** without shipping raw data.

### 4. PR Deliverables
When submitting assignments via Pull Request:
- Commit only **code + metadata samples (≤10 rows)**.
- Never commit raw FASTQ/FASTA/NGS output.
- Example PR contents:
```markdown
labs/02_alignment/my_alignment.py
labs/02_alignment/sample.csv
NOTES.md (explain source IDs + method)
```

### 5. Reuse Across Labs
- Data collected in **Weeks 1–3** (metadata, small samples) will be reused in:
- Alignment,
- NGS validation,
- Phylogenetics,
- Co-expression networks,
- ML classification,
- Graph analysis (diseasome/GNN),
- Federated learning,
- Digital Twin / Quantum demos.
- Students work with **the same core records** throughout, exploring new methods on familiar data.