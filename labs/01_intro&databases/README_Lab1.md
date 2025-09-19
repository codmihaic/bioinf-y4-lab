# Week 1 — Databases & GitHub (Onboarding)

This week also includes a short **environment onboarding** to ensure everyone can run notebooks the same way.

## Part 0 — Environment check (5–10 min)

Pick **one**:

### A) Codespaces
1. Open a Codespace on this repo (main branch).
2. In the terminal:  
   ```bash
   python labs/00_smoke/smoke.py
   ```
   Expect `ok`.

### B) Local Docker
```powershell
docker pull ghcr.io/bozdogalex/bioinf-y4:base
docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4:base `
  bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
```
Open http://localhost:8890/lab → run a cell: `print("ok")`

## Part 1 — Databases competency

- Open an NCBI/TCGA/GEO/Ensembl record.
- Extract **two fields** into a CSV.
- Commit via a **pull request**.

Deliverables (PR contents):
- `labs/01_databases/<yourname>/sample.csv`
- A 2–3 line note in the PR describing the source record and fields.
