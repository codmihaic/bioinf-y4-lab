# Săptămâna 1 — Databases & GitHub 

Scopuri:
1) Toți studentii rulează același mediu (Codespaces sau Docker + Jupyter).
2) Toți studentii trec o dată prin fluxul GitHub (fork → branch → PR) cu un task scurt.

> În săptămânile 2–12 vom **reutiliza datele** colectate devreme (GEO/TCGA/NCBI/Ensembl) pentru aliniere, NGS, filogenie, co-expresie, ML, etc. (vezi calendarul). 

---

## Partea 0 — Verificare mediu 

Alege **una**:

### A) Codespaces
1. Deschide un Codespace pe branch `main`.
2. Terminal:
   ```bash
   python labs/00_smoke/smoke.py
   ```
Așteptat: ok.
B) Docker local (Windows PowerShell)
```powershell
docker pull ghcr.io/bozdogalex/bioinf-y4-lab:base
docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base `
  bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
```
Deschide http://localhost:8890/lab și rulează print("ok").

**Notă (căi cu & pe Windows)**: când rulezi scripturi din acest folder, pune calea între ghilimele:
python "labs/01_intro&databases/demo_entrez_brca1.py"

## Partea 1 — Task PR 
Adaugă handle-ul tău GitHub în labs/01_intro&databases/roster/handles.csv (format: Nume Prenume,github_handle) prin Pull Request.

Pași:

- fork repo → creează branch feat/roster-<handle>.
- edit labs/01_intro&databases/roster/handles.csv → adaugă o singură linie.
- git commit -m "Add <handle> to roster" → git push → creează Pull Request.
- Completează checklist-ul PR.

**Checklist PR:**
 - Am rulat Partea 0 (Codespaces sau Docker) cu succes.
 - Am adăugat exact un rând în handles.csv, fără spații în github_handle.
 - PR trece verificările CI.
Criterii de acceptare: format corect CSV + PR verde + PR îmbinat.

## Partea 2 - Demo Biopython (NCBI Entrez)
- Rulează:
```bash
python "labs/01_intro&databases/demo_entrez_brca1.py"
```
Efect: Descarcă un GenBank BRCA1 (H. sapiens), scrie data/brca1.gb și afișează rezumat. Aliniat cu instrumentele oficiale recomandate (NCBI/Entrez, Biopython).

**Competențe:** 
Studentul poate rula mediul reproducibil și poate deschide un PR corect (fork → branch → PR) — fundament pentru proiectul și temele ulterioare.

**Lectură** : [labs_01_intro&databases](docs/lab_onepagers/01_intro&databases.md) 
[labs/01_intro&databases](labs/01_intro&databases)