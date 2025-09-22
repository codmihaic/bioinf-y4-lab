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
### B) Docker local (Windows PowerShell)
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
- Completează checklist-ul PR (șablonul week1_roster.md).

## Partea 2 - Demo / Exercitii
**Rulati**
- `demo01_entrez_brca1.py` — căutare + descărcare BRCA1 (GenBank) și sumar GC.
- `demo02_seq_ops.py` — operații de bază pe secvență (transcriere, traducere, GC, motif).
- `dem03_dbsnp.py` — dbSNP: interogare rapidă și sumar.
**Completati si rulati**
- `ex01_multifasta_gc.py` — **schelet** pentru multi-FASTA + GC (de completat în laborator).

## Deliverables
În PR trebuie să apară:

1. O linie nouă în `labs/01_intro&databases/roster/handles.csv`.
2. Un fișier `labs/01_intro&databases/<github_handle>_notes.md` care conține:
   - Confirmarea că ai rulat toate scripturile demo și exercițiul.
   - Un rezultat simplu observat (ex.: “GC fraction = 0.47”).
3. Completarea checklist-ului din șablonul PR (`week1_roster.md`).


## **Competențe:** 
- Rularea mediului reproducibil (Codespaces/Docker).
- Deschiderea și completarea corectă a unui PR (fork → branch → PR).
- Primele interogări și operații de bază pe secvențe biologice.

## Resurse și lecturi recomandate : 
- [Fișa scurtă (One-Pager)](../../docs/lab_onepagers/01_intro&databases.md)  


