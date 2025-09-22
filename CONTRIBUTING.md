# Contributing Guidelines – BIOINF-Y4 Lab

## Pair Work
- Laboratoarele se rezolvă individual sau în perechi (Rol A și Rol B).
- Rolurile se rotesc între laboratoare:
  - Rol A = „driver” (scrie codul, rulează comenzi).
  - Rol B = „navigator” (verifică, întreabă, documentează).

## Workflow
1. Fork repo-ul principal pe contul vostru GitHub.
2. Creați un branch:
``` bash
git checkout -b feat/labNN-handle
```
(unde NN = numărul laboratorului, handle = username GitHub).
3. Adăugați livrabilele:
- cod (scripturi / notebook-uri completate),
- rezultate mici (CSV ≤10 linii),
- *NOTES.md* scurt (≤10 linii).
Livrabilele merg în folderul:
```bash
labs/NN_topic/submissions/<handle>/
```
4. Faceți commit cu un mesaj clar (ex.: "lab02: adaugare aliniere globală").
5. Deschideți un Pull Request spre branch-ul main al cursului.
- Completați template-ul implicit de PR (Săptămâna, livrabile bifate).
- CI trebuie să fie verde (syntax + smoke + mlflow).

### Acceptarea PR
- Un PR este acceptat dacă:
  - respectă structura de fișiere,
  - trece verificările CI,
  - are PR template completat.
- PR-urile se îmbină în branch-ul principal de către asistent.

## Style Guide
- **Python**: urmați PEP8 (CI verifică doar erorile de sintaxă).
- **Notebook-uri**: ștergeți output-urile înainte de commit.
- **Rapoarte**: maxim 2 pagini PDF per assignment.
- **Politica de date** : respectați [GDPR_and_DataPolicy](docs/GDPR_and_DataPolicy.md) si [policies](docs/policies.md)— nu încărcați fișiere mari sau date sensibile.

### Nota
- Toate livrabilele vor fi evaluate la **finalul semestrului**.
- Datele colectate în primele săptămâni (GEO/TCGA/NCBI/Ensembl) vor fi reutilizate în toate laboratoarele ulterioare (aliniere, NGS, filogenie, co-expresie, ML etc.).



