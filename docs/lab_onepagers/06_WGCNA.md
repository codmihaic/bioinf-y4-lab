# One-Pager: Lab 6 — Gene Co-Expression Networks (Building Modules)

## Obiective
- Construirea unei rețele de co-expresie genică pe baza datelor RNA-Seq.
- Preprocesarea și normalizarea datelor de expresie.
- Calcularea matricilor de corelație și adiacență.
- Detectarea modulelor funcționale prin algoritmi de tip community detection (ex: Louvain).
- Înțelegerea limitărilor și capcanelor interpretative.

---

## Competențe
- Încărcarea și filtrarea datelor RNA-Seq.
- Aplicarea transformărilor log și a filtrării pe bază de varianță.
- Construirea rețelelor de co-expresie cu NetworkX.
- Detectarea modulelor și exportul mapping-ului gene → module.
- Evaluarea critică a rezultatelor.

---

## Pași principali
1. Încărcați dataset-ul RNA-Seq (GSE115469).
2. Normalizați și log-transformați datele.
3. Filtrați genele cu varianță scăzută.
4. Calculați corelația (Pearson/Spearman).
5. Aplicați un prag pentru a construi matricea de adiacență.
6. Construiți graful și detectați module cu Louvain.

---

## Resurse
- [Fișa laborator](../../docs/lab_onepagers/06_WGCNA.md)
- [NetworkX Documentation](https://networkx.org/documentation/stable/)
- [GEO Accession GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469)
- [van Dam et al., Brief Bioinform 2018](https://doi.org/10.1093/bib/bbw139)
- [Langfelder & Horvath, BMC Bioinformatics 2008](https://doi.org/10.1186/1471-2105-9-559)
