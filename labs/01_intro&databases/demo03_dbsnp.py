# Interogare rapidă dbSNP pentru BRCA1 (Python only)
from Bio import Entrez
Entrez.email = "student@example.com"

res = Entrez.esearch(db="snp", term="BRCA1[gene] AND Homo sapiens[organism]", retmax=5)
ids = Entrez.read(res)["IdList"]
print(f"Am găsit {len(ids)} SNP IDs.")
if ids:
    summ = Entrez.esummary(db="snp", id=",".join(ids), retmode="xml")
    docsums = Entrez.read(summ)
    for d in docsums:
        print("SNP_ID:", d.get("SNP_ID"), "| DOCSUM:", d.get("DOCSUM"))
