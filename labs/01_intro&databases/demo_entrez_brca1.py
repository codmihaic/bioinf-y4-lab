from pathlib import Path
from Bio import Entrez, SeqIO

# NCBI cere un email; treceti email-ul real
Entrez.email = "student@example.com"
# opțional: Entrez.api_key = "NCBI_API_KEY"

Path("data").mkdir(exist_ok=True)

search = Entrez.esearch(db="nucleotide",
                        term="BRCA1[Gene] AND Homo sapiens[Organism]",
                        retmax=1)
ids = Entrez.read(search)["IdList"]
if not ids:
    raise SystemExit("No BRCA1 nucleotide records found")

acc = ids[0]
gb = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
out = Path("data/brca1.gb")
out.write_text(gb.read(), encoding="utf-8")

rec = SeqIO.read(out, "genbank")
gc_fraction = (rec.seq.count("G") + rec.seq.count("C")) / max(1, len(rec.seq))

print("Accession:", rec.id)
print("Length:", len(rec.seq))
print("GC fraction:", round(gc_fraction, 3))
print("First 50 nt:", rec.seq[:50])
