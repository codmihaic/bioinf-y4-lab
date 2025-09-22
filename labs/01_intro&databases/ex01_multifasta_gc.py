"""
Scop: parcurgeți un fișier multi-FASTA și calculați GC pentru fiecare secvență.

Instrucțiuni:
1) Citiți dintr-un fișier FASTA (ideal: data/sample/mitochondrial_sequences.fasta).
2) Pentru fiecare înregistrare, calculați fracția GC și afișați <id>, <GC>.
3) (Bonus) sortați după GC descrescător.

Completați TODO-urile de mai jos consultând documentația Biopython (SeqIO).
"""

from pathlib import Path
from Bio import SeqIO

# TODO: setați calea către un fișier FASTA mic din repo (data/sample/...)
fasta_path = Path("data/sample/mitochondrial_sequences.fasta")  # exemplu sugerat

# TODO: verificați că fișierul există; dacă nu, afișați un mesaj prietenos și ieșiți
if not fasta_path.exists():
    print("Fișierul FASTA nu există:", fasta_path)
    print("Sugestie: creați un fișier mic în data/sample/ sau folosiți un fișier disponibil la curs.")
    raise SystemExit(1)

# TODO: parcurgeți înregistrările din FASTA (SeqIO.parse)
records = list(SeqIO.parse(str(fasta_path), "fasta"))

# TODO: pentru fiecare înregistrare, calculați fracția GC și afișați rezultatul
def gc_fraction(seq):
    return (seq.count("G") + seq.count("C")) / max(1, len(seq))

results = []
for rec in records:
    gc = gc_fraction(rec.seq)
    results.append((rec.id, gc))

# TODO (bonus): sortați după GC descrescător
results.sort(key=lambda x: x[1], reverse=True)

# Afișare finală
for rid, gc in results:
    print(f"{rid}\tGC={gc:.3f}")
