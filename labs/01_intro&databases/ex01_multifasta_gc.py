"""
Scop: parcurgeți un fișier multi-FASTA și calculați GC pentru fiecare secvență.

Instrucțiuni:
1) Citiți dintr-un fișier FASTA (ideal: data/sample/mitochondrial_sequences.fasta).
2) Pentru fiecare înregistrare, calculați fracția GC și afișați <id>\t<GC>.
3) (Bonus) sortați descrescător după GC.

Completați TODO-urile consultând documentația Biopython (SeqIO).
"""

from pathlib import Path
from Bio import SeqIO


def gc_fraction(seq) -> float:
    """Fracție GC pentru o secvență; protejat la împărțire la zero."""
    length = max(1, len(seq))
    return (seq.count("G") + seq.count("C")) / length


def main() -> None:
    # TODO: setați calea către un FASTA mic din repo (data/sample/...)
    fasta_path = Path("data/sample/mitochondrial_sequences.fasta")

    # TODO: verificați că fișierul există; altfel, mesaj prietenos și exit
    if not fasta_path.exists():
        print("Fișierul FASTA nu există:", fasta_path)
        print(
            "Sugestie: creați un fișier mic în data/sample/ sau "
            "folosiți un fișier oferit la curs."
        )
        raise SystemExit(1)

    # TODO: parcurgeți înregistrările din FASTA (SeqIO.parse)
    records = list(SeqIO.parse(str(fasta_path), "fasta"))

    # TODO: calculați GC pentru fiecare înregistrare
    results = []
    for rec in records:
        results.append((rec.id, gc_fraction(rec.seq)))

    # TODO (bonus): sortați descrescător după GC
    results.sort(key=lambda x: x[1], reverse=True)

    # Afișare finală
    for rec_id, gc in results:
        print(f"{rec_id}\tGC={gc:.3f}")


if __name__ == "__main__":
    main()
