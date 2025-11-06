"""
Task02-a4.py
Comanda:
    
"""

import sys, matplotlib.pyplot as plt
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def detect_seq_type(msa):
    """Întoarce 'identity' (ADN/ARN) sau 'blosum62' (proteine) pentru DistanceCalculator."""
    dna_letters = set("ACGTUN-")
    all_letters = set("".join(str(rec.seq).upper() for rec in msa))
    proteinish = any(ch.isalpha() and ch not in dna_letters for ch in all_letters)
    return "blosum62" if proteinish else "identity"

if __name__ == "__main__":
   fasta = Path("data/work/codmihaic/lab04/assign_seq.fasta")
   if not fasta.exists():
      sys.exit(f"[EROARE] Nu găsesc {fasta}. Creează un multi-FASTA cu ≥3 secvențe.")

   alignment = AlignIO.read(str(fasta), "fasta")
   lengths = {len(rec.seq) for rec in alignment}
   if len(lengths) != 1:
      print("[AVERTISMENT] Secvențele nu au aceeași lungime; pare că FASTA-ul nu e aliniat.")
      print("            Vă rog furnizați un multi-FASTA ALINIAT (cu '-').")
   print(f"[INFO] {len(alignment)} secvențe încărcate din {fasta.name}")

   model = detect_seq_type(alignment)
   calculator = DistanceCalculator(model)
   dm = calculator.get_distance(alignment)
   print("[INFO] Matricea de distanțe:")
   print(dm)

   constructor = DistanceTreeConstructor()
   nj_tree = constructor.nj(dm)
   out_dir = Path(f"labs/04_phylogenetics/submissions/codmihaic")
   out_dir.mkdir(parents=True, exist_ok=True)
   out_newick = out_dir / f"tree_codmihaic.nwk"
   Phylo.write(nj_tree, str(out_newick), "newick")
   print(f"[OK] Arbore salvat în Newick: {out_newick}")

   print("\n[ASCII tree]")
   Phylo.draw_ascii(nj_tree)
   out_png = out_dir / f"tree_codmihaic.png"
   fig = plt.figure(figsize=(6, 6), dpi=150)
   ax = fig.add_subplot(1, 1, 1)
   Phylo.draw(nj_tree, do_show=False, axes=ax)
   plt.tight_layout()
   fig.savefig(out_png)
   print(f"[OK] PNG salvat: {out_png}")