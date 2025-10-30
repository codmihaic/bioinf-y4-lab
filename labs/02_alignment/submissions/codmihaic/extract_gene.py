# extract_subset.py
# Usage:
#   python labs/02_alignment/submissions/codmihaic/extract_gene.py

from Bio import SeqIO
from pathlib import Path

handle = "codmihaic"
input_path = Path(f"data/work/{handle}/lab01/my_tp53.fa")
output_path = Path(f"data/work/{handle}/lab03/msa_subset.fasta")
num_bases = 1500
if not input_path.exists():
    raise FileNotFoundError(f"Nu există fișierul {input_path}")

output_path.parent.mkdir(parents=True, exist_ok=True)
records_out = []
for record in SeqIO.parse(input_path, "fasta"):
    trimmed_seq = record.seq[:num_bases]
    new_record = record[:num_bases]
    new_record.id = record.id + "_subset"
    new_record.description = f"{record.description} (primele {num_bases} baze)"
    records_out.append(new_record)


SeqIO.write(records_out, output_path, "fasta")
print(f"[OK] S-au extras primele {num_bases} baze din fiecare genă.")
print(f"Fișierul rezultat: {output_path}")
