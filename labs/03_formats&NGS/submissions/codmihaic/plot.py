# plot.py
# Comanda:
#   python labs/03_formats\&NGS/submissions/codmihaic/plot.py 

import os
import gzip
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


handle = "codmihaic" # handle
in_fastq_plain = Path(f"data/work/{handle}/lab03/your_reads.fastq")
in_fastq_gz    = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
if in_fastq_plain.exists():
    in_fastq = in_fastq_plain
elif in_fastq_gz.exists():
    in_fastq = in_fastq_gz
else:
    raise FileNotFoundError("No file found!")

root = Path(os.path.dirname(os.path.abspath(__file__)))
out_png = root / f"qc_plot_{handle}.png"

length_counts = Counter()           
phred_max = 50                      
phred_counts = np.zeros(phred_max + 1, dtype=np.int64)  


def parse_fastq_stream(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            yield from SeqIO.parse(fh, "fastq")
    else:
        with open(path, "rt") as fh:
            yield from SeqIO.parse(fh, "fastq")

total_reads = 0
for rec in parse_fastq_stream(in_fastq):
    L = len(rec)
    if L == 0:
        continue
    length_counts[L] += 1

    q = rec.letter_annotations.get("phred_quality", [])
    if q:
        mean_q = sum(q) / L
        bin_q = int(round(mean_q))
        bin_q = max(0, min(phred_max, bin_q))
        phred_counts[bin_q] += 1
    total_reads += 1

plt.figure(figsize=(11, 4.5))
plt.subplot(1, 2, 1) #lungime
if length_counts:
    xs = np.array(sorted(length_counts.keys()))
    ys = np.array([length_counts[x] for x in xs], dtype=np.int64)
    plt.bar(xs, ys, width=1)  
plt.title("Distribuția lungimilor citirilor")
plt.xlabel("Lungime (baze)")
plt.ylabel("Număr de citiri")

plt.subplot(1, 2, 2) #scoruri medii
bins_q = np.arange(phred_max + 1)
plt.bar(bins_q, phred_counts, width=0.9)
plt.title("Distribuția scorurilor Phred medii per citire")
plt.xlabel("Scor Phred mediu (rotunjit)")
plt.ylabel("Număr de citiri")

plt.tight_layout()
plt.savefig(out_png, dpi=150)
print(f"Done!")
