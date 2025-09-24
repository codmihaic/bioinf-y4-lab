#!/usr/bin/env python
"""
Calculează p-distance/Hamming pentru toate perechile dintr-un FASTA.
- Refolosim data/sample/* ; pentru lungimi diferite trunchiem la min(len).
Rulare:
  python labs/02_alignment/demo02_distance_matrix.py --fasta data/sample/tp53_dna_multi.fasta
"""
import argparse
from itertools import combinations
from Bio import SeqIO

def hamming_equal(a, b):
    return sum(x != y for x, y in zip(a, b))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    args = ap.parse_args()

    recs = list(SeqIO.parse(args.fasta, "fasta"))
    ids = [r.id for r in recs]
    seqs = [str(r.seq) for r in recs]

    print("pair,hamming,p_distance,len_used")
    for (i, j) in combinations(range(len(seqs)), 2):
        a, b = seqs[i], seqs[j]
        L = min(len(a), len(b))
        a2, b2 = a[:L], b[:L]  # trunchiere pentru comparație simplă
        d_h = hamming_equal(a2, b2)
        d_p = d_h / float(L) if L > 0 else 0.0
        print(f"{ids[i]}-{ids[j]},{d_h},{d_p:.4f},{L}")

if __name__ == "__main__":
    main()
