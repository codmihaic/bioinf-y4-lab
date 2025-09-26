"""
Demo 1 — Distance Matrix
Calculate pairwise distances (Hamming / p-distance) from a multi-FASTA file.
"""

from Bio import SeqIO
import numpy as np

def hamming_distance(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))

if __name__ == "__main__":
    fasta = "data/sample/tp53_dna_multi.fasta"
    records = list(SeqIO.parse(fasta, "fasta"))
    n = len(records)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            d = hamming_distance(str(records[i].seq), str(records[j].seq))
            p_dist = d / len(records[i].seq)
            matrix[i, j] = matrix[j, i] = p_dist

    print("Sequences:", [rec.id for rec in records])
    print("Distance matrix:\n", matrix)
