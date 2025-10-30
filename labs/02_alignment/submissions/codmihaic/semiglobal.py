#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere semi-globală (ends-free)
- Gap-urile de la capete NU se penalizează.
- Intern (în mijlocul secvențelor) gap-urile se penalizează normal.

Rulare:
  python labs/02_alignment/submissions/codmihaic/semiglobal.py --fasta data/work/codmihaic/lab01/my_tp53.fa --i1 0 --i2 1
"""

from pathlib import Path
import argparse
from Bio import SeqIO

def init_score_matrix_semiglobal(m: int, n: int):
    return [[0]*(n+1) for _ in range(m+1)]

def score_cell_semiglobal(score, i: int, j: int, a: str, b: str,
                          match: int, mismatch: int, gap: int):
    diag  = score[i-1][j-1] + (match if a == b else mismatch)
    up    = score[i-1][j]   + gap
    left  = score[i][j-1]   + gap
    return max(diag, up, left)

def semiglobal_align(seq1: str, seq2: str, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    score = init_score_matrix_semiglobal(m, n)
    for i in range(1, m+1):
        ai = seq1[i-1]
        for j in range(1, n+1):
            bj = seq2[j-1]
            score[i][j] = score_cell_semiglobal(score, i, j, ai, bj, match, mismatch, gap)

    i_best, j_best, best = m, n, score[m][n]
    for j in range(n+1):
        if score[m][j] > best:
            i_best, j_best, best = m, j, score[m][j]
    for i in range(m+1):
        if score[i][n] > best:
            i_best, j_best, best = i, n, score[i][n]

    align1, align2 = "", ""
    i, j = i_best, j_best
    while i > 0 and j > 0:
        current = score[i][j]
        diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
        up   = score[i-1][j] + gap
        left = score[i][j-1] + gap

        if current == diag:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1; j -= 1
        elif current == up:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        if i == 0 or j == 0:
            break

    return align1, align2, best

def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[eroare] Fișierul trebuie să conțină cel puțin 2 secvențe.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[eroare] Indici invalizi (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Cale către FASTA-ul propriu din data/work/<handle>/lab01/")
    ap.add_argument("--i1", type=int, default=0, help="Index prima secvență (implicit 0)")
    ap.add_argument("--i2", type=int, default=1, help="Index a doua secvență (implicit 1)")
    ap.add_argument("--limit", type=int, default=500, help="Trunchiază fiecare secvență la primele L baze (default 500)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"File not found: {fasta_path}!!")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)

    # Trunchiere pentru viteză/memorie (datele tale sunt uriașe)
    L = min(len(s1), len(s2), args.limit)
    s1, s2 = s1[:L], s2[:L]

    a1, a2, sc = semiglobal_align(s1, s2, match=1, mismatch=-1, gap=-2)

    print("=== Aliniere semi-globală (ends-free) ===")
    print(f"{id1}  vs  {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)

if __name__ == "__main__":
    main()
