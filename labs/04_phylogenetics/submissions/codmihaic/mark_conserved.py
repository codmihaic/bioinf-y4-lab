# mark_conserved.py
# Usage:
#   python labs/04_phylogenetics/submissions/codmihaic/mark_conserved.py

import argparse
from pathlib import Path
from collections import Counter
from Bio import AlignIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv
import textwrap

GAP_CHARS = set("-.")

def consensus_symbols(col_chars, thr_strong=0.80, thr_mod=0.60):
    chars = [c.upper() for c in col_chars if c not in GAP_CHARS]
    if not chars:
        return (' ', 0.0, '-')
    counts = Counter(chars)
    maj_base, maj_cnt = counts.most_common(1)[0]
    n = len(chars)
    frac = maj_cnt / n
    if frac == 1.0 and n > 0:
        return ('*', frac, maj_base)
    elif frac >= thr_strong:
        return (':', frac, maj_base)
    elif frac >= thr_mod:
        return ('.', frac, maj_base)
    else:
        return (' ', frac, maj_base)

def find_conserved_runs(symbols):
    runs = []
    i = 0
    while i < len(symbols):
        if symbols[i] == '*':
            j = i
            while j < len(symbols) and symbols[j] == '*':
                j += 1
            runs.append((i+1, j))  # 1-based inclusive
            i = j
        else:
            i += 1
    return runs

def write_alignment_with_consensus(aln, symbols, out_path, block_width=60):
    names = [rec.id for rec in aln]
    seqs  = [str(rec.seq) for rec in aln]
    L = aln.get_alignment_length()

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("CLUSTAL W (with consensus)\n\n")
        for start in range(0, L, block_width):
            end = min(L, start+block_width)
            for name, seq in zip(names, seqs):
                f.write(f"{name:<15} {seq[start:end]}\n")
            f.write(" " * 15 + "".join(symbols[start:end]) + "\n\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", default="labs/04_phylogenetics/submissions/codmihaic/msa/clustalo.aln-clustal_num")
    ap.add_argument("--outdir", default="labs/04_phylogenetics/submissions/codmihaic/msa")
    ap.add_argument("--format", default="clustal", help="clustal / fasta etc.")
    args = ap.parse_args()

    in_path = Path(args.in_path)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    aln = AlignIO.read(str(in_path), args.format)
    L = aln.get_alignment_length()
    print(f"[INFO] alignment: {len(aln)} secvențe, lungime {L} coloane")

    symbols = []
    percol = []
    for pos in range(L):
        col = aln[:, pos]
        sym, frac, maj = consensus_symbols(col)
        symbols.append(sym)
        percol.append((pos+1, frac, maj, sym))


    runs = find_conserved_runs(symbols)
    out_clw = outdir / "alignment_with_consensus.clw"
    write_alignment_with_consensus(aln, symbols, out_clw)
    print(f"[OK] {out_clw}")
    out_tsv = outdir / "conserved_regions.tsv"

    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["start", "end", "length"])
        for s, e in runs:
            w.writerow([s, e, e - s + 1])
    print(f"[OK] {out_tsv} (intervale 100% conservate)")

    out_csv = outdir / "conservation_per_column.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["position_1based", "fraction_major_base", "major_base", "symbol"])
        w.writerows(percol)
    print(f"[OK] {out_csv}")

    xs = [p for p, _, _, _ in percol]
    ys = [frac for _, frac, _, _ in percol]
    plt.figure(figsize=(10, 3), dpi=150)
    plt.plot(xs, ys)
    plt.xlabel("Poziție în aliniere (1-based)")
    plt.ylabel("Fracție identică (fără gap-uri)")
    plt.title("Conservarea pe coloană")
    plt.tight_layout()
    out_png = outdir / "conservation_plot.png"
    plt.savefig(out_png)
    plt.close()
    print(f"[OK] {out_png}")

    preview = "\n".join([f"  {s}-{e} (len={e-s+1})" for s, e in runs[:5]])
    if runs:
        print("[SUMAR] Primele intervale 100% conservate:\n" + preview)
    else:
        print("[SUMAR] Nu s-au găsit coloane 100% conservate (poate există multe mutații/gap-uri).")

if __name__ == "__main__":
    main()
