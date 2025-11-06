# make_ten_genes.py
# Usage:
#   python labs/04_phylogenics/submissions/codmihaic/make_genes.py --k 10 --len 800 --seed 42

import argparse, random
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def windows(seq, win):
    for start in range(0, len(seq) - win + 1, win):
        yield start, seq[start:start+win]

def clean_acgt(s):
    return ''.join(ch if ch in 'ACGT' else 'A' for ch in s.upper())

def mutate_light(s, rate=0.01, rng=None):
    if rate <= 0: return s
    rng = rng or random
    bases = ('A','C','G','T')
    out = []
    for ch in s:
        if rng.random() < rate:
            out.append(rng.choice([b for b in bases if b != ch]))
        else:
            out.append(ch)
    return ''.join(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in',  dest='in_path', default="data/work/codmihaic/lab04/your_sequences.fasta")
    ap.add_argument('--out', dest='out_path', default="data/work/codmihaic/lab04/assign_seq.fasta")
    ap.add_argument('--k',   type=int, default=10, help='număr secvențe')
    ap.add_argument('--len', dest='win_len', type=int, default=800, help='lungimea fiecărei secvențe')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--mut',  type=float, default=0.01, help='rata de mutații ușoare (0-0.05)')
    args = ap.parse_args()

    random.seed(args.seed)
    in_path  = Path(args.in_path)
    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(str(in_path), 'fasta'))
    if not records:
        raise SystemExit("Nu am găsit nicio secvență în fișierul sursă.")

    pool = []
    for rec in records:
        for start, chunk in windows(str(rec.seq), args.win_len):
            chunk = clean_acgt(chunk)
            pool.append((rec.id, start, chunk))
    if not pool:
        raise SystemExit("Nu există ferestre suficiente de lungimea cerută; micșorează --len.")
    picks = [random.choice(pool) for _ in range(args.k)] if len(pool) < args.k else random.sample(pool, args.k)

    out_recs = []
    for i, (src_id, start, chunk) in enumerate(picks, 1):
        chunk_mut = mutate_light(chunk, rate=args.mut, rng=random)
        rec = SeqRecord(
            Seq(chunk_mut),
            id=f"gene_{i:02d}",
            description=f"src={src_id} start={start} len={args.win_len}"
        )
        out_recs.append(rec)
    SeqIO.write(out_recs, str(out_path), 'fasta')
    print(f"[OK] Scris: {out_path} cu {len(out_recs)} secvențe de {args.win_len} bp.")

if __name__ == '__main__':
    main()
