# ex04_fastq_stats.py
# Comanda:
#   python labs/03_formats\&NGS/submissions/codmihaic/ex04_fastq_stats.py

import sys, gzip
from pathlib import Path

def phred_sanger(qual_line: str):
    return [ord(c) - 33 for c in qual_line.strip()]

def main():
    handle = "codmihaic" #handle
    in_fastq = Path(f"data/work/{handle}/lab03/your_reads.fastq.gz")
    out_report = Path(f"labs/03_formats&NGS/submissions/{handle}/qc_report_{handle}.txt")
    out_report.parent.mkdir(parents=True, exist_ok=True)
    if not in_fastq.exists():
        print(f"Error! {in_fastq} doesn't exist!")
        return

    total_reads = total_bases = total_N = total_phred = 0
    with gzip.open(in_fastq, "rt", encoding="utf-8", errors="replace") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            _ = fh.readline()  # plus
            qual = fh.readline().strip()
            if not seq or not qual:
                break

            total_reads += 1
            L = len(seq)
            q = phred_sanger(qual)
            if len(q) != L:
                Lmin = min(L, len(q))
                seq, q, L = seq[:Lmin], q[:Lmin], Lmin
            total_bases += L
            total_N += seq.count('N') + seq.count('n')
            total_phred += sum(q)

    avg_len = (total_bases / total_reads) if total_reads else 0.0
    prop_N = (total_N / total_bases) if total_bases else 0.0
    avg_phred = (total_phred / total_bases) if total_bases else 0.0
    
    lines = [
        f"FASTQ QC report for: {in_fastq}",
        f"Total reads: {total_reads}",
        f"Average read length: {avg_len:.2f}",
        f"Proportion of 'N' bases: {prop_N:.6f}",
        f"Average Phred score (per base): {avg_phred:.2f}",
        "",
        "Assumptions: Sanger FASTQ (Phred+33).",
    ]
    out_report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Done!")

if __name__ == "__main__":
    main()
