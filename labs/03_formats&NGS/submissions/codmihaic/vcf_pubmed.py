# vcf_pubmed.py
# Comanda:
#   python labs/03_formats\&NGS/submissions/codmihaic/vcf_pubmed.py --handle codmihaic --k 2 --auto-download

from pathlib import Path
import argparse, gzip, io
from Bio import Entrez


def ensure_vcf_exists(folder: Path, filename: str = "sample.vcf") -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / filename
    urls = [
        # small public VCFs for tests 
        "https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/VCF/example.vcf",
        "https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/tests/data/variants/clinvar.vcf",
    ]

    try:
        import requests
        for url in urls:
            r = requests.get(url, timeout=15)
            if r.ok and len(r.text) > 200:
                dest.write_text(r.text, encoding="utf-8")
                return dest
    except Exception:
        pass

    toy_vcf = """##fileformat=VCFv4.2
        ##source=toy
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        17\t7676158\trs1042522\tC\tT\t.\tPASS\tGENE=TP53
        17\t7673803\t.\tG\tA\t.\tPASS\tGENE=TP53 """
    dest.write_text(toy_vcf, encoding="utf-8")
    return dest


def read_variants(vcf_path: Path, max_variants: int = 2):
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    variants = []
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            chrom, pos, vid = fields[0], fields[1], fields[2]
            variants.append({"CHROM": chrom, "POS": pos, "ID": vid})
            if len(variants) >= max_variants:
                break
    return variants


def pubmed_search(term: str, retmax: int = 3):
    Entrez.email = "codreanumihaic@gmail.com"  # email
    results = []
    with Entrez.esearch(db="pubmed", term=term, retmax=retmax) as r:
        ids = Entrez.read(r).get("IdList", [])
    if not ids:
        return results
    with Entrez.esummary(db="pubmed", id=",".join(ids), retmode="xml") as r:
        docsums = Entrez.read(r)

    for ds in docsums:
        results.append({
            "pmid": ds.get("Id", ""),
            "title": ds.get("Title", "").strip()
        })
    return results


def main():
    ap = argparse.ArgumentParser(description="VCF to PubMed")
    ap.add_argument("--vcf", help="Path to VCF (.vcf or .vcf.gz)")
    ap.add_argument("--k", type=int, default=2, help="processed samples (default 2)")
    ap.add_argument("--handle", default="student", help="Handle for github name")
    ap.add_argument("--auto-download", action="store_true", help="if no --vcf given, download/generate a VCF in data/work/<handle>/lab03/")
    args = ap.parse_args()

    if args.vcf:
        vcf_path = Path(args.vcf)
        if not vcf_path.exists():
            raise SystemExit(f"Eroare: nu găsesc fișierul VCF: {vcf_path}")
    else:
        if not args.auto_download:
            raise SystemExit("No --vcf. use --auto-download for a new VCF.")
        folder = Path(f"data/work/{args.handle}/lab03")
        vcf_path = ensure_vcf_exists(folder, filename="sample.vcf")
        print(f"VCF done in: {vcf_path}")

    variants = read_variants(vcf_path, max_variants=args.k)
    if not variants:
        raise SystemExit("No samples in the VCF.")

    out_path = Path(f"variants_{args.handle}.txt")
    lines = []
    lines.append(f"Input VCF: {vcf_path}")
    lines.append(f"Processed: {len(variants)}")
    lines.append("")

    for i, v in enumerate(variants, start=1):
        chrom, pos, vid = v["CHROM"], v["POS"], v["ID"]
        term = vid if (vid and vid != "." and vid.lower().startswith("rs")) else f"chr{chrom}:{pos} AND TP53"
        hits = pubmed_search(term, retmax=3)

        lines.append(f"=== Sample {i} ===")
        lines.append(f"CHROM: {chrom}  POS: {pos}  ID: {vid}")
        lines.append(f"PubMed: {term}")
        if not hits:
            lines.append("Results: (nothing found)")
        else:
            for j, h in enumerate(hits, start=1):
                lines.append(f"- [{j}] PMID {h['pmid']}: {h['title']}")
        lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Done!")

if __name__ == "__main__":
    main()
