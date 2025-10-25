# pubmed_query.py
# Comanda:
#   python labs/03_formats\&NGS/submissions/codmihaic/pubmed_query.py 

from Bio import Entrez
from pathlib import Path

def main():
    handle = "codmihaic" # handle name
    out_file = Path(f"pubmed_{handle}.txt")

    Entrez.email = "codreanumihaic@gmail.com" #mail
    query = "TP53 AND cancer"
    with Entrez.esearch(db="pubmed", term=query, retmax=5) as r:
        ids = Entrez.read(r)["IdList"]

    if not ids:
        print("Articles weren't found!")
        return
    with Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="xml") as r:
        records = Entrez.read(r)["PubmedArticle"]
        
    lines = []
    for i, rec in enumerate(records, start=1):
        art = rec["MedlineCitation"]["Article"]
        title = str(art.get("ArticleTitle", "")).strip()
        authors = []
        for a in art.get("AuthorList", []):
            last = a.get("LastName", "")
            initials = a.get("Initials", "")
            if last or initials:
                authors.append(f"{last} {initials}".strip())
        authors_str = ", ".join(authors) if authors else "N/A"

        abstract_texts = art.get("Abstract", {}).get("AbstractText", [])
        if isinstance(abstract_texts, list):
            abstract = " ".join(str(t) for t in abstract_texts).strip()
        else:
            abstract = str(abstract_texts).strip()
        if not abstract:
            abstract = "N/A"

        block = [
            f"=== Articol {i} ===",
            f"Titlu: {title}",
            f"Autori: {authors_str}",
            "Rezumat:",
            abstract,
            ""
        ]
        lines.extend(block)

    out_file.write_text("\n".join(lines), encoding="utf-8")
    print(f"Done!")

if __name__ == "__main__":
    main()