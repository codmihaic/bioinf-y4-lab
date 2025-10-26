# fasta_copy.py
# Comanda:
#   python labs/03_formats\&NGS/submissions/codmihaic/copy.py

from pathlib import Path
import shutil

def main():
    handle = "codmihaic"
    src = Path(f"data/work/{handle}/lab01/my_tp53.fa")
    dst = Path(f"fasta_{handle}.fasta")

    if not src.exists():
        print(f"Eroare: fișierul sursă nu există -> {src}")
        return
    shutil.copy(src, dst)
    print(f"Fișierul a fost copiat cu succes în: {dst}")

if __name__ == "__main__":
    main()