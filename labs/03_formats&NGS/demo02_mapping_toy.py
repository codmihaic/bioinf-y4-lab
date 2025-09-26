"""
Demo 02 — Toy Mapping

Verificăm dacă citirile (reads) se aliniază exact pe o secvență de referință prin potrivire de șiruri.
Acest exemplu NU este un mapper real (nu suportă gap-uri, scoruri etc.).
"""

reference = "ATGCTAGCTAGGCTAATCGGATCGATCGTACGATCG"
reads = [
    "ATGCTAGC",   # match la început
    "GATCGATC",   # match la mijloc
    "TACGATCG",   # match la sfârșit
    "GGGGGGGG"    # nu există match
]

for read in reads:
    pos = reference.find(read)
    if pos != -1:
        print(f"Read {read} mapped at position {pos}")
    else:
        print(f"Read {read} did not map")
