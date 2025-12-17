from __future__ import annotations
from pathlib import Path
import random
import pandas as pd

HANDLE = "codmihaic"  # <-- pune handle-ul tău
OUT = Path(f"data/work/{HANDLE}/lab09")
OUT.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUT / f"drug_gene_{HANDLE}.csv"

random.seed(42)

# --- Gene pool (mix: real + synthetic-like) ---
core_genes = [
    "PTGS1","PTGS2","IL6","TNF","NFKB1","JUN","FOS","MAPK1","MAPK3","STAT3","CXCL8","CRP",
    "HMGCR","LDLR","PCSK9","APOB","SREBF2",
    "VKORC1","CYP2C9","CYP2C19","CYP3A4","ABCB1",
    "P2RY12","F2","F5","SERPINC1",
    "ACE","AGTR1","REN","NOS3",
    "EGFR","PIK3CA","MTOR","AMPK","AKT1",
    "SLC22A1","SLC47A1"
]

# extra "noise" genes
noise_genes = [f"G{idx:03d}" for idx in range(1, 161)]  # G001..G160
all_genes = core_genes + noise_genes

# --- Drug clusters with shared hubs ---
clusters = {
    "NSAID": {
        "drugs": ["Aspirin","Ibuprofen","Diclofenac","Naproxen","Indomethacin","Ketorolac","Meloxicam","Celecoxib"],
        "genes": ["PTGS1","PTGS2","IL6","TNF","NFKB1","MAPK1","MAPK3","CXCL8"]
    },
    "STATIN": {
        "drugs": ["Atorvastatin","Simvastatin","Rosuvastatin","Pravastatin","Lovastatin","Fluvastatin"],
        "genes": ["HMGCR","LDLR","PCSK9","APOB","SREBF2","NFKB1","CRP"]
    },
    "ANTICOAG": {
        "drugs": ["Warfarin","Heparin","Enoxaparin","Apixaban","Rivaroxaban","Dabigatran"],
        "genes": ["VKORC1","F2","F5","SERPINC1","ABCB1","CYP3A4"]
    },
    "ANTIPLATELET": {
        "drugs": ["Clopidogrel","Prasugrel","Ticagrelor","Aspirin_lowdose"],
        "genes": ["P2RY12","CYP2C19","ABCB1","NFKB1","IL6"]
    },
    "RAAS": {
        "drugs": ["Lisinopril","Enalapril","Losartan","Valsartan","Ramipril","Captopril"],
        "genes": ["ACE","AGTR1","REN","NOS3","STAT3","CRP"]
    },
    "METABOLIC": {
        "drugs": ["Metformin","Pioglitazone","Sitagliptin","Empagliflozin","Glibenclamide","Liraglutide"],
        "genes": ["AMPK","MTOR","AKT1","SLC22A1","SLC47A1","IL6","STAT3"]
    },
    "ANTIBIOTIC": {
        "drugs": ["Amoxicillin","Azithromycin","Ciprofloxacin","Doxycycline","Ceftriaxone","Vancomycin"],
        "genes": ["NFKB1","IL6","TNF","CRP","ABCB1"]
    },
    "ONCO": {
        "drugs": ["Erlotinib","Gefitinib","Imatinib","Sorafenib","Sunitinib","Trastuzumab"],
        "genes": ["EGFR","PIK3CA","MTOR","AKT1","STAT3","JUN","FOS"]
    }
}

rows = []

def add_pairs(drug: str, genes: set[str]):
    for g in genes:
        rows.append({"drug": drug, "gene": g})

for cname, c in clusters.items():
    for d in c["drugs"]:
        base = set(c["genes"])

        # each drug gets some extra core genes (cross-talk)
        base |= set(random.sample(core_genes, k=random.randint(1, 3)))

        # add noise: a few random genes (makes it "stufos" + imperfect)
        base |= set(random.sample(noise_genes, k=random.randint(2, 7)))

        add_pairs(d, base)

# Add a few "lonely" drugs (almost disconnected) to make graph more interesting
lonely_drugs = ["VitaminC","Melatonin","Biotin","ProbioticMix"]
for d in lonely_drugs:
    genes = set(random.sample(noise_genes, k=random.randint(3, 6)))
    add_pairs(d, genes)

# Intentionally add some duplicates (common in real noisy CSVs)
rows += random.sample(rows, k=30)

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)
print(f"[OK] Wrote {len(df)} rows -> {OUT_CSV}")
print(df.head(10).to_string(index=False))
