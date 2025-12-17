"""
Exercise 9.2 — Disease Proximity and Drug Ranking

Scop:
- să calculați distanța medie dintre fiecare medicament și un set de gene asociate unei boli
- să ordonați medicamentele în funcție de proximitate (network-based prioritization)

TODO-uri principale:
- încărcați graful bipartit drug–gene (din exercițiul 9.1) sau reconstruiți-l
- încărcați lista de disease genes
- pentru fiecare medicament, calculați distanța minimă / medie până la genele bolii
- exportați un tabel cu medicamente și scorul lor de proximitate
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple

import networkx as nx
import pandas as pd
import pickle

# --------------------------
# Config
# --------------------------
HANDLE = "codmihaic"

# Input: graful bipartit (salvat anterior) SAU tabelul drug-gene
GRAPH_DRUG_GENE = Path(f"labs/09_repurposing/submissions/{HANDLE}/network_drug_gene_{HANDLE}.gpickle")
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Input: lista genelor bolii
DISEASE_GENES_TXT = Path(f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")

# Output directory & file
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"
OUT_REPURPOSING = OUT_DIR / f"REPURPOSING_{HANDLE}.csv"



# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError
    """
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] Missing required file: {path}")

def load_drug_gene_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.lower().strip() for c in df.columns]
    if "drug" not in df.columns or "gene" not in df.columns:
        raise ValueError(f"[ERROR] CSV must contain columns drug,gene. Found: {list(df.columns)}")

    df["drug"] = df["drug"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[(df["drug"] != "") & (df["gene"] != "")]
    df = df.drop_duplicates(subset=["drug", "gene"]).reset_index(drop=True)
    if df.empty:
        raise ValueError("[ERROR] Drug-gene table is empty after cleaning.")
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    return df.groupby("drug")["gene"].apply(lambda s: set(s.dropna().astype(str))).to_dict()

def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    B = nx.Graph()
    drugs = list(drug2genes.keys())
    genes = sorted({g for gs in drug2genes.values() for g in gs})

    B.add_nodes_from(drugs, bipartite="drug", node_type="drug")
    B.add_nodes_from(genes, bipartite="gene", node_type="gene")

    for d, gs in drug2genes.items():
        for g in gs:
            B.add_edge(d, g)
    return B


def load_bipartite_graph_or_build() -> nx.Graph:
    """
    TODO:
    - dacă GRAPH_DRUG_GENE există, încărcați-l direct
    - altfel, reconstruiți graful plecând de la DRUG_GENE_CSV
      (puteți reutiliza logica din ex09_drug_similarity_network.py)
    """
    if GRAPH_DRUG_GENE.exists():
        with open(GRAPH_DRUG_GENE, "rb") as f:
            B = pickle.load(f)
        if not isinstance(B, nx.Graph):
            raise TypeError("[ERROR] Loaded object is not a NetworkX Graph.")
        return B

    # 2) altfel reconstruiește din CSV
    ensure_exists(DRUG_GENE_CSV)
    df = load_drug_gene_table(DRUG_GENE_CSV)
    drug2genes = build_drug2genes(df)
    return build_bipartite_graph(drug2genes)


def load_disease_genes(path: Path) -> Set[str]:
    """
    TODO:
    - încărcați fișierul text cu gene (una pe linie)
    - returnați un set de gene (string)
    """
    ensure_exists(path)
    genes: Set[str] = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.add(g)
    if not genes:
        raise ValueError("[ERROR] Disease genes file is empty.")
    return genes


def get_drug_nodes(B: nx.Graph) -> List[str]:
    """
    TODO:
    - extrageți lista nodurilor de tip 'drug'
    - presupunem atributul bipartite="drug"
    """
    return [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 5,
) -> float:
    """
    TODO:
    - pentru un medicament:
      - calculați distanța (de ex. shortest_path_length) până la fiecare genă din disease_genes
      - ignorați genele care nu sunt în graf
      - dacă nu există niciun drum, puteți seta o distanță penalizantă (ex. max_dist + 1)
    - returnați media (sau minimul) distanțelor; controlați cu parametrul 'mode'
    """
    genes_in_graph = [g for g in disease_genes if g in B]
    if not genes_in_graph:
        return float(max_dist + 1)

    dists: List[int] = []
    for g in genes_in_graph:
        try:
            d = nx.shortest_path_length(B, source=drug, target=g)
            dists.append(min(d, max_dist + 1))
        except nx.NetworkXNoPath:
            dists.append(max_dist + 1)

    if mode == "min":
        return float(min(dists))
    return float(sum(dists) / len(dists))



def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    """
    TODO:
    - pentru fiecare medicament din graf:
      - calculați scorul de distanță (ex. media distanțelor către genele bolii)
    - construiți un DataFrame cu:
      drug, distance
    - sortați crescător după distance (distanță mai mică = proximitate mai mare)
    """
    drugs = get_drug_nodes(B)
    genes_in_graph = sorted([g for g in disease_genes if g in B])
    coverage = len(genes_in_graph) / max(1, len(disease_genes))
    rows = []

    for d in drugs:
        dist = compute_drug_disease_distance(B, d, disease_genes, mode=mode, max_dist=5)
        rows.append({"drug": d, "distance": dist})
    out = pd.DataFrame(rows).sort_values(["distance", "drug"], ascending=[True, True]).reset_index(drop=True)

    out["mode"] = mode
    out["disease_genes_total"] = len(disease_genes)
    out["disease_genes_in_graph"] = len(genes_in_graph)
    out["disease_gene_coverage"] = coverage
    return out


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați input-urile
    ensure_exists(DISEASE_GENES_TXT)

    # TODO 2: încărcați / construiți graful bipartit
    B = load_bipartite_graph_or_build()

    # TODO 3: încărcați setul de disease genes
    disease_genes = load_disease_genes(DISEASE_GENES_TXT)

    # TODO 4: calculați ranking-ul medicamentelor după proximitate
    df_rank = rank_drugs_by_proximity(B, disease_genes, mode="mean")

    # TODO 5: salvați rezultatele
    df_rank.to_csv(OUT_DRUG_PRIORITY, index=False)
    df_rank.to_csv(OUT_REPURPOSING, index=False)

    print(f"[OUT]  {OUT_REPURPOSING}")
    print("[INFO] Done.")
    print(f"[INFO] Drugs in graph: {len(get_drug_nodes(B))}")
    print(f"[INFO] Disease genes: {len(disease_genes)} | Present in graph: {df_rank.loc[0,'disease_genes_in_graph']}")
    print(f"[OUT]  {OUT_DRUG_PRIORITY}")
    print("[INFO] Top 10 drugs:")
    print(df_rank.head(10).to_string(index=False))
