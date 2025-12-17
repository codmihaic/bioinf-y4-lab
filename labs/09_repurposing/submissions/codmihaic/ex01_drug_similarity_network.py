"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network

Scop:
- să construiți o rețea bipartită drug–gene plecând de la un CSV
- să proiectați layer-ul de medicamente folosind similaritatea dintre seturile de gene
- să exportați un fișier cu muchiile de similaritate între medicamente

TODO:
- încărcați datele drug–gene
- construiți dict-ul drug -> set de gene țintă
- construiți graful bipartit drug–gene (NetworkX)
- calculați similaritatea dintre medicamente (ex. Jaccard)
- construiți graful drug similarity
- exportați tabelul cu muchii: drug1, drug2, weight
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List

import itertools
import networkx as nx
import pickle
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------
# Config — adaptați pentru handle-ul vostru
# --------------------------
HANDLE = "codmihaic"

# Input: fișier cu coloane cel puțin: drug, gene
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Output directory & files
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"
OUT_GRAPH_DRUG_SIM = OUT_DIR / f"network_drug_similarity_{HANDLE}.gpickle"
OUT_GRAPH_DRUG_GENE_PNG = OUT_DIR / f"network_drug_gene_{HANDLE}.png"


def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError cu un mesaj clar
    """
    if not path.exists():
        raise FileNotFoundError(
            f"[ERROR] Input file not found: {path}\n"
            f"Check that you placed drug_gene_{HANDLE}.csv under: data/work/{HANDLE}/lab09/"
        )
    if not path.is_file():
        raise FileNotFoundError(f"[ERROR] Path exists but is not a file: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - validați că există cel puțin coloanele: 'drug', 'gene'
    - returnați DataFrame-ul
    """
    df = pd.read_csv(path)
    required = {"drug", "gene"}
    missing = required - set(df.columns.str.lower())
    df.columns = [c.lower().strip() for c in df.columns]

    if missing:
        raise ValueError(
            f"[ERROR] Missing required columns in CSV: {sorted(missing)}\n"
            f"Found columns: {list(df.columns)}"
        )

    df["drug"] = df["drug"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[(df["drug"] != "") & (df["gene"] != "")]
    df = df.drop_duplicates(subset=["drug", "gene"]).reset_index(drop=True)

    if df.empty:
        raise ValueError("[ERROR] After cleaning, the drug-gene table is empty.")
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    TODO:
    - construiți un dict: drug -> set de gene țintă
    - sugestie: folosiți groupby("drug") și aplicați set() pe coloana gene
    """
    drug2genes = (
        df.groupby("drug")["gene"]
        .apply(lambda s: set(s.dropna().astype(str)))
        .to_dict()
    )
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    TODO:
    - construiți graful bipartit:
      - nodurile 'drug' cu atribut bipartite="drug"
      - nodurile 'gene' cu atribut bipartite="gene"
      - muchii drug-gene
    """
    G = nx.Graph()
    drugs = list(drug2genes.keys())
    genes = sorted({g for gs in drug2genes.values() for g in gs})
    G.add_nodes_from(drugs, bipartite="drug", node_type="drug")
    G.add_nodes_from(genes, bipartite="gene", node_type="gene")

    for d, gs in drug2genes.items():
        for g in gs:
            G.add_edge(d, g)
    return G


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    TODO:
    - construiți un DataFrame cu:
        drug, num_targets (numărul de gene țintă)
    - returnați DataFrame-ul
    """
    rows = [{"drug": d, "num_targets": len(gs)} for d, gs in drug2genes.items()]
    out = pd.DataFrame(rows).sort_values(["num_targets", "drug"], ascending=[False, True]).reset_index(drop=True)
    return out


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """
    Calculați similaritatea Jaccard între două seturi de gene:
    J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    """
    TODO:
    - pentru toate perechile de medicamente (combinații de câte 2),
      calculați similaritatea Jaccard între seturile de gene
    - rețineți doar muchiile cu similaritate >= min_sim
    - returnați o listă de tuple (drug1, drug2, weight)
    """
    drugs = sorted(drug2genes.keys())
    edges: List[Tuple[str, str, float]] = []

    for d1, d2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[d1], drug2genes[d2])
        if sim >= min_sim:
            edges.append((d1, d2, float(sim)))

    edges.sort(key=lambda x: x[2], reverse=True)
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """
    TODO:
    - transformați lista de muchii (drug1, drug2, weight) într-un DataFrame
      cu coloanele: drug1, drug2, similarity
    """
    return pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])

def build_drug_similarity_graph(edges: List[Tuple[str, str, float]]) -> nx.Graph:
    G = nx.Graph()
    for d1, d2, w in edges:
        G.add_edge(d1, d2, weight=w)
    return G


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați că fișierul de input există
    ensure_exists(DRUG_GENE_CSV)

    # TODO 2: încărcați tabelul drug-gene
    df = load_drug_gene_table(DRUG_GENE_CSV)

    # TODO 3: construiți mapping-ul drug -> set de gene
    drug2genes = build_drug2genes(df)

    # TODO 4: construiți graful bipartit și salvați-l (opțional)
    G_bip = build_bipartite_graph(drug2genes)
    plt.figure(figsize=(12, 9))
    pos = nx.spring_layout(G_bip, seed=42)

    node_colors = [
        "skyblue" if G_bip.nodes[n].get("bipartite") == "drug" else "salmon"
        for n in G_bip.nodes()
    ]
    node_sizes = [
        220 if G_bip.nodes[n].get("bipartite") == "drug" else 60
        for n in G_bip.nodes()
    ]

    nx.draw(
        G_bip, pos,
        node_color=node_colors,
        node_size=node_sizes,
        with_labels=False,
        width=0.3,
        alpha=0.85
    )

    # etichete doar pentru drugs (altfel e prea aglomerat)
    drug_nodes = [n for n, d in G_bip.nodes(data=True) if d.get("bipartite") == "drug"]
    nx.draw_networkx_labels(G_bip, pos, labels={d: d for d in drug_nodes}, font_size=7)

    plt.title("Drug–Gene Bipartite Network")
    plt.tight_layout()
    plt.savefig(OUT_GRAPH_DRUG_GENE_PNG, dpi=250)
    plt.close()
    print(f"[OUT]  {OUT_GRAPH_DRUG_GENE_PNG}")


    # TODO 5: generați și salvați sumarul pe medicamente
    df_summary = summarize_drugs(drug2genes)
    df_summary.to_csv(OUT_DRUG_SUMMARY, index=False)

    # TODO 6: calculați similaritatea între medicamente    
    edges = compute_drug_similarity_edges(drug2genes, min_sim=0.05)
    df_edges = edges_to_dataframe(edges)
    df_edges.to_csv(OUT_DRUG_SIMILARITY, index=False)

    
    G_sim = build_drug_similarity_graph(edges)
    with open(OUT_GRAPH_DRUG_SIM, "wb") as f:
        pickle.dump(G_sim, f)

    print("[INFO] Done.")
    print(f"[INFO] Loaded pairs: {len(df)} (unique drug-gene)")
    print(f"[INFO] Drugs: {len(drug2genes)} | Genes: {len({g for s in drug2genes.values() for g in s})}")
    print(f"[INFO] Bipartite graph: nodes={G_bip.number_of_nodes()} edges={G_bip.number_of_edges()}")
    print(f"[INFO] Similarity edges (min_sim=0.05): {len(edges)}")
    print(f"[OUT]  {OUT_DRUG_SUMMARY}")
    print(f"[OUT]  {OUT_DRUG_SIMILARITY}")
    print(f"[OUT]  {OUT_GRAPH_DRUG_GENE}")
    print(f"[OUT]  {OUT_GRAPH_DRUG_SIM}")