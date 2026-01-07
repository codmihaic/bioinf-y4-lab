from pathlib import Path
import pandas as pd
import networkx as nx
import numpy as np

HANDLE = "codmihaic"

IN_CSV = Path(f"data/work/codmihaic/lab09/drug_gene_{HANDLE}.csv")
DISEASE_TXT = Path(f"data/work/codmihaic/lab09/disease_genes_{HANDLE}.txt")
OUT_DIR = Path("labs/09_repurposing/submissions/codmihaic")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_RWR = OUT_DIR / f"drug_rwr_{HANDLE}.csv"
OUT_COMPARE = OUT_DIR / f"rwr_vs_proximity_{HANDLE}.csv"
PRIORITY_CSV = OUT_DIR / f"drug_priority_{HANDLE}.csv"

RESTART = 0.5 
MAX_ITERS = 200
TOL = 1e-10

def load_disease_genes(path: Path) -> list[str]:
    genes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g)
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out

def build_bipartite(df: pd.DataFrame) -> tuple[nx.Graph, list[str], list[str]]:
    df = df[["drug", "gene"]].copy()
    df["drug"] = df["drug"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df.dropna(subset=["drug", "gene"])
    df = df[(df["drug"] != "") & (df["gene"] != "")]
    df = df.drop_duplicates(subset=["drug", "gene"])

    drugs = sorted(df["drug"].unique())
    genes = sorted(df["gene"].unique())

    B = nx.Graph()
    B.add_nodes_from(drugs, bipartite=0, node_type="drug")
    B.add_nodes_from(genes, bipartite=1, node_type="gene")
    B.add_edges_from(df.itertuples(index=False, name=None))
    return B, drugs, genes

def rwr(B: nx.Graph, seeds: list[str], restart: float, max_iters: int, tol: float):
    nodes = list(B.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    seeds_in_graph = [s for s in seeds if s in idx]
    if not seeds_in_graph:
        raise ValueError("Niciun seed (disease gene) nu există în graf.")

    p0 = np.zeros(n, dtype=float)
    for s in seeds_in_graph:
        p0[idx[s]] = 1.0
    p0 /= p0.sum()

    p = p0.copy()
    neighbors = []
    invdeg = np.zeros(n, dtype=float)
    for i, node in enumerate(nodes):
        neigh = list(B.neighbors(node))
        neighbors.append([idx[v] for v in neigh])
        deg = len(neigh)
        invdeg[i] = 1.0 / deg if deg > 0 else 0.0

    for it in range(max_iters):
        p_next = np.zeros(n, dtype=float)
        for u in range(n):
            if invdeg[u] == 0.0:
                continue
            share = (1.0 - restart) * p[u] * invdeg[u]
            for v in neighbors[u]:
                p_next[v] += share

        p_next += restart * p0
        diff = np.abs(p_next - p).sum()
        p = p_next
        if diff < tol:
            print(f"[OK] Converged at iter {it+1} (L1 diff={diff:.2e})")
            break
    else:
        print(f"[WARN] Not fully converged after {max_iters} iters (last L1 diff={diff:.2e})")
    return nodes, p, seeds_in_graph

def main():
    df = pd.read_csv(IN_CSV)
    if not {"drug", "gene"}.issubset(df.columns):
        raise ValueError(f"CSV must contain columns drug, gene. Found: {list(df.columns)}")

    B, drugs, genes = build_bipartite(df)
    disease_genes = load_disease_genes(DISEASE_TXT)
    nodes, p, seeds_used = rwr(B, disease_genes, RESTART, MAX_ITERS, TOL)
    score = {node: float(p[i]) for i, node in enumerate(nodes)}

    rows = []
    for d in drugs:
        rows.append({
            "drug": d,
            "rwr_score": score.get(d, 0.0),
            "n_targets": int(B.degree(d))
        })

    rwr_df = pd.DataFrame(rows).sort_values(["rwr_score", "n_targets", "drug"], ascending=[False, False, True])
    rwr_df["rwr_rank"] = range(1, len(rwr_df) + 1)
    rwr_df.to_csv(OUT_RWR, index=False)
    print(f"[OK] Seeds used in graph: {seeds_used}")
    print(f"[OK] Wrote: {OUT_RWR}")

    if PRIORITY_CSV.exists():
        prox = pd.read_csv(PRIORITY_CSV)
        prox = prox.copy()
        prox["prox_rank"] = range(1, len(prox) + 1)
        comp = prox.merge(rwr_df[["drug", "rwr_score", "rwr_rank"]], on="drug", how="left")
        comp.to_csv(OUT_COMPARE, index=False)

        k = 20
        top_rwr = set(rwr_df.head(k)["drug"])
        top_prox = set(prox.head(k)["drug"])
        overlap = len(top_rwr & top_prox)

        print(f"[OK] Wrote: {OUT_COMPARE}")
        print(f"[INFO] Top-{k} overlap RWR vs Proximity: {overlap}/{k}")
    else:
        print(f"[INFO] {PRIORITY_CSV} not found, skipped comparison. Run Task 3 first for proximity ranking.")

if __name__ == "__main__":
    main()
