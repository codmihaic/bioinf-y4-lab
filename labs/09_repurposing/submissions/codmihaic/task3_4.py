from pathlib import Path
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

HANDLE = "codmihaic"

IN_CSV = Path(f"data/work/codmihaic/lab09/drug_gene_{HANDLE}.csv")
DISEASE_TXT = Path(f"data/work/codmihaic/lab09/disease_genes_{HANDLE}.txt")
OUT_DIR = Path("labs/09_repurposing/submissions/codmihaic")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"
OUT_PNG = OUT_DIR / f"network_drug_gene_{HANDLE}.png"
TOPK_BY_PROX = 30      
TOPK_BY_DEG = 30       

def load_disease_genes(path: Path) -> list[str]:
    genes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g)
    seen = set()
    uniq = []
    for g in genes:
        if g not in seen:
            uniq.append(g)
            seen.add(g)
    return uniq

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

def compute_proximity(B: nx.Graph, drugs: list[str], disease_genes: list[str]) -> pd.DataFrame:
    disease_in_graph = [g for g in disease_genes if B.has_node(g)]
    missing = [g for g in disease_genes if not B.has_node(g)]

    if len(disease_in_graph) == 0:
        raise ValueError(
            "Niciuna dintre disease genes nu există în graful drug-gene. "
            "Alege gene care apar în drug_gene_<handle>.csv."
        )

    n_targets = {d: int(B.degree(d)) for d in drugs}
    distances_by_drug = {d: [] for d in drugs}
    for g in disease_in_graph:
        dist = nx.single_source_shortest_path_length(B, g)
        for d in drugs:
            if d in dist:
                distances_by_drug[d].append(dist[d])

    rows = []
    total = len(disease_in_graph)
    for d in drugs:
        dists = distances_by_drug[d]
        reachable = len(dists)
        if reachable == 0:
            mean_dist = float("inf")
        else:
            mean_dist = sum(dists) / reachable

        rows.append({
            "drug": d,
            "mean_distance": mean_dist,
            "n_reachable": reachable,
            "n_disease_genes_total": total,
            "n_targets": n_targets[d],
        })

    out = pd.DataFrame(rows)
    out["is_inf"] = out["mean_distance"].apply(lambda x: 1 if x == float("inf") else 0)
    out = out.sort_values(
        ["is_inf", "mean_distance", "n_reachable", "n_targets", "drug"],
        ascending=[True, True, False, False, True]
    ).drop(columns=["is_inf"])

    if missing:
        print(f"[WARN] Disease genes missing from graph ({len(missing)}): {missing}")
    print(f"[OK] Disease genes used (in graph): {disease_in_graph}")
    return out, disease_in_graph, missing

def make_visualization(B: nx.Graph, priority_df: pd.DataFrame, drugs: list[str], genes: list[str]):
    top_prox = priority_df.head(TOPK_BY_PROX)["drug"].tolist()
    top_deg = (
        priority_df.sort_values("n_targets", ascending=False)
        .head(TOPK_BY_DEG)["drug"].tolist()
    )
    selected_drugs = sorted(set(top_prox + top_deg))
    selected_genes = set()
    for d in selected_drugs:
        selected_genes.update(B.neighbors(d))

    nodes = set(selected_drugs) | set(selected_genes)
    H = B.subgraph(nodes).copy()
    pos = nx.spring_layout(H, seed=42, k=None)
    drug_nodes = [n for n, data in H.nodes(data=True) if data.get("node_type") == "drug"]
    gene_nodes = [n for n, data in H.nodes(data=True) if data.get("node_type") == "gene"]

    deg_full = {d: int(B.degree(d)) for d in drug_nodes}
    drug_sizes = [max(80, deg_full[d] * 30) for d in drug_nodes]
    gene_sizes = [30 for _ in gene_nodes]
    plt.figure(figsize=(14, 10))

    nx.draw_networkx_nodes(
        H, pos,
        nodelist=gene_nodes,
        node_color="red",
        node_size=gene_sizes,
        alpha=0.7
    )

    nx.draw_networkx_nodes(
        H, pos,
        nodelist=drug_nodes,
        node_color="blue",
        node_size=drug_sizes,
        alpha=0.8
    )

    nx.draw_networkx_edges(H, pos, alpha=0.25, width=1.0)
    nx.draw_networkx_labels(H, pos, labels={d: d for d in drug_nodes}, font_size=8)
    plt.title(f"Drug–Gene bipartite network (subset) — {HANDLE}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    plt.close()

    print(f"[OK] Wrote: {OUT_PNG}")
    print(f"[OK] Visualized nodes: {H.number_of_nodes()} | edges: {H.number_of_edges()}")

def main():
    df = pd.read_csv(IN_CSV)
    if not {"drug", "gene"}.issubset(df.columns):
        raise ValueError(f"CSV must contain columns drug, gene. Found: {list(df.columns)}")
    
    B, drugs, genes = build_bipartite(df)
    print(f"[OK] Graph: nodes={B.number_of_nodes()} edges={B.number_of_edges()} drugs={len(drugs)} genes={len(genes)}")
    disease_genes = load_disease_genes(DISEASE_TXT)
    print(f"[OK] Loaded disease genes from {DISEASE_TXT}: {disease_genes}")

    # Task 3: proximity + export
    priority_df, disease_in_graph, missing = compute_proximity(B, drugs, disease_genes)
    priority_df.to_csv(OUT_PRIORITY, index=False)
    print(f"[OK] Wrote: {OUT_PRIORITY}")

    # Task 4: visualization
    make_visualization(B, priority_df, drugs, genes)

if __name__ == "__main__":
    main()
