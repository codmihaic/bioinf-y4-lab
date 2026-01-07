from pathlib import Path
import pandas as pd
import networkx as nx

HANDLE = "codmihaic"
IN_CSV = Path(f"data/work/codmihaic/lab09/drug_gene_{HANDLE}.csv")

OUT_DIR = Path("labs/09_repurposing/submissions/codmihaic")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_GRAPHML = OUT_DIR / f"drug_gene_network_{HANDLE}.graphml"  # util pt task 3-4 (opțional)

def main():
    df = pd.read_csv(IN_CSV)
    expected_cols = {"drug", "gene"}
    if not expected_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns {expected_cols}. Found: {list(df.columns)}")

    df = df[["drug", "gene"]].copy()
    df["drug"] = df["drug"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df.dropna(subset=["drug", "gene"])
    df = df[(df["drug"] != "") & (df["gene"] != "")]
    df = df.drop_duplicates(subset=["drug", "gene"])

    B = nx.Graph()
    drugs = sorted(df["drug"].unique())
    genes = sorted(df["gene"].unique())

    B.add_nodes_from(drugs, bipartite=0, node_type="drug")
    B.add_nodes_from(genes, bipartite=1, node_type="gene")
    edges = list(df.itertuples(index=False, name=None))  # (drug, gene)
    B.add_edges_from(edges)
    rows = []
    for d in drugs:
        n_targets = B.degree(d)
        rows.append({"drug": d, "n_targets": int(n_targets)})

    summary = pd.DataFrame(rows).sort_values(["n_targets", "drug"], ascending=[False, True])
    summary.to_csv(OUT_SUMMARY, index=False)
    nx.write_graphml(B, OUT_GRAPHML)

    print(f"[OK] Loaded edges: {len(edges)}")
    print(f"[OK] Drugs: {len(drugs)} | Genes: {len(genes)} | Nodes total: {B.number_of_nodes()} | Edges: {B.number_of_edges()}")
    print(f"[OK] Wrote: {OUT_SUMMARY}")
    print(f"[OK] Wrote: {OUT_GRAPHML}")

if __name__ == "__main__":
    main()
