from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
from community import community_louvain


HANDLE = "codmihaic"
FASTA_PATH = Path(f"data/work/{HANDLE}/lab01/my_tp53.fa")  
MAX_NT = 1000                 
K = 4                         
WINDOW = 200                 
STEP = 50                     
VAR_DROP_PERCENTILE = 25     
CORR_METHOD = "pearson"      
THRESH = 0.7                
HUB_PERCENTILE = 95           

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

MODULES_CSV = OUT / f"modules_tp53_{HANDLE}.csv"
HUBS_CSV = OUT / f"hubs_tp53_{HANDLE}.csv"
NETWORK_PNG = OUT / f"network_tp53_{HANDLE}.png"
KMERS_MATRIX_CSV = OUT / f"kmers_matrix_tp53_{HANDLE}.csv"        
KMERS_MATRIX_LOG_CSV = OUT / f"kmers_matrix_log_tp53_{HANDLE}.csv" 


def read_fasta_sequence(path: Path) -> str:
    """Citește un FASTA simplu și întoarce secvența concatenată (A/C/G/T/N)."""
    seq = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)

def windows_from_seq(seq: str, window: int, step: int):
    """Generează ferestre (subsecvențe) dintr-o secvență."""
    for start in range(0, max(1, len(seq) - window + 1), step):
        yield start, seq[start:start + window]

def kmer_counts(s: str, k: int) -> Counter:
    """Numără k-mer-urile (doar cele fără N)."""
    c = Counter()
    for i in range(0, len(s) - k + 1):
        kmer = s[i:i+k]
        if "N" in kmer:
            continue
        c[kmer] += 1
    return c

def build_kmer_matrix(seq: str, k: int, window: int, step: int) -> pd.DataFrame:
    """
    Construiește matrice windows × k-mers (frecvențe relative).
    """
    rows = []
    idx = []

    for start, w in windows_from_seq(seq, window, step):
        counts = kmer_counts(w, k)
        total = sum(counts.values()) if sum(counts.values()) > 0 else 1
        freqs = {km: v / total for km, v in counts.items()}
        rows.append(freqs)
        idx.append(f"win_{start}_{start+len(w)}")

    df = pd.DataFrame(rows, index=idx).fillna(0.0)
    return df


seq = read_fasta_sequence(FASTA_PATH)
seq = seq[:MAX_NT]
seq = "".join([ch for ch in seq if ch in {"A", "C", "G", "T", "N"}])
if len(seq) < WINDOW:
    raise ValueError(f"Secvența după subset/curățare are {len(seq)} nt, < WINDOW={WINDOW}. Micșorează WINDOW.")

X = build_kmer_matrix(seq, K, WINDOW, STEP)
X.to_csv(KMERS_MATRIX_CSV)
print("k-mer matrix saved:", KMERS_MATRIX_CSV, "shape=", X.shape)

# =========================
# 1) PREPROCESS: log2(x+1) + variance filter
# =========================
X_log = np.log2(X + 1.0)
X_log.to_csv(KMERS_MATRIX_LOG_CSV)
vars_ = X_log.var(axis=0)  # varianță pe coloane (k-mers)
thr = np.percentile(vars_.values, VAR_DROP_PERCENTILE)
keep_cols = vars_[vars_ > thr].index
X_filt = X_log[keep_cols]
print("After variance filter:", X_filt.shape)


# =========================
# 2) CORRELATION + GRAPH + LOUVAIN
# =========================
corr = X_filt.corr(method=CORR_METHOD)
adj = corr.abs() > THRESH
np.fill_diagonal(adj.values, 0)
G = nx.from_pandas_adjacency(adj)
print("Graph:", "nodes=", G.number_of_nodes(), "edges=", G.number_of_edges())

partition = community_louvain.best_partition(G)
modules_df = pd.DataFrame({
    "kmer": list(partition.keys()),
    "module": [partition[k] for k in partition.keys()]
}).sort_values(["module", "kmer"])
modules_df.to_csv(MODULES_CSV, index=False)
print("Modules saved:", MODULES_CSV)


# =========================
# 3) HUBS + VISUALIZATION
# =========================
deg = dict(G.degree())
hub_thr = np.percentile(list(deg.values()), HUB_PERCENTILE)
hubs = [n for n, d in deg.items() if d >= hub_thr]

hubs_df = pd.DataFrame({
    "kmer": hubs,
    "degree": [deg[h] for h in hubs],
    "module": [partition[h] for h in hubs]
}).sort_values(["degree", "kmer"], ascending=[False, True])

hubs_df.to_csv(HUBS_CSV, index=False)
print("Hubs saved:", HUBS_CSV)

plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, seed=42)
node_colors = [partition[n] for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_size=35, node_color=node_colors, cmap=plt.cm.tab20)
nx.draw_networkx_edges(G, pos, alpha=0.15)

nx.draw_networkx_nodes(G, pos, nodelist=hubs, node_size=140, node_color="red")
plt.title(f"TP53 k-mer co-occurrence network (K={K}, window={WINDOW}, step={STEP})\nHub k-mers highlighted")
plt.axis("off")
plt.savefig(NETWORK_PNG, dpi=300, bbox_inches="tight")
plt.close()
print("Network image saved:", NETWORK_PNG)

# =========================
# BONUS: Export pentru Cytoscape / Gephi
# =========================

# 1) Node table (pentru Cytoscape/Gephi)
nodes_df = pd.DataFrame({
    "id": list(G.nodes()),
    "module": [partition[n] for n in G.nodes()],
    "degree": [deg[n] for n in G.nodes()],
    "is_hub": [1 if n in set(hubs) else 0 for n in G.nodes()]
}).sort_values(["module", "degree"], ascending=[True, False])

NODES_CSV = OUT / f"nodes_tp53_{HANDLE}.csv"
nodes_df.to_csv(NODES_CSV, index=False)
print("Nodes table saved:", NODES_CSV)

# 2) Edge list (pentru Cytoscape/Gephi)
edges_df = nx.to_pandas_edgelist(G)
edges_df["weight"] = edges_df.apply(
    lambda r: float(abs(corr.loc[r["source"], r["target"]])),
    axis=1
)

EDGES_CSV = OUT / f"edges_tp53_{HANDLE}.csv"
edges_df.to_csv(EDGES_CSV, index=False)
print("Edges table saved:", EDGES_CSV)


G_attr = G.copy()
nx.set_node_attributes(G_attr, partition, "module")
nx.set_node_attributes(G_attr, deg, "degree")
nx.set_node_attributes(G_attr, {n: (n in set(hubs)) for n in G.nodes()}, "is_hub")

GRAPHML = OUT / f"tp53_{HANDLE}.graphml"
nx.write_graphml(G_attr, GRAPHML)
print("GraphML saved:", GRAPHML)
