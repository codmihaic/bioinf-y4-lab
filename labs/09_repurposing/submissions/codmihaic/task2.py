from pathlib import Path
import itertools
import pandas as pd

HANDLE = "codmihaic"

IN_CSV = Path(f"data/work/codmihaic/lab09/drug_gene_{HANDLE}.csv")
OUT_DIR = Path("labs/09_repurposing/submissions/codmihaic")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
MIN_JACCARD = 0.0  
TOPK_PER_DRUG = None   


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    if inter == 0:
        return 0.0
    union = len(a | b)
    return inter / union


def main():
    df = pd.read_csv(IN_CSV)
    if not {"drug", "gene"}.issubset(df.columns):
        raise ValueError(f"CSV must have columns: drug, gene. Found: {list(df.columns)}")

    df = df[["drug", "gene"]].copy()
    df["drug"] = df["drug"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df.dropna(subset=["drug", "gene"])
    df = df[(df["drug"] != "") & (df["gene"] != "")]
    df = df.drop_duplicates(subset=["drug", "gene"])

    drug_to_genes = (
        df.groupby("drug")["gene"]
        .apply(lambda s: set(s.tolist()))
        .to_dict()
    )
    drugs = sorted(drug_to_genes.keys())

    rows = []
    for d1, d2 in itertools.combinations(drugs, 2):
        s1 = drug_to_genes[d1]
        s2 = drug_to_genes[d2]

        inter = len(s1 & s2)
        if inter == 0:
            continue 
        union = len(s1 | s2)
        jac = inter / union

        if jac < MIN_JACCARD:
            continue

        rows.append({
            "drug1": d1,
            "drug2": d2,
            "jaccard": float(jac),
            "intersection_size": int(inter),
            "union_size": int(union),
            "n_targets_drug1": int(len(s1)),
            "n_targets_drug2": int(len(s2)),
        })

    sim = pd.DataFrame(rows)
    if sim.empty:
        sim = pd.DataFrame(columns=[
            "drug1", "drug2", "jaccard",
            "intersection_size", "union_size",
            "n_targets_drug1", "n_targets_drug2"
        ])
        sim.to_csv(OUT_CSV, index=False)
        print("[WARN] No overlapping target sets found. Exported empty similarity file.")
        print(f"[OK] Wrote: {OUT_CSV}")
        return

    if TOPK_PER_DRUG is not None:
        a = sim.rename(columns={"drug1": "drug", "drug2": "other"})
        b = sim.rename(columns={"drug2": "drug", "drug1": "other"})
        bi = pd.concat([a.assign(direction=1), b.assign(direction=2)], ignore_index=True)

        bi = bi.sort_values(["drug", "jaccard"], ascending=[True, False])
        bi = bi.groupby("drug", as_index=False).head(TOPK_PER_DRUG)
        bi["u"] = bi.apply(lambda r: min(r["drug"], r["other"]), axis=1)
        bi["v"] = bi.apply(lambda r: max(r["drug"], r["other"]), axis=1)

        sim = (
            bi.groupby(["u", "v"], as_index=False)
            .agg({
                "jaccard": "max",
                "intersection_size": "max",
                "union_size": "max",
                "n_targets_drug1": "max",
                "n_targets_drug2": "max",
            })
            .rename(columns={"u": "drug1", "v": "drug2"})
        )

    sim = sim.sort_values(["jaccard", "intersection_size"], ascending=[False, False])
    sim.to_csv(OUT_CSV, index=False)

    print(f"[OK] Drugs: {len(drugs)} | Edges exported: {len(sim)}")
    print(f"[OK] Wrote: {OUT_CSV}")


if __name__ == "__main__":
    main()
