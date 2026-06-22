###this is from cp

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path

subset_dir = Path("./subsetted_h5ads")
out_dir    = Path("./markers_0.1_reclustered")
out_dir.mkdir(exist_ok=True)

files = {
    "Myeloid": "Myeloid_scgpt.h5ad",
    "NK":      "NK_scgpt.h5ad",
    "Plasma":  "Plasma_scgpt.h5ad",
    "B":       "B_scgpt.h5ad",
    "T":       "T_scgpt.h5ad",
}

for name, fname in files.items():
    print(f"\n{'─'*50}")
    print(f"── {name} ──")

    adata = sc.read_h5ad(subset_dir / fname)
    print(f"  cells: {adata.n_obs}  genes: {adata.n_vars}")

    # ── coerce X_scGPT to numpy array if it came in as a DataFrame ──────────
    emb = adata.obsm["X_scGPT"]
    if isinstance(emb, pd.DataFrame):
        emb = emb.values
    adata.obsm["X_scGPT"] = emb
    print(f"  X_scGPT shape: {emb.shape}")

    # ── set logcounts as adata.X for marker computation ──────────────────────
    adata.X = adata.layers["logcounts"].copy()
    print(f"  ✓ adata.X set to logcounts (max={adata.X.data.max() if sp.issparse(adata.X) else adata.X.max():.4f})")

    # ── neighbors using scGPT embedding ─────────────────────────────────────
    sc.pp.neighbors(adata, use_rep="X_scGPT", n_neighbors=15)
    print(f"  ✓ neighbors computed from X_scGPT")

    # ── UMAP ─────────────────────────────────────────────────────────────────
    sc.tl.umap(adata)
    print(f"  ✓ UMAP done")

    # ── clustering ───────────────────────────────────────────────────────────
    sc.tl.leiden(adata, resolution=0.1, key_added="leiden_res_0.1_reclustered")
    adata.obs["cluster"] = adata.obs["leiden_res_0.1_reclustered"].astype("category")
    clusters = sorted(adata.obs["cluster"].cat.categories.tolist())
    print(f"  ✓ Leiden done — clusters: {clusters}")
    print(f"  cluster sizes:\n{adata.obs['cluster'].value_counts().sort_index().to_string()}")

    # ── marker genes (on logcounts via adata.X) ──────────────────────────────
    print(f"  computing markers...")
    sc.tl.rank_genes_groups(
        adata,
        groupby = "cluster",
        method  = "wilcoxon",
        pts     = True,
        use_raw = False,
    )

    results = sc.get.rank_genes_groups_df(
        adata,
        group       = None,
        pval_cutoff = 0.05,
        log2fc_min  = 1,
    )
    results = results[results["pct_nz_group"] >= 0.5]
    results = results.sort_values(["group", "logfoldchanges"], ascending=[True, False])

    out = out_dir / f"markers_{name}.csv"
    results.to_csv(out, index=False)
    print(f"  ✓ saved {out}  ({len(results)} marker genes total)")
    print(results.groupby("group")["names"].count().to_string())

    # ── top 10 per cluster by % expressed ────────────────────────────────────
    top10 = (
        results
        .sort_values(["group", "pct_nz_group"], ascending=[True, False])
        .groupby("group", observed=True)
        .head(10)
    )
    print(f"\n  Top 10 markers per cluster (by % expressed):")
    for group, grp_df in top10.groupby("group", observed=True):
        genes = grp_df["names"].tolist()
        print(f"    cluster {group}: {', '.join(genes)}")

    # ── save ─────────────────────────────────────────────────────────────────
    out_h5ad = out_dir / fname
    adata.write_h5ad(out_h5ad)
    print(f"  ✓ saved reclustered h5ad → {out_h5ad}")
