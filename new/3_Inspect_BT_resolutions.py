###crid

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from pathlib import Path

recluster_dir = Path("./markers_0.1_reclustered")
out_dir       = Path("./markers_BT_multiresolution")
plot_dir      = Path("./plots_BT_multiresolution")
out_dir.mkdir(exist_ok=True)
plot_dir.mkdir(exist_ok=True)

files = {
    "B": "B_scgpt.h5ad",
    "T": "T_scgpt.h5ad",
}

lineage_markers = {
    "B": [
        "MS4A1", "CD19", "CD79A", "CD79B", "PAX5",
        "TCL1A", "FCER2", "IGHD",
        "CD27", "AIM2",
        "MEF2B", "RGS13", "AICDA",
        "MZB1", "JCHAIN", "IGHG1", "IGHA1",
        "VPREB3", "MME",
        "CD3D", "LYZ", "GNLY",
    ],
    "T": [
        "CD3D", "CD3E",
        "CD4", "CD8A", "CD8B",
        "CCR7", "SELL", "LEF1", "TCF7", "IL7R",
        "GZMK", "GZMB", "PRF1", "GNLY",
        "FOXP3", "IL2RA", "CTLA4",
        "LAG3", "PDCD1", "HAVCR2", "TIGIT",
        "CXCR5", "ICOS", "BCL6",
        "TRDV2", "TRDC",
        "MX1", "STAT1", "ISG15",
        "SRGN", "ITM2A",
        "NPM1", "HSP90AB1",
        "MS4A1", "LYZ",
    ],
}

resolutions = {
    "B": [0.2, 0.3, 0.5],
    "T": [0.2, 0.3, 0.5, 0.8],
}

MAX_CELLS_FOR_MARKERS = 50000

for name, fname in files.items():
    print(f"\n{'─'*60}")
    print(f"── {name} cells ──")

    adata = sc.read_h5ad(recluster_dir / fname)
    print(f"  cells: {adata.n_obs}  genes: {adata.n_vars}")

    # ── fix embedding ────────────────────────────────────────────
    emb = adata.obsm["X_scGPT"]
    if hasattr(emb, "values"):
        emb = emb.values
    adata.obsm["X_scGPT"] = emb

    # ── set expression layer ─────────────────────────────────────
    adata.X = adata.layers["logcounts"].copy()

    # ── neighbors + UMAP (once, shared across resolutions) ───────
    print("  computing neighbors + UMAP...")
    sc.pp.neighbors(adata, use_rep="X_scGPT", n_neighbors=15)
    sc.tl.umap(adata)
    print("  ✓ UMAP done")

    # ── cluster at each resolution ────────────────────────────────
    for res in resolutions[name]:
        key = f"leiden_res_{res}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        adata.obs[key] = adata.obs[key].astype("category")
        n = adata.obs[key].nunique()
        print(f"\n  Resolution {res}: {n} clusters")
        print(adata.obs[key].value_counts().sort_index().to_string())

    # ── save h5ad with all clusterings ───────────────────────────
    out_h5ad = out_dir / fname
    adata.write_h5ad(out_h5ad)
    print(f"\n  ✓ saved h5ad → {out_h5ad}")

    # ── UMAP plots ────────────────────────────────────────────────
    for res in resolutions[name]:
        key = f"leiden_res_{res}"
        n   = adata.obs[key].nunique()
        sc.pl.umap(
            adata,
            color=key,
            legend_loc="on data",
            legend_fontsize=8,
            title=f"{name} cells — Leiden {res} ({n} clusters)",
            show=False,
            frameon=False,
        )
        plt.savefig(plot_dir / f"umap_{name}_res{res}.pdf", bbox_inches="tight", dpi=150)
        plt.close()
        print(f"  ✓ saved umap_{name}_res{res}.pdf")

    # ── dot plots — lineage markers ───────────────────────────────
    valid_markers = [g for g in lineage_markers[name] if g in adata.var_names]
    missing = [g for g in lineage_markers[name] if g not in adata.var_names]
    if missing:
        print(f"  markers not in object: {missing}")

    for res in resolutions[name]:
        key = f"leiden_res_{res}"
        n   = adata.obs[key].nunique()
        sc.pl.dotplot(
            adata,
            var_names=valid_markers,
            groupby=key,
            use_raw=False,
            standard_scale="var",
            show=False,
            figsize=(18, max(3, n * 0.45 + 1.5)),
            title=f"{name} cells — lineage markers at resolution {res}",
        )
        plt.savefig(plot_dir / f"dotplot_{name}_res{res}.pdf", bbox_inches="tight", dpi=150)
        plt.close()
        print(f"  ✓ saved dotplot_{name}_res{res}.pdf")

    # ── marker genes — subsample for speed ───────────────────────
    print(f"\n  subsampling to {MAX_CELLS_FOR_MARKERS} cells for marker computation...")
    adata_sub = sc.pp.subsample(
        adata, n_obs=min(MAX_CELLS_FOR_MARKERS, adata.n_obs),
        random_state=42, copy=True,
    )
    print(f"  ✓ subsampled: {adata_sub.n_obs} cells")

    for res in resolutions[name]:
        key = f"leiden_res_{res}"

        if key not in adata_sub.obs.columns:
            print(f"  skipping markers for res {res} — key missing after subsample")
            continue

        adata_sub.obs["cluster"] = adata_sub.obs[key]
        print(f"\n  computing markers at resolution {res}...")

        sc.tl.rank_genes_groups(
            adata_sub,
            groupby="cluster",
            method="wilcoxon",
            pts=True,
            use_raw=False,
        )

        results = sc.get.rank_genes_groups_df(
            adata_sub,
            group=None,
            pval_cutoff=0.05,
            log2fc_min=1,
        )
        results = results[results["pct_nz_group"] >= 0.5]
        results = results.sort_values(
            ["group", "logfoldchanges"], ascending=[True, False]
        )

        csv_path = out_dir / f"markers_{name}_res{res}.csv"
        results.to_csv(csv_path, index=False)
        print(f"  ✓ saved {csv_path.name}  ({len(results)} marker genes)")
        print(results.groupby("group")["names"].count().to_string())

        top10 = (
            results
            .sort_values(["group", "pct_nz_group"], ascending=[True, False])
            .groupby("group", observed=True)
            .head(10)
        )
        print(f"\n  Top 10 markers per cluster (by % expressed) — res {res}:")
        for group, grp_df in top10.groupby("group", observed=True):
            print(f"    cluster {group}: {', '.join(grp_df['names'].tolist())}")

# ── T cell final annotation at res 0.8 ───────────────────────────
print(f"\n{'─'*60}")
print("── T cell final annotation ──")

adata_T = sc.read_h5ad(out_dir / "T_scgpt.h5ad")
adata_T.X = adata_T.layers["logcounts"].copy()
adata_T.obs["leiden_res_0.8"] = adata_T.obs["leiden_res_0.8"].astype(str)

print("\nCluster sizes at res 0.8:")
print(adata_T.obs["leiden_res_0.8"].value_counts().sort_index().to_string())

# meaningful clusters with known identity
meaningful = {
    "0":  "CD8+ effector T cell",
    "1":  "Early activated T cell",
    "3":  "CD4+ effector memory T cell",
    "4":  "IFN-stimulated T cell",
    "6":  "Naive T cell",
    "9":  "Proliferating T cell",
    "10": "Quiescent / central memory T cell",
}

# artifact clusters to remove
artifacts = {"11", "12"}

all_clusters = adata_T.obs["leiden_res_0.8"].unique().tolist()
naive_clusters = [c for c in all_clusters
                  if c not in meaningful and c not in artifacts]
print(f"\nClusters merged into Naive T cell: {sorted(naive_clusters)}")

annotation_map = {}
for c in all_clusters:
    if c in meaningful:
        annotation_map[c] = meaningful[c]
    elif c in artifacts:
        annotation_map[c] = None
    else:
        annotation_map[c] = "Naive T cell"

adata_T.obs["cell_type_T"] = adata_T.obs["leiden_res_0.8"].map(annotation_map)
adata_T = adata_T[adata_T.obs["cell_type_T"].notna()].copy()

print(f"\nAfter removing artifacts: {adata_T.n_obs} cells")
print("\nFinal cell type counts:")
print(adata_T.obs["cell_type_T"].value_counts().to_string())

# UMAP final annotation
sc.pl.umap(
    adata_T,
    color="cell_type_T",
    legend_loc="right margin",
    legend_fontsize=9,
    show=False,
    frameon=False,
    title="T cells — final annotation",
)
plt.savefig(plot_dir / "umap_T_final_annotation.pdf", bbox_inches="tight", dpi=150)
plt.close()
print("saved umap_T_final_annotation.pdf")

# final markers on merged annotation
adata_T_sub = sc.pp.subsample(
    adata_T, n_obs=min(MAX_CELLS_FOR_MARKERS, adata_T.n_obs),
    random_state=42, copy=True,
)

sc.tl.rank_genes_groups(
    adata_T_sub,
    groupby="cell_type_T",
    method="wilcoxon",
    pts=True,
    use_raw=False,
)

results = sc.get.rank_genes_groups_df(
    adata_T_sub,
    group=None,
    pval_cutoff=0.05,
    log2fc_min=1,
)
results = results[results["pct_nz_group"] >= 0.5]
results = results.sort_values(["group", "logfoldchanges"], ascending=[True, False])

results.to_csv(out_dir / "markers_T_final.csv", index=False)
print(f"saved markers_T_final.csv ({len(results)} markers)")

top10 = (
    results
    .sort_values(["group", "pct_nz_group"], ascending=[True, False])
    .groupby("group", observed=True)
    .head(10)
)
print("\nTop 10 markers per cell type (by % expressed):")
for group, grp_df in top10.groupby("group", observed=True):
    print(f"  {group}: {', '.join(grp_df['names'].tolist())}")

adata_T.write_h5ad(out_dir / "T_scgpt_annotated.h5ad")
print("saved T_scgpt_annotated.h5ad")

print("\n── all done ──")
