#!/usr/bin/env python3
"""
02_run_scgpt.py
Zero-shot cell-type annotation with scGPT on subsetted h5ad files.

Usage:
    python 02_run_scgpt.py \
        --h5ad_dir   results/h5ad \
        --model_dir  /path/to/scgpt_human_model \
        --out_dir    results/scgpt \
        [--n_hvg      3000] \
        [--batch_size 64]

Requirements:
    pip install scgpt scanpy torch anndata pandas
    scGPT model weights: https://github.com/bowang-lab/scGPT
"""

import argparse
import os
import glob
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch

# scGPT imports — requires: pip install scgpt
try:
    import scgpt
    from scgpt.tasks import GeneEmbedding
    from scgpt.tokenizer.gene_tokenizer import GeneVocab
    from scgpt.model import TransformerModel
    from scgpt.preprocess import Preprocessor
    from scgpt.tasks.cell_emb import embed_data
except ImportError:
    raise ImportError(
        "scGPT not found. Install with:\n"
        "  pip install scgpt\n"
        "and download model weights from https://github.com/bowang-lab/scGPT"
    )


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description="scGPT cell-type annotation")
    p.add_argument("--h5ad_dir",   required=True,
                   help="Directory containing *.h5ad files")
    p.add_argument("--model_dir",  required=True,
                   help="Path to scGPT pretrained model directory "
                        "(contains best_model.pt, vocab.json, args.json)")
    p.add_argument("--out_dir",    default="results/scgpt",
                   help="Output directory [default: results/scgpt]")
    p.add_argument("--n_hvg",      type=int, default=3000,
                   help="Number of highly variable genes [default: 3000]")
    p.add_argument("--batch_size", type=int, default=64,
                   help="Inference batch size [default: 64]")
    p.add_argument("--n_bins",     type=int, default=51,
                   help="Number of expression bins for scGPT [default: 51]")
    p.add_argument("--device",     default="auto",
                   help="'auto', 'cpu', 'cuda', or 'mps' [default: auto]")
    return p.parse_args()


# =============================================================================
# Helpers
# =============================================================================
def resolve_device(device_arg: str) -> torch.device:
    if device_arg == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        elif torch.backends.mps.is_available():
            return torch.device("mps")
        else:
            return torch.device("cpu")
    return torch.device(device_arg)


def preprocess_adata(adata: ad.AnnData, n_hvg: int, n_bins: int,
                     vocab: GeneVocab) -> ad.AnnData:
    """
    Standard scGPT preprocessing:
      1. Keep raw counts in .layers["counts"]
      2. Normalise + log1p
      3. Select HVGs that overlap the model vocabulary
      4. Bin expression values
    """
    print(f"  Input: {adata.n_obs} cells × {adata.n_vars} genes")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Normalise
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=False)

    # Keep only genes in the scGPT vocabulary
    vocab_genes  = set(vocab.get_stoi().keys())
    overlap_mask = adata.var_names.isin(vocab_genes)
    hvg_mask     = adata.var["highly_variable"]
    keep_mask    = overlap_mask & hvg_mask

    print(f"  HVGs overlapping vocab: {keep_mask.sum()} / {n_hvg}")
    adata = adata[:, keep_mask].copy()

    # Bin expression (scGPT uses discrete tokens)
    preprocessor = Preprocessor(
        use_key        = "X",
        filter_gene_by_counts   = False,
        filter_cell_by_counts   = False,
        normalize_total          = False,   # already done above
        result_normed_key        = "X_normed",
        log1p                    = False,   # already done above
        result_log1p_key         = "X_log1p",
        subset_hvg               = False,   # already done above
        binning                  = n_bins,
        result_binned_key        = "X_binned",
    )
    preprocessor(adata, batch_key=None)

    return adata


def run_scgpt_annotation(adata: ad.AnnData, model_dir: str,
                         device: torch.device, batch_size: int,
                         n_bins: int, vocab: GeneVocab) -> ad.AnnData:
    """
    Embed cells with scGPT and store embeddings in adata.obsm['X_scGPT'].
    Then run Leiden clustering on the embedding for label-free annotation.
    """
    print("  Generating scGPT cell embeddings...")

    # embed_data returns an np.ndarray of shape (n_cells, embed_dim)
    embeddings = embed_data(
        adata,
        model_dir     = model_dir,
        gene_col       = "index",          # use adata.var_names
        max_length     = 1200,
        batch_size     = batch_size,
        device         = device,
        use_fast_transformer = True,
        return_new_adata     = False,
    )

    adata.obsm["X_scGPT"] = embeddings

    # Cluster on scGPT embedding
    print("  Clustering on scGPT embedding (Leiden)...")
    sc.pp.neighbors(adata, use_rep="X_scGPT", n_neighbors=15)
    sc.tl.leiden(adata, resolution=0.5, key_added="scgpt_leiden")
    sc.tl.umap(adata)

    return adata


# =============================================================================
# Main per-file pipeline
# =============================================================================
def process_file(h5ad_path: str, model_dir: str, out_dir: str,
                 n_hvg: int, n_bins: int, batch_size: int,
                 device: torch.device, vocab: GeneVocab) -> None:

    lineage = os.path.basename(h5ad_path).replace(".h5ad", "")
    print(f"\n{'='*60}")
    print(f"  Lineage: {lineage}")
    print(f"{'='*60}")

    # ── Load ──────────────────────────────────────────────────────────────────
    print("  Loading h5ad...")
    adata = sc.read_h5ad(h5ad_path)

    # ── Preprocess ────────────────────────────────────────────────────────────
    adata = preprocess_adata(adata, n_hvg=n_hvg, n_bins=n_bins, vocab=vocab)

    # ── Embed + cluster ───────────────────────────────────────────────────────
    adata = run_scgpt_annotation(adata, model_dir=model_dir, device=device,
                                 batch_size=batch_size, n_bins=n_bins, vocab=vocab)

    # ── Save outputs ──────────────────────────────────────────────────────────
    # Annotated h5ad
    out_h5ad = os.path.join(out_dir, f"{lineage}_scgpt.h5ad")
    adata.write_h5ad(out_h5ad)
    print(f"  Saved h5ad: {out_h5ad}")

    # Embeddings as numpy
    emb_path = os.path.join(out_dir, f"{lineage}_scgpt_embeddings.npy")
    np.save(emb_path, adata.obsm["X_scGPT"])
    print(f"  Saved embeddings: {emb_path}")

    # Cell-level metadata + cluster labels as CSV
    meta_cols = ["scgpt_leiden"] + \
                [c for c in adata.obs.columns if c not in ["scgpt_leiden"]]
    meta_path = os.path.join(out_dir, f"{lineage}_scgpt_metadata.csv")
    adata.obs[meta_cols].to_csv(meta_path)
    print(f"  Saved metadata: {meta_path}")

    # UMAP coordinates
    umap_df   = pd.DataFrame(adata.obsm["X_umap"],
                              index   = adata.obs_names,
                              columns = ["UMAP_1", "UMAP_2"])
    umap_df["scgpt_leiden"] = adata.obs["scgpt_leiden"].values
    umap_path = os.path.join(out_dir, f"{lineage}_scgpt_umap.csv")
    umap_df.to_csv(umap_path)
    print(f"  Saved UMAP coords: {umap_path}")


# =============================================================================
# Entry point
# =============================================================================
def main():
    args   = parse_args()
    device = resolve_device(args.device)
    print(f"Device: {device}")

    os.makedirs(args.out_dir, exist_ok=True)

    # Load vocab once — shared across all lineages
    vocab_path = os.path.join(args.model_dir, "vocab.json")
    if not os.path.exists(vocab_path):
        raise FileNotFoundError(
            f"vocab.json not found in model_dir: {args.model_dir}\n"
            "Download scGPT weights from https://github.com/bowang-lab/scGPT"
        )
    print(f"Loading vocab from: {vocab_path}")
    vocab = GeneVocab.from_file(vocab_path)

    # Find h5ad files
    h5ad_files = sorted(glob.glob(os.path.join(args.h5ad_dir, "*.h5ad")))
    if not h5ad_files:
        raise FileNotFoundError(f"No .h5ad files found in: {args.h5ad_dir}")

    print(f"\nFound {len(h5ad_files)} h5ad file(s):")
    for f in h5ad_files:
        print(f"  - {os.path.basename(f)}")

    # Process each lineage
    for h5ad_path in h5ad_files:
        process_file(
            h5ad_path  = h5ad_path,
            model_dir  = args.model_dir,
            out_dir    = args.out_dir,
            n_hvg      = args.n_hvg,
            n_bins     = args.n_bins,
            batch_size = args.batch_size,
            device     = device,
            vocab      = vocab,
        )

    print(f"\n=== scGPT annotation complete ===")
    print(f"Results written to: {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
