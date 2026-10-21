"""Merge QC'd samples -> normalize -> HVG -> PCA -> integrate -> neighbors/Leiden/UMAP."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def merge(adatas: dict[str, ad.AnnData], metadata: dict[str, dict], cfg: dict) -> ad.AnnData:
    for name, a in adatas.items():
        for k, v in metadata[name].items():
            a.obs[k] = v
        a.obs_names = [f"{name}_{bc}" for bc in a.obs_names]  # barcodes collide across samples
    # Outer join keeps genes seen in any sample (missing = 0), like Seurat's merge.
    adata = ad.concat(adatas, join="outer", merge="same", fill_value=0)
    for col in ["sample", "study", "treatment", "mouse", "bacSpecies", "Gram"]:
        adata.obs[col] = adata.obs[col].astype("category")
    sc.pp.filter_genes(adata, min_cells=cfg["qc"]["min_cells_per_gene"])
    return adata


def normalize_and_reduce(adata: ad.AnnData, cfg: dict) -> None:
    p = cfg["processing"]
    adata.X = adata.layers["counts"].copy()
    # HVGs chosen within each sample then combined, so one lab can't dominate.
    sc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"], flavor="seurat_v3",
                                layer="counts", batch_key=p["batch_key"])
    if p["normalization"] == "pearson_residuals":
        # Analytic Pearson residuals: the Scanpy analogue of SCTransform v2.
        hv = adata[:, adata.var["highly_variable"]].copy()
        sc.experimental.pp.normalize_pearson_residuals(hv)
        sc.pp.pca(hv, n_comps=p["n_pcs"], random_state=p["seed"])
        adata.obsm["X_pca"] = hv.obsm["X_pca"]
        adata.uns["pca"] = hv.uns["pca"]
        sc.pp.normalize_total(adata, target_sum=1e4)  # log-normalised X kept for plots/DE
        sc.pp.log1p(adata)
    elif p["normalization"] == "lognorm":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=p["n_pcs"], mask_var="highly_variable", random_state=p["seed"])
    else:
        raise ValueError(f"unknown normalization {p['normalization']!r}")


def integrate(adata: ad.AnnData, cfg: dict) -> str:
    """Returns the obsm key to build neighbors on."""
    p = cfg["processing"]
    method = p["integration"]
    if method == "none":
        return "X_pca"
    if method == "harmony":
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, p["batch_key"],
                                   random_state=p["seed"], verbose=False)
        Z = np.asarray(ho.Z_corr)
        if Z.shape[0] != adata.n_obs:  # older harmonypy returns PCs x cells
            Z = Z.T
        adata.obsm["X_pca_harmony"] = Z
        return "X_pca_harmony"
    raise ValueError(f"unknown integration {method!r}")


def cluster(adata: ad.AnnData, rep: str, cfg: dict) -> None:
    p = cfg["processing"]
    n_pcs = min(p["neighbor_pcs"], adata.obsm[rep].shape[1])
    sc.pp.neighbors(adata, n_neighbors=p["n_neighbors"], n_pcs=n_pcs, use_rep=rep,
                    random_state=p["seed"])
    sc.tl.leiden(adata, resolution=p["leiden_resolution"], key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False, random_state=p["seed"])
    sc.tl.umap(adata, random_state=p["seed"])


def sample_summary(adata: ad.AnnData, cells: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """One row per sample: depth, retention and a sex call. Feeds the Streamlit overview."""
    g = cells.groupby("sample")
    kept = cells[cells["keep"]]
    out = pd.DataFrame({
        "cells_raw": g.size(),
        "cells_kept": kept.groupby("sample").size(),
        "median_counts_raw": g["total_counts"].median(),
        "median_genes_raw": g["n_genes_by_counts"].median(),
        "median_counts_kept": kept.groupby("sample")["total_counts"].median(),
        "pct_doublets": 100 * g["predicted_doublet"].mean(),
        "female_cpm": 1e6 * g["female_counts"].sum() / g["total_counts"].sum(),
        "male_cpm": 1e6 * g["male_counts"].sum() / g["total_counts"].sum(),
    }).fillna(0)
    out["pct_kept"] = 100 * out["cells_kept"] / out["cells_raw"]
    ratio = (out["female_cpm"] + 1) / (out["male_cpm"] + 1)
    # Samples pool several mice (~50:50 male:female per Griffin), so "mixed" is the expected call.
    out["sex_call"] = np.select([ratio > 10, ratio < 0.1], ["female only", "male only"], "mixed")
    meta = adata.obs.groupby("sample", observed=True)[["study", "treatment", "mouse", "Gram"]].first()
    return out.join(meta).reset_index(names="sample")
