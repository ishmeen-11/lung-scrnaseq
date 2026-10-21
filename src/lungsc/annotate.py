"""Marker-based cluster annotation using Griffin's curated gene lists.

For every cell type, score each cell with sc.tl.score_genes (mean expression of the
markers minus a matched random background). Each Leiden cluster is then labelled
with the cell type that has the highest mean score in that cluster. The margin to
the runner-up is kept so weak calls can be reviewed by hand instead of trusted.
"""
from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

log = logging.getLogger("lungsc")


def score_cell_types(adata: ad.AnnData, markers: dict[str, list[str]], seed: int) -> list[str]:
    scored = []
    for ct, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            log.warning("annotation: skipping %s (only %d/%d markers in data)", ct, len(present), len(genes))
            continue
        sc.tl.score_genes(adata, present, score_name=f"score_{ct}", random_state=seed, use_raw=False)
        scored.append(ct)
    return scored


def label_clusters(adata: ad.AnnData, cell_types: list[str], min_margin: float,
                   coarse: dict[str, str], cluster_key: str = "leiden") -> pd.DataFrame:
    cols = [f"score_{ct}" for ct in cell_types]
    means = adata.obs.groupby(cluster_key, observed=True)[cols].mean()
    means.columns = cell_types
    ranked = np.argsort(-means.to_numpy(), axis=1)
    top = means.columns[ranked[:, 0]]
    second = means.columns[ranked[:, 1]] if len(cell_types) > 1 else top
    top_score = means.to_numpy()[np.arange(len(means)), ranked[:, 0]]
    second_score = means.to_numpy()[np.arange(len(means)), ranked[:, 1]] if len(cell_types) > 1 else top_score
    table = pd.DataFrame({
        "cluster": means.index.astype(str),
        "n_cells": adata.obs[cluster_key].value_counts().reindex(means.index).to_numpy(),
        "cell_type": top, "score": top_score,
        "runner_up": second, "runner_up_score": second_score,
    })
    table["margin"] = table["score"] - table["runner_up_score"]
    table["confident"] = table["margin"] >= min_margin
    table["cell_type_coarse"] = table["cell_type"].map(lambda c: coarse.get(c, c))

    lab = dict(zip(table["cluster"], table["cell_type"]))
    lab_c = dict(zip(table["cluster"], table["cell_type_coarse"]))
    conf = dict(zip(table["cluster"], table["confident"]))
    cl = adata.obs[cluster_key].astype(str)
    adata.obs["cell_type"] = pd.Categorical(cl.map(lab))
    adata.obs["cell_type_coarse"] = pd.Categorical(cl.map(lab_c))
    adata.obs["annotation_confident"] = cl.map(conf).astype(bool)
    return table.sort_values("n_cells", ascending=False)


def composition(adata: ad.AnnData, by: str = "cell_type") -> pd.DataFrame:
    """Cells per (sample, cell type) with within-sample fractions, plus sample metadata."""
    obs = adata.obs
    counts = obs.groupby(["sample", by], observed=True).size().rename("n_cells").reset_index()
    counts["fraction"] = counts["n_cells"] / counts.groupby("sample")["n_cells"].transform("sum")
    meta = obs.groupby("sample", observed=True)[["treatment", "mouse", "Gram", "study"]].first().reset_index()
    return counts.merge(meta, on="sample").rename(columns={by: "cell_type"})


def marker_dotplot(adata: ad.AnnData, markers: dict[str, list[str]], cell_types: list[str], path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    present_types = [ct for ct in cell_types if ct in set(adata.obs["cell_type"].astype(str))]
    var_names = {ct: [g for g in markers[ct] if g in adata.var_names] for ct in present_types}
    var_names = {k: v for k, v in var_names.items() if v}
    if not var_names:
        return
    adata.obs["cell_type"] = adata.obs["cell_type"].cat.set_categories(
        [c for c in var_names if c in adata.obs["cell_type"].cat.categories])
    dp = sc.pl.dotplot(adata, var_names, groupby="cell_type", standard_scale="var",
                       show=False, return_fig=True)
    dp.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


def annotate(adata: ad.AnnData, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = cfg["annotation"]
    cell_types = score_cell_types(adata, a["markers"], cfg["processing"]["seed"])
    table = label_clusters(adata, cell_types, a["min_margin"], a.get("coarse", {}))
    return table, composition(adata)
