"""Per-sample QC: metrics, filters, diagnostics, and a funnel of what each step removed.

Design rule: every cell gets a pass/fail flag per filter and is written to a QC
table, so the Streamlit app can show exactly what each filter removes. Cells are
only dropped from the AnnData at the end.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# ---------- metrics -----------------------------------------------------------

def add_metrics(adata, cfg: dict) -> None:
    q = cfg["qc"]
    adata.var["mt"] = adata.var_names.str.startswith(q["mt_prefix"])
    adata.var["hb"] = adata.var_names.str.match(q["hb_regex"])
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], percent_top=None,
                               log1p=True, inplace=True)


def _gene_mask(adata, genes):
    present = [g for g in genes if g in adata.var_names]
    return present


def _sum_genes(adata, genes) -> np.ndarray:
    present = _gene_mask(adata, genes)
    if not present:
        return np.zeros(adata.n_obs)
    X = adata[:, present].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return X.sum(axis=1)


def add_diagnostics(adata, cfg: dict) -> None:
    """Neutrophil score + sex genes. Diagnostics only, never used to filter."""
    c = cfg["checks"]
    # Neutrophil score: fraction of a cell's counts coming from neutrophil markers.
    # Normalising by library size matters here because neutrophils are low-count cells.
    neu = _sum_genes(adata, c["neutrophil_markers"])
    lib = np.asarray(adata.obs["total_counts"]).clip(min=1)
    adata.obs["neutrophil_score"] = neu / lib
    thr = np.quantile(adata.obs["neutrophil_score"], c["neutrophil_score_quantile"])
    adata.obs["putative_neutrophil"] = (adata.obs["neutrophil_score"] > thr) & (neu > 0)

    adata.obs["female_counts"] = _sum_genes(adata, c["female_genes"])
    adata.obs["male_counts"] = _sum_genes(adata, c["male_genes"])


# ---------- filters -----------------------------------------------------------

def _mad_bounds(x: pd.Series, nmads: float) -> tuple[float, float]:
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826  # consistent with SD for normal data
    return med - nmads * mad, med + nmads * mad


def compute_thresholds(obs: pd.DataFrame, cfg: dict) -> dict:
    q = cfg["qc"]
    if q["mode"] == "fixed":
        f = q["fixed"]
        return {
            "min_counts": f["min_counts"], "max_counts": np.inf,
            "min_genes": max(f["min_genes"], q["hard_min_genes"]), "max_genes": np.inf,
            "max_pct_mt": f["max_pct_mt"],
        }
    if q["mode"] == "mad":
        m = q["mad"]
        lo_c, hi_c = _mad_bounds(obs["log1p_total_counts"], m["nmads_counts"])
        lo_g, hi_g = _mad_bounds(obs["log1p_n_genes_by_counts"], m["nmads_genes"])
        _, hi_mt = _mad_bounds(obs["pct_counts_mt"], m["nmads_mt"])
        return {
            "min_counts": float(np.expm1(lo_c)), "max_counts": float(np.expm1(hi_c)),
            "min_genes": max(float(np.expm1(lo_g)), q["hard_min_genes"]),
            "max_genes": float(np.expm1(hi_g)),
            "max_pct_mt": float(min(hi_mt, m["max_pct_mt_ceiling"])),
        }
    raise ValueError(f"qc.mode must be 'fixed' or 'mad', got {q['mode']!r}")


def flag_cells(adata, cfg: dict) -> dict:
    """Add one boolean column per filter (True = cell passes). Returns thresholds used."""
    t = compute_thresholds(adata.obs, cfg)
    o = adata.obs
    o["pass_counts"] = (o["total_counts"] > t["min_counts"]) & (o["total_counts"] < t["max_counts"])
    o["pass_genes"] = (o["n_genes_by_counts"] > t["min_genes"]) & (o["n_genes_by_counts"] < t["max_genes"])
    o["pass_mt"] = o["pct_counts_mt"] < t["max_pct_mt"]
    max_hb = cfg["qc"].get("max_pct_hb")
    o["pass_hb"] = (o["pct_counts_hb"] < max_hb) if max_hb is not None else True
    o["pass_qc"] = o["pass_counts"] & o["pass_genes"] & o["pass_mt"] & o["pass_hb"]
    t["max_pct_hb"] = max_hb if max_hb is not None else np.inf
    return t


# ---------- doublets ----------------------------------------------------------

def flag_doublets(adata, cfg: dict) -> None:
    """Run Scrublet on QC-passing cells only; failing cells get NaN / False."""
    d = cfg["doublets"]
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    keep = adata.obs["pass_qc"].to_numpy()
    if keep.sum() < 50:  # too few cells for a meaningful simulation
        return
    sub = adata[keep].copy()
    sc.pp.scrublet(sub, expected_doublet_rate=d["expected_rate"],
                   random_state=cfg["processing"]["seed"], verbose=False)
    adata.obs.loc[sub.obs_names, "doublet_score"] = sub.obs["doublet_score"].to_numpy()
    adata.obs.loc[sub.obs_names, "predicted_doublet"] = sub.obs["predicted_doublet"].astype(bool).to_numpy()


# ---------- funnel ------------------------------------------------------------

FUNNEL_STEPS = [
    ("start", None),
    ("counts", "pass_counts"),
    ("genes", "pass_genes"),
    ("mito", "pass_mt"),
    ("hemoglobin", "pass_hb"),
    ("doublets", "singlet"),
]


def funnel(obs: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Cells left after each filter, applied in order, overall and for putative neutrophils."""
    obs = obs.assign(singlet=~obs["predicted_doublet"].astype(bool))
    keep = pd.Series(True, index=obs.index)
    rows = []
    for step, col in FUNNEL_STEPS:
        if col is not None:
            keep &= obs[col].astype(bool)
        rows.append({
            "sample": sample, "step": step,
            "cells_remaining": int(keep.sum()),
            "neutrophils_remaining": int((keep & obs["putative_neutrophil"]).sum()),
        })
    df = pd.DataFrame(rows)
    df["cells_removed"] = -df["cells_remaining"].diff().fillna(0).astype(int)
    df["pct_neutrophils_lost"] = 100 * (1 - df["neutrophils_remaining"] / max(df["neutrophils_remaining"].iloc[0], 1))
    return df


# ---------- driver ------------------------------------------------------------

QC_COLUMNS = [
    "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb",
    "neutrophil_score", "putative_neutrophil", "female_counts", "male_counts",
    "pass_counts", "pass_genes", "pass_mt", "pass_hb", "pass_qc",
    "doublet_score", "predicted_doublet", "keep",
]


def run_sample_qc(adata, cfg: dict, sample: str):
    """Return (filtered AnnData, per-cell QC table, funnel table, thresholds)."""
    adata.layers["counts"] = adata.X.copy()
    add_metrics(adata, cfg)
    add_diagnostics(adata, cfg)
    thresholds = flag_cells(adata, cfg)
    flag_doublets(adata, cfg)

    remove_dbl = cfg["doublets"]["remove"]
    adata.obs["keep"] = adata.obs["pass_qc"] & (~adata.obs["predicted_doublet"] if remove_dbl else True)

    cells = adata.obs[QC_COLUMNS].copy()
    cells.insert(0, "sample", sample)
    cells.index.name = "barcode"
    cells = cells.reset_index()
    fun = funnel(adata.obs, sample)

    filtered = adata[adata.obs["keep"].to_numpy()].copy()
    thr = pd.DataFrame([{"sample": sample, "mode": cfg["qc"]["mode"], **thresholds}])
    return filtered, cells, fun, thr
