"""Reading count matrices in the formats we expect to meet."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import scanpy as sc

SUPPORTED = ("cellranger", "cellbender", "soupx_mtx", "h5ad")


def read_counts(path: str | Path, input_format: str) -> ad.AnnData:
    """Return raw integer-ish counts, cells x genes, unique gene symbols.

    cellranger : filtered_feature_bc_matrix.h5, or the matching mtx directory
    cellbender : CellBender *_filtered.h5 (ambient RNA already removed)
    soupx_mtx  : Griffin's SoupX output folders (DropletUtils::write10xCounts)
    h5ad       : anything already converted
    """
    path = Path(path)
    if input_format not in SUPPORTED:
        raise ValueError(f"input_format must be one of {SUPPORTED}, got {input_format!r}")
    if not path.exists():
        raise FileNotFoundError(path)

    if input_format == "h5ad":
        adata = ad.read_h5ad(path)
    elif path.is_dir():
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    else:
        adata = sc.read_10x_h5(path)  # works for Cell Ranger and CellBender h5

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata
