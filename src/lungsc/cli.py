"""Command-line entry points. Snakemake calls these; you can also run them by hand.

    python -m lungsc.cli qc --config config/config.yaml --sample Sa_WT
    python -m lungsc.cli integrate --config config/config.yaml
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from . import config as C
from . import io, process, qc

log = logging.getLogger("lungsc")


def _out(cfg, *parts) -> Path:
    p = C.resolve(cfg, cfg["results_dir"]).joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def cmd_qc(args):
    cfg = C.load_config(args.config)
    samples = C.load_samples(cfg)
    row = samples.loc[args.sample]
    sc.settings.verbosity = 1

    adata = io.read_counts(C.resolve(cfg, row["counts_path"]), row["input_format"])
    log.info("%s: %d cells x %d genes read", args.sample, adata.n_obs, adata.n_vars)
    filtered, cells, fun, thr = qc.run_sample_qc(adata, cfg, args.sample)

    filtered.write_h5ad(_out(cfg, "per_sample", f"{args.sample}.h5ad"), compression="gzip")
    cells.to_parquet(_out(cfg, "qc", "cells", f"{args.sample}.parquet"), index=False)
    fun.to_parquet(_out(cfg, "qc", "funnel", f"{args.sample}.parquet"), index=False)
    thr.to_parquet(_out(cfg, "qc", "thresholds", f"{args.sample}.parquet"), index=False)
    log.info("%s: kept %d / %d cells", args.sample, filtered.n_obs, adata.n_obs)


def cmd_integrate(args):
    cfg = C.load_config(args.config)
    samples = C.load_samples(cfg)
    names = list(samples.index)

    adatas = {s: ad.read_h5ad(_out(cfg, "per_sample", f"{s}.h5ad")) for s in names}
    meta = {s: C.sample_metadata(cfg, samples.loc[s]) for s in names}
    adata = process.merge(adatas, meta, cfg)
    process.normalize_and_reduce(adata, cfg)
    rep = process.integrate(adata, cfg)
    process.cluster(adata, rep, cfg)

    cells = pd.concat([pd.read_parquet(_out(cfg, "qc", "cells", f"{s}.parquet")) for s in names])
    summary = process.sample_summary(adata, cells, cfg)
    summary.to_parquet(_out(cfg, "qc", "sample_summary.parquet"), index=False)

    emb = pd.DataFrame(adata.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adata.obs_names)
    cols = ["sample", "treatment", "mouse", "Gram", "leiden", "total_counts",
            "n_genes_by_counts", "pct_counts_mt", "neutrophil_score"]
    emb = emb.join(adata.obs[cols]).reset_index(names="cell")
    emb.to_parquet(_out(cfg, "qc", "umap.parquet"), index=False)

    adata.write_h5ad(_out(cfg, "integrated.h5ad"), compression="gzip")
    log.info("integrated: %d cells, %d genes, %d clusters",
             adata.n_obs, adata.n_vars, adata.obs["leiden"].nunique())


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    ap = argparse.ArgumentParser(prog="lungsc")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("qc", help="QC + doublet calling for one sample")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("integrate", help="merge all samples, normalize, integrate, cluster")
    p.add_argument("--config", required=True)
    p.set_defaults(func=cmd_integrate)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
