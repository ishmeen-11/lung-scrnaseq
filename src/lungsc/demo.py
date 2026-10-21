"""Generate a small synthetic dataset + config so the pipeline and the Streamlit app
can be tried without the real data.

    python -m lungsc.demo --out demo
    snakemake --cores 2 --configfile demo/config/config.yaml
    streamlit run app/streamlit_app.py -- --results demo/results
"""
from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import io as sio
from scipy import sparse

NEU = ["Csf3r", "S100a8", "S100a9", "Retnlg", "Cxcr2", "Mmp9", "Lcn2"]
MAC = ["Siglecf", "Pparg", "Ear1", "Krt79", "Car4"]  # alveolar macrophage markers
BCELL = ["Cd79a", "Cd79b", "Cd19", "Pax5", "Ly6d"]
EC = ["Cldn5", "Ramp2", "Pecam1", "Lyve1"]
MT = [f"mt-Gene{i}" for i in range(13)]
HB = ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"]
SEX = ["Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty"]

SAMPLES = [  # name, treatment, mouse, study, sex
    ("Sa_WT", "Sa", "WT", "Tomlinson2023_Sa", "female"),
    ("Sa_IRG1", "Sa", "Irg1", "Tomlinson2023_Sa", "female"),
    ("Kp_WT", "Kp", "WT", "WongFokLung2022_Kp", "male"),
    ("PBS_WT_Tom", "PBS", "WT", "Tomlinson2023_shared", "female"),
]


def _genes(n_other=1500):
    marker = NEU + MAC + BCELL + EC + MT + HB + SEX
    return marker + [f"Gene{i}" for i in range(n_other)]


def _simulate(rng, genes, n_cells, sex, depth_scale, shared_base):
    g = len(genes)
    idx = {x: i for i, x in enumerate(genes)}
    # same tissue in every sample, plus a mild per-sample (batch) perturbation
    base = shared_base * rng.lognormal(0, 0.15, g)
    for x in MT:
        base[idx[x]] = 3.0
    for x in HB + SEX:
        base[idx[x]] = 0.0
    types = rng.choice(["neu", "mac", "b", "ec"], size=n_cells, p=[0.45, 0.2, 0.2, 0.15])
    rows = []
    for t in types:
        mu = base.copy()
        genes_t = {"neu": NEU, "mac": MAC, "b": BCELL, "ec": EC}[t]
        for x in genes_t:
            mu[idx[x]] += 40
        mu[idx["Xist" if sex == "female" else "Ddx3y"]] = 2.0
        if sex == "male":
            for x in ["Eif2s3y", "Kdm5d", "Uty"]:
                mu[idx[x]] = 1.0
        # neutrophils carry much less RNA: this is what the fixed 1,000-UMI cutoff hits
        depth = rng.lognormal(np.log(700 if t == "neu" else 3500), 0.45) * depth_scale
        rows.append(rng.poisson(mu / mu.sum() * depth))
    X = np.vstack(rows)
    n = n_cells
    dmg = rng.choice(n, n // 25, replace=False)          # damaged cells: high mito
    for i in dmg:
        X[i, [idx[x] for x in MT]] += rng.poisson(80, len(MT))
    rbc = rng.choice(n, n // 40, replace=False)          # RBC contamination
    for i in rbc:
        X[i, [idx[x] for x in HB]] += rng.poisson(60, len(HB))
    return sparse.csr_matrix(X)


def _write_10x(outdir: Path, X, genes):
    outdir.mkdir(parents=True, exist_ok=True)
    sio.mmwrite(str(outdir / "matrix.mtx"), X.T.tocoo().astype(np.int32))
    bcs = [f"AAAC{i:012d}-1" for i in range(X.shape[0])]
    pd.Series(bcs).to_csv(outdir / "barcodes.tsv", index=False, header=False)
    pd.DataFrame({"id": [f"ENSMUSG{i:011d}" for i in range(len(genes))], "name": genes,
                  "type": "Gene Expression"}).to_csv(outdir / "features.tsv", sep="\t",
                                                      index=False, header=False)
    for f in ["matrix.mtx", "barcodes.tsv", "features.tsv"]:
        with open(outdir / f, "rb") as a, gzip.open(outdir / f"{f}.gz", "wb") as b:
            shutil.copyfileobj(a, b)
        (outdir / f).unlink()


def make_demo(out: str | Path, n_cells: int = 900, seed: int = 0) -> Path:
    out = Path(out).resolve()
    rng = np.random.default_rng(seed)
    genes = _genes()
    shared_base = rng.gamma(0.3, 1.0, len(genes))
    rows = []
    for k, (name, trt, mouse, study, sex) in enumerate(SAMPLES):
        X = _simulate(rng, genes, n_cells, sex, depth_scale=1.0 + 0.3 * k, shared_base=shared_base)
        d = out / "data" / "counts" / name
        _write_10x(d, X, genes)
        rows.append({"sample": name, "study": study, "treatment": trt, "mouse": mouse,
                     "run_accessions": "demo", "input_format": "soupx_mtx",
                     "counts_path": f"data/counts/{name}"})
    (out / "config").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "config" / "samples.tsv", sep="\t", index=False)

    repo = Path(__file__).resolve().parents[2]
    cfg = yaml.safe_load(open(repo / "config" / "config.yaml"))
    cfg["samples"] = "config/samples.tsv"
    cfg["processing"].update({"n_top_genes": 500, "n_pcs": 20, "neighbor_pcs": 15, "n_neighbors": 15})
    with open(out / "config" / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    return out / "config" / "config.yaml"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="demo")
    ap.add_argument("--n-cells", type=int, default=900)
    a = ap.parse_args()
    print(make_demo(a.out, a.n_cells))
