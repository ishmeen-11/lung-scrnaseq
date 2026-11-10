# lung-scrnaseq

Reproducible preprocessing for the Prince lab 3-species lung scRNA-seq data
(*S. aureus*, *P. aeruginosa* PAO1/CS17, *K. pneumoniae*, PBS; WT mice. The *Irg1*-/- samples
are listed in `config/samples_with_irg1.tsv` but excluded from the current run).
Python / Scanpy, run with Snakemake (laptop or Ginsburg/Slurm), with a Streamlit
QC viewer.

```
raw FASTQ ─► Cell Ranger ─► per-sample QC + doublets ─► merge ─► normalize ─► Harmony ─► Leiden/UMAP
                              │                                                          │
                              └──────── QC tables (parquet) ──► Streamlit viewer ◄───────┘
```

## Quick start (no real data needed)

```bash
conda env create -f environment.yml && conda activate lungsc && pip install -e .

python -m lungsc.demo --out demo                        # synthetic 4-sample dataset
snakemake --cores 2 --configfile demo/config/config.yaml
streamlit run app/streamlit_app.py -- --results demo/results
pytest -q
```

## Real data

0. Unzip Griffin's `SoupX_Outs.zip` into `data/` so you have `data/SoupX_Outs/<sample>_strained/`.
   `config/samples.tsv` already points there, with GEO/SRA accessions, dose and timepoint filled in.
1. To add new samples, edit `config/samples.tsv`: one row per sample. `input_format` is one of
   `cellranger` (filtered `.h5` or mtx folder), `cellbender`, `soupx_mtx`
   (Griffin's SoupX folders) or `h5ad`. Add SRA/ENA run IDs to `run_accessions` to
   download FASTQs and run Cell Ranger automatically.
2. Set paths in `config/config.yaml` (Cell Ranger binary, reference, Slurm account).
3. Run:
   ```bash
   snakemake --profile workflow/profiles/laptop     # 8 GB laptop: 1 job at a time, 4,000 cells/sample
   snakemake --profile workflow/profiles/insomnia   # Columbia Insomnia (set slurm_account first)
   snakemake --profile workflow/profiles/slurm   # Ginsburg
   snakemake --cores 8                           # local, if counts already exist
   snakemake --cores 8 qc_all                    # QC only, to tune thresholds first
   ```
4. Copy `results/qc/` to your laptop and open the viewer.

## Full run on AWS EC2

An 8 GB laptop can only merge a 4,000-cell subsample per sample. For the full run, one
r6i.xlarge (32 GB, 4 vCPU, roughly $0.25/hr) is enough for the 6 WT samples (~30k cells, ~5 min):

1. Launch an Ubuntu 24.04 instance with an SSH key and a security group open to your IP on port 22.
2. `rsync` the repo with `data/SoupX_Outs/` to it (exclude `results/` and `.snakemake/`).
3. On the instance: install Miniforge, `conda env create -f environment.yml`, `conda activate lungsc`,
   `pip install -e .`, then `snakemake --cores 4 --rerun-incomplete` inside `tmux`.
4. `rsync` `results/` back to `results_full/` on your laptop, then **terminate the instance**.

## Streamlit viewer

```bash
streamlit run app/streamlit_app.py --server.fileWatcherType none -- --results results_full
```

With no `--results` folder present (e.g. on Streamlit Community Cloud), the app shows the QC
tables from the full run bundled in `app/bundled_results/` (parquet only, no expression data).
Deploy on Streamlit Cloud with main file `app/streamlit_app.py`; it installs `app/requirements.txt`.

## Outputs (`results/`)

| File | What |
|---|---|
| `per_sample/{sample}.h5ad` | QC'd, singlet cells, raw counts |
| `integrated.h5ad` | merged object: `layers["counts"]`, log-normalized `X`, `X_pca_harmony`, `X_umap`, `obs.leiden` + metadata (`sample, study, treatment, mouse, bacSpecies, Gram`) |
| `qc/cells/*.parquet` | every cell before filtering, with metrics and pass/fail per filter |
| `qc/funnel/*.parquet` | cells (and putative neutrophils) left after each filter |
| `qc/thresholds/*.parquet` | thresholds actually applied per sample |
| `qc/sample_summary.parquet` | depth, retention, doublet rate, sex call per sample |
| `annotation/cluster_labels.parquet` | cell type per Leiden cluster (Griffin's markers), with confidence margin |
| `annotation/composition.parquet` | cells and fractions per sample × cell type |
| `annotation/marker_dotplot.png` | marker expression by assigned cell type |

To open the result in R/Seurat: `zellkonverter::readH5AD("results/integrated.h5ad")`.

## Layout

```
config/      config.yaml (all paths + thresholds), samples.tsv
workflow/    Snakefile, Slurm profile
src/lungsc/  io, qc, process, cli, demo
app/         Streamlit QC viewer
tests/       end-to-end test on synthetic data
docs/        decisions vs. the original R pipeline, open items
```

See `docs/methods.md` for how each step maps to the original R/Seurat pipeline.
