"""Config + sample sheet loading. Single source of truth for paths and thresholds."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

REQUIRED_SAMPLE_COLS = ["sample", "study", "treatment", "mouse", "input_format", "counts_path"]


def load_config(path: str | Path) -> dict:
    path = Path(path)
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_root"] = str(path.resolve().parent.parent)  # repo root = parent of config/
    return cfg


def resolve(cfg: dict, p: str | Path) -> Path:
    """Resolve a config path relative to the repo root."""
    p = Path(p)
    return p if p.is_absolute() else Path(cfg["_root"]) / p


def load_samples(cfg: dict) -> pd.DataFrame:
    df = pd.read_csv(resolve(cfg, cfg["samples"]), sep="\t", dtype=str).fillna("")
    missing = [c for c in REQUIRED_SAMPLE_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"samples.tsv is missing columns: {missing}")
    if df["sample"].duplicated().any():
        raise ValueError("samples.tsv has duplicate sample names")
    unknown = set(df["treatment"]) - set(cfg["treatment_metadata"])
    if unknown:
        raise ValueError(f"treatments not in config.treatment_metadata: {unknown}")
    return df.set_index("sample", drop=False)


def sample_metadata(cfg: dict, row: pd.Series) -> dict:
    """All per-sample labels that get attached to every cell."""
    meta = {
        "sample": row["sample"],
        "study": row["study"],
        "treatment": row["treatment"],
        "mouse": row["mouse"],
    }
    meta.update(cfg["treatment_metadata"][row["treatment"]])
    return meta
