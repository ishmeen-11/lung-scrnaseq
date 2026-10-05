"""End-to-end tests on the synthetic demo data. Run: pytest -q"""
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest
import yaml

from lungsc import cli, qc
from lungsc.config import load_config, load_samples
from lungsc.demo import make_demo


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    root = tmp_path_factory.mktemp("demo")
    cfg_path = make_demo(root, n_cells=400, seed=1)
    for s in load_samples(load_config(cfg_path)).index:
        cli.main(["qc", "--config", str(cfg_path), "--sample", s])
    cli.main(["integrate", "--config", str(cfg_path)])
    return root, cfg_path


def test_outputs_exist(demo):
    root, _ = demo
    res = root / "results"
    assert (res / "integrated.h5ad").exists()
    for t in ["sample_summary.parquet", "umap.parquet"]:
        assert (res / "qc" / t).exists()


def test_funnel_is_monotonic(demo):
    root, _ = demo
    for f in (root / "results" / "qc" / "funnel").glob("*.parquet"):
        df = pd.read_parquet(f)
        assert df["cells_remaining"].is_monotonic_decreasing
        assert df["neutrophils_remaining"].is_monotonic_decreasing


def test_integrated_object(demo):
    root, _ = demo
    a = ad.read_h5ad(root / "results" / "integrated.h5ad")
    for col in ["sample", "study", "treatment", "mouse", "bacSpecies", "Gram", "leiden"]:
        assert col in a.obs
    assert "counts" in a.layers and "X_umap" in a.obsm
    assert a.obs_names.is_unique
    # QC actually removed the damaged / RBC cells we simulated
    assert a.obs["pct_counts_mt"].max() < 15
    assert a.obs["pct_counts_hb"].max() < 5


def test_sex_check(demo):
    root, _ = demo
    s = pd.read_parquet(root / "results" / "qc" / "sample_summary.parquet").set_index("sample")
    assert s.loc["Kp_WT", "sex_call"] == "male"
    assert s.loc["Sa_WT", "sex_call"] == "female"


def test_fixed_mode_matches_griffin(demo):
    _, cfg_path = demo
    cfg = load_config(cfg_path)
    cfg["qc"]["mode"] = "fixed"
    obs = pd.DataFrame({"total_counts": [500, 1500], "n_genes_by_counts": [300, 300],
                        "pct_counts_mt": [1, 1], "log1p_total_counts": [0, 0],
                        "log1p_n_genes_by_counts": [0, 0]})
    t = qc.compute_thresholds(obs, cfg)
    assert t["min_counts"] == 1000 and t["max_pct_mt"] == 7.5 and t["min_genes"] == 200
