"""Streamlit QC viewer for the lung scRNA-seq preprocessing.

    streamlit run app/streamlit_app.py -- --results results
    streamlit run app/streamlit_app.py -- --results demo/results

Reads only the parquet tables written by the pipeline, so it runs on a laptop
even when the h5ad files live on the HPC (just copy results/qc/).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

alt.data_transformers.disable_max_rows()
st.set_page_config(page_title="Lung scRNA-seq QC", layout="wide")

MAX_POINTS = 30_000  # downsample for plotting only


def _results_dir() -> Path:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.environ.get("LUNGSC_RESULTS", "results"))
    args, _ = ap.parse_known_args(sys.argv[1:])
    return Path(args.results)


@st.cache_data
def load(results: str):
    qc = Path(results) / "qc"
    read_dir = lambda d: pd.concat([pd.read_parquet(f) for f in sorted((qc / d).glob("*.parquet"))],
                                   ignore_index=True)
    cells = read_dir("cells")
    funnel = read_dir("funnel")
    thr = read_dir("thresholds")
    summary = pd.read_parquet(qc / "sample_summary.parquet") if (qc / "sample_summary.parquet").exists() else None
    umap = pd.read_parquet(qc / "umap.parquet") if (qc / "umap.parquet").exists() else None
    ann = Path(results) / "annotation"
    labels = pd.read_parquet(ann / "cluster_labels.parquet") if (ann / "cluster_labels.parquet").exists() else None
    comp = pd.read_parquet(ann / "composition.parquet") if (ann / "composition.parquet").exists() else None
    return cells, funnel, thr, summary, umap, labels, comp


RESULTS = _results_dir()
if not (RESULTS / "qc" / "cells").exists():
    st.error(f"No QC tables found under {RESULTS}/qc. Run the pipeline first (see README).")
    st.stop()

cells, funnel, thr, summary, umap, labels, comp = load(str(RESULTS))
samples = sorted(cells["sample"].unique())

st.title("Lung scRNA-seq · preprocessing QC")
st.caption(f"Results: `{RESULTS}` · {len(samples)} samples · {len(cells):,} cells before QC · "
           f"QC mode: **{thr['mode'].iloc[0]}**")

tab_over, tab_funnel, tab_dist, tab_explore, tab_sex, tab_umap, tab_ct = st.tabs(
    ["Overview", "Filter funnel", "Distributions", "Threshold explorer", "Sex check", "UMAP", "Cell types"])

# ---------------------------------------------------------------- overview
with tab_over:
    if summary is not None:
        c1, c2, c3 = st.columns(3)
        c1.metric("Cells before QC", f"{summary['cells_raw'].sum():,}")
        c2.metric("Cells kept", f"{summary['cells_kept'].sum():,}",
                  f"{100 * summary['cells_kept'].sum() / summary['cells_raw'].sum():.1f}% kept")
        c3.metric("Median UMIs / cell (range across samples)",
                  f"{summary['median_counts_raw'].min():,.0f} – {summary['median_counts_raw'].max():,.0f}")
        st.markdown("**Per-sample summary.** A large spread in median UMIs between studies "
                    "means sequencing depth differs by lab and can look like biology.")
        st.dataframe(summary.round(1), width="stretch", hide_index=True)
        st.altair_chart(
            alt.Chart(summary).mark_bar().encode(
                x=alt.X("sample:N", sort=samples), y=alt.Y("median_counts_raw:Q", title="median UMIs/cell"),
                color="study:N", tooltip=list(summary.columns)).properties(height=260),
            width="stretch")
    else:
        st.info("Integration step not run yet: per-sample summary appears after `integrate`.")

# ---------------------------------------------------------------- funnel
with tab_funnel:
    st.markdown("Cells remaining after each filter, applied in order. "
                "**Neutrophils lost** tracks putative neutrophils (high neutrophil-marker score), "
                "since they are low-RNA cells and the project's main finding.")
    order = ["start", "counts", "genes", "mito", "hemoglobin", "doublets"]
    sel = st.multiselect("Samples", samples, default=samples, key="fun_s")
    f = funnel[funnel["sample"].isin(sel)]
    st.altair_chart(
        alt.Chart(f).mark_line(point=True).encode(
            x=alt.X("step:N", sort=order), y="cells_remaining:Q", color="sample:N",
            tooltip=list(f.columns)).properties(height=320), width="stretch")
    last = f[f["step"] == "doublets"][["sample", "cells_remaining", "pct_neutrophils_lost"]]
    st.dataframe(last.rename(columns={"cells_remaining": "cells_kept"}).round(1),
                 hide_index=True, width="stretch")
    with st.expander("Full funnel table"):
        st.dataframe(f.round(1), hide_index=True, width="stretch")

# ---------------------------------------------------------------- distributions
METRICS = {
    "total_counts": ("UMIs per cell", True, ["min_counts", "max_counts"]),
    "n_genes_by_counts": ("Genes per cell", True, ["min_genes", "max_genes"]),
    "pct_counts_mt": ("% mitochondrial", False, ["max_pct_mt"]),
    "pct_counts_hb": ("% hemoglobin", False, ["max_pct_hb"]),
}

with tab_dist:
    metric = st.selectbox("Metric", list(METRICS), format_func=lambda m: METRICS[m][0])
    label, logx, tcols = METRICS[metric]
    d = cells[["sample", metric, "putative_neutrophil"]]
    if len(d) > MAX_POINTS:
        d = d.sample(MAX_POINTS, random_state=0)
    d = d.assign(value=np.log10(d[metric] + 1) if logx else d[metric])
    base = alt.Chart(d).transform_density(
        "value", groupby=["sample"], as_=["value", "density"]).mark_area(opacity=0.45).encode(
        x=alt.X("value:Q", title=f"log10({label} + 1)" if logx else label), y="density:Q",
        color="sample:N")
    rules = thr.melt(id_vars="sample", value_vars=tcols).replace([np.inf, -np.inf], np.nan).dropna()
    rules["value"] = np.log10(rules["value"] + 1) if logx else rules["value"]
    rule_chart = alt.Chart(rules).mark_rule(strokeDash=[4, 3]).encode(
        x="value:Q", color="sample:N", tooltip=["sample", "variable", "value"])
    st.altair_chart((base + rule_chart).properties(height=360), width="stretch")
    st.caption("Dashed lines = thresholds actually used for each sample.")

# ---------------------------------------------------------------- explorer
with tab_explore:
    st.markdown("Try thresholds and see what they would remove, **before** rerunning the pipeline. "
                "Doublet calls are taken from the pipeline run.")
    c1, c2, c3, c4 = st.columns(4)
    min_counts = c1.number_input("min UMIs", 0, 20000, 1000, step=100)
    min_genes = c2.number_input("min genes", 0, 5000, 200, step=50)
    max_mt = c3.number_input("max % mito", 0.0, 100.0, 7.5, step=0.5)
    max_hb = c4.number_input("max % hemoglobin", 0.0, 100.0, 5.0, step=0.5)
    st.caption("Defaults = Griffin's original fixed cutoffs (+ hemoglobin).")

    keep = ((cells["total_counts"] > min_counts) & (cells["n_genes_by_counts"] > min_genes)
            & (cells["pct_counts_mt"] < max_mt) & (cells["pct_counts_hb"] < max_hb)
            & ~cells["predicted_doublet"].astype(bool))
    e = cells.assign(keep_try=keep)
    g = e.groupby("sample")
    res = pd.DataFrame({
        "cells_raw": g.size(),
        "kept (these settings)": g["keep_try"].sum(),
        "kept (pipeline run)": g["keep"].sum(),
        "neutrophils_raw": g["putative_neutrophil"].sum(),
        "neutrophils kept (these settings)": e[e["putative_neutrophil"]].groupby("sample")["keep_try"].sum(),
    }).fillna(0).astype(int)
    res["% neutrophils lost (these settings)"] = (
        100 * (1 - res["neutrophils kept (these settings)"] / res["neutrophils_raw"].clip(lower=1))).round(1)
    st.dataframe(res.reset_index(), hide_index=True, width="stretch")

    sub = e if len(e) <= MAX_POINTS else e.sample(MAX_POINTS, random_state=0)
    sub = sub.assign(log_counts=np.log10(sub["total_counts"] + 1))
    st.altair_chart(
        alt.Chart(sub).mark_circle(size=6, opacity=0.5).encode(
            x=alt.X("log_counts:Q", title="log10(UMIs + 1)"), y=alt.Y("pct_counts_mt:Q", title="% mito"),
            color=alt.Color("keep_try:N", title="kept",
                            scale=alt.Scale(domain=[True, False], range=["#4c78a8", "#e45756"])), tooltip=["sample", "total_counts", "pct_counts_mt"],
            facet=alt.Facet("sample:N", columns=4)).properties(width=220, height=180))

# ---------------------------------------------------------------- sex check
with tab_sex:
    st.markdown("Xist = female; Ddx3y / Eif2s3y / Kdm5d / Uty = male (counts per million). "
                "If studies used different sexes, sex is confounded with pathogen.")
    if summary is not None:
        s = summary[["sample", "study", "female_cpm", "male_cpm", "sex_call"]]
        st.dataframe(s.round(1), hide_index=True, width="stretch")
        m = s.melt(id_vars=["sample", "sex_call"], value_vars=["female_cpm", "male_cpm"])
        st.altair_chart(alt.Chart(m).mark_bar().encode(
            x=alt.X("sample:N", sort=samples), y="value:Q", color="variable:N", xOffset="variable:N")
            .properties(height=280), width="stretch")
    else:
        st.info("Available after the integrate step.")

# ---------------------------------------------------------------- umap
with tab_umap:
    if umap is None:
        st.info("Available after the integrate step.")
    else:
        opts = [c for c in ["cell_type", "cell_type_coarse", "leiden", "sample", "treatment", "mouse", "Gram",
                            "neutrophil_score", "pct_counts_mt", "total_counts"] if c in umap.columns]
        color = st.selectbox("Color by", opts)
        u = umap if len(umap) <= MAX_POINTS else umap.sample(MAX_POINTS, random_state=0)
        ctype = "Q" if pd.api.types.is_numeric_dtype(u[color]) else "N"
        sort = sorted(u[color].astype(str).unique(), key=lambda v: (len(v), v)) if ctype == "N" else None
        st.altair_chart(alt.Chart(u).mark_circle(size=5, opacity=0.6).encode(
            x=alt.X("umap_1:Q", axis=None), y=alt.Y("umap_2:Q", axis=None),
            color=alt.Color(f"{color}:{ctype}", sort=sort), tooltip=[c for c in ["sample", "leiden", "cell_type", "treatment", "mouse"] if c in u.columns])
            .properties(height=560), width="stretch")

# ---------------------------------------------------------------- cell types
with tab_ct:
    if labels is None or comp is None:
        st.info("Available after the integrate step with annotation enabled.")
    else:
        st.markdown("Clusters labelled with Griffin's marker lists. **Margin** = top score minus runner-up; "
                    "low-margin clusters need a manual look before anyone trusts the label.")
        level = st.radio("Level", ["cell_type", "cell_type_coarse"], horizontal=True,
                         format_func=lambda x: "fine" if x == "cell_type" else "coarse")
        c = comp.copy()
        if level == "cell_type_coarse":
            fine_to_coarse = dict(zip(labels["cell_type"], labels["cell_type_coarse"]))
            c["cell_type"] = c["cell_type"].map(lambda x: fine_to_coarse.get(x, x))
            c = c.groupby(["sample", "cell_type", "treatment", "mouse", "Gram", "study"],
                          observed=True, as_index=False)[["n_cells", "fraction"]].sum()
        st.markdown("**Composition per sample** (fraction of each sample's cells)")
        st.altair_chart(alt.Chart(c).mark_bar().encode(
            x=alt.X("sample:N", sort=samples), y=alt.Y("fraction:Q", stack="normalize", title="fraction"),
            color=alt.Color("cell_type:N", sort="-y"),
            tooltip=["sample", "cell_type", "n_cells", alt.Tooltip("fraction:Q", format=".1%")])
            .properties(height=380), width="stretch")
        pick = st.selectbox("Compare one cell type across samples", sorted(c["cell_type"].unique()),
                            index=sorted(c["cell_type"].unique()).index("Neu") if "Neu" in set(c["cell_type"]) else 0)
        one = c[c["cell_type"] == pick]
        st.altair_chart(alt.Chart(one).mark_bar().encode(
            x=alt.X("sample:N", sort=samples), y=alt.Y("fraction:Q", axis=alt.Axis(format="%")),
            color="mouse:N", tooltip=["sample", "n_cells", alt.Tooltip("fraction:Q", format=".1%")])
            .properties(height=260), width="stretch")
        st.caption("One sample per condition, so these are descriptive, not tested differences.")
        st.markdown("**Cluster labels**")
        st.dataframe(labels.round(3), hide_index=True, width="stretch")
        dot = RESULTS / "annotation" / "marker_dotplot.png"
        if dot.exists():
            st.image(str(dot), caption="Griffin's markers by assigned cell type")
