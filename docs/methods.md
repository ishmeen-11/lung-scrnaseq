# Methods: original R pipeline → this repo

| Step | Griffin (R / Seurat) | Here (Python / Scanpy) | Why |
|---|---|---|---|
| Counts | Cell Ranger 9.0.1, GRCm39-2024-A | Same (Snakemake rule, `--create-bam=true`) | Same reference for every study; BAMs kept for RNA velocity |
| Ambient RNA | SoupX per sample | Reads SoupX output (`soupx_mtx`) or CellBender output (`cellbender`) | No maintained SoupX in Python. CellBender is the planned comparison |
| QC metrics | nCount, nFeature, % mito, ribo | Same, plus **% hemoglobin**; ribo dropped (not used later) | Griffin's slides: RBCs are the main ambient source |
| Cell filters | Fixed: mito < 7.5%, UMIs > 1,000, genes > 200 | `qc.mode: mad` = per-sample, 5 MADs on log counts/genes, 3 MADs on mito (ceiling 15%); `qc.mode: fixed` reproduces the original | Studies differ in depth; fixed UMI cutoff may drop low-RNA neutrophils |
| Gene filter | Expressed in ≥ 3 cells, per sample | ≥ 3 cells, after merging | Avoids different gene sets per sample |
| Doublets | scDblFinder (clusters = T) | Scrublet (`sc.pp.scrublet`), per sample | Python-native; scDblFinder via R is a possible later comparison |
| Normalization | SCTransform, regress % mito | `lognorm` (default) or `pearson_residuals` (≈ SCTransform v2); no mito regression | Mito can be real infection biology; regression is optional to revisit |
| HVGs | SCT variable features | 3,000 seurat_v3 HVGs, selected per sample then combined | One lab can't dominate gene selection |
| Integration | RPCA (Seurat), 50 dims | Harmony on 50 PCs (default); scVI planned | RPCA is Seurat-only. Harmony is the closest well-tested "loose" integration |
| Clustering | Leiden, res 0.4, 30 dims | Leiden, res 0.4, 30 dims | Same |
| Annotation | Manual, guided by marker lists + UCell | Automatic first pass: `sc.tl.score_genes` with Griffin's final marker lists; each cluster gets the top-scoring type, low-margin calls flagged for review | Reproducible starting point; final labels still need a human check |

## New diagnostics (not filters)

- **Neutrophil retention**: putative neutrophils = top 10% by share of counts from
  *Csf3r, S100a8, S100a9, Retnlg, Cxcr2, Mmp9, Lcn2*. The funnel reports how many each
  filter removes.
- **Sex check**: *Xist* vs *Ddx3y, Eif2s3y, Kdm5d, Uty* per sample, to catch sex
  confounded with study or pathogen.
- **Depth by study**: median UMIs/genes per sample in the summary table.

## Validation plan

1. Run with `qc.mode: fixed` on Griffin's SoupX outputs and compare cell counts
   per sample and cluster labels with `currentObject.rds`.
2. Switch to `mad` and report the difference, especially neutrophils retained.
