# Export per-sample counts + metadata from Griffin's currentObject.rds so the
# Python pipeline can read them (input_format: soupx_mtx).
# Note: these cells are already QC-filtered, so use this to validate the pipeline
# against Griffin's results, not to compare QC filters.
#
#   Rscript scripts/export_from_seurat.R path/to/currentObject.rds data/griffin_export
# Needs ~16 GB RAM, Seurat v5 and DropletUtils.

suppressPackageStartupMessages({ library(Seurat); library(DropletUtils) })
args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "currentObject.rds"
out <- if (length(args) >= 2) args[2] else "data/griffin_export"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(rds)
counts <- LayerData(JoinLayers(obj[["RNA"]]), "counts")
for (s in levels(factor(obj$orig.ident))) {
  cells <- colnames(obj)[obj$orig.ident == s]
  write10xCounts(file.path(out, s), counts[, cells], version = "3", overwrite = TRUE)
  message(s, ": ", length(cells), " cells")
}
write.csv(obj@meta.data, file.path(out, "metadata.csv"))
message("done -> ", out)
