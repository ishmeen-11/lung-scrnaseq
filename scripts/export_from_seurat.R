# Export per-sample counts + metadata from Griffin's currentObject.rds so the
# Python pipeline can read them (input_format: soupx_mtx).
# Note: these cells are already QC-filtered, so use this to validate the pipeline
# against Griffin's results, not to compare QC filters.
#
#   Rscript scripts/export_from_seurat.R path/to/currentObject.rds data/griffin_export
# Needs ~16 GB RAM and Seurat v5 (Matrix comes with it). No Bioconductor packages.

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "currentObject.rds"
out <- if (length(args) >= 2) args[2] else "data/griffin_export"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

write_10x <- function(m, dir) {
  # Same layout Cell Ranger writes: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  mtx <- file.path(dir, "matrix.mtx")
  writeMM(m, mtx)
  system2("gzip", c("-f", shQuote(mtx)))
  feats <- data.frame(id = rownames(m), name = rownames(m), type = "Gene Expression")
  gz <- gzfile(file.path(dir, "features.tsv.gz"), "w")
  write.table(feats, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(gz)
  gz <- gzfile(file.path(dir, "barcodes.tsv.gz"), "w")
  writeLines(colnames(m), gz); close(gz)
}

message("reading ", rds, " (this can take a few minutes)")
obj <- readRDS(rds)
counts <- LayerData(JoinLayers(obj[["RNA"]]), "counts")
for (s in levels(factor(obj$orig.ident))) {
  cells <- colnames(obj)[obj$orig.ident == s]
  write_10x(counts[, cells, drop = FALSE], file.path(out, s))
  message(s, ": ", length(cells), " cells")
}
write.csv(obj@meta.data, file.path(out, "metadata.csv"))
message("done -> ", out)
