filter_cells <- function(seurat_obj, basic_qc, ml_scores, ml_threshold = NULL) {
  if (is.null(ml_threshold)) {
    ml_threshold <- quantile(ml_scores, 0.95)  # top 5% outliers
  }
  keep <- ml_scores < ml_threshold
  return(subset(seurat_obj, cells = colnames(seurat_obj)[keep]))
}