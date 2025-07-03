run_basic_qc <- function(seurat_obj, mito_pattern = "^MT-") {
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = mito_pattern)
  qc_metrics <- Seurat::FetchData(seurat_obj, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  return(qc_metrics)
}