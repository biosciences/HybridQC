#' Run Isolation Forest QC using Python pyod
#'
#' @param seurat_obj A Seurat object with 'nFeature_RNA', 'nCount_RNA', 'percent.mt'
#' @param python_path Optional Python path (e.g., from conda or venv)
#' @param add_umap_plot If TRUE, display UMAP colored by outlier scores
#'
#' @return Seurat object with metadata column 'isoforest_score'
#' @export
run_isolation_forest_qc <- function(seurat_obj, python_path = NULL, add_umap_plot = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Please install the 'reticulate' package.")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Please install the 'Seurat' package.")

  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  # Ensure pyod is available
  if (!reticulate::py_module_available("pyod")) {
    stop("The 'pyod' module is not available in the selected Python environment.")
  }

  py_iforest <- reticulate::import("pyod.models.iforest", convert = TRUE)

  # Check required features exist
  required_vars <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  missing_vars <- setdiff(required_vars, colnames(seurat_obj@meta.data))
  if (length(missing_vars) > 0) {
    stop("Missing QC features in Seurat object: ", paste(missing_vars, collapse = ", "))
  }

  # Extract and clean QC matrix
  qc_data <- Seurat::FetchData(seurat_obj, vars = required_vars)
  qc_data <- as.matrix(na.omit(qc_data))

  clf <- py_iforest$IForest()
  clf$fit(qc_data)
  scores <- clf$decision_scores_

  seurat_obj$isoforest_score <- NA
  seurat_obj$isoforest_score[rownames(qc_data)] <- scores

  if (add_umap_plot && "umap" %in% names(seurat_obj@reductions)) {
    p <- Seurat::FeaturePlot(seurat_obj, features = "isoforest_score", reduction = "umap") +
      ggplot2::ggtitle("Isolation Forest Outlier Score")
    print(p)
  }

  return(seurat_obj)
}