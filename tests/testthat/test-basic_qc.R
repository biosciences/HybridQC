library(testthat)

test_that("HybridQC workflow runs end-to-end", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("pyod"), "pyod not available")

  library(Seurat)

  seurat_obj <- LoadPBMC2k()
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

  qc_basic <- run_basic_qc(seurat_obj)

  # Test Isolation Forest scoring
  seurat_obj <- run_isolation_forest_qc(
    seurat_obj,
    python_path = "/opt/homebrew/Caskroom/miniconda/base/envs/hybridqc_py311/bin/python",
    add_umap_plot = FALSE
  )

  expect_true("isoforest_score" %in% colnames(seurat_obj@meta.data))
  expect_equal(length(seurat_obj$isoforest_score), ncol(seurat_obj))

  # Test filtering
  ml_scores <- seurat_obj$isoforest_score
  filtered_obj <- filter_cells(seurat_obj, qc_basic, ml_scores)

  expect_s4_class(filtered_obj, "Seurat")
  expect_true(ncol(filtered_obj) < ncol(seurat_obj))
})