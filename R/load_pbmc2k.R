#' Load PBMC 2k example data
#'
#' Load Seurat object from pbmc2k.csv (used in HybridQC examples)
#' @param dev_path Optional fallback path to pbmc2k.csv when running locally
#' @return A Seurat object
#' @export
LoadPBMC2k <- function(dev_path = NULL) {
  if (!is.null(dev_path) && file.exists(dev_path)) {
    message("✅ Using dev_path: ", dev_path)
    counts <- read.csv(dev_path, row.names = 1, check.names = FALSE)
  } else {
    data_path <- system.file("extdata", "pbmc2k.csv", package = "HybridQC")
    if (data_path == "") stop("pbmc2k.csv not found. Provide dev_path = 'inst/extdata/pbmc2k.csv'")
    counts <- read.csv(data_path, row.names = 1, check.names = FALSE)
  }
  Seurat::CreateSeuratObject(counts = as.matrix(counts))
}