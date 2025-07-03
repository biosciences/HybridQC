#' Load PBMC2k Example Dataset
#'
#' This function loads a synthetic 2,000-cell × 1,000-gene count matrix from
#' the `inst/extdata/` directory and returns it as a Seurat object.
#'
#' @return A Seurat object containing the pbmc2k data
#' @export
#'
#' @examples
#' pbmc <- LoadPBMC2k()
LoadPBMC2k <- function(dev_path = NULL) {
  if (!is.null(dev_path) && file.exists(dev_path)) {
    counts <- read.csv(dev_path, row.names = 1, check.names = FALSE)
  } else {
    data_path <- system.file("extdata", "pbmc2k.csv", package = "HybridQC")
    if (data_path == "") stop("pbmc2k.csv not found. Please provide dev_path.")
    counts <- read.csv(data_path, row.names = 1, check.names = FALSE)
  }
  Seurat::CreateSeuratObject(counts = as.matrix(counts))
}