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
LoadPBMC2k <- function() {
  # Directly point to inst/extdata/ path (for development use)
  data_path <- file.path("inst", "extdata", "pbmc2k.csv")

  if (!file.exists(data_path)) {
    stop("pbmc2k.csv not found in inst/extdata/. Please ensure the file exists.")
  }

  counts <- read.csv(data_path, row.names = 1, check.names = FALSE)
  seurat_obj <- Seurat::CreateSeuratObject(counts = as.matrix(counts))
  return(seurat_obj)
}