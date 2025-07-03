
# data-raw/generate_pbmc2k.R
# Script to generate synthetic PBMC-like scRNA-seq count matrix

set.seed(42)
n_genes <- 1000
n_cells <- 2000

# Simulate sparse count matrix using negative binomial distribution
counts <- matrix(rnbinom(n_genes * n_cells, size = 2, prob = 0.98), nrow = n_genes, ncol = n_cells)
rownames(counts) <- paste0("Gene", seq_len(n_genes))
colnames(counts) <- paste0("Cell", seq_len(n_cells))

# Convert to data frame for export
df_counts <- as.data.frame(counts)

# Save as CSV to inst/extdata/
out_dir <- file.path("inst", "extdata")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(df_counts, file = file.path(out_dir, "pbmc2k.csv"), row.names = TRUE)

message("✅ Synthetic PBMC dataset saved to inst/extdata/pbmc2k.csv")
