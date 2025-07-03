# HybridQC <img src="https://img.shields.io/badge/R-4.0+-brightgreen" align="right"/>

**HybridQC** is a lightweight R package that improves quality control (QC) for single-cell RNA-seq data by combining traditional QC metrics (gene/cell counts, mitochondrial % content) with machine learning–based outlier detection (e.g., Isolation Forest).

---

## 📄 Project Links
- 📂 [Source Code](https://github.com/biosciences/HybridQC): Explore the full repository
- 🔗 [Live Report](https://biosciences.github.io/HybridQC/HybridQC.html): View the interactive HTML output

## 🚀 Features

- Rule-based filtering (nFeatures, nCounts, percent.mito)
- Machine learning QC using Isolation Forest
- Adaptive, data-driven filtering
- Works with Seurat objects
- Supports integration with Python ML tools via `reticulate`

---

## 🧪 Installation

```r
# Install devtools if not already
install.packages("devtools")

# Install HybridQC from GitHub
devtools::install_github("biosciences/HybridQC")

# If you are using zsh
echo ". /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh" >> ~/.zshrc
source ~/.zshrc

# If you are uising bash
echo ". /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh" >> ~/.bash_profile
source ~/.bash_profile

conda create -n hybridqc_py311 python=3.11 -y
conda activate hybridqc_py311
conda install pyod scikit-learn numpy numba -c conda-forge -y

# Test using hybridqc_py311
reticulate::use_condaenv("hybridqc_py311", required = TRUE)
reticulate::py_config()  # Confirm that the path points to Conda’s hybridqc_py311 environment
reticulate::import("pyod.models.iforest")  # ✅ If successful, it means the environment is ready
```

## 📦 Dependencies
• R (>= 4.0)
• Seurat
• dplyr
• reticulate
• Python: pyod (will be installed automatically)

## 📖 Usage Example
```r
library(Seurat)
library(HybridQC)
library(reticulate)
library(ggplot2)

# Set working directory (replace this with your path)
setwd("path/to/output/folder")

# Load Seurat object
seurat_obj <- LoadPBMC2k()  # Replace with your own data loader

# Step 1: Basic QC
qc_basic <- run_basic_qc(seurat_obj)
write.csv(qc_basic, file = "results/qc/qc_basic.csv", row.names = TRUE)

# Step 2: Configure Python environment
reticulate::use_condaenv("hybridqc_py311", required = TRUE)
reticulate::py_config()

# Step 3: Preprocessing for ML-based QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Step 4: Run Isolation Forest QC and generate UMAP
seurat_obj <- run_isolation_forest_qc(
  seurat_obj,
  python_path = "/opt/homebrew/Caskroom/miniconda/base/envs/hybridqc_py311/bin/python",
  add_umap_plot = FALSE  # turn off interactive plot so we can save it below
)

# Save Isolation Forest scores
iso_scores <- FetchData(seurat_obj, vars = "isoforest_score")
write.csv(iso_scores, file = "results/qc/isoforest_scores.csv", row.names = TRUE)

# Save UMAP plot
p <- FeaturePlot(seurat_obj, features = "isoforest_score", reduction = "umap") +
  ggtitle("Isolation Forest Outlier Score (UMAP)")
ggsave("results/plots/isoforest_umap.png", plot = p, width = 7, height = 5)

# Step 5: Hybrid Filtering using vector of scores
ml_scores <- seurat_obj$isoforest_score
filtered_obj <- filter_cells(seurat_obj, qc_basic, ml_scores)

# Save filtered Seurat metadata
write.csv(filtered_obj@meta.data, file = "results/qc/filtered_cells_metadata.csv")
```
💡 See vignettes/HybridQC.Rmd for a full walkthrough.

## 📄 Documentation
• Vignette: vignettes/HybridQC.Rmd
• Package Manual: See man/ folder
• Example dataset: inst/extdata/pbmc2k.csv


## 📝 Citation

This package is accompanied by a JOSS submission:

Lai K. (2025). HybridQC: Machine Learning-Augmented Quality Control for Single-Cell RNA-seq Data. Journal of Open Source Software (under review). [link pending]

## 🛠️ License

MIT © Kaitao Lai

⸻

## 🌐 Links
• GitHub: https://github.com/biosciences/HybridQC
• JOSS: https://joss.theoj.org

---

Developed by Kaitao Lai, University of Sydney.
