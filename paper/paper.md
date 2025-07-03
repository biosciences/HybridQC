---
title: 'HybridQC: Machine Learning-Augmented Quality Control for Single-Cell RNA-seq Data'
tags:
  - R
  - single-cell RNA-seq
  - quality control
  - machine learning
  - Seurat
  - bioinformatics
authors:
  - name: Kaitao Lai
    orcid: 0000-0002-9420-9352
    affiliation: "1"
affiliations:
  - name: The University of Sydney
    index: 1
date: 2025-09-03
bibliography: paper.bib
---

# Summary

HybridQC is an R package that streamlines quality control (QC) of single-cell RNA sequencing (scRNA-seq) data by combining traditional threshold-based filtering with machine learning–based outlier detection. It provides an efficient and adaptive framework to identify low-quality cells in noisy or shallow-depth datasets using techniques such as Isolation Forest, while remaining compatible with widely adopted formats such as Seurat objects.

The package is lightweight, easy to install, and suitable for small-to-medium scRNA-seq datasets in research settings. HybridQC is especially useful for projects involving non-model organisms, rare samples, or pilot studies, where automated and flexible QC is critical for reproducibility and downstream analysis.

# Statement of Need

scRNA-seq experiments often suffer from technical noise, dropout events, and variability in sequencing depth. Traditional quality control relies on static cutoffs for metrics such as gene count, UMI count, and mitochondrial content, which may be suboptimal for non-standard datasets. 

HybridQC fills this gap by integrating machine learning methods—specifically unsupervised outlier detection—with traditional QC to improve filtering precision and robustness. This dual-level approach can better preserve informative but unconventional cell types and adapt dynamically to diverse datasets. No existing R packages provide this hybrid QC strategy as a standalone tool with seamless integration into Seurat-based pipelines.

# Features

- Computes standard QC metrics: `nFeature_RNA`, `nCount_RNA`, `percent.mt`
- Supports Isolation Forest outlier detection via `reticulate` and `pyod`
- Filters cells using a hybrid decision rule
- Works on Seurat objects
- Lightweight and suitable for quick prototyping or small studies

# Example Usage

```r
library(Seurat)
library(HybridQC)

pbmc <- LoadPBMC2k()  # Load example data
qc_basic <- run_basic_qc(pbmc)
ml_scores <- run_isolation_forest_qc(pbmc)
filtered <- filter_cells(pbmc, qc_basic, ml_scores)
```

# Acknowledgements

The author thanks collaborators at University of Sydney for feedback on early concepts.

# References

• Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A., Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M., Colomé-Tatché, M., & Theis, F. J. (2022). Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19(1), 41–50. https://doi.org/10.1038/s41592-021-01336-8
• McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv preprint arXiv:1802.03426. https://arxiv.org/abs/1802.03426
• Zhao, Y., Nasrullah, Z., & Li, Z. (2018). PyOD: A Python Toolbox for Scalable Outlier Detection. arXiv preprint arXiv:1901.01588. https://arxiv.org/abs/1901.01588