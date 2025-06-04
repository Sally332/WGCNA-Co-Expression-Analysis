# **WGCNA Co-Expression Analysis Pipeline**

## **Analysis Purposes**

* **Bulk RNA-seq**: Identify co-expression modules and correlate with sample traits (e.g., disease status, clinical measures).
* **Single-Cell RNA-seq**: Adapt WGCNA methods to pseudo-bulk or cell-type–specific analyses to uncover cell subpopulation networks.
* **Spatial Transcriptomics**: Integrate spatial coordinates with co-expression modules to reveal tissue architecture–driven gene networks.
* **Hypothesis Generation**: Discover key regulatory modules and candidate driver genes for further validation.

---

## **Key Methodological Steps**

1. **Data Import & QC**

   * Load expression matrix (bulk, single-cell pseudo-bulk, or spatial spot matrix) and sample/spot metadata.
   * Filter low-expression features and (optionally) normalize (e.g., variance-stabilizing transformation, SCTransform).

2. **Soft-Threshold Power Selection**

   * Determine optimal power for scale-free network using `pickSoftThreshold()`.
   * Plot scale-free topology fit and mean connectivity to guide selection.

3. **Network Construction**

   * Compute adjacency matrix from feature-feature correlations (Pearson/Spearman).
   * Calculate Topological Overlap Matrix (TOM) to measure network interconnectedness.

4. **Module Detection**

   * Perform hierarchical clustering on TOM-based dissimilarity.
   * Use dynamic tree cutting to define co-expression modules and merge similar modules.

5. **Module Eigengene & Trait Analysis**

   * Compute module eigengenes (MEs) representing module expression profiles.
   * Correlate MEs with traits (bulk phenotype, cell-type proportions, or spatial features) to identify biologically relevant modules.

6. **Functional Enrichment**

   * Perform GO/KEGG enrichment on genes within significant modules using `clusterProfiler`.
   * Visualize enriched pathways to interpret module functionality.

7. **Reporting & Visualization**

   * Generate plots: sample/spot dendrogram, power selection plots, module dendrogram, ME–trait heatmap.
   * Compile an HTML report summarizing results and figures.

---

## **Quick Start**

```bash
git clone https://github.com/<YOUR-USERNAME>/WGCNA-Coexpression-Analysis.git
cd WGCNA-Coexpression-Analysis

# Install R dependencies
Rscript -e "install.packages(c('data.table','WGCNA','clusterProfiler','org.Hs.eg.db','ggplot2','pheatmap','optparse')); BiocManager::install(c('DESeq2','limma','GSVA'))"

# Run pipeline (bulk, single-cell pseudo-bulk, or spatial)
Rscript scripts/wgcna_pipeline.R \
  --expr data/expression_matrix.tsv \
  --traits data/sample_traits.tsv \
  --out results/wgcna_results
```

## **Folder Structure**

```text
WGCNA-Coexpression-Analysis/
├── README.md             
├── data/                 # Input (expression, metadata)
├── scripts/
│   └── wgcna_pipeline.R  # Main pipeline script
└── results/              # Auto-generated outputs
```
For full documentation, see docs/README.pdf.

**Contact**: Sally Yepes (sallyepes233@gmail.com)| **License**: MIT | **Last Updated**: June 2025
