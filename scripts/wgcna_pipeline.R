# Co-Expression Network Analysis for Multi-Omic Integration and Functional Characterization

# ===================================
# 1. Load Required Libraries
# ===================================
suppressPackageStartupMessages({
  library(WGCNA)  # Main package for co-expression network analysis
  library(ggplot2)  # For data visualization
  library(dplyr)  # For data manipulation
  library(tidyr)  # For data tidying
  library(data.table)  # For fast data loading
  library(reshape2)  # For matrix reshaping
  library(pheatmap)  # For heatmap visualization
})

# ===================================
# Project Overview and Analysis Purposes
# ===================================
# This pipeline is designed for gene co-expression network analysis using WGCNA.
# It includes several key steps:
#
#   - Data preprocessing and quality control
#   - Network construction and module detection
#   - Module significance assessment
#   - Functional enrichment analysis
#   - Correlation with clinical traits or custom gene sets
#   - Exporting data for downstream analyses
#
# Analysis Purposes:
#
#   1. Gene Prioritization for Functional Studies
#   2. Network-Based Biomarker Discovery
#   3. Pathway and Functional Enrichment
#   4. Data Integration Across Omic Layers
#   5. Single-Cell or Spatial Transcriptomics
#   7. Regulatory Network Reconstruction (Future Project)

# ===================================
# 2. Load and Preprocess Expression Data
# ===================================
expr_file <- "data/raw/example_expression_data.tsv"
expr_data <- fread(expr_file, header = TRUE, row.names = 1)
cat("Expression data loaded. Dimensions:", dim(expr_data), "\n")

# Data Quality Check and Preprocessing
expr_data <- expr_data[complete.cases(expr_data), ]  # Remove rows with missing values
expr_data <- expr_data[rowSums(expr_data) > 0, ]  # Remove zero-sum rows
expr_data <- expr_data[, apply(expr_data, 2, var) > 0.01]  # Remove low-variance genes
cat("After filtering low-variance genes, dimensions:", dim(expr_data), "\n")

# Transpose the expression data for WGCNA
expr_data <- t(expr_data)
cat("Data transposed for WGCNA. Final dimensions:", dim(expr_data), "\n")

# ===================================
# 3. Load and Match Trait Data
# ===================================
trait_file <- "data/raw/example_traits.tsv"
trait_data <- fread(trait_file, header = TRUE, row.names = 1)
cat("Trait data loaded. Dimensions:", dim(trait_data), "\n")

# Match Trait Data to Expression Data
matched_traits <- trait_data[rownames(expr_data), , drop = FALSE]
cat("Matched trait data dimensions:", dim(matched_traits), "\n")

# ===================================
# 4. Network Construction and Module Detection
# ===================================
# Build the Network
softPower <- 6  # Update this based on your data analysis
adjacency_matrix <- adjacency(expr_data, power = softPower)
dissTOM <- 1 - TOMsimilarity(adjacency_matrix)
geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf("figures/gene_clustering_tree.pdf")
plot(geneTree, main = "Gene Clustering Dendrogram", xlab = "Genes", sub = "", cex.lab = 1.2, cex.axis = 0.8, cex.main = 1.5)
dev.off()

# ===================================
# 5. Module Characterization and Enrichment
# ===================================

# Module Characterization and Gene Prioritization
# -----------------------------------
# This section focuses on extracting biologically meaningful insights from the detected modules.
# It directly supports the original analysis purposes, including:
#
#   1. Gene Prioritization for Functional Studies:
#      - Identify hub genes with high Module Membership (MM) and high Gene Significance (GS) values.
#      - Prioritize genes that are central to key modules and strongly correlated with important clinical traits or experimental conditions.
#      - Ideal for selecting candidates for targeted validation studies like CRISPR knockout or siRNA knockdown.
#
#   2. Correlation with Clinical Traits or Custom Gene Sets:
#      - Use Gene Significance (GS) to measure the correlation between individual gene expression profiles and traits.
#      - Highly flexible, as it can incorporate clinical traits, pathway scores, or any custom gene sets depending on the research question.
#      - Can also be extended to include pathway activity scores or other functional annotations for a more comprehensive analysis.
#
# The calculations below include both MM and GS for each gene, supporting downstream prioritization and functional studies.
# ===================================

minModuleSize <- 30
moduleColors <- labels2colors(cutreeDynamic(dendro = geneTree, distM = as.dist(dissTOM), deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize))

# Calculate Module Eigengenes and Statistics
MEs <- moduleEigengenes(expr_data, colors = moduleColors)$eigengenes
moduleTraitCor <- cor(MEs, matched_traits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(expr_data))

# Calculate Module Membership (MM) and Gene Significance (GS)
module_stats <- data.frame(Gene = rownames(expr_data), Module = moduleColors)

# Flexible GS Calculation (can use clinical traits or gene sets)
for (trait in colnames(matched_traits)) {
  gs_values <- cor(expr_data, matched_traits[[trait]], use = "p")
  module_stats[[paste0("GS_", trait)]] <- gs_values
}

# Calculate MM for each Module
for (module in unique(moduleColors)) {
  module_genes <- module_stats$Gene[module_stats$Module == module]
  module_mm <- cor(expr_data[, module_genes], MEs[, paste0("ME", module)], use = "p")
  module_stats[module_stats$Module == module, paste0("MM_", module)] <- module_mm
}

# Save Module-Level Statistics
write.csv(module_stats, "results/module_statistics.csv", row.names = FALSE)
cat("Module statistics, including MM and GS, saved.
")

# Functional Enrichment Analysis
gmt_file <- "data/raw/example_gene_sets.gmt"
gene_sets <- readLines(gmt_file)
gene_set_list <- lapply(gene_sets, function(line) strsplit(line, "\t")[[1]][-c(1,2)])
names(gene_set_list) <- sapply(gene_sets, function(line) strsplit(line, "\t")[[1]][1])

module_enrichment <- lapply(unique(moduleColors), function(mod) {
  mod_genes <- module_stats$Gene[module_stats$Module == mod]
  enriched_sets <- sapply(gene_set_list, function(genes) {
    overlap <- length(intersect(mod_genes, genes))
    total_genes <- length(mod_genes)
    total_pathway <- length(genes)
    p_value <- phyper(overlap - 1, total_pathway, length(module_stats$Gene) - total_pathway, total_genes, lower.tail = FALSE)
    return(p_value)
  })
  return(enriched_sets)
})
names(module_enrichment) <- unique(moduleColors)

# Save Enrichment Results
# Note: Enrichment results can guide the selection of specific modules for further study,
# focusing on those enriched for pathways relevant to your research question or experimental context.
# This can also help prioritize modules for targeted validation or integration with other omic layers.
enrichment_df <- do.call(rbind, lapply(names(module_enrichment), function(mod) {
  data.frame(Module = mod, Gene_Set = names(module_enrichment[[mod]]),
             P_Value = module_enrichment[[mod]])
}))
enrichment_df$Adjusted_P_Value <- p.adjust(enrichment_df$P_Value, method = "BH")  # Multiple testing correction
write.csv(enrichment_df, "results/module_enrichment.csv", row.names = FALSE)
cat("Gene set enrichment results saved.
")

# ===================================
# Module Prioritization Visualizations
# ===================================

# Bar plot of top enriched pathways
library(ggplot2)

top_enriched <- enrichment_df[enrichment_df$Adjusted_P_Value < 0.05, ]
top_enriched$-log10P <- -log10(top_enriched$Adjusted_P_Value)

ggplot(top_enriched, aes(x = reorder(Gene_Set, -log10P), y = -log10P, fill = Module)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Enriched Pathways by Module",
       x = "Gene Set",
       y = "-log10(Adjusted P-Value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
ggsave("figures/top_enriched_pathways.png", width = 12, height = 8)
cat("Top enriched pathways bar plot saved as figures/top_enriched_pathways.png.
")
cat("Gene set enrichment results saved.\n")

cat("WGCNA Pipeline completed successfully with enrichment analysis.\n")

# ===================================
# 6. Multi-Omic Integration
# ===================================

# Multi-Omic Integration via Correlation Analysis
# This section explores the association between module eigengenes (MEs) and other omic layers (e.g., genomic, proteomic, methylation).
# Here, we focus on simple Pearson correlations, but this can be extended to partial correlations or mutual information for non-linear relationships.

# Load Genomic Data
genomic_file <- "data/raw/example_genomic_data.tsv"
genomic_data <- fread(genomic_file, header = TRUE, row.names = 1)
cat("Genomic data loaded. Dimensions:", dim(genomic_data), "
")

# Match Genomic Data to Expression Data
matched_genomic <- genomic_data[rownames(expr_data), , drop = FALSE]

# Calculate Correlation
multiomic_correlation_matrix <- cor(MEs, matched_genomic, use = "p")

# Apply Multiple Testing Correction
multiomic_pvals <- apply(multiomic_correlation_matrix, 2, function(col) corPvalueStudent(col, nrow(expr_data)))
multiomic_adj_pvals <- p.adjust(multiomic_pvals, method = "BH")

# Save Multi-Omic Correlations
multiomic_results <- data.frame(Feature = colnames(matched_genomic), Correlation = multiomic_correlation_matrix, Adjusted_P_Value = multiomic_adj_pvals)
write.csv(multiomic_results, "results/multiomic_module_correlations.csv", row.names = FALSE)
cat("Multi-omic module correlations saved.
")

# Visualize Correlation Heatmap
pheatmap(multiomic_correlation_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Multi-Omic Module Correlations",
         filename = "figures/multiomic_correlation_heatmap.png")
cat("Multi-omic correlation heatmap saved.
")


# ===================================
# 7. Single-Cell / Spatial Transcriptomics Integration
# ===================================

# Single-Cell or Spatial Data Integration via Correlation Analysis
# This section focuses on correlating module eigengenes (MEs) with single-cell or spatial transcriptomics data.
# As with the multi-omic section, this can be extended to more complex metrics if needed.

# Load Single-Cell Data
single_cell_file <- "data/raw/example_single_cell_data.tsv"
single_cell_data <- fread(single_cell_file, header = TRUE, row.names = 1)
cat("Single-cell data loaded. Dimensions:", dim(single_cell_data), "
")

# Match Single-Cell Data to Expression Data
matched_single_cell <- single_cell_data[rownames(expr_data), , drop = FALSE]

# Calculate Correlation
single_cell_correlation_matrix <- cor(MEs, matched_single_cell, use = "p")

# Apply Multiple Testing Correction
single_cell_pvals <- apply(single_cell_correlation_matrix, 2, function(col) corPvalueStudent(col, nrow(expr_data)))
single_cell_adj_pvals <- p.adjust(single_cell_pvals, method = "BH")

# Save Single-Cell Correlations
single_cell_results <- data.frame(Feature = colnames(matched_single_cell), Correlation = single_cell_correlation_matrix, Adjusted_P_Value = single_cell_adj_pvals)
write.csv(single_cell_results, "results/single_cell_module_correlations.csv", row.names = FALSE)
cat("Single-cell module correlations saved.
")

# Visualize Correlation Heatmap
pheatmap(single_cell_correlation_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Single-Cell Module Correlations",
         filename = "figures/single_cell_correlation_heatmap.png")
cat("Single-cell correlation heatmap saved.
")


