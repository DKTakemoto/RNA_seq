library(tximport)
library(DESeq2)
library(readr)
library(pheatmap)
library(ggplot2)

# === 1. Read sample metadata ===
base_dir <- "/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo"
samples <- read.table(file.path(base_dir, "salmon", "salmon_metadata.txt"), header=TRUE, stringsAsFactors=FALSE)

# === 2. Prepare Salmon quantification file paths ===
files <- file.path(
  "/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/salmon",
  samples$sample,
  "quant.sf"
)
names(files) <- samples$sample

# === 3. Import transcript-level quantification (no tx2gene) ===
txi <- tximport(files, type = "salmon", txOut = TRUE)

# === 4. Create DESeq2 dataset ===
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

# === 5. Prefiltering: keep transcripts with enough counts ===
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# === 6. Run DESeq2 differential analysis ===
dds <- DESeq(dds)

# === 7. Extract results ===
res <- results(dds)

# === 8. Order results by adjusted p-value and save ===
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), file="/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/salmon/dte_results_transcripts_C3xN.csv")

# === 9. Volcano plot (ggplot2, manual) ===
pCutoff <- 0.05
FCcutoff <- 1.0

res_df <- as.data.frame(res)
res_df$log10padj <- -log10(res_df$padj)
res_df$significance <- "NS"
res_df$significance[res_df$padj < pCutoff & abs(res_df$log2FoldChange) >= FCcutoff] <- "FDR < 0.05 & |log2FC| ≥ 1"
res_df$significance[res_df$padj < pCutoff & abs(res_df$log2FoldChange) < FCcutoff] <- "FDR < 0.05"
res_df$significance[res_df$padj >= pCutoff & abs(res_df$log2FoldChange) >= FCcutoff] <- "|log2FC| ≥ 1"

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = log10padj, color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("gray", "blue", "red", "purple")) +
  geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pCutoff), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    subtitle = "Treated vs. Control",
    x = expression(Log[2]~fold~change),
    y = expression(-Log[10]~adjusted~italic(P)),
    caption = paste0("log2 FC cutoff: ", FCcutoff,
                     "; FDR cutoff: ", pCutoff,
                     "\nTotal transcripts = ", nrow(res_df))
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/salmon/DTE_VolcanoPlot_ggplot2.png", volcano_plot, width = 7, height = 7, dpi = 300)

# === 10. MA plot ===
png("/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/salmon/DTE_MAplot.Salmon.png", width = 7, height = 7, units = "in", res = 300)
plotMA(res, main = "DTE results C3xN", ylim = c(-5, 5))
dev.off()

# === 11. PCA plot ===
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  ggtitle("PCA of Transcript Expression (VST-transformed)")

ggsave("/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/salmon/DTE_PCAplot.Salmon.png", pca_plot, width = 7, height = 7, dpi = 300)

# Heatmap
# Use vst-transformed data
vsd_mat <- assay(vsd)

# Select top 30 most variable transcripts
topVarGenes <- head(order(rowVars(vsd_mat), decreasing = TRUE), 30)

# Optional: center each gene (mean 0)
mat_scaled <- t(scale(t(vsd_mat[topVarGenes, ])))

# Annotate columns by condition
annotation_col <- data.frame(condition = colData(dds)$condition)
rownames(annotation_col) <- colnames(vsd_mat)

# Save heatmap
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  fontsize_col = 10,
  main = "Top 30 Variable Transcripts",
  filename = "DTE_heatmap_top30.png",
  width = 8,
  height = 10
)

# === 12. Summary ===
summary(res)
