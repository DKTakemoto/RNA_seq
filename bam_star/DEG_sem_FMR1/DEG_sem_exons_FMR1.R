#!/usr/bin/env Rscript

# Carregar bibliotecas
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)

# Diretórios
root_dir <- "/work/pancreas/takemoto/RNA_seq"
figures_dir <- file.path(root_dir, "figures", "DEG_sem_FMR1")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Ler os metadados
colData <- read.csv(file.path(root_dir, "data", "reference", "metadata_DEG.csv"), header = TRUE, row.names = 1)

# Ler matriz de contagens
countData <- read.csv(file.path(root_dir, "results", "feature_counts", "counts_FMR1_exons_removidos.csv"), header = TRUE, row.names = 1)

# Ajustar os nomes das colunas para combinar com colData
colnames(countData) <- gsub("trimmed_(.*?)_.*", "\\1", colnames(countData))
countData <- countData[, rownames(colData)]

# Ajustar o nível do fator
colData$Condition <- factor(colData$Condition)
colData$Condition <- relevel(colData$Condition, ref = "Control")

# Criar o objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ Condition)
dds <- DESeq(dds)

# Comparações
res_1_vs_control <- results(dds, contrast = c("Condition", "circRNAx13x14x15x16", "Control"))
res_2_vs_control <- results(dds, contrast = c("Condition", "circRNAx13x16", "Control"))

# Importar anotação
anotacao <- read.delim(file.path(root_dir, "data", "reference",  "mart_export.txt"), header = FALSE, col.names = c("Gene_ID", "Gene_Type", "Gene_Name"))
anotacao <- as.data.frame(apply(anotacao, 2, function(x) gsub('\\s+', '', x)))

# Adicionar coluna de nomes dos genes
res_tb_1 <- as.data.frame(res_1_vs_control)
res_tb_2 <- as.data.frame(res_2_vs_control)
res_tb_1$Genes <- rownames(res_tb_1)
res_tb_2$Genes <- rownames(res_tb_2)

# Adicionar anotação
res_tb_1 <- merge(res_tb_1, anotacao, by.x = "Genes", by.y = "Gene_ID", all.x = TRUE)
res_tb_2 <- merge(res_tb_2, anotacao, by.x = "Genes", by.y = "Gene_ID", all.x = TRUE)

# Adicionar coluna de regulação
res_tb_1$Regulation <- "No sig"
res_tb_1$Regulation[res_tb_1$log2FoldChange >= 0.58 & res_tb_1$padj <= 0.05] <- "Upregulated"
res_tb_1$Regulation[res_tb_1$log2FoldChange <= -0.58 & res_tb_1$padj <= 0.05] <- "Downregulated"

res_tb_2$Regulation <- "No sig"
res_tb_2$Regulation[res_tb_2$log2FoldChange >= 0.58 & res_tb_2$padj <= 0.05] <- "Upregulated"
res_tb_2$Regulation[res_tb_2$log2FoldChange <= -0.58 & res_tb_2$padj <= 0.05] <- "Downregulated"

# Filtrar DEGs significativos
sig_DEG_1_vs_control <- res_tb_1 %>% filter(Regulation %in% c("Upregulated", "Downregulated"))
sig_DEG_2_vs_control <- res_tb_2 %>% filter(Regulation %in% c("Upregulated", "Downregulated"))

# Salvar resultados
write.csv(sig_DEG_1_vs_control, file = file.path(root_dir, "results", "DEG_sem_FMR1", "DEGs_C3_vs_N.csv"), row.names = FALSE)
write.csv(sig_DEG_2_vs_control, file = file.path(root_dir, "results", "DEG_sem_FMR1", "DEGs_C2_vs_N.csv"), row.names = FALSE)
write.csv(rbind(sig_DEG_1_vs_control, sig_DEG_2_vs_control), file = file.path(root_dir, "results", "DEG_sem_FMR1", "DEGs_combined_C3x_C2_vs_N.csv"), row.names = FALSE)

# PCA
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = c("Condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca <- ggplot(pcaData, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 0.6),
    legend.title = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(filename = file.path(figures_dir, "PCA_sem_FMR1.tiff"), plot = pca, units = "in", width = 9, height = 6, dpi = 500)

# Volcano plots
mycolors <- c("Upregulated" = "firebrick", "Downregulated" = "royalblue", "No sig" = "gray70")

volcano1 <- ggplot(res_tb_1, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  labs(title = "circx13x14x15x16 vs control", x = "log2 Fold Change", y = "-log10 Adjusted p-value")

ggsave(filename = file.path(figures_dir, "Volcano_C3_vs_N.tiff"), plot = volcano1, units = "in", width = 7, height = 6, dpi = 500)

volcano2 <- ggplot(res_tb_2, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  labs(title = "circx13x16 vs control", x = "log2 Fold Change", y = "-log10 Adjusted p-value")

ggsave(filename = file.path(figures_dir, "Volcano_C2_vs_N.tiff"), plot = volcano2, units = "in", width = 7, height = 6, dpi = 500)

# Heatmaps
genes_deg_1 <- sig_DEG_1_vs_control$Genes
heatmap_data_1 <- assay(vsd)[genes_deg_1, ]
samples_1 <- rownames(colData)[colData$Condition %in% c("Control", "circRNAx13x14x15x16")]

tiff(file.path(figures_dir, "Heatmap_C3_vs_N.tiff"), width = 9, height = 8, units = "in", res = 300)
pheatmap(heatmap_data_1[, samples_1],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = colData[samples_1, , drop = FALSE],
         fontsize_col = 12,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

genes_deg_2 <- sig_DEG_2_vs_control$Genes
heatmap_data_2 <- assay(vsd)[genes_deg_2, ]
samples_2 <- rownames(colData)[colData$Condition %in% c("Control", "circRNAx13x16")]

tiff(file.path(figures_dir, "Heatmap_C2_vs_N.tiff"), width = 9, height = 8, units = "in", res = 300)
pheatmap(heatmap_data_2[, samples_2],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = colData[samples_2, , drop = FALSE],
         fontsize_col = 12,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
