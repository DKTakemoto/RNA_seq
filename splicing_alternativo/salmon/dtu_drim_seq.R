library(tximport)
library(DRIMSeq)
library(stageR)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(GenomicFeatures)
library(pheatmap)
library(txdbmaker)
library(matrixStats)
library(reshape2)
library(data.table)
library(rtracklayer)
library(ggrepel)


# --- 1. Configuration ---
base_dir <- "/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo"
metadata_path <- file.path(base_dir, "salmon", "salmon_metadata.txt")
samples <- read.table(metadata_path, header = TRUE, stringsAsFactors = FALSE)
gtf_path <- "/work/pancreas/takemoto/RNA_seq/data/reference/database/gencode.v47.annotation.gtf"

# --- 2. Create TxDb and tx2gene (KEEP VERSIONS) ---
txdb <- makeTxDbFromGFF(file = gtf_path, format = "gtf")
tx2gene_df <- AnnotationDbi::select(txdb,
                                    keys = keys(txdb, keytype = "TXNAME"),
                                    columns = c("TXNAME", "GENEID"),
                                    keytype = "TXNAME")  # No version stripping here

# Add transcript count per gene
tab <- table(tx2gene_df$GENEID)
tx2gene_df$ntx <- tab[match(tx2gene_df$GENEID, names(tab))]

# --- 3. tximport (KEEP VERSIONS) ---
files <- file.path(base_dir, "salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "scaledTPM")
txi$counts <- txi$counts[rowSums(txi$counts) > 0, ]  # No version stripping of rownames

# --- 4. Prepare data for DRIMSeq (KEEP VERSIONS) ---
sample_info <- data.frame(sample_id = samples$sample, group = samples$condition)
cts <- as.data.frame(txi$counts)
cts$feature_id <- rownames(cts)  # Preserve versions

# Merge with tx2gene pairs (maintain versions)
cts <- merge(cts, tx2gene_df[, c("TXNAME", "GENEID")], by.x = "feature_id", by.y = "TXNAME")
count_data <- data.frame(gene_id = cts$GENEID,
                         feature_id = cts$feature_id,
                         cts[, samples$sample])

# --- 5. DRIMSeq object ---
d <- dmDSdata(counts = count_data, samples = sample_info)

# --- 6. Filtering ---
n <- nrow(sample_info)
n.small <- min(table(sample_info$group))

d <- dmFilter(d,
              min_samps_feature_expr = n.small, min_feature_expr = 10,
              min_samps_feature_prop = n.small, min_feature_prop = 0.1,
              min_samps_gene_expr = n, min_gene_expr = 10)

# --- 7. Model fitting ---
design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = design_full)
d <- dmFit(d, design = design_full)
d <- dmTest(d, coef = 2)

# --- 8. Results with NA correction ---
res_gene <- DRIMSeq::results(d, level = "gene") %>%
  mutate(pvalue = if_else(is.na(pvalue), 1, pvalue),
         adj_pvalue = if_else(is.na(adj_pvalue), 1, adj_pvalue))

res_feature <- DRIMSeq::results(d, level = "feature") %>%
  mutate(pvalue = if_else(is.na(pvalue), 1, pvalue),
         adj_pvalue = if_else(is.na(adj_pvalue), 1, adj_pvalue))

# Low proportion variance filter
# Low proportion variance filter
smallProportionSD <- function(d, filter = 0.1) {
  cts <- as.matrix(subset(counts(d), select = -c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)), ]
  props <- cts / total.cts
  propSD <- sqrt(matrixStats::rowVars(props))
  propSD < filter
}

filt <- smallProportionSD(d)
res_feature$pvalue[filt] <- 1
res_feature$adj_pvalue[filt] <- 1

# --- 9. stageR analysis (MAINTAIN VERSIONS) ---
pScreen <- res_gene$pvalue
names(pScreen) <- res_gene$gene_id  # With versions

pConfirmation <- matrix(res_feature$pvalue, ncol = 1)
rownames(pConfirmation) <- res_feature$feature_id  # With versions

tx2gene_map <- res_feature[, c("feature_id", "gene_id")]  # With versions

stageR_obj <- stageRTx(pScreen = pScreen,
                       pConfirmation = pConfirmation,
                       pScreenAdjusted = FALSE,
                       tx2gene = tx2gene_map)

stageR_obj <- stageWiseAdjustment(stageR_obj, method = "dtu", alpha = 0.05)

results_stageR_all <- getAdjustedPValues(stageR_obj, order = FALSE, onlySignificantGenes = TRUE) %>%
  rename(adjP = transcript) %>%
  mutate(log10p = -log10(adjP + 1e-10),
         significant = adjP < 0.05)

# Write results (with versions)
write.csv(results_stageR_all, file = "results/dtu/results_stageR_significant_genes.csv", row.names = FALSE)

# Resultados em nível de transcrito (DTU)
results_stageR_tx <- results_stageR_all %>%
  filter(!is.na(txID)) %>%
  rename(
    feature_id = txID,
    gene_id = geneID
  )

results_stageR_tx <- results_stageR_tx %>%
  mutate(log10p = -log10(adjP + 1e-10),
         significant = adjP < 0.05)

volcano_plot <- ggplot(results_stageR_tx, aes(x = log10p, fill = significant)) +
  geom_histogram(bins = 50, color = "black") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), name = "Significativo") +
  labs(title = "Distribuição dos -log10(p-valor) ajustados (Transcritos)",
       x = "-log10(p-valor ajustado)", y = "Número de transcritos") +
  theme_minimal(base_size = 14)

ggsave("results/dtu/hist_adjP_transcripts.png", volcano_plot, width = 8, height = 6, dpi = 300)



sig_tx <- results_stageR_tx %>% filter(adjP < 0.05)

tx_by_gene <- sig_tx %>%
  count(gene_id, name = "Significant_Transcripts") %>%
  arrange(desc(Significant_Transcripts)) %>%
  head(30)

barplot_tx_gene <- ggplot(tx_by_gene, aes(x = reorder(gene_id, -Significant_Transcripts), y = Significant_Transcripts)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top 30 Genes com Mais Transcritos Significativos",
       x = "Gene", y = "Nº de Transcritos Significativos") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank())

ggsave("results/dtu/barplot_top30_genes_transcripts.png", barplot_tx_gene, width = 10, height = 7, dpi = 300)


#  ================= Segunda etapa =================
# Seleciona apenas os transcritos filtrados no objeto d
cts_filtered <- counts(d)

# Usa melt para transformar em formato longo
drim.prop <- reshape2::melt(cts_filtered, id = c("gene_id", "feature_id"),
                            variable.name = "sample", value.name = "count")

# Ordena para garantir consistência
drim.prop <- drim.prop[order(drim.prop$gene_id, drim.prop$sample, drim.prop$feature_id), ]

# Calcula as proporções de expressão de cada transcrito por gene em cada amostra
drim.prop <- drim.prop %>%
  group_by(gene_id, sample) %>%
  mutate(total = sum(count),
         prop = count / total) %>%
  ungroup()

write.csv(drim.prop, file = "results/dtu/drim_prop.csv", row.names = FALSE)

# Escolha um gene de interesse
gene_alvo <- "ENSG00000000003.16"

# Filtra para esse gene
drim_gene <- drim.prop %>% filter(gene_id == gene_alvo)

ggplot(drim_gene, aes(x = sample, y = prop, fill = feature_id)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ gene_id, scales = "free_x") +
  labs(title = paste("Proporção dos Transcritos do Gene", gene_alvo),
       x = "Amostra", y = "Proporção") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom")
# Exemplo: genes com transcritos significativos
genes_top <- unique(sig_tx$gene_id)[1:20] 

drim_top <- drim.prop %>%
  filter(gene_id %in% genes_top) %>%
  dplyr::select(gene_id, feature_id, sample, prop) %>%
  tidyr::unite("gene_transcript", gene_id, feature_id, sep = ":") %>%
  tidyr::pivot_wider(names_from = sample, values_from = prop, values_fill = 0)

mat <- as.matrix(drim_top[,-1])
rownames(mat) <- drim_top$gene_transcript

# Histograma dos -log10(p-valor ajustados) dos transcritos
volcano_plot <- ggplot(results_stageR_tx, aes(x = log10p, fill = significant)) +
  geom_histogram(bins = 50, color = "black") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), name = "Significativo") +
  labs(title = "Distribuição dos -log10(p-valor) ajustados (Transcritos)",
       x = "-log10(p-valor ajustado)", y = "Número de transcritos") +
  theme_minimal(base_size = 14)

ggsave("results/dtu/hist_adjP_transcripts.png", volcano_plot, width = 8, height = 6, dpi = 300)


# Barplot top 30 genes com mais transcritos significativos
barplot_tx_gene <- ggplot(tx_by_gene, aes(x = reorder(gene_id, -Significant_Transcripts), y = Significant_Transcripts)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top 30 Genes com Mais Transcritos Significativos",
       x = "Gene", y = "Nº de Transcritos Significativos") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank())

ggsave("results/dtu/barplot_top30_genes_transcripts.png", barplot_tx_gene, width = 10, height = 7, dpi = 300)


# Heatmap das proporções dos transcritos (top genes)
pheatmap(mat, scale = "row",
         main = "Heatmap - Proporção de transcritos (Top genes)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         fontsize_row = 7,
         filename = "results/dtu/heatmap_top_genes_transcripts.png")


# Proporção média dos transcritos por condição (barras)
barplot_gene_cond <- ggplot(drim_gene_group, aes(x = feature_id, y = mean_prop, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = paste("Proporção média dos transcritos por condição - gene", gene_alvo),
       x = "Transcrito", y = "Proporção média") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/dtu/barplot_mean_prop_by_condition.png", barplot_gene_cond, width = 10, height = 7, dpi = 300)

# Proporção acumulada dos transcritos ao longo das amostras (area plot)
area_prop_gene <- ggplot(drim_gene, aes(x = sample, y = prop, fill = feature_id)) +
  geom_area(alpha = 0.8 , position = 'stack') +
  labs(title = paste("Proporção acumulada dos transcritos - gene", gene_alvo),
       x = "Amostra", y = "Proporção acumulada") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("results/dtu/area_prop_transcripts_gene.png", area_prop_gene, width = 10, height = 6, dpi = 300)


gene_alvo <- "ENSG00000102081.16"  # FMR1 com versão;

# 2. Filtrar as proporções do gene
fmr1_prop <- drim.prop %>%
  filter(gene_id == gene_alvo)

# 3. Adicionar a condição de cada amostra
fmr1_prop <- fmr1_prop %>%
  left_join(sample_info, by = c("sample" = "sample_id"))

# 4. Criar boxplot
boxplot_fmr1 <- ggplot(fmr1_prop, aes(x = feature_id, y = prop, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.2), alpha = 0.7) +
  labs(title = "FMR1 isoforms proportions per condition",
       x = "Transcrito", y = "Expression proportion") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

# 5. Salvar a figura
ggsave("results/dtu/boxplot_FMR1_isoforms_by_condition.png", boxplot_fmr1, width = 10, height = 6, dpi = 300)