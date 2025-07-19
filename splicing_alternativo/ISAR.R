library(IsoformSwitchAnalyzeR)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(gprofiler2)
library(VennDiagram)
library(grid)

# --- 1. Import samples and metadata ---
samples <- read.table(
  "/work/pancreas/takemoto/RNA_seq/data/reference/database/salmon_metadata.txt",
  header = TRUE,
  sep = ","
)
isoformData <- importIsoformExpression(parentDir = "/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/salmon")

samples$sample_id <- basename(dirname(samples$path))
designMatrix <- data.frame(
  sampleID = samples$sample_id,
  condition = samples$condition
)

# --- 2. Import to switchAnalyzeRlist ---
aSwitchList <- importRdata(
  isoformCountMatrix = isoformData$counts,
  isoformRepExpression = isoformData$abundance,
  designMatrix = designMatrix,
  isoformExonAnno = "/work/pancreas/takemoto/RNA_seq/data/reference/database/Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf.gz",
  removeNonConvensionalChr = TRUE,
  ignoreAfterPeriod = TRUE,
  showProgress = TRUE
)

# --- 2.1 Remove specific transcripts (optional cleanup) ---
transcripts_to_remove <- c("ENST00000478848", "ENST00000611273")
aSwitchList <- subsetSwitchAnalyzeRlist(
  switchAnalyzeRlist = aSwitchList,
  !aSwitchList$isoformFeatures$isoform_id %in% transcripts_to_remove
)

# --- 3. Filtering ---
aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 1,
  removeSingleIsoformGenes = TRUE
)

# --- 4. Isoform switch test ---
switchTestResults <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered,
  reduceToSwitchingGenes = TRUE
)

switchTestResults <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = switchTestResults,
  onlySwitchingGenes = TRUE,
  alpha = 0.05,
  dIFcutoff = 0.1
)

# --- 5. Subset for a condition ---
SwitchListAnalyzedSubset <- switchTestResults

# --- 6. Extract top switching genes ---
topQgenes <- extractTopSwitches(
  SwitchListAnalyzedSubset, 
  n = 10,
  sortByQvals = TRUE,
)
write.csv(topQgenes, file = "/work/pancreas/takemoto/RNA_seq/results/ISAR/top_Q_genes_switches.csv", row.names = FALSE)

topdIFgenes <- extractTopSwitches(
  SwitchListAnalyzedSubset, 
  n = 10, 
  sortByQvals = FALSE,
)
write.csv(topdIFgenes, file = "/work/pancreas/takemoto/RNA_seq/results/ISAR/top_dIF_genes_switches.csv", row.names = FALSE)

# --- 7. All switching isoforms ---
switchingIso <- extractTopSwitches( 
  SwitchListAnalyzedSubset, 
  n = NA,                  
  extractGenes = FALSE,    
  sortByQvals = TRUE
)

write.csv(
  switchingIso,
  file = "/work/pancreas/takemoto/RNA_seq/results/ISAR/switchingIso_all.csv",
  row.names = FALSE
)


# --- 8. Splicing summary and enrichment ---
# --- Splicing Summary ---
pdf("/work/pancreas/takemoto/RNA_seq/results/ISAR/splicing_summary.pdf")
extractSplicingSummary(
  switchTestResults,
  asFractionTotal = FALSE,
  plotGenes = TRUE
)
dev.off()

# ---  Splicing Enrichment ---
pdf("/work/pancreas/takemoto/RNA_seq/results/ISAR/splicing_enrichment.pdf")
splicingEnrichment <- extractSplicingEnrichment(
  switchTestResults,
  splicingToAnalyze = 'all',
  returnResult = TRUE,
  returnSummary = TRUE
)
dev.off()

# --- Splicing Genome Wide ---
pdf("/work/pancreas/takemoto/RNA_seq/results/ISAR/splicing_genome_wide.pdf")
extractSplicingGenomeWide(
  switchTestResults,
  featureToExtract = 'all',
  splicingToAnalyze = c('A3','MES','ATSS'),
  plot = TRUE,
  returnResult = FALSE
)
dev.off()

# --- 9. Mechanism behind isoform switches ---
consequenceResults <- analyzeSwitchConsequences(
  SwitchListAnalyzedSubset,
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)

myConsequences <- consequenceResults$switchConsequence
myConsequences <- myConsequences[myConsequences$isoformsDifferent == TRUE, ]
myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated)

bioMechanismeAnalysis <- bioMechanismeAnalysis[bioMechanismeAnalysis$isoformsDifferent == TRUE, ]
bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)

bioMechanismeAnalysis <- bioMechanismeAnalysis[bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair, ]

AS   <- bioMechanismeAnalysis$isoPair[bioMechanismeAnalysis$featureCompared == 'intron_structure']
aTSS <- bioMechanismeAnalysis$isoPair[bioMechanismeAnalysis$featureCompared == 'tss']
aTTS <- bioMechanismeAnalysis$isoPair[bioMechanismeAnalysis$featureCompared == 'tts']

mechList <- list(
  AS = AS,
  aTSS = aTSS,
  aTTS = aTTS
)

# --- 10. Venn diagram ---
myVenn <- venn.diagram(
  x = mechList,
  col = 'transparent',
  alpha = 0.4,
  fill = RColorBrewer::brewer.pal(n = 3, name = 'Dark2'),
  filename = NULL
)

output_folder <- "/work/pancreas/takemoto/RNA_seq/results/ISAR"
grid.newpage()
grid.draw(myVenn)
venn_file <- file.path(output_folder, "venn_diagram_switch_mechanisms.pdf")
pdf(venn_file, width = 6, height = 6)
grid.draw(myVenn)
dev.off()


#---------------------------------------------------------------#
# --- VISUALIZAÇÃO E EXPORTAÇÃO DE RESULTADOS ---

# Volcano plot
volcano <- ggplot(data = switchTestResults$isoformFeatures, aes(x = dIF, y = -log10(isoform_switch_q_value))) +
  geom_point(
    aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
    size = 1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dashed') +
  facet_wrap('condition_2') + 
  scale_color_manual('Significant\nIsoform Switch', values = c('black', 'red')) +
  labs(x = 'dIF', y = '-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

# Salve o gráfico como PDF
ggsave("/work/pancreas/takemoto/RNA_seq/results/ISAR/isoform_switch_volcano_plot.pdf", plot = volcano, width = 8, height = 6)