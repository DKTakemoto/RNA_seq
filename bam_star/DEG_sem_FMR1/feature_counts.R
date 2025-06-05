#!/usr/bin/env Rscript

library(Rsubread)
# Diretório raiz
root_dir <- "/work/pancreas/takemoto/RNA_seq/"


# Caminhos
bam_dir <- file.path(root_dir, "data", "bam", "bam_star")
gtf_path <- file.path(root_dir, "data", "reference", "database", "gencode.v47.FMR1_exon13-16_removed.gtf")

results_dir <- file.path(root_dir, "results", "feature_counts")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
output_file <- file.path(results_dir, "counts_FMR1_exons_removidos.csv")

# 1. Verificar arquivo GTF
if(!file.exists(gtf_path)) {
  stop("ERRO: Arquivo GTF não encontrado em:\n", gtf_path)
}

cat("\n=== Verificação do GTF ===\n")
cat("Tamanho do arquivo:", file.size(gtf_path), "bytes\n")
cat("Primeiras linhas:\n")
system(paste("head -n 2", gtf_path))

# 2. Listar arquivos BAM
bamFiles <- list.files(path = bam_dir, 
                       pattern = ".sortedByCoord.out.bam$", 
                       full.names = TRUE)

if(length(bamFiles) == 0) {
  stop("Nenhum arquivo BAM encontrado em: ", bam_dir)
}

# 3. Executar featureCounts
fcResults <- featureCounts(
  files = bamFiles,
  annot.ext = gtf_path,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  strandSpecific = 0,
  countReadPairs = TRUE,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  countMultiMappingReads = FALSE,
  countChimericFragments = FALSE,
  nthreads = 8  # Número mais seguro de threads
)

# 4. Salvar resultados
write.csv(fcResults$counts, output_file)

cat("\nAnálise concluída com sucesso!\n")
cat("Resultados salvos em:", output_file, "\n")