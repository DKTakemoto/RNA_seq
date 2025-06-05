#!/usr/bin/env Rscript

library(Rsubread)

# Caminhos
bam_dir <- "/work/pancreas/takemoto/alinhamento_bam_novo"
gtf_path <- "/work/pancreas/takemoto/bam_star/DEG_sem_FMR1/gencode.v47.FMR1_exon13-16_removed.gtf"

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
output_file <- "counts_FMR1_exons_removidos2.csv"
write.csv(fcResults$counts, output_file)

cat("\nAnálise concluída com sucesso!\n")
cat("Resultados salvos em:", output_file, "\n")