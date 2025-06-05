library(Rsubread)
bamFiles <- list.files(path = "/work/pancreas/takemoto/bam_star/", pattern = ".sortedByCoord.out.bam$", full.names = T)

fcResultsBowtie <- featureCounts(#vetor que contém os caminhos para os arquivos bam
  bamFiles,
  annot.ext = "/work/pancreas/takemoto/database/gencode.v47.basic.annotation.gtf", 
  isGTFAnnotationFile = T,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  strandSpecific = 0,
  countReadPairs = T,
  isPairedEnd = T,
  requireBothEndsMapped = T,
  countMultiMappingReads = F,
  countChimericFragments = F,
  nthreads = 30)

counts <- fcResultsBowtie[["counts"]]
write.csv(counts, "counts_primário.csv")