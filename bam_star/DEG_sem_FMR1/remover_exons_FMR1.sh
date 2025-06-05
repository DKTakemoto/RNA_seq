#!/bin/bash

# Caminho para o GTF original
GTF_ORIGINAL="/work/pancreas/takemoto/RNA_seq/data/reference/database/gencode.v47.annotation.gtf"

# Nome do GTF de saída
GTF_FILTRADO="/work/pancreas/takemoto/RNA_seq/data/reference/database/gencode.v47.FMR1_exon13-16_removed.gtf"

awk '!/gene_name "FMR1";/ || ($0 !~ /exon_number 13[;]/ && $0 !~ /exon_number 14[;]/ && $0 !~ /exon_number 15[;]/ && $0 !~ /exon_number 16[;]/)' "$GTF_ORIGINAL" > "$GTF_FILTRADO"

echo "terminou"