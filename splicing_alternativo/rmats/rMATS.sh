#!/bin/bash
# Defina o caminho do GTF e o read length aqui:
GTF="gencode.v47.basic.annotation.gtf"
READ_LENGTH=101 #samtools view C21_S747_S16_L001_STAR_HsAligned.sortedByCoord.out.bam | head -n 1 | awk '{print length($10)}'

# Número de threads:
NTHREADS=8

# Lib type: fr-unstranded, fr-firststrand, fr-secondstrand
LIBTYPE="fr-unstranded"

# Diretório atual com os BAMs
BAM_DIR="$PWD/bam_files"

# Criar arquivos com listas de BAMs
ls ${BAM_DIR}/N*_STAR_HsAligned.sortedByCoord.out.bam | paste -sd "," - > grupo_N.txt
ls ${BAM_DIR}/C2*_STAR_HsAligned.sortedByCoord.out.bam | paste -sd "," - > grupo_C2.txt
ls ${BAM_DIR}/C3*_STAR_HsAligned.sortedByCoord.out.bam | paste -sd "," - > grupo_C3.txt

# Rodar C2 vs N
echo "Iniciando rMATS: C2 vs N"
mkdir -p output_C2_vs_N tmp_C2_vs_N
rmats.py \
--b1 grupo_N.txt \
--b2 grupo_C2.txt \
--gtf "$GTF" \
--od output_C2_vs_N \
--tmp tmp_C2_vs_N \
--readLength $READ_LENGTH \
--nthread $NTHREADS \
--libType $LIBTYPE

# Rodar C3 vs N
echo "Iniciando rMATS: C3 vs N"
mkdir -p output_C3_vs_N tmp_C3_vs_N
rmats.py \
--b1 grupo_N.txt \
--b2 grupo_C3.txt \
--gtf "$GTF" \
--od output_C3_vs_N \
--tmp tmp_C3_vs_N \
--readLength $READ_LENGTH \
--nthread $NTHREADS \
--libType $LIBTYPE

echo "Análises rMATS concluídas!"
