#!/bin/bash
#PBS -N rmats_analysis
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -l walltime=24:00:00
#PBS -q batch
#PBS -j oe
#PBS -o /work/pancreas/takemoto/RNA_seq/logs/rmats_output.log
#PBS -V

# Ir para o diretório de onde o script foi submetido
cd $PBS_O_WORKDIR

# Ativar o conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rmats_env

# Diretórios e arquivos
BASE_DIR="/work/pancreas/takemoto/RNA_seq"
GTF="$BASE_DIR/data/reference/database/gencode.v47.annotation.gtf"
BAM_DIR="$BASE_DIR/data/bam/bam_star"
OUTPUT_DIR="$BASE_DIR/results/splicing_alternativo/rmats"

READ_LENGTH=101
NTHREADS=8
LIBTYPE="fr-unstranded"

# Criar arquivos com listas de BAMs
ls ${BAM_DIR}/*N*_STAR_Aligned.sortedByCoord.out.bam | paste -sd "," - > grupo_N.txt
ls ${BAM_DIR}/*C3*_STAR_Aligned.sortedByCoord.out.bam | paste -sd "," - > grupo_C3.txt


# Rodar comparação C3 vs N
echo "Iniciando rMATS: C3 vs N"
mkdir -p "$OUTPUT_DIR/output_C3_vs_N" "$OUTPUT_DIR/tmp_C3_vs_N"
rmats.py \
  --b1 grupo_N.txt \
  --b2 grupo_C3.txt \
  --gtf "$GTF" \
  --od "$OUTPUT_DIR/output_C3_vs_N" \
  --tmp "$OUTPUT_DIR/tmp_C3_vs_N" \
  --readLength $READ_LENGTH \
  --nthread $NTHREADS \
  --libType $LIBTYPE

echo "Análises rMATS concluídas!"
