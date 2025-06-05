#!/bin/bash
#PBS -N STAR_Align
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=12,walltime=48:00:00,vmem=100gb
#PBS -M daniel.takemoto@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/RNA_seq/
#PBS -o /work/pancreas/takemoto/RNA_seq/logs/bam/log_stdout/
#PBS -e /work/pancreas/takemoto/RNA_seq/logs/bam/log_stderr/
#PBS -W group_list=pancreas
#PBS -V

set -euo pipefail
ROOT_DIR=$(pwd)

GENOME_DIR="/data/reference/star_index/anotacao_primaria"
FASTQ_DIR="/data/trimmed_fastq/fastp_trim/"
OUTPUT_DIR="/data/bam/bam_star"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/log_reports" "$OUTPUT_DIR/log_reports_mapeamento"

SAMPLES_FILE="$OUTPUT_DIR/samples_to_align.txt"

# Criar arquivo samples_to_align.txt automaticamente
ls "$ROOT_DIR/data/trimmed_fastq/fastp_trim/"*"_R1.fastq" | sed 's/_R1.fastq//' | xargs -n 1 basename | sort | uniq > "$SAMPLES_FILE"


THREADS=12

while read -r sample; do
  echo "Iniciando alinhamento da amostra: $sample"

  STAR --runThreadN $THREADS \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$FASTQ_DIR/${sample}_R1.fastq" "$FASTQ_DIR/${sample}_R2.fastq" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$OUTPUT_DIR/${sample}_STAR_" \
    --outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA

  echo "Finalizado alinhamento da amostra: $sample"
done < "$SAMPLES_FILE"
