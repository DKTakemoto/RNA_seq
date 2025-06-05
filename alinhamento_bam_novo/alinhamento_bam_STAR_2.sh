#!/bin/bash
#PBS -N STAR_Align
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=12,walltime=48:00:00,vmem=100gb
#PBS -M daniel.takemoto@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/bam_star
#PBS -o /work/pancreas/takemoto/bam_star/log_reports/
#PBS -e /work/pancreas/takemoto/bam_star/log_reports_mapeamento/
#PBS -W group_list=pancreas
#PBS -V

mkdir -p /work/pancreas/takemoto/bam_star

# Criar arquivo samples_to_align.txt automaticamente
ls /work/pancreas/takemoto/trimmed/trimmed2/*_R1.fastq | sed 's/_R1.fastq//' | xargs -n 1 basename | sort | uniq > /work/pancreas/takemoto/bam_star/samples_to_align.txt

set -euo pipefail

GENOME_DIR="/work/pancreas/takemoto/star_index/anotacao_primaria"
FASTQ_DIR="/work/pancreas/takemoto/trimmed/trimmed2"
OUTPUT_DIR="/work/pancreas/takemoto/alinhamento_bam_novo"

mkdir -p "$OUTPUT_DIR/log_reports" "$OUTPUT_DIR/log_reports_mapeamento"

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
done < /work/pancreas/takemoto/bam_star/samples_to_align.txt
