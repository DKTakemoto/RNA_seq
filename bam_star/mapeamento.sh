#!/bin/bash
#PBS -N STAR
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=12,walltime=120:00:00,vmem=200gb
#PBS -M raunimarques@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/bam_star
#PBS -o /work/pancreas/takemoto/bam_star/log_reports
#PBS -e /work/pancreas/takemoto/bam_star/log_reports_mapeamento
#PBS -W group_list=pancreas
#PBS -V

for i in `cat files4.txt`; do
    STAR --runThreadN 12 \
         --limitGenomeGenerateRAM 20000000000 \
         --genomeDir /work/pancreas/takemoto/star_index/anotacao_primaria \
         --readFilesIn /work/pancreas/takemoto/bbtols/trimmed/${i}_R1.fastq /work/pancreas/takemoto/bbtols/trimmed/${i}_R2.fastq \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${i}_STAR_Hs
done
