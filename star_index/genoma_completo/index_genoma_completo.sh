#!/bin/bash
#PBS -N STAR
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=30,walltime=100:00:00,vmem=200gb
#PBS -M raunimarques@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/star_index/genoma_completo
#PBS -o /work/pancreas/takemoto/star_index/genoma_completo/log_reports
#PBS -e /work/pancreas/takemoto/star_index/genoma_completo/log_reports_mapeamento
#PBS -W group_list=pancreas
#PBS -V


STAR --runThreadN 30 \
     --runMode genomeGenerate \
     --limitGenomeGenerateRAM 20000000000 \
     --genomeDir /work/pancreas/takemoto/star_index/genoma_completo \
     --genomeFastaFiles /work/pancreas/takemoto/database/GRCh38.p14.genome.fa \
     --sjdbGTFfile /work/pancreas/takemoto/database/gencode.v47.annotation.gtf \
     --sjdbOverhang 149
