#!/bin/bash
#PBS -N STAR_Index
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=70,walltime=100:00:00,vmem=300gb
#PBS -M daniel.takemoto@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/star_index/anotacao_primaria/
#PBS -o /work/pancreas/takemoto/star_index/anotacao_primaria/log_reports/
#PBS -e /work/pancreas/takemoto/star_index/anotacao_primaria/log_reports_mapeamento/
#PBS -W group_list=pancreas
#PBS -V


STAR --runThreadN 70 \
     --runMode genomeGenerate \
     --limitGenomeGenerateRAM 20000000000 \
     --genomeDir /work/pancreas/takemoto/star_index/anotacao_primaria \
     --genomeFastaFiles /work/pancreas/takemoto/database/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile /work/pancreas/takemoto/database/gencode.v47.basic.annotation.gtf \
     --sjdbOverhang 149
