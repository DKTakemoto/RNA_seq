#!/bin/bash
#PBS -N Bowtie2_index
#PBS -l nodes=1:ppn=30,walltime=200:00:00,vmem=200gb
#PBS -M raunimarques@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/bowtie_index/
#PBS -o /work/pancreas/takemoto/bowtie_index/log_reports/
#PBS -e /work/pancreas/takemoto/bowtie_index/log_reports_mapeamento/
#PBS -W group_list=pancreas
#PBS -V

bowtie2-build /work/pancreas/takemoto/database/GRCh38.primary_assembly.genome.fa /work/pancreas/takemoto/bowtie_index
