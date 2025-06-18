#!/bin/bash
#PBS -N DRIM_seq           # Nome do job
#PBS -l walltime=10:00:00     # Tempo máximo (4 horas)
#PBS -l vmem=100gb              # Memória RAM
#PBS -l nodes=1:ppn=4         # 1 nó com 4 CPUs
#PBS -o /work/pancreas/takemoto/RNA_seq/logs/logs_dtu/drimseq_out.log       # Arquivo de saída
#PBS -e /work/pancreas/takemoto/RNA_seq/logs/logs_dtu/drimseq_err.log       # Arquivo de erro
#PBS -M daniel.takemoto@usp.br
#PBS -m abe
#PBS -V                       # Exporta variáveis de ambiente

cd /work/pancreas/takemoto/RNA_seq/scripts/splicing_alternativo/salmon
# Executa o script R
Rscript dtu_drim_seq.R