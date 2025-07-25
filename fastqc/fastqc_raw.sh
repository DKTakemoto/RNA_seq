#!/bin/bash
#PBS -N fastqc                         
#PBS -l nodes=1:ppn=8                      
#PBS -l walltime=08:00:00                  
#PBS -l mem=20gb                           
#PBS -l vmem=25gb                          
#PBS -j oe                                 
#PBS -V                                    
#PBS -M daniel.takemoto@usp.br             
#PBS -m abe                                # Envia email em: (a) início, (b) fim, (e) erro

cd /work/pancreas/takemoto/RNA_seq/data/raw_fastq
fastqc *.fastq -o /work/pancreas/takemoto/RNA_seq/results/quality_check
echo "FastQC finalizado em $(date)"