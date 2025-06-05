#!/bin/bash

# Define o diretório com os arquivos FASTQ e o diretório de saída
fastq_dir="/work/pancreas/takemoto/rna_circ/S747_DanielTakamoto_mRNALigation-423566592/BCLConvert_07_05_2024_12_54_28Z-746216483"  # Altere para o caminho onde estão os arquivos .fastq

# Cria o diretório de saída se ainda não existir
mkdir -p "$fastq_dir/trimmed"

# Loop para processar cada arquivo FASTQ com "_R2_001.fastq" no nome
for fastq_file_2 in "$fastq_dir"/*_R2_001.fastq; do
    # Extrai o nome base do arquivo para nomear o arquivo de saída
    dir=$(basename "$fastq_file_2" "_R2_001.fastq")

    # Comando fastp para processar apenas o arquivo R2
    fastp -w 16 \
          -h "$fastq_dir/trimmed/${dir}_fastp_R2.html" \
          -j "$fastq_dir/trimmed/${dir}_fastp_R2.json" \
          -i "$fastq_file_2" \
          -o "$fastq_dir/trimmed/trimmed_${dir}_R2.fastq"
done

