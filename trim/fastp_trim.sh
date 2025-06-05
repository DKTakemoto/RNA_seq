#!/bin/bash

# Define o diretório com os arquivos FASTQ
fastq_dir="/data/raw_fastq"

# Define o diretório de saída
output_dir="/data/trimmed_fastq/fastp_trim"

# Cria o diretório de saída, se não existir
mkdir -p "$output_dir"

# Loop para processar cada par de arquivos FASTQ emparelhados
for fastq_file_1 in "$fastq_dir"/*_R1.fastq; do
    # Define o arquivo R2 correspondente
    fastq_file_2="${fastq_file_1/_R1.fastq/_R2.fastq}"

    # Extrai o nome base (com sufixo _R1) para manter no nome de saída
    base_R1=$(basename "$fastq_file_1" ".fastq")   # ex: C21_S747_S16_L001_R1
    base_R2=$(basename "$fastq_file_2" ".fastq")   # ex: C21_S747_S16_L001_R2

    echo "Processando: $base_R1 e $base_R2"

    # Comando fastp
    fastp -w 16 \
          -h "$output_dir/${base_R1%_R1}_fastp.html" \
          -j "$output_dir/${base_R1%_R1}_fastp.json" \
          -i "$fastq_file_1" \
          -I "$fastq_file_2" \
          -o "$output_dir/trimmed_${base_R1}.fastq" \
          -O "$output_dir/trimmed_${base_R2}.fastq"
done
