#!/bin/bash
#PBS -N BBTools
#PBS -l nodes=1:ppn=20,walltime=200:00:00,vmem=140gb
#PBS -M raunimarques@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/bbtols/
#PBS -o /work/pancreas/takemoto/bbtols/log_reports/
#PBS -e /work/pancreas/takemoto/bbtols/log_reports_mapeamento/
#PBS -W group_list=pancreas
#PBS -V


# Diretórios de origem e destino
input_dir="/work/pancreas/takemoto/fastq"
output_dir="/work/pancreas/takemoto/bbtols"

# Arquivo com os nomes dos arquivos fastq
files="files.txt"

# Ler os nomes de arquivos do files.txt em pares (R1 e R2)
while read -r r1 && read -r r2; do
    # Gerar caminhos completos para os arquivos de entrada
    in1="${input_dir}/${r1}"
    in2="${input_dir}/${r2}"

    # Gerar nomes para os arquivos de saída
    out1="${output_dir}/${r1}"
    out2="${output_dir}/${r2}"
    outs="${output_dir}/${r1%%_R1_001.fastq}_singletons.fastq"

    # Executar o comando bbtools
    /usr/lib/bbtools/repair.sh \
        in1="${in1}" \
        in2="${in2}" \
        out1="${out1}" \
        out2="${out2}" \
        outs="${outs}"
done < "$files"
