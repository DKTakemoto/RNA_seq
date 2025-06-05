#!/bin/bash
#PBS -N Bowtie2_mapeamento
#PBS -l nodes=2:ppn=30,walltime=200:00:00,vmem=200gb
#PBS -M raunimarques@usp.br
#PBS -m bea
#PBS -d /work/pancreas/takemoto/bam_bowtie
#PBS -o /work/pancreas/takemoto/bam_bowtie/log_reports
#PBS -e /work/pancreas/takemoto/bam_bowtie/log_reports_mapeamento
#PBS -W group_list=pancreas
#PBS -V

# Diretórios
input_dir="/work/pancreas/takemoto/rna_circ/S747_DanielTakamoto_mRNALigation-423566592/BCLConvert_07_05_2024_12_54_28Z-746216483/trimmed"
output_dir="/work/pancreas/takemoto/bam_bowtie"
bwt2IdxPath="/work/pancreas/takemoto/bowtie_index"  # Atualize para o caminho correto do índice Bowtie2

# Número de threads
threads=30  # Ajuste o número de threads conforme necessário

# Loop pelos arquivos _R1.fastq
for file_r1 in "$input_dir"/*_L001_R1.fastq; do
    # Nome base sem a extensão e sem o sufixo _L001_R1.fastq
    base_name=$(basename "$file_r1" _L001_R1.fastq)

    # Arquivo _R2 correspondente
    file_r2="${input_dir}/${base_name}_L001_R2.fastq"

    # Verificar se o arquivo _R2 existe
    if [[ -f "$file_r2" ]]; then
        echo "Processando: $base_name"
        # Saída SAM temporária
        sam_output="${output_dir}/${base_name}.sam"

        # Alinhamento com Bowtie2
        bowtie2 -p "$threads" -x "$bwt2IdxPath" -1 "$file_r1" -2 "$file_r2" -S "$sam_output" 2> "${output_dir}/${base_name}_bowtie2.log"

        # Converter SAM para BAM
        bam_output="${output_dir}/${base_name}.bam"
        samtools view -Sb "$sam_output" > "$bam_output"

        # Ordenar BAM
        sorted_bam="${output_dir}/${base_name}_sorted.bam"
        samtools sort "$bam_output" -o "$sorted_bam"

        # Indexar BAM
        samtools index "$sorted_bam"

        # Remover o SAM e o BAM intermediário
        rm -f "$sam_output" "$bam_output"
    else
        echo "Arquivo correspondente _R2 para $file_r1 não encontrado, ignorando."
    fi
done

# Mensagem final correta
echo "Processamento concluído!"
