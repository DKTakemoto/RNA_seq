#!/bin/bash
#PBS -N salmon_quant
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l vmem=70gb
#PBS -o /work/pancreas/takemoto/RNA_seq/logs/salmon_quant.out
#PBS -e /work/pancreas/takemoto/RNA_seq/logs/salmon_quant.err
#PBS -V
#PBS -M daniel.takemoto@usp.br
#PBS -m abe

# Caminhos e parâmetros
BASE_DIR="/work/pancreas/takemoto/RNA_seq"
TRANSCRIPTOME_DIR="$BASE_DIR/data/reference/database"
FASTQ_DIR="$BASE_DIR/data/trimmed_fastq/fastp_trim"
SALMON_INDEX="$BASE_DIR/results/splicing_alternativo/salmon/salmon_index"
OUTPUT_DIR="$BASE_DIR/results/splicing_alternativo/salmon"
THREADS=8

# Arquivo de referência
CDNA_FASTA=$(ls "$TRANSCRIPTOME_DIR"/*.fa.gz | head -n 1)

# Criar índice do Salmon (assume que será criado sempre)
echo "Criando índice para Salmon..."
salmon index -t "$CDNA_FASTA" -i "$SALMON_INDEX"

# Criar diretório de saída
mkdir -p "$OUTPUT_DIR"

# Processamento de amostras
echo "Iniciando quantificação com Salmon em $(date)"
for fq1 in "$FASTQ_DIR"/*_R1.fastq; do
fq2="${fq1/_R1.fastq/_R2.fastq}"
sample=$(basename "$fq1" | cut -d'_' -f2-4)
sample_output_dir="$OUTPUT_DIR/$sample"

echo "Processando $sample..."
mkdir -p "$sample_output_dir"

# -l A tells automatically if unstranded or stranded
# -p ammount of threads
# -o output dir
# --gcBias to learn and correct for fragment-level GC biases in the input data
salmon quant -i "$SALMON_INDEX" -l A \
-1 "$fq1" -2 "$fq2" \
-p $THREADS -o "$sample_output_dir" \
--gcBias --validateMappings

echo "Concluído $sample"
done

echo "Quantificação finalizada em $(date)"
