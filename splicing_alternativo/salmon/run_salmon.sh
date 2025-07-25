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

# Caminhos
BASE_DIR="/work/pancreas/takemoto/RNA_seq"
TRANSCRIPTOME_DIR="$BASE_DIR/data/reference/database"
FASTQ_DIR="$BASE_DIR/data/trimmed_fastq/fastp_trim"
SALMON_INDEX="$BASE_DIR/results/splicing_alternativo/salmon/salmon_index"
OUTPUT_DIR="$BASE_DIR/results/splicing_alternativo/salmon"
THREADS=8

# Verificar diretórios de entrada
if [ ! -d "$TRANSCRIPTOME_DIR" ]; then
    echo "Erro: Diretório $TRANSCRIPTOME_DIR não encontrado" >&2
    exit 1
fi

if [ ! -d "$FASTQ_DIR" ]; then
    echo "Erro: Diretório $FASTQ_DIR não encontrado" >&2
    exit 1
fi

# Arquivo cDNA compactado
CDNA_FASTA=$(ls "$TRANSCRIPTOME_DIR"/gencode.v47.transcripts.fa.gz 2>/dev/null | head -n 1)

if [ -z "$CDNA_FASTA" ]; then
    echo "Erro: Nenhum arquivo .fa.gz encontrado em $TRANSCRIPTOME_DIR" >&2
    exit 1
fi

# Criação de índice:
if [ ! -d "$SALMON_INDEX" ]; then
    echo "Criando índice para Salmon..."
    salmon index -t "$CDNA_FASTA" -i "$SALMON_INDEX"
else
    echo "Índice para Salmon já existe. Pulando criação de índice."
fi

# Criar diretório de saída
mkdir -p "$OUTPUT_DIR"

# Função para verificar se a saída da amostra está completa
check_salmon_complete() {
    local sample_dir="$1"
    [[ -s "$sample_dir/quant.sf" && -s "$sample_dir/cmd_info.json" && -s "$sample_dir/lib_format_counts.json" ]]
}

# Rodar Salmon para todos os pares de arquivos
echo "Iniciando processamento com Salmon em $(date)"
for fq1 in "$FASTQ_DIR"/*_R1.fastq; do
    fq2="${fq1/_R1.fastq/_R2.fastq}"
    sample=$(basename "$fq1" | cut -d'_' -f2-4)
    sample_output_dir="$OUTPUT_DIR/$sample"

    if [[ "$sample" =~ ^(C|N)[0-9]+ ]]; then
        if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
            echo "Erro: Arquivos FASTQ para amostra $sample não encontrados" >&2
            continue
        fi

        if check_salmon_complete "$sample_output_dir"; then
            echo "Amostra $sample já processada. Pulando..."
            continue
        else
            echo "Processando $sample com Salmon em $(date)..."
            rm -rf "$sample_output_dir"
            mkdir -p "$sample_output_dir"

            if ! salmon quant -i "$SALMON_INDEX" -l A -1 "$fq1" -2 "$fq2" \
                 -p $THREADS -o "$sample_output_dir" \
                 --gcBias --validateMappings; then
                echo "Erro ao processar $sample. Pulando para a próxima." >&2
                continue
            fi
            echo "Concluído processamento de $sample em $(date)"
        fi
    else
        echo "Arquivo $fq1 não segue o padrão esperado (C* ou N*). Pulando..."
    fi
done

echo "Processamento concluído em $(date)"
