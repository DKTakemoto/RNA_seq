#!/bin/bash
set -e

# Caminhos
TRANSCRIPTOME_DIR="./transcriptoma_database"
FASTQ_DIR="./trimmed2"
KALLISTO_INDEX="$TRANSCRIPTOME_DIR/kallisto_index.idx"
SALMON_INDEX="$TRANSCRIPTOME_DIR/salmon_index"

# Arquivo cDNA compactado
CDNA_FASTA=$(ls $TRANSCRIPTOME_DIR/*.fa.gz | head -n 1)

# Criar índice para Kallisto
echo "Criando índice para Kallisto..."
kallisto index -i "$KALLISTO_INDEX" "$CDNA_FASTA"

# Rodar Kallisto para todos os pares de arquivos
echo "Rodando Kallisto..."
for fq1 in $FASTQ_DIR/*_R1.fastq; do
    	fq2="${fq1/_R1.fastq/_R2.fastq}"
    	sample=$(basename "$fq1" | cut -d'_' -f2-4)
	if [[ "$sample" == C* || "$sample" == N* ]]; then
		mkdir -p "kallisto_out/$sample"
    		kallisto quant -i "$KALLISTO_INDEX" -o "kallisto_out/$sample" -b 100 "$fq1" "$fq2"
	else
		echo "Arquivo $fq1 não segue padrão"
	fi
done

echo "Pipeline concluída!"

