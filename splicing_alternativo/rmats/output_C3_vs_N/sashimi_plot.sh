#!/bin/bash

# Diretórios configuráveis
BASE_DIR="/work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/rmats"
GRUPO_C3="$BASE_DIR/grupo_C3.txt"
GRUPO_N="$BASE_DIR/grupo_N.txt"
SASHIMI_DIR="$BASE_DIR/output_C3_vs_N/sashimi_plots"
INPUT_DIR="$BASE_DIR/output_C3_vs_N"

# Lista de eventos do rMATS
event_types=("SE" "RI" "MXE" "A3SS" "A5SS")

# Ativar ambiente Conda
source /home/takemoto/miniconda3/etc/profile.d/conda.sh
conda activate sashimi_env

# Criar diretório de saída
mkdir -p "$SASHIMI_DIR"

# Ler os caminhos dos BAMs e gerar strings separadas por vírgula
B1=$(cat "$BASE_DIR/grupo_C3.txt")
B2=$(cat "$BASE_DIR/grupo_N.txt")

# Rótulos
L1="C3"
L2="Controle"

# Loop por tipo de evento
for event in "${event_types[@]}"; do
    infile="${INPUT_DIR}/${event}.MATS.JC.txt"

    if [[ ! -f "$infile" ]]; then
        echo "Arquivo ${infile} não encontrado. Pulando ${event}..."
        continue
    fi

    echo "Processando evento: $event"

    filtered_out="${BASE_DIR}/${event}_filtered.txt"
    top10_out="${BASE_DIR}/${event}_top10.txt"

    # Cabeçalho
    head -n 1 "$infile" > "$filtered_out"

    # Filtro de eventos
    awk -F'\t' '
    NR > 1 {
        split($22, psi1, ","); split($23, psi2, ",");

        c1=0; s1=0;
        for(i in psi1) { if(psi1[i] != "NA") { s1+=psi1[i]; c1++ } }
        avg1 = (c1>0) ? s1/c1 : 0

        c2=0; s2=0;
        for(i in psi2) { if(psi2[i] != "NA") { s2+=psi2[i]; c2++ } }
        avg2 = (c2>0) ? s2/c2 : 0

        dpsi = (avg1 - avg2); if (dpsi < 0) dpsi = -dpsi

        fdr = $20 + 0
        r1 = $14 + 0; r2 = $15 + 0

        if (r1 >= 10 && r2 >= 10 &&
            avg1 >= 0.05 && avg1 <= 0.95 &&
            avg2 >= 0.05 && avg2 <= 0.95 &&
            fdr <= 0.01 && dpsi >= 0.2)
        {
            print $0
        }
    }' "$infile" >> "$filtered_out"

    # Top 10 por FDR
    (head -n 1 "$infile" && tail -n +2 "$filtered_out" | sort -k20,20g | head -n 10) > "$top10_out"

    # Diretórios de saída
    out_fdr_dir="${SASHIMI_DIR}/${event}_FDR"
    out_top10_dir="${SASHIMI_DIR}/${event}_TOP10"
    mkdir -p "$out_fdr_dir" "$out_top10_dir"

    # Sashimi plots filtrados
    echo "Gerando SashimiPlot (eventos filtrados) para $event..."
    rmats2sashimiplot \
        --b1 "$B1" \
        --b2 "$B2" \
        -t "$event" \
        -e "$filtered_out" \
        --l1 "$L1" \
        --l2 "$L2" \
        -o "$out_fdr_dir"

    # Sashimi plots top 10
    echo "Gerando SashimiPlot (top 10 por FDR) para $event..."
    rmats2sashimiplot \
        --b1 "$B1" \
        --b2 "$B2" \
        -t "$event" \
        -e "$top10_out" \
        --l1 "$L1" \
        --l2 "$L2" \
        -o "$out_top10_dir"
done

rmats2summary \
  -i /work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/rmats/output_C3_vs_N \
  -o /work/pancreas/takemoto/RNA_seq/results/splicing_alternativo/rmats/output_C3_vs_N
echo "Todos os plots foram gerados com sucesso."
