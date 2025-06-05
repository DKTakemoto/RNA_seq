#!/bin/bash
# Lista de tipos de eventos do rMATS
event_types=("SE" "RI" "MXE" "A3SS" "A5SS")

# BAMs e rótulos dos grupos
B1=$(paste -sd, /work/pancreas/takemoto/splicing_alternativo/rmats/grupo_C3.txt)
B2=$(paste -sd, /work/pancreas/takemoto/splicing_alternativo/rmats/grupo_N.txt)
L1="C3"
L2="Controle"

# Cria diretório principal de saída
mkdir -p sashimi_plots

for event in "${event_types[@]}"; do
    infile="${event}.MATS.JC.txt"

    # Verifica se o arquivo existe
    if [[ ! -f "$infile" ]]; then
        echo "Arquivo ${infile} não encontrado. Pulando ${event}..."
        continue
    fi

    # Filtra FDR < 0.05
    fdr_out="${event}_fdr05.txt"
    (head -n 1 "$infile" && awk 'NR > 1 && $20 < 0.05' "$infile") > "$fdr_out"

    # Top 10 eventos mais significativos
    top10_out="${event}_top10.txt"
    (head -n 1 "$infile" && tail -n +2 "$infile" | sort -k20,20g | head -n 10) > "$top10_out"

    # Diretórios de saída
    out_fdr_dir="sashimi_plots/${event}_FDR"
    out_top10_dir="sashimi_plots/${event}_TOP10"
    mkdir -p "$out_fdr_dir" "$out_top10_dir"

    # Plots FDR < 0.05
    rmats2sashimiplot \
        --b1 "$B1" \
        --b2 "$B2" \
        -t "$event" \
        -e "$fdr_out" \
        --l1 "$L1" \
        --l2 "$L2" \
        -o "$out_fdr_dir"

    # Plots top 10
    rmats2sashimiplot \
        --b1 "$B1" \
        --b2 "$B2" \
        -t "$event" \
        -e "$top10_out" \
        --l1 "$L1" \
        --l2 "$L2" \
        -o "$out_top10_dir"
done