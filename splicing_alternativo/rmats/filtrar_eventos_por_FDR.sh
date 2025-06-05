#!/bin/bash

# Caminho para os arquivos JCEC
EVENTOS=("SE" "MXE" "A3SS" "A5SS" "RI")

for evento in "${EVENTOS[@]}"; do
arquivo="${evento}.MATS.JCEC.txt"
saida="${evento}_fdr05.txt"

if [[ -f "$arquivo" ]]; then
echo "Filtrando $arquivo..."
awk 'NR==1 || $20 < 0.05' "$arquivo" > "$saida"
else
  echo "Arquivo $arquivo não encontrado."
fi
done
