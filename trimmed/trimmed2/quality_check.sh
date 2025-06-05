#!/bin/bash

trimmed_dir="/work/pancreas/takemoto/trimmed/trimmed2"
multiqc_dir="$trimmed_dir/multiqc"

mkdir -p "$multiqc_dir"

echo "Rodando MultiQC no diretório: $trimmed_dir (buscando JSONs do fastp)"
multiqc "$trimmed_dir" -o "$multiqc_dir"

echo "Relatório MultiQC pronto em: $multiqc_dir/multiqc_report.html"
