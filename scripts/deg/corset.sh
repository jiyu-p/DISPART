#!/bin/bash
# Corset clustering of contigs into gene-level clusters from Salmon eq. classes
# Usage: bash corset.sh <QUANT_BASE_DIR> <SAMPLES_FILE> <OUT_DIR>
set -euo pipefail

QUANT_BASE="$1" SAMPLES_FILE="$2" OUT="$3"

mkdir -p "$OUT"
[ -f "$OUT/corset-clusters.txt" ] && { echo "[corset] skip"; exit 0; }

GROUPS_ARG="" NAMES_ARG="" EQ_FILES=()
while IFS=$'\t' read -r condition sample r1 r2; do
    [[ "$condition" =~ ^# ]] && continue
    eq="$QUANT_BASE/$sample/aux_info/eq_classes.txt"
    # Salmon >= 1.10 writes eq_classes.txt.gz
    [ ! -f "$eq" ] && [ -f "$eq.gz" ] && gunzip -k "$eq.gz"
    GROUPS_ARG+="${condition},"
    NAMES_ARG+="${sample},"
    EQ_FILES+=("$eq")
done < "$SAMPLES_FILE"

cd "$OUT"
corset -i salmon_eq_classes -g "${GROUPS_ARG%,}" -n "${NAMES_ARG%,}" -f TRUE -p corset "${EQ_FILES[@]}"
echo "[corset] done -> $OUT"
