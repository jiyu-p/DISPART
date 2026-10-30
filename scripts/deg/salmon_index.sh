#!/bin/bash
# Salmon index of a transcript FASTA (reference rna.fna or an assembly)
# Usage: bash salmon_index.sh <FASTA> <INDEX_DIR> [THREADS]
set -euo pipefail

FASTA="$1" INDEX_DIR="$2"
THREADS="${3:-${THREADS:-32}}"

[ -f "$INDEX_DIR/versionInfo.json" ] && { echo "[salmon_index] skip"; exit 0; }
mkdir -p "$INDEX_DIR"
salmon index -t "$FASTA" -i "$INDEX_DIR" -p "$THREADS" > "$INDEX_DIR/index.log" 2>&1
echo "[salmon_index] done -> $INDEX_DIR"
