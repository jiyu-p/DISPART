#!/bin/bash
# Keep sequences >= MIN_LEN bp
# Usage: bash length_filter.sh <FASTA> <OUT_DIR> [MIN_LEN]   (default: $LENGTH_MIN or 300)
set -euo pipefail

FASTA="$1" OUT="$2"
MIN_LEN="${3:-${LENGTH_MIN:-300}}"

mkdir -p "$OUT"
[ -f "$OUT/filtered.fasta" ] && { echo "[length_filter] skip"; exit 0; }

seqkit seq -m "$MIN_LEN" "$FASTA" > "$OUT/filtered.fasta"
echo "[length_filter] done -> $OUT/filtered.fasta"
