#!/bin/bash
# CD-HIT-EST redundancy removal
# Usage: bash cdhit.sh <FASTA> <OUT_DIR>
set -euo pipefail

FASTA="$1" OUT="$2"
THREADS="${THREADS:-32}"
CDHIT_ID="${CDHIT_ID:-0.95}"
CDHIT_WORDLEN="${CDHIT_WORDLEN:-10}"

mkdir -p "$OUT"
[ -f "$OUT/cdhit.fasta" ] && { echo "[cdhit] skip"; exit 0; }

cd-hit-est -i "$FASTA" -o "$OUT/cdhit.fasta" \
    -c "$CDHIT_ID" -n "$CDHIT_WORDLEN" -T "$THREADS" -M 32000 -d 0 \
    > "$OUT/cdhit.log" 2>&1
echo "[cdhit] done -> $OUT/cdhit.fasta"
