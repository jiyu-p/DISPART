#!/bin/bash
# CAP3 overlap-based merging of contigs (contigs + singlets)
# Usage: bash cap3.sh <FASTA> <OUT_DIR>
# -p 90 (overlap identity >= 90%), -o 40 (overlap >= 40 bp)
set -euo pipefail

FASTA="$1" OUT="$2"
mkdir -p "$OUT"
[ -f "$OUT/cap3.fasta" ] && { echo "[cap3] skip"; exit 0; }

cp "$FASTA" "$OUT/input.fasta"
cd "$OUT"
cap3 input.fasta -p 90 -o 40 > cap3.log 2>&1
cat input.fasta.cap.contigs input.fasta.cap.singlets > cap3.fasta
echo "[cap3] done -> $OUT/cap3.fasta"
