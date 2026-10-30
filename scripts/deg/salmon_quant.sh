#!/bin/bash
# Salmon quantification of one paired-end sample
# Usage: bash salmon_quant.sh <INDEX> <R1> <R2> <OUT_DIR> [THREADS] [LIB_TYPE]
# --dumpEq writes equivalence classes for Corset.
set -euo pipefail

INDEX="$1" R1="$2" R2="$3" OUT="$4" THREADS="${5:-32}" LIB_TYPE="${6:-A}"

[ -f "$OUT/quant.sf" ] && { echo "[salmon_quant] skip: $(basename "$OUT")"; exit 0; }
mkdir -p "$OUT"
salmon quant -i "$INDEX" -l "$LIB_TYPE" -1 "$R1" -2 "$R2" -p "$THREADS" \
    --validateMappings --gcBias --dumpEq -o "$OUT" \
    > "$OUT/quant.log" 2>&1
echo "[salmon_quant] done -> $OUT/quant.sf"
