#!/bin/bash
# Full de novo assembly of one paired-end sample (Trinity v2.15.2, unstranded).
# Paper run: THREADS=52, MAX_MEM=200G.
# Usage: [THREADS=8] [MAX_MEM=16G] [TMP_BASE=/tmp] bash 04b_assemble_full.sh R1.fastq R2.fastq OUT_DIR
set -euo pipefail

R1=$1
R2=$2
OUT_DIR=$3
THREADS=${THREADS:-8}
MAX_MEM=${MAX_MEM:-16G}
TMP_BASE=${TMP_BASE:-/tmp}

FASTA_OUT="$OUT_DIR/Trinity.fasta"
[ -f "$FASTA_OUT" ] && { echo "exists, skipping: $FASTA_OUT"; exit 0; }
mkdir -p "$OUT_DIR"

T=$(mktemp -d "$TMP_BASE/trinity_full_XXXX")
trap 'rm -rf "$T"' EXIT

Trinity --seqType fq --left "$R1" --right "$R2" \
        --CPU "$THREADS" --max_memory "$MAX_MEM" \
        --output "$T/trinity_out" > "$OUT_DIR/trinity.log" 2>&1

mv "$T/trinity_out.Trinity.fasta" "$FASTA_OUT"
echo "done: $FASTA_OUT"
