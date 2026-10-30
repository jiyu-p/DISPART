#!/bin/bash
# Model (classifier-partitioned) assembly: one Trinity run per cluster, then merge all clusters.
# Input: cluster_manifest.csv from 03_split_reads.py (cluster,n_reads_top1,n_reads_total,fq1,fq2).
# Trinity 2.15.2, unstranded, default normalization.
# Paper run: THREADS=8 per cluster, MAX_JOBS=7, MAX_MEM=100G.
#
# Usage: [THREADS=4] [MAX_JOBS=2] [MAX_MEM=8G] [TMP_BASE=/tmp] bash 04a_assemble_partitioned.sh MANIFEST OUT_DIR
set -euo pipefail

MANIFEST=$1
OUT_DIR=$2
MAX_JOBS=${MAX_JOBS:-2}
export THREADS=${THREADS:-4} MAX_MEM=${MAX_MEM:-8G} TMP_BASE=${TMP_BASE:-/tmp} OUT_DIR

assemble() {   # args: cluster fq1 fq2
    set -euo pipefail
    local pad out
    pad=$(printf "%04d" "$1")
    out="$OUT_DIR/clusters/trinity_cluster_$pad"
    [ -f "$out/Trinity.fasta" ] && return 0
    mkdir -p "$out"
    T=$(mktemp -d "$TMP_BASE/trinity_${pad}_XXXX")
    trap 'rm -rf "$T"' EXIT
    Trinity --seqType fq --left "$2" --right "$3" \
            --CPU "$THREADS" --max_memory "$MAX_MEM" \
            --output "$T/trinity_out" > "$out/trinity.log" 2>&1 \
        || { echo "cluster_$pad failed, see $out/trinity.log" >&2; exit 1; }
    mv "$T/trinity_out.Trinity.fasta" "$out/Trinity.fasta"
    echo "cluster_$pad done"
}
export -f assemble

# Largest clusters first (n_reads_total, descending) to shorten total wall time.
tail -n +2 "$MANIFEST" | sort -t, -k3,3nr | cut -d, -f1,4,5 | tr ',' '\n' \
    | xargs -d '\n' -n 3 -P "$MAX_JOBS" bash -c 'assemble "$@"' _

# Contig names get a cluster prefix: >cluster_0000_TRINITY_DN...
for f in "$OUT_DIR"/clusters/trinity_cluster_*/Trinity.fasta; do
    cid=$(basename "$(dirname "$f")"); cid=${cid#trinity_cluster_}
    awk -v cid="$cid" '/^>/{print ">cluster_"cid"_"substr($0,2); next}{print}' "$f"
done > "$OUT_DIR/merged_assembly.fasta"
echo "done: $OUT_DIR/merged_assembly.fasta ($(grep -c '^>' "$OUT_DIR/merged_assembly.fasta") contigs)"
