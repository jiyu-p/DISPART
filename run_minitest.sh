#!/bin/bash
# Mini test, end to end: tokenize -> classify -> split -> assemble (partitioned, full)
# -> post-process (cd-hit-est, CAP3, length filter) -> reference and assembly DEG -> summary.
# Usage: bash run_minitest.sh [CONFIG]   (default: config.sh)
set -euo pipefail
cd "$(dirname "$0")"
source "${1:-config.sh}"

for v in DATA READS SAMPLES REF_RNA REF_GTF WEIGHTS OUT; do
    printf -v "$v" '%s' "$(realpath -m "${!v}")"
done
S=$PWD/scripts
read -ra RUN <<< "${CONDA_RUN:-conda run --no-capture-output}"
inf() { "${RUN[@]}" -n "$ENV_INFER" "$@"; }
asm() { "${RUN[@]}" -n "$ENV_ASSEMBLY" "$@"; }
ev()  { "${RUN[@]}" -n "$ENV_EVAL" "$@"; }
crs() { "${RUN[@]}" -n "$ENV_CORSET" "$@"; }
export THREADS
mkdir -p "$OUT"

# Sample sheet with read paths: condition, sample, read1, read2
SHEET=$OUT/samples_reads.tsv
awk -v d="$READS" 'BEGIN{OFS="\t"} !/^#/ && NF >= 2 {print $1, $2, d"/"$2"_1.fastq", d"/"$2"_2.fastq"}' \
    "$SAMPLES" > "$SHEET"
R1=$READS/${ASSEMBLY_SAMPLE}_1.fastq
R2=$READS/${ASSEMBLY_SAMPLE}_2.fastq

quant_all() {   # args: salmon_index out_dir
    while IFS=$'\t' read -r _ sample r1 r2 <&3; do
        ev bash "$S/deg/salmon_quant.sh" "$1" "$r1" "$r2" "$2/$sample" "$THREADS" A
    done 3< "$SHEET"
}

echo "== classify and split"
inf python "$S/01_tokenize.py" --r1 "$R1" --out "$OUT/tokens" --tokenizer "$WEIGHTS" --workers "$THREADS"
inf python "$S/02_classify.py" --tokenized "$OUT/tokens" --model "$WEIGHTS" --out "$OUT/classify" \
    ${CLASSIFY_DEVICE:+--device "$CLASSIFY_DEVICE"}
inf python "$S/03_split_reads.py" --r1 "$R1" --r2 "$R2" \
    --predictions "$OUT/classify/predictions.npz" --out "$OUT/split"

echo "== assemble"
export MAX_MEM=$TRINITY_MAX_MEM MAX_JOBS=$PARTITION_JOBS
asm bash "$S/04a_assemble_partitioned.sh" "$OUT/split/cluster_manifest.csv" "$OUT/partitioned/assembly"
# Full de novo gets the same total cores as the partitioned run (THREADS x PARTITION_JOBS).
THREADS=$((THREADS * PARTITION_JOBS)) asm bash "$S/04b_assemble_full.sh" "$R1" "$R2" "$OUT/full/assembly"

echo "== reference DEG"
REFD=$OUT/reference
mkdir -p "$REFD"
ev python "$S/deg/make_tx2gene.py" "$REF_GTF" "$REFD/tx2gene.tsv"
ev bash "$S/deg/salmon_index.sh" "$REF_RNA" "$REFD/salmon_index" "$THREADS"
quant_all "$REFD/salmon_index" "$REFD/quant"
ev Rscript "$S/deg/deseq2_reference.R" "$SHEET" "$REFD/quant" "$REFD/deg" 0.05 1 \
    "$REFD/tx2gene.tsv" "$REF_CONDITION"

declare -A RAW=([partitioned]=$OUT/partitioned/assembly/merged_assembly.fasta
                [full]=$OUT/full/assembly/Trinity.fasta)
for A in partitioned full; do
    echo "== post-process and DEG: $A"
    D=$OUT/$A
    F=$D/final/filtered.fasta
    ev bash "$S/postprocess/cdhit.sh" "${RAW[$A]}" "$D/cdhit"
    ev bash "$S/postprocess/cap3.sh" "$D/cdhit/cdhit.fasta" "$D/cap3"
    ev bash "$S/postprocess/length_filter.sh" "$D/cap3/cap3.fasta" "$D/final" 300
    ev bash "$S/postprocess/assembly_stats.sh" "$F" "$D/stats"

    ev bash "$S/deg/salmon_index.sh" "$F" "$D/deg/salmon_index" "$THREADS"
    quant_all "$D/deg/salmon_index" "$D/deg/quant"
    ev bash "$S/deg/blast_best_hit.sh" "$F" "$REF_RNA" "$D/deg/blast" "$REFD/tx2gene.tsv" 0.001 98 "$THREADS"
    crs bash "$S/deg/corset.sh" "$D/deg/quant" "$SHEET" "$D/deg/corset"
    ev Rscript "$S/deg/count_matrix.R" "$SHEET" "$D/deg/corset" "$D/deg/counts"
    ev Rscript "$S/deg/deseq2_clusters.R" "$D/deg/counts/count_matrix.rds" "$D/deg/counts/coldata.rds" \
        "$D/deg/deseq2" "$REF_CONDITION"
    ev Rscript "$S/deg/deg_recall.R" "$REFD/deg/deg_reference_gene_deseq2_ids.txt" "$D/deg/deseq2" \
        "$D/deg/corset/corset-clusters.txt" "$D/deg/blast/best_hit.tsv" "$D/deg/recall" "$A"
done

echo "== summary"
ev python "$S/summarize.py" "$OUT/summary.tsv" partitioned="$OUT/partitioned" full="$OUT/full"
