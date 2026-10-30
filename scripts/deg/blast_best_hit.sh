#!/bin/bash
# BLASTn de novo transcripts -> reference mRNA, gene-level best hit
# Usage: bash blast_best_hit.sh <DENOVO_FASTA> <REF_MRNA> <OUT_DIR> <TX2GENE_TSV> [EVALUE] [IDENTITY] [THREADS]
set -euo pipefail

DENOVO="$1" REF="$2" OUT="$3" TX2GENE="$4"
EVALUE="${5:-0.001}" IDENTITY="${6:-98}" THREADS="${7:-32}"
ALN_MIN=100

mkdir -p "$OUT"
DONE="$OUT/best_hit.tsv"
[ -f "$DONE" ] && { echo "[blast_best_hit] skip"; exit 0; }

makeblastdb -in "$REF" -dbtype nucl -out "$OUT/ref_db" > /dev/null

blastn \
    -query "$DENOVO" \
    -db "$OUT/ref_db" \
    -out "$OUT/blast_raw.tsv" \
    -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue "$EVALUE" \
    -perc_identity "$IDENTITY" \
    -num_threads "$THREADS" \
    -max_target_seqs 20

python3 - "$OUT/blast_raw.tsv" "$OUT/best_hit.tsv" "$TX2GENE" "$ALN_MIN" << 'PY'
import sys, re, csv
from collections import defaultdict

blast_file   = sys.argv[1]
out_file     = sys.argv[2]
tx2gene_file = sys.argv[3]
aln_min      = int(sys.argv[4])

tx2gene = {}
with open(tx2gene_file) as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        tx_id, gene_id = parts[0], parts[1]
        tx_base = re.sub(r'\.\d+$', '', tx_id)
        tx2gene[tx_base] = gene_id

def strip_version(s):
    return re.sub(r'\.\d+$', '', s)

def union_length(intervals):
    covered, cs, ce = 0, *intervals[0]
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            covered += ce - cs + 1
            cs, ce = s, e
    return covered + ce - cs + 1

hits = defaultdict(lambda: {"qlen": 0, "intervals": [], "bitscore": 0.0, "pident": 0.0})

with open(blast_file) as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) < 13:
            continue
        q, s     = row[0], row[1]
        pident   = float(row[2])
        q_len    = int(row[4])
        qs, qe   = int(row[7]), int(row[8])
        bitscore = float(row[12])
        if qs > qe:
            qs, qe = qe, qs
        d = hits[(q, s)]
        d["qlen"] = q_len
        d["intervals"].append((qs, qe))
        if bitscore > d["bitscore"]:
            d["bitscore"] = bitscore
            d["pident"]   = pident

best_gene = defaultdict(dict)

for (q, s), d in hits.items():
    if d["qlen"] == 0:
        continue
    aln_len = union_length(sorted(d["intervals"]))
    if aln_len < aln_min:
        continue
    ref_tx_base = strip_version(s)
    gene_id = tx2gene.get(ref_tx_base)
    if gene_id is None:
        continue

    cur = best_gene[q].get(gene_id)
    if cur is None or d["bitscore"] > cur[3]:
        best_gene[q][gene_id] = (ref_tx_base, d["pident"], aln_len, d["bitscore"])

best = {}
for q, gene_hits in best_gene.items():
    best_gid = max(gene_hits, key=lambda g: gene_hits[g][3])
    tx_base, pident, aln_len, bitscore = gene_hits[best_gid]
    best[q] = (best_gid, tx_base, pident, aln_len, bitscore)

with open(out_file, "w") as f:
    f.write("trinity_id\tref_gene_id\tref_tx_base\tpident\taln_len\tbitscore\n")
    for q, (gene_id, tx_base, pident, aln_len, bitscore) in sorted(best.items()):
        f.write(f"{q}\t{gene_id}\t{tx_base}\t{pident:.1f}\t{aln_len}\t{bitscore:.1f}\n")
PY
echo "[blast_best_hit] done -> $DONE"
