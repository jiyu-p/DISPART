#!/bin/bash
# Basic assembly statistics (count, length, Nx, GC)
# Usage: bash assembly_stats.sh <FASTA> <OUT_DIR>
set -euo pipefail

FASTA="$1" OUT="$2"
mkdir -p "$OUT"
[ -f "$OUT/assembly_stats.txt" ] && { echo "[assembly_stats] skip"; exit 0; }

python3 - "$FASTA" "$OUT" << 'PY'
import sys, json
from pathlib import Path
import numpy as np

fasta, out = sys.argv[1], Path(sys.argv[2])

lengths, gc_sum, seq = [], 0, []
def flush():
    global gc_sum
    if seq:
        s = "".join(seq).upper()
        lengths.append(len(s))
        gc_sum += s.count("G") + s.count("C")
with open(fasta) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            flush(); seq = []
        else:
            seq.append(line)
flush()

lengths.sort(reverse=True)
total = sum(lengths)

def nx(pct):
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total * pct / 100:
            return l

stats = {
    "n_transcripts": len(lengths),
    "total_bases":   total,
    "mean_length":   round(float(np.mean(lengths)), 1),
    "median_length": int(np.median(lengths)),
    "max_length":    lengths[0],
    "min_length":    lengths[-1],
    "n25": nx(25), "n50": nx(50),
    "n75": nx(75), "n90": nx(90),
    "gc_pct": round(gc_sum / total * 100, 2),
}
json.dump(stats, open(out / "assembly_stats.json", "w"), indent=2)
(out / "assembly_stats.txt").write_text(
    f"Assembly stats: {Path(fasta).name}\n" + "".join(f"  {k}: {v}\n" for k, v in stats.items()))
PY
echo "[assembly_stats] done -> $OUT"
