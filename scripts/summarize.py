#!/usr/bin/env python3
"""One table comparing assemblies: final contig statistics and DEG recovery.
Usage: python summarize.py <OUT_TSV> <label>=<assembly_dir> [<label>=<assembly_dir> ...]
Each assembly_dir holds stats/assembly_stats.json and deg/recall/gene_recall_summary.tsv.
"""
import csv
import json
import sys
from pathlib import Path

COLS = ["assembly", "n_contigs", "n50", "mean_length", "ref_DEGs", "deg_genes",
        "TP", "precision_pct", "recall_pct", "f1"]

rows = []
for arg in sys.argv[2:]:
    label, d = arg.split("=", 1)
    d = Path(d)
    stats = json.loads((d / "stats/assembly_stats.json").read_text())
    with open(d / "deg/recall/gene_recall_summary.tsv") as f:
        rec = next(csv.DictReader(f, delimiter="\t"))
    rows.append([label, stats["n_transcripts"], stats["n50"], stats["mean_length"], rec["ref_DEGs"],
                 rec["denovo_DEGs"], rec["n_TP"], rec["precision_pct"], rec["recall_pct"], rec["f1"]])

with open(sys.argv[1], "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(COLS)
    w.writerows(rows)
print("\t".join(COLS))
for r in rows:
    print("\t".join(map(str, r)))
