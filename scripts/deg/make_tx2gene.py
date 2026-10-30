#!/usr/bin/env python3
"""Transcript-to-gene table from a GTF (transcript lines): tx_id<TAB>gene_id, with header.
Usage: python make_tx2gene.py <GTF> <OUT_TSV>
"""
import re
import sys

gtf, out = sys.argv[1], sys.argv[2]
seen = set()
with open(gtf) as f, open(out, "w") as o:
    o.write("tx_id\tgene_id\n")
    for line in f:
        if "\ttranscript\t" not in line:
            continue
        tx = re.search(r'transcript_id "([^"]+)"', line)
        gn = re.search(r'gene_id "([^"]+)"', line)
        if tx and gn and tx.group(1) not in seen:
            seen.add(tx.group(1))
            o.write(f"{tx.group(1)}\t{gn.group(1)}\n")
