#!/usr/bin/env python3
"""Split read pairs by predicted class (top-1). Input order must match; stops on length mismatch.
Paper run: --min-reads 100 (all 46 classes kept).
Output: trinity_cluster_XXXX/reads_{1,2}.fq.gz, cluster_manifest.csv
"""
import argparse
import gzip
from itertools import zip_longest
from pathlib import Path

import numpy as np
import pandas as pd


def open_fastq(path):
    with open(path, "rb") as f:
        is_gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if is_gz else open(path, "rb")


def records(f):
    while True:
        rec = f.readline()
        if not rec:
            return
        yield rec + f.readline() + f.readline() + f.readline()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, required=True)
    ap.add_argument("--predictions", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-reads", type=int, default=1)
    a = ap.parse_args()

    preds = np.load(a.predictions)["preds"]
    counts = np.bincount(preds)
    keep = [c for c in range(len(counts)) if counts[c] >= a.min_reads]

    a.out.mkdir(parents=True, exist_ok=True)
    paths = {c: a.out / f"trinity_cluster_{c:04d}" for c in keep}
    handles = {}
    for c, d in paths.items():
        d.mkdir(exist_ok=True)
        handles[c] = (gzip.open(d / "reads_1.fq.gz", "wb", compresslevel=4),
                      gzip.open(d / "reads_2.fq.gz", "wb", compresslevel=4))

    n = 0
    with open_fastq(a.r1) as f1, open_fastq(a.r2) as f2:
        for c, rec1, rec2 in zip_longest(preds, records(f1), records(f2)):
            if c is None or rec1 is None or rec2 is None:
                raise ValueError(f"R1, R2 and predictions differ in length (stopped at read {n + 1:,}; "
                                 f"{len(preds):,} predictions)")
            if c in handles:
                handles[c][0].write(rec1)
                handles[c][1].write(rec2)
            n += 1
    for h1, h2 in handles.values():
        h1.close()
        h2.close()

    pd.DataFrame([{"cluster": c, "n_reads_top1": int(counts[c]), "n_reads_total": int(counts[c]),
                   "fq1": str(paths[c] / "reads_1.fq.gz"), "fq2": str(paths[c] / "reads_2.fq.gz")}
                  for c in keep]).to_csv(a.out / "cluster_manifest.csv", index=False)
    print(f"{n:,} read pairs -> {len(keep)} clusters in {a.out}")


if __name__ == "__main__":
    main()
