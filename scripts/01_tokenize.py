#!/usr/bin/env python3
"""Tokenize R1 reads with the DNABERT-2 tokenizer for 02_classify.py.

Output (<out>): input_ids.mmap (int32), attention_mask.mmap (uint8), both n_reads x max_len
in FASTQ order, and tokenize_meta.json (n_reads, max_len).
max_len 50 tokens (150-bp read ~33 tokens).
"""
import argparse
import gzip
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from transformers import AutoTokenizer

_TOK = None
_MAX_LEN = None


def open_fastq(path):
    with open(path, "rb") as f:
        is_gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gz else open(path)


def iter_seqs(path):
    with open_fastq(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.strip()


def iter_chunks(path, size):
    buf = []
    for seq in iter_seqs(path):
        buf.append(seq)
        if len(buf) == size:
            yield buf
            buf = []
    if buf:
        yield buf


def init_tokenizer(tokenizer, max_len):
    global _TOK, _MAX_LEN
    _TOK = AutoTokenizer.from_pretrained(tokenizer, trust_remote_code=True, use_fast=True)
    if _TOK.pad_token_id is None:
        _TOK.pad_token = _TOK.eos_token
    _MAX_LEN = max_len


def tokenize(seqs):
    enc = _TOK(seqs, add_special_tokens=True, truncation=True, max_length=_MAX_LEN,
               padding="max_length", return_attention_mask=True)
    return np.array(enc["input_ids"], dtype=np.int32), np.array(enc["attention_mask"], dtype=np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tokenizer", default="zhihan1996/DNABERT-2-117M",
                    help="Hugging Face id or local directory with the tokenizer files")
    ap.add_argument("--max-len", type=int, default=50)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunk-size", type=int, default=100000)
    a = ap.parse_args()

    n = sum(1 for _ in iter_seqs(a.r1))
    a.out.mkdir(parents=True, exist_ok=True)
    ids = np.memmap(a.out / "input_ids.mmap", dtype=np.int32, mode="w+", shape=(n, a.max_len))
    mask = np.memmap(a.out / "attention_mask.mmap", dtype=np.uint8, mode="w+", shape=(n, a.max_len))

    pos = 0
    with Pool(a.workers, initializer=init_tokenizer, initargs=(a.tokenizer, a.max_len)) as pool:
        for chunk_ids, chunk_mask in pool.imap(tokenize, iter_chunks(a.r1, a.chunk_size)):
            k = len(chunk_ids)
            ids[pos:pos + k] = chunk_ids
            mask[pos:pos + k] = chunk_mask
            pos += k
    ids.flush()
    mask.flush()

    (a.out / "tokenize_meta.json").write_text(json.dumps({"n_reads": n, "max_len": a.max_len}, indent=2))
    print(f"{n:,} reads -> {a.out}")


if __name__ == "__main__":
    main()
