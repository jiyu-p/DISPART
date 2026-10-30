#!/usr/bin/env python3
"""Classify tokenized R1 reads into 45 transcript clusters + NC (class 45).
Output: <out>/predictions.npz with preds (top-1 class, int32), in the tokenized read order.

Original (paper) run: 4 GPUs, batch size 16384, fp16 weights + CUDA autocast.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModelForSequenceClassification


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tokenized", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--batch-size", type=int, default=1024)
    a = ap.parse_args()

    meta = json.loads((a.tokenized / "tokenize_meta.json").read_text())
    n, L = meta["n_reads"], meta["max_len"]
    ids  = np.memmap(a.tokenized / "input_ids.mmap", dtype=np.int32, mode="r", shape=(n, L))
    mask = np.memmap(a.tokenized / "attention_mask.mmap", dtype=np.uint8, mode="r", shape=(n, L))

    device = torch.device(a.device)
    model = AutoModelForSequenceClassification.from_pretrained(a.model, trust_remote_code=True)
    if device.type == "cuda":
        model = model.half()
    else:
        # CPU: disable flash-attn/Triton (CUDA only)
        for mod in model.modules():
            if hasattr(mod, "_flash_attn_varlen"):
                mod._flash_attn_varlen = None
                sys.modules[type(mod).__module__].flash_attn_qkvpacked_func = None
    model = model.to(device).eval()
    assert model.config.num_labels == 46, f"expected 46 classes (45 clusters + NC), got {model.config.num_labels}"
    ctx = torch.autocast("cuda", enabled=(device.type == "cuda"))

    preds = np.empty(n, np.int32)
    with torch.no_grad(), ctx:
        for s in range(0, n, a.batch_size):
            e = min(s + a.batch_size, n)
            x = torch.from_numpy(ids[s:e].astype(np.int64)).to(device)
            m = torch.from_numpy(mask[s:e].astype(np.int64)).to(device)
            logits = model(input_ids=x, attention_mask=m).logits
            preds[s:e] = logits.float().argmax(-1).cpu().numpy()

    a.out.mkdir(parents=True, exist_ok=True)
    np.savez(a.out / "predictions.npz", preds=preds)
    print(f"{n:,} reads -> {a.out / 'predictions.npz'}")


if __name__ == "__main__":
    main()
