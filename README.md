# DISPART

Classifier-partitioned *de novo* transcriptome assembly.
RNA-seq reads are classified into 45 transcript clusters plus an unassigned class with a
DNABERT-2-based classifier, each cluster is assembled with Trinity, and the result is compared
with a full *de novo* assembly by DEG recovery. This repository contains the pipeline and a
small *C. elegans* test dataset.

## Installation

Linux with conda (or mamba/micromamba).

```bash
conda env create -f envs/infer.yml       # tokenize, classify, split
conda env create -f envs/assembly.yml    # Trinity
conda env create -f envs/eval.yml        # post-processing, DEG
conda env create -f envs/corset.yml      # Corset
```

With micromamba, use `micromamba create -f ...` and set `CONDA_RUN="micromamba run"` in `config.sh`.

GPU is optional. The paper run used flash-attn on GPU:

```bash
conda run -n dispart-infer pip install flash-attn==2.8.3.post1 --no-build-isolation
```

Without it (e.g. on CPU), a few reads may be assigned to a different class.

## Classifier weights

Download `dispart_classifier_k45_v1.tar.gz` from Zenodo
([10.5281/zenodo.23719397](https://doi.org/10.5281/zenodo.23719397)) and extract it into `weights/`:

```bash
tar xzf dispart_classifier_k45_v1.tar.gz -C weights/
```

This adds `weights/classifier/model.safetensors`; the other files in that directory are
already in this repository.

## Running the mini test

Check `config.sh` (threads, memory, conda command), then:

```bash
bash run_minitest.sh            # or: bash run_minitest.sh my_config.sh
```

Steps: tokenize, classify, split, assemble (partitioned and full *de novo*), post-process
(CD-HIT-EST, CAP3, length >= 300 bp), DEG (Salmon, Corset, DESeq2), summary.
Results are written to `out/minitest/`; `out/minitest/summary.tsv` compares the two assemblies.

At this scale the partitioned assembly is slower than the full assembly because each Trinity run has a fixed
start-up cost.

## Expected output

`expected/summary.tsv` (GPU with flash-attn, `THREADS=2`, `PARTITION_JOBS=2`):

```
assembly      n_contigs  n50   mean_length  ref_DEGs  deg_genes  TP  precision_pct  recall_pct  f1
partitioned   219        1407  1092.8       30        26         26  100            86.67       92.8591
full          257        1385  1046.9       30        25         24  96             80          87.2727
```

Small differences are possible with other hardware or thread settings.

## Data

RNA-seq reads: NCBI BioProject [PRJNA1006044](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1006044)
/ GEO [GSE241011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241011)
(*C. elegans*, control vs. heat shock, 3 replicates each).

Reference genome: *C. elegans* WBcel235 (GCF_000002985.6, NCBI RefSeq).
GMAP database name used in the paper: `celegans_PRJNA13758`.

The mini dataset was sampled from these data: 150 genes were selected
(30 reference DEGs and 120 non-DEGs covering all 46 classes; `data/minitest/genes_selected.tsv`),
read pairs assigned by Salmon only to these genes were kept together with 5% Salmon-unmapped pairs,
and each sample was downsampled to at most 70,000 pairs.

## Third-party code

The model code in `weights/classifier/` is from
[DNABERT-2](https://huggingface.co/zhihan1996/DNABERT-2-117M)
(revision `7bce263b15377fc15361f52cfab88f8b586abda0`, Apache License 2.0, see
`weights/classifier/LICENSE_DNABERT2`). `bert_padding.py` and `configuration_bert.py` are
unmodified; `bert_layers.py` and `flash_attn_triton.py` are modified (noted at the top of each file).

Zhou Z, et al. DNABERT-2: Efficient Foundation Model and Benchmark for Multi-Species Genome.
ICLR 2024. arXiv:2306.15006.
