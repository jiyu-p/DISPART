# Settings for run_minitest.sh. Relative paths are resolved from the repository root.

DATA=data/minitest
READS=$DATA/reads              # <sample>_1.fastq and <sample>_2.fastq
SAMPLES=$DATA/samples.tsv      # condition<TAB>sample (lines starting with # are ignored)
ASSEMBLY_SAMPLE=SRR25658444    # sample that is classified and assembled
REF_CONDITION=CTR
REF_RNA=$DATA/ref/rna.fna      # reference transcripts (DEG reference and BLAST target)
REF_GTF=$DATA/ref/genomic.gtf  # transcript-to-gene table
WEIGHTS=weights/classifier     # classifier directory (weights, config, tokenizer)
OUT=out/minitest

THREADS=2                      # per Trinity run; Full de novo uses THREADS x PARTITION_JOBS
TRINITY_MAX_MEM=8G
PARTITION_JOBS=8               # Trinity runs in parallel for the partitioned assembly
CLASSIFY_DEVICE=               # cuda or cpu; empty = cuda if available

CONDA_RUN="conda run --no-capture-output"
ENV_INFER=dispart-infer
ENV_ASSEMBLY=dispart-assembly
ENV_EVAL=dispart-eval
ENV_CORSET=dispart-corset
