#!/usr/bin/env Rscript
# Corset cluster count matrix (columns in samples_file order) for DESeq2
# Usage: Rscript count_matrix.R <samples_file> <corset_dir> <out_dir>
args         <- commandArgs(trailingOnly = TRUE)
samples_file <- args[1]
corset_dir   <- args[2]
out_dir      <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.table(samples_file, header = FALSE, sep = "\t",
                   col.names = c("condition", "sample", "read1", "read2"),
                   comment.char = "#", stringsAsFactors = FALSE)
counts <- read.table(file.path(corset_dir, "corset-counts.txt"), header = TRUE, row.names = 1,
                     sep = "\t", check.names = FALSE)
counts <- counts[, meta$sample, drop = FALSE]

coldata <- data.frame(sample = meta$sample, condition = factor(meta$condition),
                      row.names = meta$sample)

saveRDS(counts,  file.path(out_dir, "count_matrix.rds"))
saveRDS(coldata, file.path(out_dir, "coldata.rds"))
write.csv(counts,  file.path(out_dir, "count_matrix.csv"))
write.csv(coldata, file.path(out_dir, "coldata.csv"))
cat(sprintf("[count_matrix] %d clusters x %d samples -> %s\n", nrow(counts), ncol(counts), out_dir))
