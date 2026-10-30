#!/usr/bin/env Rscript
# DESeq2 on Corset cluster counts
# DEG: padj < 0.05 AND |log2FC| >= 1 (log2FC not shrunk)
# Usage: Rscript deseq2_clusters.R <count_matrix.rds> <coldata.rds> <out_dir> [ref_condition]
suppressPackageStartupMessages(library(DESeq2))

args    <- commandArgs(trailingOnly = TRUE)
counts  <- readRDS(args[1])
coldata <- readRDS(args[2])
out_dir <- args[3]
if (length(args) >= 4) coldata$condition <- relevel(coldata$condition, ref = args[4])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds, parallel = FALSE)

ref <- levels(coldata$condition)[1]
tgt <- levels(coldata$condition)[2]
res_df <- as.data.frame(results(dds, contrast = c("condition", tgt, ref),
                                alpha = 0.05, independentFiltering = TRUE))
res_df$gene_id <- rownames(res_df)
deg_std <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 &
                  abs(res_df$log2FoldChange) >= 1, ]

label <- paste0(tgt, "_vs_", ref)
write.csv(res_df,  file.path(out_dir, sprintf("deseq2_all_%s.csv",     label)))
write.csv(deg_std, file.path(out_dir, sprintf("deseq2_deg_std_%s.csv", label)))
writeLines(deg_std$gene_id, file.path(out_dir, sprintf("deg_std_ids_%s.txt", label)))
cat(sprintf("[deseq2_clusters] %s: %d DEG clusters -> %s\n", label, nrow(deg_std), out_dir))
