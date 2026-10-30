#!/usr/bin/env Rscript
# Reference gene-level DEG (tximport + DESeq2)
# DEG: padj < fdr AND |log2FC| >= lfc (log2FC not shrunk).
# Usage: Rscript deseq2_reference.R <samples_file> <quant_dir> <out_dir> <fdr> <lfc> <tx2gene_tsv> [ref_condition]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    stop("Usage: Rscript deseq2_reference.R ",
         "<samples_file> <quant_dir> <out_dir> <fdr> <lfc> <tx2gene_tsv> [ref_condition]")
}
samples_file  <- args[1]
quant_dir     <- args[2]
out_dir       <- args[3]
fdr_thr       <- as.numeric(args[4])
lfc_thr       <- as.numeric(args[5])
tx2gene_file  <- args[6]
ref_condition <- if (length(args) >= 7 && nchar(args[7]) > 0) args[7] else "CTR"

done_file <- file.path(out_dir, "deg_reference_gene_deseq2_ids.txt")
if (file.exists(done_file)) { cat("[deseq2_reference] skip\n"); quit(save = "no", status = 0) }
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
    library(tximport)
    library(DESeq2)
})

tx2gene <- read.table(tx2gene_file, header = TRUE, sep = "\t",
                      col.names = c("tx_id", "gene_id"), stringsAsFactors = FALSE)
tx2gene <- tx2gene[!duplicated(tx2gene$tx_id), ]

meta <- read.table(samples_file, header = FALSE, sep = "\t",
                   col.names = c("condition", "sample", "read1", "read2"),
                   comment.char = "#", stringsAsFactors = FALSE)
files <- setNames(file.path(quant_dir, meta$sample, "quant.sf"), meta$sample)

txi <- tximport(files, type = "salmon",
                tx2gene = tx2gene[, c("tx_id", "gene_id")],
                ignoreTxVersion = FALSE)

cond <- relevel(factor(meta$condition), ref = ref_condition)
col_data <- data.frame(condition = cond, row.names = meta$sample)
dds <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~condition)
dds <- DESeq(dds, quiet = TRUE)

treat_lv <- levels(cond)[2]
res <- as.data.frame(results(dds, contrast = c("condition", treat_lv, ref_condition),
                             alpha = fdr_thr))
res <- res[order(res$pvalue, na.last = TRUE), ]
res$DEG <- !is.na(res$padj) & res$padj < fdr_thr & abs(res$log2FoldChange) >= lfc_thr
res$direction <- ifelse(res$log2FoldChange > 0,
                        paste0("up_", treat_lv), paste0("up_", ref_condition))
res$direction[!res$DEG | is.na(res$padj)] <- "not_sig"
res$gene_id <- rownames(res)

deg_flag <- res$DEG & !is.na(res$DEG)
write.csv(res,             file.path(out_dir, "deg_reference_gene_deseq2_all.csv"))
write.csv(res[deg_flag, ], file.path(out_dir, "deg_reference_gene_deseq2_sig.csv"))
writeLines(rownames(res)[deg_flag], done_file)
cat(sprintf("[deseq2_reference] reference DEGs: %d -> %s\n", sum(deg_flag), out_dir))
