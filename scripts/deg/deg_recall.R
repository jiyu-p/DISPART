#!/usr/bin/env Rscript
# Gene-level DEG precision/recall with direction check
# Usage: Rscript deg_recall.R <ref_ids.txt> <deseq2_dir> <corset-clusters.txt> <best_hit.tsv> <out_dir> [label]

args  <- commandArgs(trailingOnly = TRUE)
out   <- args[5]
label <- if (length(args) >= 6) args[6] else "assembly"
done  <- file.path(out, "gene_recall_summary.tsv")
if (file.exists(done)) { cat("[deg_recall] skip\n"); quit(status = 0) }
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Reference DEGs and log2FC (_sig.csv next to the id list)
ref_genes <- setdiff(trimws(readLines(args[1])), "")
ref_sig   <- read.csv(sub("_ids\\.txt$", "_sig.csv", args[1]))
ref_lfc   <- setNames(ref_sig$log2FoldChange, ref_sig$gene_id)

# Assembly DEG clusters and log2FC
deg      <- read.csv(list.files(args[2], "^deseq2_deg_std_.*\\.csv$", full.names = TRUE)[1])
clusters <- as.character(deg$gene_id)
cl_lfc   <- setNames(deg$log2FoldChange, clusters)

members   <- read.table(args[3], sep = "\t", col.names = c("contig", "cluster"))
members   <- split(members$contig, members$cluster)
hit       <- read.table(args[4], header = TRUE, sep = "\t")
hit_gene  <- setNames(hit$ref_gene_id, hit$trinity_id)
hit_score <- setNames(hit$bitscore, hit$trinity_id)

cluster_gene <- vapply(clusters, function(cl) {
    g <- hit_gene[members[[cl]]]
    g <- g[!is.na(g)]
    if (length(g) == 0) return(NA_character_)
    n <- table(g)
    tied <- names(n)[n == max(n)]
    if (length(tied) == 1) return(tied)
    tied[which.max(sapply(tied, function(x) max(hit_score[names(g)[g == x]])))]
}, character(1))

classify <- function(g) {
    if (!g %in% ref_genes) return("FP_nonDEG")
    r <- sign(ref_lfc[g])
    d <- sign(cl_lfc[names(cluster_gene)[which(cluster_gene == g)]])
    d <- d[!is.na(d) & d != 0]
    if (is.na(r) || r == 0 || length(d) == 0) return("unknown")
    if (all(d == r)) "TP" else if (all(d != r)) "FP_direction" else "FP_conflict"
}
genes      <- unique(unname(cluster_gene[!is.na(cluster_gene)]))
gene_class <- vapply(genes, classify, character(1))

n_tp      <- sum(gene_class == "TP")
precision <- if (length(genes) > 0) round(n_tp / length(genes) * 100, 2) else 0
recall    <- if (length(ref_genes) > 0) round(n_tp / length(ref_genes) * 100, 2) else 0
f1        <- if (precision + recall > 0) round(2 * precision * recall / (precision + recall), 4) else 0
n_class   <- function(x) sum(gene_class == x)

write.table(data.frame(
    assembly = label, ref_DEGs = length(ref_genes), denovo_DEGs = length(genes), n_TP = n_tp,
    n_FP_nonDEG = n_class("FP_nonDEG"), n_FP_direction = n_class("FP_direction"),
    n_FP_conflict = n_class("FP_conflict"), n_unknown = n_class("unknown"),
    n_unmapped_clusters = sum(is.na(cluster_gene)),
    precision_pct = precision, recall_pct = recall, f1 = f1
), done, sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(
    cluster_id  = clusters,
    gene_id     = unname(cluster_gene),
    gene_class  = ifelse(is.na(cluster_gene), "unmapped", unname(gene_class[cluster_gene])),
    cluster_lfc = unname(cl_lfc[clusters])
), file.path(out, "cluster_classification.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("[deg_recall] %s: TP=%d  Precision=%.1f%%  Recall=%.1f%%  F1=%.4f -> %s\n",
            label, n_tp, precision, recall, f1, out))