#### 07 - RNA-seq: differential expression (edgeR) and PCA ####################
# Gene-level Salmon quantifications (tximport) -> edgeR quasi-likelihood models;
# PCA on DESeq2 variance-stabilised counts (Figure 6).
# Original script: Kevin Mulder (Feb 2025); adapted by Andrea Barbi.
#
# Comparisons
#   four_groups   AA vs AP vs CA vs CP (QL F-test over all group coefficients)
#   treatment     disinfected (axenic) vs control
#   water         autoclaved pond water vs pond water
#   left_vs_rest  exploratory: samples on the left of the PCA (AA1, AA3, AA4, AP3, CP1)
#                 vs the remaining samples. The grouping is data-driven, so its p-values
#                 are descriptive only.
#
# For every comparison the QL F-test is reported (logFC, logCPM, F, p, FDR).
# For two-group comparisons, edgeR's glmTreat is also run: it tests against a minimum
# fold change (H0: |logFC| <= lfc, default 1), with FDR < 0.05.
#
# Input : AxenicTadpole_RNA_metadata.csv, tx2gene_17709.txt, quants_17709/<Sample>_qt/quant.sf.gz
# Output: results/DE_<comparison>_QLF.csv, DE_<comparison>_TREAT.csv, DE_summary.csv,
#         SupplementaryTable1_top20_left_vs_rest.csv, PCA_groups.png, sessionInfo.txt
# Run from the project folder (relative paths).
################################################################################

suppressPackageStartupMessages({
  library(tximport)
  library(edgeR)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

## Analysis settings --------------------------------------------------------------
settings <- list(
  gene_filter        = "none",   # "none", "rowSums>=5", "rowSums>=50" or "filterByExpr"
  use_length_offsets = TRUE,     # TRUE: gene-length offsets from tximport; FALSE: plain TMM
  robust_fit         = TRUE,     # robust empirical Bayes in glmQLFit
  lfc                = 1         # log2 fold-change threshold (glmTreat; counts of large changes)
)
left_side <- c("AA1", "AA3", "AA4", "AP3", "CP1")     # left side of the PCA (Figure 6)

meta_file    <- "AxenicTadpole_RNA_metadata.csv"
tx2gene_file <- "tx2gene_17709.txt"
quant_dir    <- "quants_17709"
out_dir      <- "results"
dir.create(out_dir, showWarnings = FALSE)

## Metadata and Salmon import -----------------------------------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta$Group  <- factor(meta$Group,  levels = c("AA", "AP", "CA", "CP"))
meta$Status <- factor(meta$Status, levels = c("Control", "Axenic"))
meta$Water  <- factor(meta$Water,  levels = c("Pond", "Autoclaved"))
# reference level = Left, so logFC > 0 means higher in the remaining samples
meta$Side   <- factor(ifelse(meta$Sample %in% left_side, "Left", "Rest"), levels = c("Left", "Rest"))
stopifnot(sum(meta$Side == "Left") == length(left_side))

files <- file.path(quant_dir, paste0(meta$Sample, "_qt"), "quant.sf.gz")
names(files) <- meta$Sample
stopifnot(all(file.exists(files)))

tx2gene <- read.table(tx2gene_file, header = FALSE)
colnames(tx2gene) <- c("transcript", "gene")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)     # gene level
cts <- txi$counts
len <- txi$length
stopifnot(identical(colnames(cts), meta$Sample))

## Differential expression (edgeR) ------------------------------------------------
run_edger <- function(cts, len, meta, var, s = settings) {
  design <- model.matrix(reformulate(var), data = meta)
  keep <- switch(s$gene_filter,
    "none"         = rep(TRUE, nrow(cts)),
    "rowSums>=5"   = rowSums(cts) >= 5,
    "rowSums>=50"  = rowSums(cts) >= 50,
    "filterByExpr" = filterByExpr(DGEList(cts), design),
    stop("unknown gene_filter"))
  cts <- cts[keep, , drop = FALSE]
  len <- len[keep, , drop = FALSE]

  y <- DGEList(counts = cts)
  if (s$use_length_offsets) {
    normMat <- len / exp(rowMeans(log(len)))
    o <- log(calcNormFactors(cts / normMat)) + log(colSums(cts / normMat))
    y$offset <- t(t(log(normMat)) + o)
  } else {
    y <- calcNormFactors(y)
  }
  y   <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = s$robust_fit)

  out <- list(QLF = topTags(glmQLFTest(fit, coef = 2:ncol(design)), n = Inf)$table,
              TREAT = NULL)
  if (ncol(design) == 2)                      # two-group comparisons only
    out$TREAT <- topTags(glmTreat(fit, coef = 2, lfc = s$lfc), n = Inf)$table
  out
}

comparisons <- c(four_groups = "Group", treatment = "Status", water = "Water", left_vs_rest = "Side")
results <- lapply(comparisons, function(v) run_edger(cts, len, meta, v))

for (a in names(results)) {
  write.csv(results[[a]]$QLF, file.path(out_dir, paste0("DE_", a, "_QLF.csv")))
  if (!is.null(results[[a]]$TREAT))
    write.csv(results[[a]]$TREAT, file.path(out_dir, paste0("DE_", a, "_TREAT.csv")))
}

summary_tab <- do.call(rbind, lapply(names(results), function(a) {
  q <- results[[a]]$QLF
  lfc_cols <- q[, grep("^logFC", names(q)), drop = FALSE]
  data.frame(comparison = a, genes_tested = nrow(q),
             QLF_FDR_below_0.05 = sum(q$FDR < 0.05),
             QLF_FDR_below_0.05_and_abs_logFC_above_lfc =
               sum(q$FDR < 0.05 & apply(abs(lfc_cols), 1, max) > settings$lfc),
             TREAT_FDR_below_0.05 = if (is.null(results[[a]]$TREAT)) NA
                                    else sum(results[[a]]$TREAT$FDR < 0.05))
}))
print(summary_tab)
write.csv(summary_tab, file.path(out_dir, "DE_summary.csv"), row.names = FALSE)

# Supplementary Table 1: top 20 genes (QL F-test), left vs rest
write.csv(head(results$left_vs_rest$QLF[, c("logFC", "logCPM", "F", "PValue", "FDR")], 20),
          file.path(out_dir, "SupplementaryTable1_top20_left_vs_rest.csv"))

## PCA (Figure 6) -----------------------------------------------------------------
# DESeq2 variance-stabilised counts, genes with >= 5 reads in total
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Group)
dds <- dds[rowSums(counts(dds)) >= 5, ]
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
cols <- c(AA = "#1F77B4", AP = "#FFD700", CA = "#D62728", CP = "#7F7F7F")

p <- ggplot(pca, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = name), size = 4, show.legend = FALSE) +
  scale_color_manual(values = cols, name = "Group") +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  coord_fixed() + theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(out_dir, "PCA_groups.png"), p, width = 8, height = 6, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
