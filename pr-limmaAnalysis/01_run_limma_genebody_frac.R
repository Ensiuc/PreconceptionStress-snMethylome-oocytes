#!/usr/bin/env Rscript

# ============================================================
# limma-trend differential methylation on FRACTIONS (0–1)
#   - Works with MCDS genebody_da_frac (no count_type)
#   - Runs: Stressed vs Control
#   - Prints BOTH nominal (P.Value) and FDR-significant
#   - Saves full + nominal + FDR tables
#
# INPUTS expected from Python export (see notes below):
#   1) genebody_frac_CGN.csv  (genes x cells; fractions)
#   2) genebody_frac_CHN.csv  (genes x cells; fractions)
#   3) sample_metadata.csv    (columns: sample, Group)
#   4) GeneMetadata_withGeneName_corrected.csv (optional mapping)
# ============================================================

suppressPackageStartupMessages({
  library(limma)
})

# ---------------------------
# USER SETTINGS (EDIT THESE)
# ---------------------------

frac_file_cgn <- "genebody_frac_CGN.csv"   # genes x samples (cells)
frac_file_chn <- "genebody_frac_CHN.csv"   # genes x samples (cells)
meta_file     <- "sample_metadata.csv"     # sample, Group
gene_meta_file <- "GeneMetadata_withGeneName_corrected.csv"  # optional

# Cutoffs
p_cutoff   <- 0.05
fdr_cutoff <- 0.10

# Make group baseline first so coef = Stressed - Control
group_levels <- c("Control", "Stressed")

# Output prefix
out_prefix <- "limma_genebody_frac"

# ---------------------------
# Helpers
# ---------------------------
read_frac <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "double"
  return(mat)
}

attach_gene_names <- function(res, gene_meta_file) {
  if (!file.exists(gene_meta_file)) {
    res$gene_name <- res$gene_id
    return(res)
  }
  gm <- read.csv(gene_meta_file, row.names = 1, stringsAsFactors = FALSE)
  if (!("gene_name" %in% colnames(gm))) {
    res$gene_name <- res$gene_id
    return(res)
  }
  # gm rownames are gene_id
  res$gene_name <- gm[res$gene_id, "gene_name"]
  res$gene_name[is.na(res$gene_name) | res$gene_name == ""] <- res$gene_id[is.na(res$gene_name) | res$gene_name == ""]
  return(res)
}

run_limma_one <- function(frac_mat, meta, label, out_prefix,
                          p_cutoff=0.05, fdr_cutoff=0.10, gene_meta_file=NULL) {

  # Ensure sample order matches
  stopifnot(all(meta$sample %in% colnames(frac_mat)))
  frac_mat <- frac_mat[, meta$sample, drop=FALSE]

  # Optional but recommended for methylation fractions:
  # logit transform helps with 0/1 bounds; add small offset to avoid Inf
 # eps <- 1e-66
  #frac_mat <- pmin(pmax(frac_mat, eps), 1 - eps)
#  M <- log_toggle(frac_mat)  # logit transform
  M <- log2(frac_mat + 1e-6)

  # Design matrix (intercept + GroupStressed)
  design <- model.matrix(~ meta$Group)
  colnames(design) <- make.names(colnames(design))

  cat("\n--------------------------------------------\n")
  cat("Running limma-trend on:", label, "\n")
  cat("--------------------------------------------\n")
  cat("Design columns:\n")
  print(colnames(design))

  # Fit
  fit <- lmFit(M, design)
  fit <- eBayes(fit, trend=TRUE)

  # coef 2 corresponds to GroupStressed if levels are Control, Stressed
  tt <- topTable(fit, coef=2, number=Inf, sort.by="P")
  tt$gene_id <- rownames(tt)

  # Rename for clarity
  # limma columns: logFC, AveExpr, t, P.Value, adj.P.Val, B
  res <- tt[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  names(res)[names(res) == "adj.P.Val"] <- "FDR"

  res$Direction <- ifelse(res$logFC > 0, "higher_in_Stressed", "higher_in_Control")

  # Attach gene names (optional)
  if (!is.null(gene_meta_file)) {
    res <- attach_gene_names(res, gene_meta_file)
  } else {
    res$gene_name <- res$gene_id
  }

  # Reorder columns
  res <- res[, c("gene_id", "gene_name", "logFC", "P.Value", "FDR", "Direction", "AveExpr", "t", "B")]

  # Filter
  sig_nominal <- res[res$P.Value < p_cutoff, ]
  sig_fdr     <- res[res$FDR < fdr_cutoff, ]

  # Print summary
  cat(sprintf("Total genes tested: %d\n", nrow(res)))
  cat(sprintf("Nominally significant (P < %.3f): %d\n", p_cutoff, nrow(sig_nominal)))
  cat(sprintf("FDR significant (FDR < %.3f): %d\n", fdr_cutoff, nrow(sig_fdr)))

  cat("\nTop 10 by nominal P.Value:\n")
  print(head(res[order(res$P.Value), c("gene_id","gene_name","logFC","P.Value","FDR","Direction")], 10))

  cat("\nTop 10 by FDR:\n")
  print(head(res[order(res$FDR), c("gene_id","gene_name","logFC","P.Value","FDR","Direction")], 10))

  # Save outputs
  full_out <- paste0(out_prefix, "_", label, "_FULL.csv")
  nom_out  <- paste0(out_prefix, "_", label, "_nominal_P", p_cutoff, ".csv")
  fdr_out  <- paste0(out_prefix, "_", label, "_FDR_", fdr_cutoff, ".csv")

  write.csv(res, full_out, row.names=FALSE)
  write.csv(sig_nominal, nom_out, row.names=FALSE)
  write.csv(sig_fdr, fdr_out, row.names=FALSE)

  cat("\nSaved:\n")
  cat("  ", full_out, "\n", sep="")
  cat("  ", nom_out,  "\n", sep="")
  cat("  ", fdr_out,  "\n", sep="")

  invisible(list(full=res, nominal=sig_nominal, fdr=sig_fdr))
}

# logit transform helper
log_toggle <- function(x) log(x / (1 - x))

# ---------------------------
# Read metadata
# ---------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
stopifnot(all(c("sample","Group") %in% colnames(meta)))

meta$Group <- factor(meta$Group, levels = group_levels)
if (any(is.na(meta$Group))) stop("Some Group values are not in group_levels.")

cat("Group counts:\n")
print(table(meta$Group))

# ---------------------------
# Read fraction matrices
# ---------------------------
frac_cgn <- read_frac(frac_file_cgn)
frac_chn <- read_frac(frac_file_chn)

# ---------------------------
# Run limma for CGN and CHN
# ---------------------------
res_cgn <- run_limma_one(frac_cgn, meta, "CGN", out_prefix,
                        p_cutoff=p_cutoff, fdr_cutoff=fdr_cutoff,
                        gene_meta_file=gene_meta_file)

res_chn <- run_limma_one(frac_chn, meta, "CHN", out_prefix,
                        p_cutoff=p_cutoff, fdr_cutoff=fdr_cutoff,
                        gene_meta_file=gene_meta_file)

cat("\nDONE.\n")
