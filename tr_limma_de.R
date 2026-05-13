#!/usr/bin/env Rscript
# Limma moderated t-tests on log2(expected_count+1) for lncRNA stage / M_stage
# contrasts. Same sample/gene filters as tr_lncrna_de_analysis.py (except DE
# rule: BH-FDR on limma adj.P.Val).

suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1L) normalizePath(args[[1L]]) else normalizePath(".")
DATA <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "tr_lncrna_output", "limma")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

META_STAGE <- c("sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage")
META_META <- c("sample_id", "cancer_type", "stage", "M_stage")

LOG2_THRESH <- 1.0
EXPR_FRAC_MIN <- 0.4
MIN_SAMPLES_CANCER <- 100L
FDR_MAX <- 0.05
# Minimum samples in each arm (early vs late) for a contrast to be run.
MIN_PER_GROUP <- 10L

STAGE_TRANSITIONS <- list(
  list(name = "I_II", early = "I", late = "II"),
  list(name = "II_III", early = "II", late = "III"),
  list(name = "III_IV", early = "III", late = "IV"),
  list(name = "E_L", early = c("I", "II"), late = c("III", "IV"))
)

M_TRANSITIONS <- list(
  list(name = "M0s_M0l", early = "M0_s", late = "M0_l"),
  list(name = "M0l_M1L", early = "M0_l", late = "M1_L"),
  list(name = "ME_M1L", early = c("M0_s", "M0_l"), late = "M1_L")
)

cancer_types_over_n <- function(df, n = MIN_SAMPLES_CANCER) {
  tb <- table(df$cancer_type)
  names(tb)[tb > n]
}

in_set <- function(v, s) {
  v <- as.character(v)
  if (length(s) == 1L) v == s else v %in% s
}

assign_group <- function(stage_vec, early, late) {
  v <- as.character(stage_vec)
  ie <- in_set(v, early)
  il <- in_set(v, late)
  grp <- rep(NA_character_, length(v))
  grp[ie] <- "early"
  grp[il] <- "late"
  grp
}

gene_pass <- function(Xe, Xl) {
  nge <- if (is.null(Xe) || nrow(Xe) == 0L) 0L else ncol(Xe)
  ngl <- if (is.null(Xl) || nrow(Xl) == 0L) 0L else ncol(Xl)
  if (nge == 0L && ngl == 0L) return(logical(0L))
  stopifnot(nge == 0L || ngl == 0L || nge == ngl)
  ng <- if (nge > 0L) nge else ngl
  Xu <- if (nrow(Xe) == 0L) Xl else if (nrow(Xl) == 0L) Xe else rbind(Xe, Xl)
  if (nrow(Xu) == 0L) return(rep(FALSE, ng))
  colMeans(Xu >= LOG2_THRESH) >= EXPR_FRAC_MIN
}

run_limma_block <- function(df, meta_cols, stage_col, transitions, analysis_label) {
  gene_cols <- setdiff(names(df), meta_cols)
  cancers <- cancer_types_over_n(df[, meta_cols, drop = FALSE])
  rows <- list()

  for (ctype in sort(cancers)) {
    dfc <- df[df$cancer_type == ctype, , drop = FALSE]
    for (tr in transitions) {
      v <- dfc[[stage_col]]
      grp <- assign_group(v, tr$early, tr$late)
      ok <- !is.na(grp)
      if (!any(ok)) next
      sub <- dfc[ok, , drop = FALSE]
      grp <- grp[ok]
      ie <- grp == "early"
      il <- grp == "late"
      if (sum(ie) < 1L || sum(il) < 1L) next
      Xe <- as.matrix(sub[ie, gene_cols, drop = FALSE])
      Xl <- as.matrix(sub[il, gene_cols, drop = FALSE])
      mode(Xe) <- "numeric"
      mode(Xl) <- "numeric"
      keep <- gene_pass(Xe, Xl)
      if (!any(keep)) next
      genes <- gene_cols[keep]
      Y <- t(rbind(sub[ie, genes, drop = FALSE], sub[il, genes, drop = FALSE]))
      mode(Y) <- "numeric"
      grp_f <- factor(c(rep("early", sum(ie)), rep("late", sum(il))), levels = c("early", "late"))
      if (sum(grp_f == "early") < MIN_PER_GROUP || sum(grp_f == "late") < MIN_PER_GROUP) next
      design <- model.matrix(~ grp_f)
      fit <- lmFit(Y, design)
      fit <- eBayes(fit)
      coef_name <- colnames(design)[ncol(design)]
      tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
      tt$gene <- rownames(tt)
      tt$cancer_type <- ctype
      tt$transition <- tr$name
      tt$analysis <- analysis_label
      tt$n_early <- sum(ie)
      tt$n_late <- sum(il)
      rows[[length(rows) + 1L]] <- tt
    }
  }
  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

message("Reading stage matrix ...")
df_s <- read.csv(
  file.path(DATA, "primary_exp_stage_lncRNA.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
message("Reading metastasis matrix ...")
df_m <- read.csv(
  file.path(DATA, "primary_exp_metastasis_lncRNA.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
df_m <- df_m[as.character(df_m$M_stage) != "M1_s", , drop = FALSE]

message("Running limma (stage) ...")
res_s <- run_limma_block(df_s, META_STAGE, "stage", STAGE_TRANSITIONS, "stage")
message("Running limma (M_stage) ...")
res_m <- run_limma_block(df_m, META_META, "M_stage", M_TRANSITIONS, "metastasis")

fmt_detail <- function(res) {
  if (nrow(res) == 0L) {
    return(data.frame())
  }
  data.frame(
    analysis = res$analysis,
    cancer_type = res$cancer_type,
    transition = res$transition,
    gene = res$gene,
    logFC = res$logFC,
    t = res$t,
    P.Value = res$P.Value,
    adj.P.Val = res$adj.P.Val,
    n_early = res$n_early,
    n_late = res$n_late,
    stringsAsFactors = FALSE
  )
}

det_s <- fmt_detail(res_s)
det_m <- fmt_detail(res_m)
write.csv(det_s, file.path(OUT, "limma_stage_all_tests.csv"), row.names = FALSE)
write.csv(det_m, file.path(OUT, "limma_metastasis_all_tests.csv"), row.names = FALSE)

sig_s <- det_s[!is.na(det_s$adj.P.Val) & det_s$adj.P.Val < FDR_MAX, , drop = FALSE]
sig_m <- det_m[!is.na(det_m$adj.P.Val) & det_m$adj.P.Val < FDR_MAX, , drop = FALSE]
write.csv(sig_s, file.path(OUT, "limma_stage_FDR0.05.csv"), row.names = FALSE)
write.csv(sig_m, file.path(OUT, "limma_metastasis_FDR0.05.csv"), row.names = FALSE)

genes_s <- unique(sig_s$gene)
genes_m <- unique(sig_m$gene)
genes_union <- sort(unique(c(genes_s, genes_m)))
writeLines(genes_union, file.path(OUT, "limma_de_genes_union_FDR0.05.txt"))

z_union_path <- file.path(ROOT, "tr_lncrna_output", "tr_genes_union.txt")
z_overlap <- NA_integer_
n_z_union <- NA_integer_
limma_z_genes <- character(0L)
if (file.exists(z_union_path)) {
  zg <- readLines(z_union_path, warn = FALSE)
  zg <- zg[nzchar(zg)]
  n_z_union <- length(zg)
  limma_z_genes <- sort(intersect(zg, genes_union))
  z_overlap <- length(limma_z_genes)
  writeLines(limma_z_genes, file.path(OUT, "limma_z_intersection_genes.txt"))
  if (length(limma_z_genes) > 0L) {
    writeLines(limma_z_genes, file.path(DATA, "canonical_significant_lncRNAs.txt"))
  }
}

summary_list <- list(
  method = "limma_eBayes_on_log2_expected_count_plus_1",
  FDR_max = FDR_MAX,
  min_samples_per_cancer_gt = MIN_SAMPLES_CANCER,
  expr_frac_min_union = EXPR_FRAC_MIN,
  log2_threshold_for_count_ge_1 = LOG2_THRESH,
  min_samples_per_group = MIN_PER_GROUP,
  n_tests_stage = nrow(det_s),
  n_tests_metastasis = nrow(det_m),
  n_sig_FDR_stage = nrow(sig_s),
  n_sig_FDR_metastasis = nrow(sig_m),
  n_unique_genes_sig_stage = length(genes_s),
  n_unique_genes_sig_metastasis = length(genes_m),
  n_unique_genes_sig_union = length(genes_union),
  n_z_union_genes = n_z_union,
  n_overlap_z_union_limma_union = z_overlap,
  n_limma_z_intersection_genes = length(limma_z_genes),
  sig_by_transition_stage = if (nrow(sig_s) > 0L) as.list(table(sig_s$transition)) else list(),
  sig_by_transition_metastasis = if (nrow(sig_m) > 0L) as.list(table(sig_m$transition)) else list()
)

write_json(summary_list, file.path(OUT, "limma_summary.json"), pretty = TRUE, auto_unbox = TRUE)

## Quick diagnostic plot: BH-adjusted p-value histogram (all tests, stage)
if (nrow(det_s) > 0L) {
  png(file.path(OUT, "hist_adjP_stage.png"), width = 640, height = 480)
  hist(det_s$adj.P.Val, breaks = 50, main = "Stage contrasts: limma adj.P.Val (all tests)", xlab = "adj.P.Val")
  dev.off()
}
if (nrow(det_m) > 0L) {
  png(file.path(OUT, "hist_adjP_metastasis.png"), width = 640, height = 480)
  hist(det_m$adj.P.Val, breaks = 50, main = "M_stage contrasts: limma adj.P.Val (all tests)", xlab = "adj.P.Val")
  dev.off()
}

message("Done. Output: ", OUT)
message("Unique DE genes (FDR<", FDR_MAX, ") stage: ", length(genes_s))
message("Unique DE genes (FDR<", FDR_MAX, ") metastasis: ", length(genes_m))
message("Unique DE genes union: ", length(genes_union))
if (!is.na(z_overlap)) {
  message("Overlap with z-score union list: ", z_overlap)
}
