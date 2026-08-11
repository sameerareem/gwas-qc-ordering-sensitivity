#!/usr/bin/env Rscript
# =============================================================================
#  compute_qc_metrics.R
#  Post-SAIGE metrics for QC ordering benchmarking
#
#  Usage:
#    Rscript compute_qc_metrics.R \
#      --saige_results  <saige_result.txt> \
#      --causal_snps    <causal_snps.tsv> \
#      --qc_log         <qc_log.txt> \
#      --baseline_snps  <baseline_snps.txt>   (optional, for Jaccard)
#      --baseline_samps <baseline_samps.txt>  (optional, for Jaccard)
#      --baseline_lambda <float>              (optional, for lambda_ratio)
#      --out            <output_prefix>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--saige_results",   type="character"),
  make_option("--causal_snps",     type="character"),
  make_option("--final_bim",       type="character", help=".bim of final QC dataset"),
  make_option("--final_fam",       type="character", help=".fam of final QC dataset"),
  make_option("--baseline_snps",   type="character", default=NULL),
  make_option("--baseline_samps",  type="character", default=NULL),
  make_option("--baseline_lambda", type="double",    default=NULL),
  make_option("--out",             type="character", default="qc_metrics")
)
opt <- parse_args(OptionParser(option_list=option_list))

results <- fread(opt$saige_results)
causal  <- fread(opt$causal_snps)
bim     <- fread(opt$final_bim,
                 col.names=c("CHR","ID","CM","POS","A1","A2"))
fam     <- fread(opt$final_fam,
                 col.names=c("FID","IID","PAT","MAT","SEX","PHENO"))

cat("=== QC Ordering Metrics ===\n\n")

# ── 1. Sample and SNP counts ──────────────────────────────────────────────────
n_samples  <- nrow(fam)
n_snps     <- nrow(bim)
n_tested   <- nrow(results)
cat(sprintf("Retained samples:  %d\n", n_samples))
cat(sprintf("Retained SNPs:     %d\n", n_snps))
cat(sprintf("SNPs tested:       %d\n", n_tested))
cat("\n")

# ── 2. Genomic inflation (lambda GC) ─────────────────────────────────────────
pvals <- results$p.value
pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
chisq <- qchisq(pvals, df=1, lower.tail=FALSE)
lambda_gc <- median(chisq, na.rm=TRUE) / qchisq(0.5, df=1, lower.tail=FALSE)
cat(sprintf("Lambda GC:         %.4f\n", lambda_gc))

# Lambda ratio relative to baseline (if provided)
if (!is.null(opt$baseline_lambda)) {
  lambda_ratio <- lambda_gc / opt$baseline_lambda
  cat(sprintf("Lambda ratio:      %.4f  (vs baseline %.4f)\n",
              lambda_ratio, opt$baseline_lambda))
}
cat("\n")

# ── 3. Causal SNP recovery ────────────────────────────────────────────────────
causal_in_results <- merge(
  causal[, .(SNP, BETA, CHR, POS)],
  results[, .(MarkerID, BETA_est=BETA, SE, p.value, AF_Allele2)],
  by.x="SNP", by.y="MarkerID"
)

n_causal        <- nrow(causal)
n_causal_tested <- nrow(causal_in_results)
n_recovered_gw  <- sum(causal_in_results$p.value < 5e-8,  na.rm=TRUE)
n_recovered_sug <- sum(causal_in_results$p.value < 1e-5,  na.rm=TRUE)
n_recovered_nom <- sum(causal_in_results$p.value < 0.05,  na.rm=TRUE)

cat(sprintf("Causal SNPs total:           %d\n", n_causal))
cat(sprintf("Causal SNPs tested:          %d  (%.1f%%)\n",
            n_causal_tested, 100*n_causal_tested/n_causal))
cat(sprintf("Recovery (p<5e-8, GW sig):   %d  (%.1f%%)\n",
            n_recovered_gw,  100*n_recovered_gw/n_causal))
cat(sprintf("Recovery (p<1e-5, suggestive):%d  (%.1f%%)\n",
            n_recovered_sug, 100*n_recovered_sug/n_causal))
cat(sprintf("Recovery (p<0.05, nominal):  %d  (%.1f%%)\n",
            n_recovered_nom, 100*n_recovered_nom/n_causal))
cat("\n")

# ── 4. Effect size concordance ────────────────────────────────────────────────
if (nrow(causal_in_results) > 2) {
  # Align effect direction — flip if alleles are swapped
  beta_rank_corr <- cor(causal_in_results$BETA,
                        causal_in_results$BETA_est,
                        method="spearman", use="complete.obs")
  beta_pearson   <- cor(causal_in_results$BETA,
                        causal_in_results$BETA_est,
                        method="pearson",  use="complete.obs")
  cat(sprintf("Beta rank correlation (Spearman): %.4f\n", beta_rank_corr))
  cat(sprintf("Beta correlation (Pearson):       %.4f\n", beta_pearson))
  cat("\n")
}

# ── 5. Jaccard overlap with baseline (if provided) ───────────────────────────
if (!is.null(opt$baseline_snps)) {
  baseline_snps  <- fread(opt$baseline_snps,  header=FALSE)$V1
  current_snps   <- bim$ID
  intersection_s <- length(intersect(current_snps, baseline_snps))
  union_s        <- length(union(current_snps, baseline_snps))
  jaccard_snps   <- intersection_s / union_s
  cat(sprintf("Jaccard SNP overlap with baseline:    %.4f\n", jaccard_snps))
  cat(sprintf("  Shared: %d, Only this QC: %d, Only baseline: %d\n",
              intersection_s,
              length(setdiff(current_snps, baseline_snps)),
              length(setdiff(baseline_snps, current_snps))))
}

if (!is.null(opt$baseline_samps)) {
  baseline_samps  <- fread(opt$baseline_samps, header=FALSE)$V1
  current_samps   <- fam$IID
  intersection_sa <- length(intersect(current_samps, baseline_samps))
  union_sa        <- length(union(current_samps, baseline_samps))
  jaccard_samps   <- intersection_sa / union_sa
  cat(sprintf("Jaccard sample overlap with baseline: %.4f\n", jaccard_samps))
  cat(sprintf("  Shared: %d, Only this QC: %d, Only baseline: %d\n",
              intersection_sa,
              length(setdiff(current_samps, baseline_samps)),
              length(setdiff(baseline_samps, current_samps))))
  cat("\n")
}

# ── 6. Write summary to TSV ───────────────────────────────────────────────────
summary_dt <- data.table(
  n_samples             = n_samples,
  n_snps                = n_snps,
  n_snps_tested         = n_tested,
  lambda_gc             = round(lambda_gc, 4),
  lambda_ratio          = if (!is.null(opt$baseline_lambda))
                            round(lambda_gc / opt$baseline_lambda, 4) else NA,
  causal_tested         = n_causal_tested,
  recovery_gw_sig       = n_recovered_gw,
  recovery_gw_pct       = round(100 * n_recovered_gw  / n_causal, 2),
  recovery_sug          = n_recovered_sug,
  recovery_sug_pct      = round(100 * n_recovered_sug / n_causal, 2),
  beta_rank_corr        = if (nrow(causal_in_results) > 2)
                            round(beta_rank_corr, 4) else NA,
  jaccard_snps          = if (!is.null(opt$baseline_snps))
                            round(jaccard_snps, 4) else NA,
  jaccard_samples       = if (!is.null(opt$baseline_samps))
                            round(jaccard_samps, 4) else NA
)

out_file <- paste0(opt$out, "_metrics.tsv")
fwrite(summary_dt, out_file, sep="\t")
cat(sprintf("Metrics written -> %s\n", out_file))
