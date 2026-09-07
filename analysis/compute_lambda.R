#!/usr/bin/env Rscript
# =============================================================================
#  compute_lambda.R
#  Computes lambda GC and lambda_1000 from a SAIGE result file.
#
#  Usage:
#    Rscript compute_lambda.R <saige_result.txt> [n_samples]
#
#  Arguments:
#    saige_result.txt  : SAIGE step2 output file
#    n_samples         : post-QC sample size (read from file if not provided)
#
#  Outputs to stdout and writes <prefix>_lambda.tsv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript compute_lambda.R <saige_result.txt> [n_samples]")
}

result_file <- args[1]
n_override  <- if (length(args) >= 2) as.integer(args[2]) else NULL

cat(sprintf("Reading SAIGE results from: %s\n", result_file))
res <- read.table(result_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)

cat(sprintf("Total variants in file: %d\n", nrow(res)))

# ── Get N ──────────────────────────────────────────────────────────────────
# SAIGE reports N per variant in the N column
if (!is.null(n_override)) {
  N <- n_override
  cat(sprintf("N (from argument): %d\n", N))
} else if ("N" %in% names(res)) {
  N <- as.integer(median(res$N, na.rm = TRUE))
  cat(sprintf("N (median from file): %d\n", N))
} else {
  stop("N column not found. Please provide n_samples as second argument.")
}

# ── Clean p-values ─────────────────────────────────────────────────────────
pvals <- res$p.value
n_before <- length(pvals)
pvals <- pvals[!is.na(pvals) & !is.nan(pvals) & pvals > 0 & pvals <= 1]
n_after <- length(pvals)
n_removed <- n_before - n_after

if (n_removed > 0) {
  cat(sprintf("Removed %d invalid p-values (NA/NaN/<=0/missing).\n", n_removed))
}
cat(sprintf("Valid p-values used: %d\n", n_after))

# ── Compute lambda GC ──────────────────────────────────────────────────────
# Lambda GC = median(observed chi-sq) / median(expected chi-sq under null)
# Expected median of chi-sq(1) = qchisq(0.5, df=1) = 0.4549
chisq_obs    <- qchisq(pvals, df = 1, lower.tail = FALSE)
median_obs   <- median(chisq_obs, na.rm = TRUE)
median_exp   <- qchisq(0.5, df = 1, lower.tail = FALSE)   # = 0.4549

lambda_gc    <- median_obs / median_exp

cat(sprintf("\n=== Genomic Inflation Results ===\n"))
cat(sprintf("Median observed chi-sq:  %.4f\n", median_obs))
cat(sprintf("Median expected chi-sq:  %.4f  (theoretical under null)\n", median_exp))
cat(sprintf("Lambda GC:               %.4f\n", lambda_gc))

# ── Interpret lambda ───────────────────────────────────────────────────────
if (lambda_gc > 1.10) {
  cat("Interpretation: INFLATED — possible stratification or cryptic relatedness\n")
} else if (lambda_gc >= 0.95 & lambda_gc <= 1.10) {
  cat("Interpretation: WELL CALIBRATED — acceptable range\n")
} else if (lambda_gc < 0.95 & lambda_gc >= 0.80) {
  cat("Interpretation: MILDLY DEFLATED — check model or data structure\n")
} else {
  cat("Interpretation: SEVERELY DEFLATED — likely structural data issue\n")
}

# ── Compute lambda_1000 ────────────────────────────────────────────────────
# Rescales lambda to what it would be in an equivalent study of N=1000.
# Formula: lambda_1000 = 1 + (lambda_obs - 1) * (1000 / N)
# NOTE: This formula was designed for INFLATION from population stratification.
# For deflation, lambda_1000 will be close to 1 regardless of the true problem,
# so interpret with caution. We report it for completeness / reviewer requests.
lambda_1000 <- 1 + (lambda_gc - 1) * (1000 / N)

cat(sprintf("\nLambda_1000 (rescaled to N=1000): %.4f\n", lambda_1000))
cat(sprintf("N used for rescaling: %d\n", N))

if (lambda_gc < 1) {
  cat("\nCAUTION: lambda_1000 is designed for stratification INFLATION.\n")
  cat("For deflation (lambda < 1), lambda_1000 approaches 1.0 as N increases,\n")
  cat("masking the problem. Your deflation is structural (block reuse),\n")
  cat("not sample-size-dependent. Use lambda_ratio across QC orderings instead.\n")
}

# ── P-value distribution summary ──────────────────────────────────────────
cat(sprintf("\n=== P-value Distribution ===\n"))
cat(sprintf("Min p-value:    %.3e\n", min(pvals)))
cat(sprintf("Median p-value: %.4f  (expected 0.5000 under null)\n",
            median(pvals)))
cat(sprintf("Max p-value:    %.4f\n", max(pvals)))
cat(sprintf("\nSuggested GW significance threshold: 5e-8\n"))
cat(sprintf("Hits at p < 5e-8:  %d\n",  sum(pvals < 5e-8)))
cat(sprintf("Hits at p < 1e-6:  %d\n",  sum(pvals < 1e-6)))
cat(sprintf("Hits at p < 0.05:  %d\n",  sum(pvals < 0.05)))

# ── Write output TSV ──────────────────────────────────────────────────────
prefix   <- sub("\\.[^.]+$", "", result_file)
out_file <- paste0(prefix, "_lambda.tsv")

out <- data.frame(
  result_file      = result_file,
  N                = N,
  n_variants       = n_after,
  lambda_gc        = round(lambda_gc, 4),
  lambda_1000      = round(lambda_1000, 4),
  median_pval      = round(median(pvals), 4),
  min_pval         = formatC(min(pvals), format = "e", digits = 3),
  hits_gw_sig      = sum(pvals < 5e-8),
  hits_suggestive  = sum(pvals < 1e-6)
)

write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nResults written to: %s\n", out_file))
