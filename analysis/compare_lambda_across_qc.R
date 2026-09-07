#!/usr/bin/env Rscript
# =============================================================================
#  compare_lambda_across_qc.R
#  Compares lambda across QC orderings using three valid approaches:
#    1. Lambda_ratio (relative to baseline ordering)
#    2. Lambda_1000 (N-normalised to 1000 samples)
#    3. Summary table for dissertation reporting
#
#  Usage:
#    Rscript compare_lambda_across_qc.R \
#      --results_dir  <directory containing *_saige_results.txt files> \
#      --baseline_qc  <name of baseline QC ordering, e.g. "qc1"> \
#      --out          <output_prefix>
#
#  Expects files named: genotypes_qc1_saige_results.txt,
#                       genotypes_qc2_saige_results.txt, etc.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--results-dir",  type="character", default=".",   dest="results_dir"),
  make_option("--baseline-qc",  type="character", default="qc1", dest="baseline_qc"),
  make_option("--pattern",      type="character",
              default="*_saige_results.txt",
              help="glob pattern for result files"),
  make_option("--out",          type="character", default="lambda_comparison")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ── Helper: compute lambda from a results file ─────────────────────────────
compute_lambda <- function(filepath) {
  res   <- fread(filepath, select=c("p.value", "N"))
  pvals <- res$p.value
  N     <- as.integer(median(res$N, na.rm=TRUE))

  # Clean p-values
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]

  if (length(pvals) < 100) {
    warning(sprintf("Too few valid p-values in %s (%d)", filepath, length(pvals)))
    return(NULL)
  }

  # Lambda GC
  chisq_obs  <- qchisq(pvals, df=1, lower.tail=FALSE)
  median_exp <- qchisq(0.5,   df=1, lower.tail=FALSE)  # 0.4549
  lambda_gc  <- median(chisq_obs, na.rm=TRUE) / median_exp

  # Lambda_1000: rescale to N=1000
  lambda_1000 <- 1 + (lambda_gc - 1) * (1000 / N)

  list(
    lambda_gc   = lambda_gc,
    lambda_1000 = lambda_1000,
    N           = N,
    n_variants  = length(pvals),
    median_pval = median(pvals),
    n_gw_hits   = sum(pvals < 5e-8)
  )
}

# ── Find all result files ──────────────────────────────────────────────────
files <- list.files(
  opt$results_dir,
  pattern  = gsub("\\*", ".*", opt$pattern),
  full.names = TRUE
)

if (length(files) == 0) {
  stop(sprintf("No files matching pattern '%s' in '%s'",
               opt$pattern, opt$results_dir))
}

cat(sprintf("Found %d result files\n", length(files)))

# ── Compute lambda for each file ───────────────────────────────────────────
results <- lapply(files, function(f) {
  cat(sprintf("  Processing: %s\n", basename(f)))
  vals <- compute_lambda(f)
  if (is.null(vals)) return(NULL)

  # Extract QC ordering name from filename
  qc_name <- regmatches(basename(f),
                        regexpr("qc[0-9]+", basename(f)))
  if (length(qc_name) == 0) qc_name <- basename(f)

  data.table(
    file        = basename(f),
    qc_ordering = qc_name,
    N           = vals$N,
    n_variants  = vals$n_variants,
    lambda_gc   = round(vals$lambda_gc,   4),
    lambda_1000 = round(vals$lambda_1000, 4),
    median_pval = round(vals$median_pval, 4),
    n_gw_hits   = vals$n_gw_hits
  )
})

results <- rbindlist(Filter(Negate(is.null), results))
results <- results[order(qc_ordering)]

# ── Compute lambda_ratio relative to baseline ──────────────────────────────
baseline_row <- results[qc_ordering == opt$baseline_qc]

if (nrow(baseline_row) == 0) {
  warning(sprintf(
    "Baseline QC '%s' not found. Using first ordering as baseline.",
    opt$baseline_qc
  ))
  baseline_row <- results[1]
}

lambda_baseline      <- baseline_row$lambda_gc
lambda_1000_baseline <- baseline_row$lambda_1000

results[, lambda_ratio      := round(lambda_gc   / lambda_baseline,      4)]
results[, lambda_1000_ratio := round(lambda_1000 / lambda_1000_baseline, 4)]

# ── Print summary ──────────────────────────────────────────────────────────
cat(sprintf("\n=== Lambda comparison across QC orderings ===\n"))
cat(sprintf("Baseline QC: %s  (lambda_gc = %.4f)\n\n",
            opt$baseline_qc, lambda_baseline))

cat(sprintf("%-12s %8s %10s %12s %12s %12s %10s %10s\n",
            "QC ordering", "N", "n_variants",
            "lambda_gc", "lambda_1000",
            "lambda_ratio", "1000_ratio", "gw_hits"))
cat(strrep("-", 92), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i]
  cat(sprintf("%-12s %8d %10d %12.4f %12.4f %12.4f %10.4f %10d\n",
              r$qc_ordering, r$N, r$n_variants,
              r$lambda_gc, r$lambda_1000,
              r$lambda_ratio, r$lambda_1000_ratio,
              r$n_gw_hits))
}

cat("\n")

# ── Interpretation guidance ────────────────────────────────────────────────
lambda_range <- max(results$lambda_gc) - min(results$lambda_gc)
ratio_range  <- max(results$lambda_ratio) - min(results$lambda_ratio)

cat("=== Interpretation ===\n")
cat(sprintf("Range of lambda_gc:     %.4f  (max %.4f, min %.4f)\n",
            lambda_range, max(results$lambda_gc), min(results$lambda_gc)))
cat(sprintf("Range of lambda_ratio:  %.4f  (max %.4f, min %.4f)\n",
            ratio_range, max(results$lambda_ratio), min(results$lambda_ratio)))
cat("\n")

if (ratio_range < 0.02) {
  cat("lambda_ratio range < 0.02: QC orderings produce essentially identical\n")
  cat("calibration — ordering does not substantially affect GWAS inflation.\n")
} else if (ratio_range < 0.05) {
  cat("lambda_ratio range 0.02-0.05: Small but potentially meaningful differences\n")
  cat("in calibration across QC orderings.\n")
} else {
  cat("lambda_ratio range > 0.05: Substantial differences in calibration across\n")
  cat("QC orderings — ordering meaningfully affects GWAS test statistics.\n")
}

cat("\n")
cat("NOTE: lambda_1000 is designed for stratification inflation (lambda > 1).\n")
cat("For your deflated lambda, lambda_ratio is the more valid comparison metric.\n")
cat("Both are reported for completeness.\n")

# ── Write output ───────────────────────────────────────────────────────────
out_file <- paste0(opt$out, "_lambda_comparison.tsv")
fwrite(results, out_file, sep="\t")
cat(sprintf("\nResults written to: %s\n", out_file))
