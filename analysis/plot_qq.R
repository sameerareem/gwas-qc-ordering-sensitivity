#!/usr/bin/env Rscript
# =============================================================================
#  plot_qq.R
#  QQ plot for SAIGE association results
#
#  Usage: Rscript plot_qq.R <saige_result.txt> <output_prefix>
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_qq.R <saige_result.txt> <output_prefix>")
}

saige_file  <- args[1]
out_prefix  <- args[2]

cat("Reading SAIGE output...\n")
results <- fread(saige_file, header=TRUE, sep="\t", data.table=TRUE)

# Show column names for debugging
cat("Columns in result file:", paste(names(results), collapse=", "), "\n")
cat("Rows in result file:", nrow(results), "\n")

# ── Locate the p-value column ──────────────────────────────────────────────
# SAIGE 1.5.1 uses "p.value"; guard against other naming conventions
pcol <- grep("^p[._]?value$", names(results), ignore.case=TRUE, value=TRUE)[1]
if (is.na(pcol)) {
  # Fallback: find any column with "p" and "value" or "pval"
  pcol <- grep("pval|p_val|p[._]val", names(results), ignore.case=TRUE, value=TRUE)[1]
}
if (is.na(pcol)) {
  stop("Cannot find p-value column in SAIGE output. Columns: ",
       paste(names(results), collapse=", "))
}
cat("Using p-value column:", pcol, "\n")

# Coerce explicitly — guards against "NA" strings read as character
pvals <- suppressWarnings(as.numeric(results[[pcol]]))
pvals <- pvals[!is.na(pvals) & !is.nan(pvals) & pvals > 0 & pvals <= 1]

if (length(pvals) == 0) {
  stop("No valid p-values found after filtering NAs and out-of-range values.")
}
cat(sprintf("Valid p-values for QQ plot: %d\n", length(pvals)))

# ── Compute lambda_GC ──────────────────────────────────────────────────────
chisq      <- qchisq(pvals, df=1, lower.tail=FALSE)
lambda_gc  <- median(chisq, na.rm=TRUE) / qchisq(0.5, df=1, lower.tail=FALSE)
n          <- length(pvals)
lambda_1000 <- 1 + (lambda_gc - 1) * (1000 / n)
cat(sprintf("lambda_GC   = %.4f\n", lambda_gc))
cat(sprintf("lambda_1000 = %.4f  (N=%d)\n", lambda_1000, n))

# ── QQ plot ────────────────────────────────────────────────────────────────
expected <- -log10(ppoints(n))
observed <- sort(-log10(pvals), decreasing=TRUE)

ci_upper <- -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1))
ci_lower <- -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1))
lim      <- max(max(observed, na.rm=TRUE), max(expected, na.rm=TRUE)) * 1.05

png(paste0(out_prefix, "_qqplot.png"), width=900, height=900, res=150)
par(mar=c(5,5,4,2), bg="white", family="sans")

plot(expected, observed,
     pch=20, cex=0.4, col="#2166AC",
     xlim=c(0, lim), ylim=c(0, lim),
     xlab=expression(Expected ~ -log[10](p)),
     ylab=expression(Observed ~ -log[10](p)),
     main=paste0("QQ Plot  -  ", basename(out_prefix)),
     cex.main=1.1, col.main="#1F3864")

polygon(c(expected, rev(expected)),
        c(ci_upper, rev(ci_lower)),
        col=adjustcolor("#2166AC", 0.15), border=NA)

abline(0, 1, col="#D73027", lwd=1.5)

legend("topleft", bty="n", cex=0.85,
       legend=c(
         paste0("lambda_GC   = ", round(lambda_gc,   4)),
         paste0("lambda_1000 = ", round(lambda_1000, 4)),
         paste0("N = ", formatC(n, big.mark=",")),
         "95% CI band"
       ),
       fill=c(NA, NA, NA, adjustcolor("#2166AC", 0.15)),
       border=c(NA, NA, NA, "#2166AC"),
       text.col="#1F3864")

dev.off()
cat(sprintf("QQ plot written -> %s_qqplot.png\n", out_prefix))
