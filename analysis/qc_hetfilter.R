#!/usr/bin/env Rscript
# =============================================================================
#  qc_hetfilter.R
#  Identifies heterozygosity outliers from a PLINK --het output file.
#
#  Usage:
#    Rscript qc_hetfilter.R <het_file> <threads> [output_file]
#
#  Arguments:
#    het_file    : PLINK .het output file (from plink --het)
#    threads     : number of threads (passed for compatibility, not used in R)
#    output_file : output file for samples to remove (default: qc_test_het_fail.txt)
#
#  Removes samples whose observed heterozygosity is > 3 SD from the mean.
#  Output format: FID IID (space-separated, for plink --remove)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript qc_hetfilter.R <het_file> [threads] [output_file]")
}

het_file    <- args[1]
# args[2] is threads — accepted but not used in R
output_file <- if (length(args) >= 3) args[3] else "qc_test_het_fail.txt"

cat(sprintf("[INFO] Reading het file: %s\n", het_file))
cat(sprintf("[INFO] Output file: %s\n", output_file))

het <- read.table(het_file, header = TRUE, stringsAsFactors = FALSE)

# Compute observed heterozygosity rate
# PLINK .het columns: FID IID O(HOM) E(HOM) N(NM) F
# Observed het = (N(NM) - O(HOM)) / N(NM)
het$het_rate <- (het$N.NM. - het$O.HOM.) / het$N.NM.

mean_het <- mean(het$het_rate, na.rm = TRUE)
sd_het   <- sd(het$het_rate,   na.rm = TRUE)

# Flag outliers > 3 SD from mean
het$outlier <- abs(het$het_rate - mean_het) > 3 * sd_het

n_outliers <- sum(het$outlier, na.rm = TRUE)
cat(sprintf("[INFO] Mean het rate: %.4f, SD: %.4f\n", mean_het, sd_het))
cat(sprintf("[INFO] Threshold: mean ± 3 SD = [%.4f, %.4f]\n",
            mean_het - 3*sd_het, mean_het + 3*sd_het))

# Print count table as PLINK-style output
print(table(het$outlier))

# Write outlier samples to output file
outliers <- het[het$outlier & !is.na(het$outlier), c("FID", "IID")]
write.table(outliers, output_file,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

cat(sprintf("[INFO] Written %d het outliers to %s\n", n_outliers, output_file))
