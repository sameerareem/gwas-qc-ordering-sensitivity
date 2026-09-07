#!/usr/bin/env Rscript
# =============================================================================
#  compute_qc_metrics.R
#  Post-SAIGE metrics for QC ordering benchmarking
#
#  Reference ordering: QC4. All ratios and paired statistics compare the
#  current ordering against QC4.
#
#  Usage:
#    Rscript compute_qc_metrics.R \
#      --saige_results       <saige_result.txt> \
#      --causal_snps         <causal_snps.tsv> \
#      --final_bim           <final_autosomes.bim> \
#      --final_fam           <final_autosomes.fam> \
#      --qc_label            <e.g. "QC2"> \
#      --ref_results         <qc4_saige_result.txt>  (QC4 output, for paired tests)
#      --ref_bim             <qc4_final_autosomes.bim>
#      --ref_fam             <qc4_final_autosomes.fam>
#      --baseline_lambda     <QC4 lambda_gc>
#      --baseline_lambda1000 <QC4 lambda_1000>
#      --baseline_snps       <qc4_snps.txt>
#      --baseline_samps      <qc4_samps.txt>
#      --out                 <output_prefix>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--saige_results",       type="character"),
  make_option("--causal_snps",         type="character"),
  make_option("--final_bim",           type="character"),
  make_option("--final_fam",           type="character"),
  make_option("--qc_label",            type="character", default="QC"),
  # QC4 reference files (for paired statistical tests)
  make_option("--ref_results",         type="character", default=NULL,
              help="QC4 SAIGE output for paired tests (optional)"),
  make_option("--ref_bim",             type="character", default=NULL),
  make_option("--ref_fam",             type="character", default=NULL),
  # QC4 scalar values (for ratio computation)
  make_option("--baseline_lambda",     type="double",    default=NULL),
  make_option("--baseline_lambda1000", type="double",    default=NULL),
  make_option("--baseline_snps",       type="character", default=NULL),
  make_option("--baseline_samps",      type="character", default=NULL),
  make_option("--seed",                type="character", default=NA_character_,
              help="Simulation seed label (e.g. seed_42). Used as blocking unit."),
  make_option("--dataset",             type="character", default=NA_character_,
              help="Dataset label (e.g. dataset_001). Used as blocking unit."),
  make_option("--out",                 type="character", default="qc_metrics")
)
opt <- parse_args(OptionParser(option_list=option_list))

results <- fread(opt$saige_results)
# Coerce p.value to numeric immediately after reading
# SAIGE may write "NA" strings which cause character inference
pcol_res <- grep("^p[._]?value$", names(results), ignore.case=TRUE, value=TRUE)[1]
if (!is.na(pcol_res) && pcol_res != "p.value") setnames(results, pcol_res, "p.value")
results[, p.value := suppressWarnings(as.numeric(p.value))]
causal  <- fread(opt$causal_snps)
bim     <- fread(opt$final_bim, col.names=c("CHR","ID","CM","POS","A1","A2"))
fam     <- fread(opt$final_fam, col.names=c("FID","IID","PAT","MAT","SEX","PHENO"))

cat("=== QC Ordering Metrics:", opt$qc_label, "(reference: QC4) ===\n\n")

# ── 1. Sample and SNP counts ─────────────────────────────────────────────────
n_samples <- nrow(fam)
n_snps    <- nrow(bim)
n_tested  <- nrow(results)
cat(sprintf("Retained samples : %d\n", n_samples))
cat(sprintf("Retained SNPs    : %d\n", n_snps))
cat(sprintf("SNPs tested      : %d\n", n_tested))
cat("\n")

# ── 2. Genomic inflation (lambda_GC and lambda_1000) ─────────────────────────
pvals    <- results$p.value
pvals    <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
chisq    <- qchisq(pvals, df=1, lower.tail=FALSE)
med_null <- qchisq(0.5, df=1, lower.tail=FALSE)   # 0.4549

lambda_gc   <- median(chisq, na.rm=TRUE) / med_null
lambda_1000 <- 1 + (lambda_gc - 1) * (1000 / n_samples)

cat(sprintf("Lambda GC        : %.4f\n", lambda_gc))
cat(sprintf("Lambda_1000      : %.4f  (N-rescaled to 1000)\n", lambda_1000))

lambda_ratio      <- NA_real_
lambda_1000_ratio <- NA_real_

if (!is.null(opt$baseline_lambda)) {
  lambda_ratio <- lambda_gc / opt$baseline_lambda
  cat(sprintf("Lambda ratio vs QC4      : %.4f\n", lambda_ratio))
}
if (!is.null(opt$baseline_lambda1000)) {
  lambda_1000_ratio <- lambda_1000 / opt$baseline_lambda1000
  cat(sprintf("Lambda_1000 ratio vs QC4 : %.4f\n", lambda_1000_ratio))
}
cat("\n")

# ── 3. Per-chromosome lambda ──────────────────────────────────────────────────
chrom_lambda <- results[!is.na(p.value) & p.value > 0 & p.value <= 1,
  .(lambda    = median(qchisq(p.value, df=1, lower.tail=FALSE), na.rm=TRUE) / med_null,
    n_snps_chr = .N), by=.(CHR)]
setorder(chrom_lambda, CHR)

# ── 4. Causal SNP recovery ────────────────────────────────────────────────────
causal_in_results <- merge(
  causal[, .(SNP, BETA, CHR, POS)],
  results[, .(MarkerID, BETA_est=BETA, SE, p.value, AF_Allele2)],
  by.x="SNP", by.y="MarkerID"
)

n_causal        <- nrow(causal)
n_causal_tested <- nrow(causal_in_results)
causal_in_results[, p.value := suppressWarnings(as.numeric(p.value))]
# Cumulative (nested): recovery_sug includes all GW, recovery_nom includes all sug
n_recovered_gw  <- sum(causal_in_results$p.value < 5e-8, na.rm=TRUE)
n_recovered_sug <- sum(causal_in_results$p.value < 1e-5, na.rm=TRUE)
n_recovered_nom <- sum(causal_in_results$p.value < 0.05, na.rm=TRUE)
# Mutually exclusive bins (sum = n_causal_tested)
n_only_gw    <- n_recovered_gw
n_only_sug   <- n_recovered_sug - n_recovered_gw
n_only_nom   <- n_recovered_nom - n_recovered_sug
n_not_sig    <- n_causal_tested  - n_recovered_nom
n_removed_qc <- n_causal - n_causal_tested

cat(sprintf("Causal SNPs total         : %d\n",   n_causal))
cat(sprintf("Causal SNPs tested        : %d (%.1f%%)\n",
            n_causal_tested, 100*n_causal_tested/n_causal))
cat(sprintf("Recovery p<5e-8 (GW)      : %d (%.1f%%)\n",
            n_recovered_gw,  100*n_recovered_gw/n_causal))
cat(sprintf("Recovery p<1e-5 (sug)     : %d (%.1f%%)\n",
            n_recovered_sug, 100*n_recovered_sug/n_causal))
cat(sprintf("Recovery p<0.05 (nom)     : %d (%.1f%%)\n",
            n_recovered_nom, 100*n_recovered_nom/n_causal))
cat("\n")

# ── 4b. False positive analysis ──────────────────────────────────────────────
# False positive = non-causal SNP that reaches significance threshold.
# null_results = all tested variants NOT in the causal SNP list.
# FP rate = FP / n_null_tested  (expected ≈ alpha under H0)
# FDP      = FP / (FP + TP)     (proportion of "hits" that are spurious)

null_results <- results[!MarkerID %in% causal$SNP]
n_null_tested <- nrow(null_results)

fp_gw  <- sum(null_results$p.value < 5e-8,  na.rm=TRUE)
fp_sug <- sum(null_results$p.value < 1e-5,  na.rm=TRUE)
fp_nom <- sum(null_results$p.value < 0.05,  na.rm=TRUE)

fp_rate_gw  <- if (n_null_tested > 0) fp_gw  / n_null_tested else NA_real_
fp_rate_sug <- if (n_null_tested > 0) fp_sug / n_null_tested else NA_real_
fp_rate_nom <- if (n_null_tested > 0) fp_nom / n_null_tested else NA_real_

# False Discovery Proportion: FP / (FP + TP) — proportion of hits that are wrong
fdp_gw  <- if ((fp_gw  + n_recovered_gw)  > 0) fp_gw  / (fp_gw  + n_recovered_gw)  else NA_real_
fdp_sug <- if ((fp_sug + n_recovered_sug) > 0) fp_sug / (fp_sug + n_recovered_sug) else NA_real_
fdp_nom <- if ((fp_nom + n_recovered_nom) > 0) fp_nom / (fp_nom + n_recovered_nom) else NA_real_

# Lambda from null variants only — inflation in truly null test statistics
# If lambda_null > lambda_gc, inflation is concentrated in non-causal variants
null_pvals   <- null_results$p.value
null_pvals   <- null_pvals[!is.na(null_pvals) & null_pvals > 0 & null_pvals <= 1]
lambda_null  <- if (length(null_pvals) > 10) median(qchisq(null_pvals, df=1, lower.tail=FALSE), na.rm=TRUE) / med_null else NA_real_

cat(sprintf("Null SNPs tested              : %d\n", n_null_tested))
cat(sprintf("False positives p<5e-8 (GW)   : %d  (rate=%.2e  FDP=%.3f)\n",
            fp_gw,  fp_rate_gw,  fdp_gw))
cat(sprintf("False positives p<1e-5 (sug)  : %d  (rate=%.2e  FDP=%.3f)\n",
            fp_sug, fp_rate_sug, fdp_sug))
cat(sprintf("False positives p<0.05 (nom)  : %d  (rate=%.2e  FDP=%.3f)\n",
            fp_nom, fp_rate_nom, fdp_nom))
cat(sprintf("Lambda (null variants only)   : %.4f\n", lambda_null))
cat("\n")


# ── 5. Effect size concordance (vs true betas) ────────────────────────────────
beta_rank_corr <- NA_real_
beta_pearson   <- NA_real_
if (nrow(causal_in_results) > 2) {
  beta_rank_corr <- cor(causal_in_results$BETA, causal_in_results$BETA_est,
                        method="spearman", use="complete.obs")
  beta_pearson   <- cor(causal_in_results$BETA, causal_in_results$BETA_est,
                        method="pearson",  use="complete.obs")
  cat(sprintf("Beta Spearman (vs truth)  : %.4f\n", beta_rank_corr))
  cat(sprintf("Beta Pearson  (vs truth)  : %.4f\n", beta_pearson))
  cat("\n")
}

# ── 6. Jaccard overlap with QC4 ──────────────────────────────────────────────
jaccard_snps  <- NA_real_
jaccard_samps <- NA_real_

if (!is.null(opt$baseline_snps)) {
  baseline_snps  <- fread(opt$baseline_snps,  header=FALSE)$V1
  current_snps   <- bim$ID
  intersection_s <- length(intersect(current_snps, baseline_snps))
  union_s        <- length(union(current_snps, baseline_snps))
  jaccard_snps   <- intersection_s / union_s
  cat(sprintf("Jaccard SNP overlap vs QC4    : %.4f\n", jaccard_snps))
  cat(sprintf("  Shared %d | Only this: %d | Only QC4: %d\n",
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
  cat(sprintf("Jaccard sample overlap vs QC4 : %.4f\n", jaccard_samps))
  cat(sprintf("  Shared %d | Only this: %d | Only QC4: %d\n",
              intersection_sa,
              length(setdiff(current_samps, baseline_samps)),
              length(setdiff(baseline_samps, current_samps))))
  cat("\n")
}

# ── 7. Paired statistical tests vs QC4 (if ref_results provided) ─────────────
# These tests require the QC4 SAIGE output and are only run when --ref_results
# is supplied. They characterise whether any observed differences are
# statistically meaningful, not merely numerical.

mcnemar_gw_p   <- NA_real_
mcnemar_fp_gw_p  <- NA_real_
mcnemar_fp_sug_p <- NA_real_
mcnemar_sug_p  <- NA_real_
mcnemar_nom_p  <- NA_real_
paired_beta_p  <- NA_real_
paired_beta_d  <- NA_real_   # Cohen's d for paired beta difference
lins_ccc       <- NA_real_   # Lin's concordance correlation vs QC4 betas
lins_ccc_p     <- NA_real_

if (!is.null(opt$ref_results)) {
  cat("=== Paired statistical tests vs QC4 ===\n\n")

  ref_res <- fread(opt$ref_results)

  # ── 7a. McNemar's test for causal SNP recovery ──────────────────────────
  # For each of the n_causal simulated SNPs, define recovery as a binary
  # outcome (1 = p < threshold, 0 = not recovered or not tested).
  # McNemar's test detects whether the two orderings differ in WHICH
  # causal SNPs they recover, accounting for the pairing.

  for (thresh_name in c("gw","sug","nom")) {
    thresh <- switch(thresh_name, gw=5e-8, sug=1e-5, nom=0.05)

    this_rec <- causal_in_results[, .(SNP, p_this = p.value)]
    ref_causal <- merge(
      causal[, .(SNP, BETA)],
      ref_res[, .(MarkerID, p_ref=p.value)],
      by.x="SNP", by.y="MarkerID"
    )
    paired_rec <- merge(this_rec, ref_causal[, .(SNP, p_ref)], by="SNP", all=TRUE)
    paired_rec[is.na(p_this), p_this := 1]
    paired_rec[is.na(p_ref),  p_ref  := 1]
    paired_rec[, rec_this := as.integer(p_this < thresh)]
    paired_rec[, rec_ref  := as.integer(p_ref  < thresh)]

    # McNemar contingency: b = this=1,ref=0; c = this=0,ref=1
    b <- sum(paired_rec$rec_this == 1 & paired_rec$rec_ref == 0)
    c <- sum(paired_rec$rec_this == 0 & paired_rec$rec_ref == 1)

    if (b + c > 0) {
      mc_p <- binom.test(min(b,c), b+c, p=0.5)$p.value
    } else {
      mc_p <- 1.0
    }

    cat(sprintf("McNemar recovery test (p<%s) : b=%d c=%d p=%.4f\n",
                switch(thresh_name,gw="5e-8",sug="1e-5",nom="0.05"),
                b, c, mc_p))

    if (thresh_name=="gw")  mcnemar_gw_p  <- mc_p
    if (thresh_name=="sug") mcnemar_sug_p <- mc_p
    if (thresh_name=="nom") mcnemar_nom_p <- mc_p
  }
  cat("\n")

  # ── 7d. McNemar's test for false positives vs QC4 ────────────────────────
  # For each non-causal SNP present in both outputs, test whether the two
  # orderings differ in which null variants they falsely declare significant.
  mcnemar_fp_gw_p  <- NA_real_
  mcnemar_fp_sug_p <- NA_real_

  for (thresh_name in c("fp_gw","fp_sug")) {
    thresh <- switch(thresh_name, fp_gw=5e-8, fp_sug=1e-5)

    # Non-causal SNPs tested in this ordering
    this_null <- null_results[, .(SNP=MarkerID, p_this=p.value)]
    # Non-causal SNPs tested in QC4
    ref_null_res  <- ref_res[!MarkerID %in% causal$SNP,
                             .(SNP=MarkerID, p_ref=p.value)]
    paired_null <- merge(this_null, ref_null_res, by="SNP", all=TRUE)
    paired_null[is.na(p_this), p_this := 1]
    paired_null[is.na(p_ref),  p_ref  := 1]
    paired_null[, fp_this := as.integer(p_this < thresh)]
    paired_null[, fp_ref  := as.integer(p_ref  < thresh)]

    b_fp <- sum(paired_null$fp_this==1 & paired_null$fp_ref==0)
    c_fp <- sum(paired_null$fp_this==0 & paired_null$fp_ref==1)

    mc_fp <- if (b_fp + c_fp > 0) binom.test(min(b_fp,c_fp), b_fp+c_fp, p=0.5)$p.value else 1.0

    cat(sprintf("McNemar FP test (p<%s)    : b=%d c=%d p=%.4f\n",
                switch(thresh_name, fp_gw="5e-8", fp_sug="1e-5"),
                b_fp, c_fp, mc_fp))

    if (thresh_name=="fp_gw")  mcnemar_fp_gw_p  <- mc_fp
    if (thresh_name=="fp_sug") mcnemar_fp_sug_p <- mc_fp
  }
  cat("\n")


  # ── 7b. Paired t-test on beta differences (this ordering vs QC4) ─────────
  # For each causal SNP present in BOTH outputs, test whether the estimated
  # beta from this ordering differs systematically from the QC4 estimate.
  # H0: mean(beta_this - beta_QC4) = 0

  ref_betas <- ref_res[MarkerID %in% causal$SNP,
                       .(SNP=MarkerID, beta_ref=BETA)]
  this_betas <- causal_in_results[, .(SNP, beta_this=BETA_est)]
  paired_betas <- merge(this_betas, ref_betas, by="SNP")
  paired_betas <- paired_betas[!is.na(beta_this) & !is.na(beta_ref)]

  if (nrow(paired_betas) > 2) {
    diffs <- paired_betas$beta_this - paired_betas$beta_ref
    tt <- t.test(diffs, mu=0)
    paired_beta_p <- tt$p.value
    paired_beta_d <- mean(diffs, na.rm=TRUE) / sd(diffs, na.rm=TRUE)  # Cohen's d

    cat(sprintf("Paired t-test (beta_this - beta_QC4):\n"))
    cat(sprintf("  mean diff = %.5f | SD = %.5f\n", mean(diffs), sd(diffs)))
    cat(sprintf("  t = %.3f | df = %d | p = %.4f\n",
                tt$statistic, round(tt$parameter), paired_beta_p))
    cat(sprintf("  Cohen's d = %.4f\n", paired_beta_d))
    cat("\n")

    # ── 7c. Lin's Concordance Correlation Coefficient (CCC) ───────────────
    # CCC = (2 * rho * sd_x * sd_y) / (var_x + var_y + (mu_x - mu_y)^2)
    # Combines Pearson correlation (precision) with location/scale accuracy.
    # CCC = 1 means perfect concordance with QC4; < 1 means either
    # different ordering of estimates OR systematic shift in magnitude.
    x   <- paired_betas$beta_ref
    y   <- paired_betas$beta_this
    rho <- cor(x, y, use="complete.obs")
    sx  <- sd(x, na.rm=TRUE)
    sy  <- sd(y, na.rm=TRUE)
    mx  <- mean(x, na.rm=TRUE)
    my  <- mean(y, na.rm=TRUE)
    lins_ccc <- (2 * rho * sx * sy) / (sx^2 + sy^2 + (mx - my)^2)

    # Approximate 95% CI via Fisher z-transformation on rho component
    n_ccc  <- sum(!is.na(x) & !is.na(y))
    z_rho  <- 0.5 * log((1 + rho) / (1 - rho))
    se_z   <- 1 / sqrt(n_ccc - 3)
    ccc_ci <- tanh(c(z_rho - 1.96*se_z, z_rho + 1.96*se_z))  # approx CI on rho

    # One-sample test: H0: CCC = 1 (perfect concordance with QC4)
    # Use Fisher z on CCC itself
    z_ccc     <- 0.5 * log((1 + lins_ccc) / (1 - lins_ccc))
    z_null    <- 0.5 * log((1 + 1.0) / (1 - 1.0 + 1e-9))  # ~Inf under H0: CCC=1
    # Practical test: CCC significantly < 1 if CI upper bound < 1
    lins_ccc_p <- 2 * pnorm(-abs((lins_ccc - 1) / (se_z * (1 - lins_ccc^2 + 1e-9))))

    cat(sprintf("Lin's CCC (this vs QC4)   : %.4f  (approx 95%% CI on rho: %.4f-%.4f)\n",
                lins_ccc, ccc_ci[1], ccc_ci[2]))
    cat(sprintf("  Pearson r component     : %.4f\n", rho))
    cat(sprintf("  Location shift (mu_diff): %.5f\n", my - mx))
    cat(sprintf("  Scale shift (sy/sx)     : %.4f\n", sy/sx))
    cat("\n")
  }
}

# ── 8. Write summary TSV ──────────────────────────────────────────────────────
summary_dt <- data.table(
  qc_label              = opt$qc_label,
  seed                  = opt$seed,
  dataset               = opt$dataset,
  block_id              = paste(opt$seed, opt$dataset, sep="_"),
  n_samples             = n_samples,
  n_snps                = n_snps,
  n_snps_tested         = n_tested,
  lambda_gc             = round(lambda_gc,   4),
  lambda_1000           = round(lambda_1000, 4),
  lambda_ratio          = round(lambda_ratio,      4),
  lambda_1000_ratio     = round(lambda_1000_ratio, 4),
  causal_tested         = n_causal_tested,
  recovery_gw_sig       = n_recovered_gw,
  recovery_gw_pct       = round(100 * n_recovered_gw  / n_causal, 2),
  recovery_sug          = n_recovered_sug,
  recovery_sug_pct      = round(100 * n_recovered_sug / n_causal, 2),
  recovery_nom          = n_recovered_nom,
  recovery_nom_pct      = round(100 * n_recovered_nom / n_causal, 2),
  bin_gw_sig         = n_only_gw,
  bin_sug_only       = n_only_sug,
  bin_nom_only       = n_only_nom,
  bin_tested_not_sig = n_not_sig,
  bin_removed_by_qc  = n_removed_qc,
  beta_spearman         = round(beta_rank_corr, 4),
  beta_pearson          = round(beta_pearson,   4),
  jaccard_snps          = round(jaccard_snps,   4),
  jaccard_samples       = round(jaccard_samps,  4),
  mcnemar_gw_p          = round(mcnemar_gw_p,  4),
  mcnemar_sug_p         = round(mcnemar_sug_p, 4),
  mcnemar_nom_p         = round(mcnemar_nom_p, 4),
  paired_beta_p         = round(paired_beta_p, 4),
  paired_beta_cohens_d  = round(paired_beta_d, 4),
  lins_ccc              = round(lins_ccc,      4),
  # False positive metrics
  n_null_tested         = n_null_tested,
  fp_gw                 = fp_gw,
  fp_rate_gw            = round(fp_rate_gw,  6),
  fdp_gw                = round(fdp_gw,      4),
  fp_sug                = fp_sug,
  fp_rate_sug           = round(fp_rate_sug, 6),
  fdp_sug               = round(fdp_sug,     4),
  fp_nom                = fp_nom,
  fp_rate_nom           = round(fp_rate_nom, 6),
  fdp_nom               = round(fdp_nom,     4),
  lambda_null           = round(lambda_null,  4),
  mcnemar_fp_gw_p       = round(mcnemar_fp_gw_p,  4),
  mcnemar_fp_sug_p      = round(mcnemar_fp_sug_p, 4)
)
out_tsv <- paste0(opt$out, "_metrics.tsv")
fwrite(summary_dt, out_tsv, sep="\t")
cat(sprintf("Metrics written -> %s\n\n", out_tsv))

# =============================================================================
#  VISUALISATIONS
# =============================================================================

chr_cols <- rep(c("#2166AC","#4DAC26"), 11)

# ── Plot 1: Manhattan plot ────────────────────────────────────────────────────
cat("Generating Manhattan plot...\n")
man_data <- results[!is.na(p.value) & p.value > 0 & p.value <= 1,
                    .(CHR, POS, p.value, MarkerID)]
chr_maxpos <- man_data[, .(maxpos=max(POS,na.rm=TRUE)), by=CHR]
setorder(chr_maxpos, CHR)
chr_maxpos[, cum_offset := cumsum(shift(maxpos,fill=0)) + (CHR-1)*5e6]
man_data <- merge(man_data, chr_maxpos[,.(CHR, cum_offset)], by="CHR")
man_data[, pos_cum   := POS + cum_offset]
man_data[, log10p    := -log10(p.value)]
man_data[, is_causal := MarkerID %in% causal$SNP]
man_data[, is_fp     := !is_causal & log10p >= -log10(5e-8)]
chr_mids <- man_data[, .(mid=mean(pos_cum)), by=CHR]
setorder(chr_mids, CHR)

png(paste0(opt$out,"_manhattan.png"), width=2400, height=900, res=150)
par(mar=c(4,5,3,2), bg="white", family="sans")
plot(man_data$pos_cum, man_data$log10p,
     col  = ifelse(man_data$is_causal, "#D73027",
            ifelse(man_data$is_fp,     "#7B2D8B",
                   chr_cols[man_data$CHR])),
     pch  = ifelse(man_data$is_causal, 18,
            ifelse(man_data$is_fp,     17, 20)),
     cex  = ifelse(man_data$is_causal, 1.0,
            ifelse(man_data$is_fp,     0.9, 0.35)),
     xaxt="n", yaxt="n",
     xlab="Chromosome", ylab=expression(-log[10](p)),
     main=paste0("Manhattan Plot  -  ", opt$qc_label,
                 "   (lGC=", round(lambda_gc,4),
                 "  l1000=", round(lambda_1000,4), ")"),
     cex.main=1.1, col.main="#1F3864",
     ylim=c(0, max(man_data$log10p,na.rm=TRUE)*1.1))
abline(h=-log10(5e-8), col="#D73027", lty=2, lwd=1.2)
abline(h=-log10(1e-5), col="#F4A582", lty=3, lwd=1.0)
axis(1, at=chr_mids$mid, labels=chr_mids$CHR, las=1, cex.axis=0.65, tick=FALSE)
axis(2, las=1, cex.axis=0.8)
legend("topright", bty="n", cex=0.72,
       legend=c("Non-causal SNP","Simulated causal SNP",
                "False positive (non-causal, p<5e-8)",
                "GW sig (p<5e-8)","Suggestive (p<1e-5)"),
       col=c("#2166AC","#D73027","#7B2D8B","#D73027","#F4A582"),
       pch=c(20,18,17,NA,NA), lty=c(NA,NA,NA,2,3),
       lwd=c(NA,NA,NA,1.2,1.0))
dev.off()
cat(sprintf("  -> %s_manhattan.png\n", opt$out))

# ── Plot 2: QQ plot ───────────────────────────────────────────────────────────
cat("Generating QQ plot...\n")
n_p      <- length(pvals)
expected <- -log10(ppoints(n_p))
observed <- sort(-log10(pvals), decreasing=TRUE)
ci_upper <- -log10(qbeta(0.025, seq_len(n_p), n_p-seq_len(n_p)+1))
ci_lower <- -log10(qbeta(0.975, seq_len(n_p), n_p-seq_len(n_p)+1))
# Axis limit: cap at expected maximum + 20% headroom.
# Using the observed maximum would compress the entire plot when a single
# extreme outlier exists (e.g. a numerically unstable SNP with -log10(p)=300+).
# The expected maximum for n_p SNPs is -log10(1/n_p).
lim_exp <- max(expected, na.rm=TRUE) * 1.20   # expected range + 20%
lim_obs <- quantile(observed[is.finite(observed)], 0.9999, na.rm=TRUE) * 1.10
lim     <- max(lim_exp, lim_obs)               # show at least the 99.99th percentile

# Count and report outliers that fall above the plot limit
n_outliers <- sum(observed > lim, na.rm=TRUE)

png(paste0(opt$out,"_qqplot.png"), width=1000, height=1000, res=150)
par(mar=c(5,5,4,9), bg="white", family="sans")
plot(expected, pmin(observed, lim),   # clip outliers to top of plot
     pch=20, cex=0.4, col="#2166AC",
     xlim=c(0, lim), ylim=c(0, lim),
     xlab=expression(Expected ~ -log[10](p)),
     ylab=expression(Observed ~ -log[10](p)),
     main=paste0("QQ Plot  -  ", opt$qc_label),
     cex.main=1.1, col.main="#1F3864")

# Mark clipped outliers as triangles at the top edge
if (n_outliers > 0) {
  out_idx <- which(observed > lim)
  points(expected[out_idx], rep(lim * 0.99, n_outliers),
         pch=17, cex=0.9, col="#D73027")
  mtext(sprintf("%d outlier(s) clipped (max = %.0f)", n_outliers,
                max(observed[is.finite(observed)], na.rm=TRUE)),
        side=3, cex=0.7, col="#D73027", line=0.1)
}
polygon(c(expected, rev(expected)), c(ci_upper, rev(ci_lower)),
        col=adjustcolor("#2166AC",0.15), border=NA)
abline(0,1, col="#D73027", lwd=1.5)
legend(x=par("usr")[2]*1.02, y=par("usr")[4],
       bty="n", cex=0.80, xpd=TRUE,
       legend=c(paste0("lambda_GC   = ",  round(lambda_gc,4)),
                paste0("lambda_1000 = ",  round(lambda_1000,4)),
                paste0("N = ", formatC(n_samples, big.mark=",")),
                "95% CI band"),
       fill=c(NA,NA,NA,adjustcolor("#2166AC",0.15)),
       border=c(NA,NA,NA,"#2166AC"),
       text.col="#1F3864")
dev.off()
cat(sprintf("  -> %s_qqplot.png\n", opt$out))

# ── Plot 3: Recovery bar chart ────────────────────────────────────────────────
cat("Generating recovery plot...\n")
rec_vals <- c(n_recovered_gw, n_recovered_sug, n_recovered_nom)
rec_pct  <- 100 * rec_vals / n_causal
rec_labs <- c("GW sig\n(p<5e-8)","Suggestive\n(p<1e-5)","Nominal\n(p<0.05)")
rec_cols <- c("#D73027","#FC8D59","#91BFDB")
mc_ps    <- c(mcnemar_gw_p, mcnemar_sug_p, mcnemar_nom_p)
mc_stars <- ifelse(is.na(mc_ps), "",
              ifelse(mc_ps < 0.001, "***",
              ifelse(mc_ps < 0.01,  "**",
              ifelse(mc_ps < 0.05,  "*", "ns"))))

png(paste0(opt$out,"_recovery.png"), width=900, height=720, res=150)
par(mar=c(6,5,4,2), bg="white", family="sans")
bp <- barplot(rec_pct, col=rec_cols, names.arg=rec_labs, ylim=c(0,125),
              ylab="% of causal SNPs recovered",
              main=paste0("Causal SNP Recovery  -  ", opt$qc_label,
                          "\n(", n_causal, " simulated causal SNPs)"),
              cex.names=0.8, cex.main=1.0, col.main="#1F3864",
              border="white", las=1)
text(bp, rec_pct + 3.5,
     labels=paste0(rec_vals,"\n(",round(rec_pct,1),"%)"),
     cex=0.75, col="#1F3864", font=2)
# McNemar significance stars
text(bp, rec_pct + 14, labels=mc_stars, cex=1.1, col="#D73027", font=2)
if (any(!is.na(mc_ps)))
  mtext("* p<0.05  ** p<0.01  *** p<0.001  (McNemar vs QC4)",
        side=1, line=5, cex=0.7, col="#595959")
abline(h=100*n_causal_tested/n_causal, lty=2, col="#595959", lwd=1.2)
dev.off()
cat(sprintf("  -> %s_recovery.png\n", opt$out))

# ── Plot 4a: Effect size vs simulated truth (PRIMARY — always produced) ────────
# This is the scientifically meaningful comparison: how close is each ordering's
# beta estimate to the ground-truth simulated effect size?
# The vs-QC4 comparison (Plot 4b) is secondary — it only tells you whether two
# pipelines agree with each other, not whether either is close to the truth.
if (nrow(causal_in_results) > 2) {
  cat("Generating effect size vs truth scatter (primary)...\n")

  # Colour by significance status at GW threshold
  sig_col <- ifelse(causal_in_results$p.value < 5e-8,  "#D73027",   # GW significant
             ifelse(causal_in_results$p.value < 1e-5,  "#FC8D59",   # suggestive
                    "#91BFDB"))                                        # not significant

  rng <- range(c(causal_in_results$BETA, causal_in_results$BETA_est), na.rm=TRUE)
  rng <- rng + c(-0.08, 0.08)*diff(rng)

  png(paste0(opt$out,"_effectsize_vs_truth.png"), width=950, height=950, res=150)
  par(mar=c(5,5,4,2), bg="white", family="sans")

  plot(causal_in_results$BETA, causal_in_results$BETA_est,
       pch  = 21,
       bg   = sig_col,
       col  = adjustcolor(sig_col, 0.8),
       cex  = 1.4,
       xlim = rng, ylim = rng,
       xlab = "True effect size (simulated \u03b2)",
       ylab = paste0("Estimated effect size (SAIGE \u03b2)  \u2014  ", opt$qc_label),
       main = paste0("Effect Size vs Ground Truth  \u2014  ", opt$qc_label),
       cex.main = 1.0, col.main = "#1F3864")

  abline(0, 1, col="#2166AC", lty=2, lwd=1.8)   # identity: perfect recovery
  abline(h=0, v=0, col="grey80", lwd=0.8)

  # OLS regression line — slope < 1 indicates attenuation bias
  fit_truth <- lm(BETA_est ~ BETA, data=causal_in_results)
  abline(fit_truth, col="#D73027", lwd=1.3, lty=3)

  # Annotate regression slope (attenuation factor)
  slope_val <- round(coef(fit_truth)[2], 3)

  legend("topleft", bty="n", cex=0.75, inset=c(0.01,0.01), bg=adjustcolor("white",0.85),
         legend=c(
           paste0("Spearman r = ", round(beta_rank_corr, 3)),
           paste0("Pearson r  = ", round(beta_pearson,   3)),
           paste0("OLS slope  = ", slope_val,
                  if (slope_val < 0.95) "  (attenuation)" else ""),
           paste0("N tested   = ", n_causal_tested, "/", n_causal),
           "Identity (slope=1)",
           "OLS fit",
           "GW significant (p<5e-8)",
           "Suggestive (p<1e-5)",
           "Not significant"
         ),
         col  = c(NA,NA,NA,NA,"#2166AC","#D73027","#D73027","#FC8D59","#91BFDB"),
         pch  = c(NA,NA,NA,NA,NA,NA,21,21,21),
         pt.bg= c(NA,NA,NA,NA,NA,NA,"#D73027","#FC8D59","#91BFDB"),
         lty  = c(NA,NA,NA,NA,2,3,NA,NA,NA),
         lwd  = c(NA,NA,NA,NA,1.8,1.3,NA,NA,NA),
         text.col="#1F3864")
  dev.off()
  cat(sprintf("  -> %s_effectsize_vs_truth.png\n", opt$out))
}

# ── Plot 4b: Effect size vs QC4 estimates (SECONDARY — only with ref_results) ─
# Compares two pipelines against each other. Useful for quantifying ordering
# differences, but should be interpreted alongside the vs-truth plot above.
if (!is.null(opt$ref_results) && exists("paired_betas") && nrow(paired_betas) > 2) {
  cat("Generating effect size vs QC4 scatter (secondary)...\n")

  rng2 <- range(c(paired_betas$beta_ref, paired_betas$beta_this), na.rm=TRUE)
  rng2 <- rng2 + c(-0.08, 0.08)*diff(rng2)

  # Colour by whether GW significance status changed vs QC4
  pt_col <- ifelse(
    (paired_betas$SNP %in% causal_in_results[p.value < 5e-8]$SNP) &
    !(paired_betas$SNP %in% ref_res[p.value < 5e-8]$MarkerID), "#D73027",
    ifelse(
      !(paired_betas$SNP %in% causal_in_results[p.value < 5e-8]$SNP) &
       (paired_betas$SNP %in% ref_res[p.value < 5e-8]$MarkerID), "#4DAC26",
      "#FC8D59"))

  png(paste0(opt$out,"_effectsize_vs_qc4.png"), width=950, height=950, res=150)
  par(mar=c(5,5,4,2), bg="white", family="sans")

  plot(paired_betas$beta_ref, paired_betas$beta_this,
       pch=21, bg=pt_col, col=adjustcolor(pt_col,0.8), cex=1.4,
       xlim=rng2, ylim=rng2,
       xlab="Beta estimate  (QC4 reference)",
       ylab=paste0("Beta estimate  (", opt$qc_label, ")"),
       main=paste0("Effect Size: ", opt$qc_label, " vs QC4  (pipeline comparison)"),
       cex.main=1.0, col.main="#1F3864")

  abline(0,1, col="#2166AC", lty=2, lwd=1.8)
  abline(h=0, v=0, col="grey80", lwd=0.8)
  fit2 <- lm(beta_this ~ beta_ref, data=paired_betas)
  abline(fit2, col="#595959", lwd=1.2, lty=3)

  legend("topleft", bty="n", cex=0.75, inset=c(0.01,0.01), bg=adjustcolor("white",0.85),
         legend=c(
           paste0("CCC = ",        round(lins_ccc,     3)),
           paste0("Pearson r = ",  round(cor(paired_betas$beta_ref,
                                             paired_betas$beta_this), 3)),
           paste0("Paired t p = ", round(paired_beta_p, 3)),
           paste0("Cohen's d = ",  round(paired_beta_d, 3)),
           "Identity (slope=1)",
           "Gained GW sig vs QC4",
           "Lost GW sig vs QC4",
           "Both / neither"
         ),
         col  = c(NA,NA,NA,NA,"#2166AC","#D73027","#4DAC26","#FC8D59"),
         pch  = c(NA,NA,NA,NA,NA,21,21,21),
         pt.bg= c(NA,NA,NA,NA,NA,"#D73027","#4DAC26","#FC8D59"),
         lty  = c(NA,NA,NA,NA,2,NA,NA,NA),
         lwd  = c(NA,NA,NA,NA,1.8,NA,NA,NA),
         text.col="#1F3864")
  dev.off()
  cat(sprintf("  -> %s_effectsize_vs_qc4.png\n", opt$out))
}

# ── Plot 5: Per-chromosome lambda ─────────────────────────────────────────────
cat("Generating per-chromosome lambda plot...\n")

# y-axis: start at 0, extend 20% above max(bars, null=1.0)
ylim_top <- max(max(chrom_lambda$lambda, na.rm=TRUE), 1.0) * 1.20
ylim_chr <- c(0, ylim_top)

n_chrs   <- nrow(chrom_lambda)
png_w    <- max(1600, n_chrs * 68)

png(paste0(opt$out,"_chr_lambda.png"), width=png_w, height=900, res=150)
# right margin=10 gives space for outside legend and reference labels
par(mar=c(7, 5.5, 4, 15), bg="white", family="sans")

bp <- barplot(chrom_lambda$lambda,
              col      = chr_cols[chrom_lambda$CHR],
              ylim     = ylim_chr,
              ylab     = expression(lambda[GC] ~ "(per chromosome)"),
              xlab     = "",
              main     = paste0("Per-Chromosome Genomic Inflation  \u2014  ", opt$qc_label),
              cex.main = 1.0, col.main = "#1F3864",
              border   = "white", las = 1, axes = TRUE,
              names.arg = rep("", n_chrs))

# X-axis labels at 45 degrees
text(x=bp, y=-0.012*ylim_top, labels=chrom_lambda$CHR,
     srt=45, adj=c(1,1), xpd=TRUE, cex=0.82, col="#1F3864")
mtext("Chromosome", side=1, line=5.2, cex=0.95, col="#1F3864")

# Reference lines
abline(h=1.0,         col="#D73027", lty=2, lwd=1.8)
abline(h=lambda_gc,   col="#4DAC26", lty=3, lwd=1.4)
abline(h=lambda_1000, col="#2166AC", lty=4, lwd=1.4)

# Single legend outside the plot — contains values, no separate margin text
usr <- par("usr")
legend(x=usr[2]+0.015*(usr[2]-usr[1]),
       y=ylim_top * 0.98,
       legend=c(
         "Null (lambda = 1.000)",
         paste0("Overall lambda_GC   = ", round(lambda_gc,   4)),
         paste0("Overall lambda_1000 = ", round(lambda_1000, 4))
       ),
       col=c("#D73027","#4DAC26","#2166AC"),
       lty=c(2,3,4), lwd=c(1.8,1.4,1.4),
       bty="n", cex=0.75, xpd=TRUE, x.intersp=0.8, y.intersp=1.6,
       seg.len=1.8)

# SNP count inside bars
text(x=bp, y=chrom_lambda$lambda/2,
     labels=chrom_lambda$n_snps_chr, cex=0.52, col="white", font=2)

dev.off()
cat(sprintf("  -> %s_chr_lambda.png\n", opt$out))


# ── Plot 6: Genome-wide significant hits — True Positives vs False Positives ───
# Shows ONLY SNPs that crossed the GW significance threshold (p < 5e-8).
# Left strip  = non-causal SNPs that are GW significant (false positives).
# Right strip = causal SNPs that are GW significant (true positives).
# Each point is one SNP. Points above the threshold by definition.
# Directly answers: of everything declared significant, what was real?
cat("Generating TP/FP significance plot...\n")

# GW-significant non-causal SNPs (false positives)
fp_gw_data <- null_results[!is.na(p.value) & p.value < 5e-8,
                            .(MarkerID, p.value, CHR)]
fp_gw_logp <- -log10(fp_gw_data$p.value)

# GW-significant causal SNPs (true positives)
tp_gw_data <- causal_in_results[!is.na(p.value) & p.value < 5e-8,
                                 .(SNP, p.value, BETA, BETA_est)]
tp_gw_logp <- -log10(tp_gw_data$p.value)

# Also include causal SNPs that were tested but missed (not GW sig) — shown faded
tp_miss_data <- causal_in_results[!is.na(p.value) & p.value >= 5e-8]
tp_miss_logp <- -log10(tp_miss_data$p.value[tp_miss_data$p.value > 0])

# y-axis top: accommodate all points plus headroom
all_logp <- c(fp_gw_logp, tp_gw_logp)
ylim_top <- max(c(all_logp[is.finite(all_logp)], -log10(5e-8)*1.3), na.rm=TRUE)
ylim_fp  <- c(0, ylim_top+50)

png(paste0(opt$out,"_fp_tp_plot.png"), width=950, height=950, res=150)
par(mar=c(7,5,4,3), bg="white", family="sans")

# Base plot — empty frame
plot(NULL,
     xlim=c(0.3, 2.7), ylim=ylim_fp,
     xaxt="n", yaxt="n",
     xlab="", ylab=expression(-log[10](p)),
     main=paste0("Genome-Wide Significant Hits: TP vs FP  \u2014  ", opt$qc_label),
     cex.main=0.95, col.main="#1F3864")

# GW threshold line (the minimum y for all plotted points)
abline(h=-log10(5e-8), col="#D73027", lty=2, lwd=1.8)
mtext(expression("GW threshold (p < 5" %*% "10"^{-8} * ")"),
      side=4, at=-log10(1e-1), las=2, cex=0.65, col="#D73027", line=-5)

# Grey background band below threshold — shows the "null" region
rect(par("usr")[1], 0, par("usr")[2], -log10(5e-8),
     col=adjustcolor("grey90", 0.4), border=NA)

# Missed causal SNPs (tested but not GW sig) — faint blue on right strip
if (length(tp_miss_logp) > 0) {
  points(jitter(rep(2, length(tp_miss_logp)), amount=0.18),
         tp_miss_logp,
         pch=21, cex=0.9,
         bg=adjustcolor("#91BFDB", 0.5),
         col=adjustcolor("#91BFDB", 0.7))
}

# False positive points — left strip (purple triangles)
if (length(fp_gw_logp) > 0) {
  points(jitter(rep(1, length(fp_gw_logp)), amount=0.18),
         fp_gw_logp,
         pch=25, cex=1.2,
         bg="#7B2D8B", col="#5A1A6B")
} else {
  # No false positives — print "0 FP" label in strip
  text(1, -log10(5e-8)*1.15, "0 FP", cex=1.0, col="#7B2D8B", font=2)
}

# True positive points — right strip (red circles)
if (length(tp_gw_logp) > 0) {
  points(jitter(rep(2, length(tp_gw_logp)), amount=0.18),
         tp_gw_logp,
         pch=21, cex=1.4,
         bg="#D73027", col="#A50026")
}

# Count labels above each strip
text(1, ylim_top * 0.97,
     sprintf("FP = %d\n(%.1f%% of hits)", fp_gw,
             ifelse((fp_gw+n_recovered_gw)>0, 100*fp_gw/(fp_gw+n_recovered_gw), 0)),
     cex=0.85, col="#7B2D8B", font=2)
text(2, ylim_top * 0.97,
     sprintf("TP = %d\n(%.1f%% of hits)", n_recovered_gw,
             ifelse((fp_gw+n_recovered_gw)>0, 100*n_recovered_gw/(fp_gw+n_recovered_gw), 0)),
     cex=0.85, col="#D73027", font=2)

# X-axis labels
axis(1, at=c(1, 2),
     labels=c(
       sprintf("Non-causal SNPs\n(FP: %d of %s tested)", fp_gw,
               formatC(n_null_tested, big.mark=",")),
       sprintf("Causal SNPs\n(TP: %d  Missed: %d  of %d)", n_recovered_gw,
               n_causal_tested - n_recovered_gw, n_causal_tested)
     ),
     cex.axis=0.78, tick=FALSE, line=1.5)
axis(2, las=1, cex.axis=0.85)

# FDP annotation
fdp_label <- ifelse(is.na(fdp_gw), "NA",
                    sprintf("%.1f%%", 100*fdp_gw))
mtext(sprintf("FDP = %s  |  lambda_null = %.4f", fdp_label, lambda_null),
      side=1, line=5.5, cex=0.78, col="#1F3864")



legend("topleft", bty="n", cex=0.82,
       legend=c(
         sprintf("False positive (non-causal, GW sig)  n=%d", fp_gw),
         sprintf("True positive  (causal, GW sig)      n=%d", n_recovered_gw),
         sprintf("Causal SNP tested but missed         n=%d", n_causal_tested - n_recovered_gw)
       ),
       col=c("#7B2D8B","#D73027","#91BFDB"),
       pch=c(25, 21, 21),
       pt.bg=c("#7B2D8B","#D73027",adjustcolor("#91BFDB",0.5)),
       pt.cex=c(1.2,1.4,0.9),
       text.col="#1F3864")

dev.off()
cat(sprintf("  -> %s_fp_tp_plot.png\n", opt$out))

cat("\nAll outputs written with prefix:", opt$out, "\n")
