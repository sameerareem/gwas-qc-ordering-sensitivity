#!/usr/bin/env Rscript
# =============================================================================
#  compare_to_qc4.R
#  Cross-ordering statistical comparison with QC4 as reference.
#
#  UNIT OF ANALYSIS: simulation seed (n=5 independent replicates).
#
#  Datasets within a seed share the same autosomal genotype matrix and
#  causal variant assignment. They differ only in sex chromosome
#  representation (% sex-distorted samples), making them technical
#  pseudo-replicates rather than independent observations. All hypothesis
#  tests therefore use seed-level means as the unit of analysis. Within-seed
#  SD across datasets is retained as a secondary measure of sensitivity to
#  sex chromosome representation.
#
#  Usage:
#    Rscript compare_to_qc4.R \
#      --metrics_dir  <directory containing all *_metrics.tsv files> \
#      --out          <output_prefix>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--metrics_dir", type="character"),
  make_option("--out",         type="character", default="qc4_comparison")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ── Load all metric TSV files ─────────────────────────────────────────────────
files <- list.files(opt$metrics_dir, pattern="_metrics\\.tsv$",
                    full.names=TRUE, recursive=TRUE)
# Exclude the comparison output directory to avoid reading our own output TSVs
files <- files[!grepl("/comparison/", files, fixed=TRUE)]
if (length(files) == 0)
  stop("No *_metrics.tsv files found in: ", opt$metrics_dir)

all_metrics <- rbindlist(lapply(files, fread), fill=TRUE)
setorder(all_metrics, qc_label, seed, dataset)

n_seeds    <- length(unique(all_metrics$seed))
n_datasets <- length(unique(all_metrics$dataset))
n_orderings <- length(unique(all_metrics$qc_label))

cat("=== Cross-ordering comparison vs QC4 ===\n")
cat(sprintf("Loaded %d rows from %d files\n", nrow(all_metrics), length(files)))
cat(sprintf("Seeds     (independent replicates) : %d  -> %s\n",
            n_seeds, paste(sort(unique(all_metrics$seed)), collapse=", ")))
cat(sprintf("Datasets  (pseudo-replicates/seed) : %d  -> %s\n",
            n_datasets, paste(sort(unique(all_metrics$dataset)), collapse=", ")))
cat(sprintf("QC labels : %s\n\n",
            paste(sort(unique(all_metrics$qc_label)), collapse=", ")))

if (n_datasets > 1) {
  cat("NOTE: Datasets within a seed share the same autosomal genotype matrix.\n")
  cat("      They are NOT independent replicates. The unit of analysis is the\n")
  cat("      seed (n=", n_seeds, " independent replicates).\n")
  cat("      Within-seed SD across datasets is reported separately as a measure\n")
  cat("      of sensitivity to sex chromosome representation only.\n\n")
}

metrics_primary <- c("lambda_1000","lambda_1000_ratio",
                     "recovery_gw_pct","recovery_sug_pct",
                     "fp_gw","fp_rate_gw","fdp_gw",
                     "beta_spearman","jaccard_snps","jaccard_samples",
                     "lambda_null")
metrics_primary <- metrics_primary[metrics_primary %in% names(all_metrics)]

stars <- function(p) ifelse(is.na(p), "  NA",
                    ifelse(p < 0.001, "***",
                    ifelse(p < 0.01,  " **",
                    ifelse(p < 0.05,  "  *", " ns"))))

# =============================================================================
#  STEP 1: Aggregate to seed level
#  For each (qc_label, seed) compute the mean across datasets.
#  This is the correct unit of analysis for all hypothesis tests.
#  Also compute within-seed SD to quantify sex-chr sensitivity.
# =============================================================================

seed_means <- all_metrics[,
  c(lapply(.SD, function(x) mean(as.numeric(x), na.rm=TRUE)),
    list(n_datasets_used = sum(!is.na(as.numeric(.SD[[1]]))))),
  by=.(qc_label, seed),
  .SDcols=metrics_primary]

seed_sd <- all_metrics[,
  lapply(.SD, function(x) sd(as.numeric(x), na.rm=TRUE)),
  by=.(qc_label, seed),
  .SDcols=metrics_primary]
setnames(seed_sd, metrics_primary, paste0(metrics_primary, "_within_seed_sd"))

seed_level <- merge(seed_means, seed_sd, by=c("qc_label","seed"))
setorder(seed_level, qc_label, seed)

# Write seed-level aggregated data
fwrite(seed_level, paste0(opt$out, "_seed_level.tsv"), sep="\t")
cat(sprintf("Seed-level means -> %s_seed_level.tsv\n\n", opt$out))

# Also write full raw metrics for reference
fwrite(all_metrics, paste0(opt$out, "_all_metrics.tsv"), sep="\t")

ref_seed  <- seed_level[grepl("qc4|QC4", qc_label, ignore.case=TRUE)]
comp_seed <- seed_level[!grepl("qc4|QC4", qc_label, ignore.case=TRUE)]
orderings <- sort(unique(comp_seed$qc_label))

# =============================================================================
#  A. Descriptive summary: mean ± SD of seed-level means per ordering
#     SD here is between-seed variation — the genuine biological variability
# =============================================================================

cat("=== A. Descriptive summary (mean ± SD of seed-level means) ===\n")
cat(sprintf("    Unit: seed mean (n=%d seeds", n_seeds))
if (n_datasets > 1)
  cat(sprintf(", each averaged over %d datasets", n_datasets))
cat(")\n\n")

summary_tbl <- seed_level[, lapply(.SD, function(x) {
  x <- as.numeric(x)
  sprintf("%.4f ± %.4f", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
}), by=qc_label, .SDcols=metrics_primary]
setorder(summary_tbl, qc_label)
print(summary_tbl)
cat("\n")

if (n_datasets > 1) {
  cat("Within-seed SD (sensitivity to sex-chr representation):\n\n")
  sd_cols <- paste0(metrics_primary, "_within_seed_sd")
  sd_cols <- sd_cols[sd_cols %in% names(seed_level)]
  within_sd_tbl <- seed_level[, lapply(.SD, function(x)
    sprintf("%.5f", mean(as.numeric(x), na.rm=TRUE))),
    by=qc_label, .SDcols=sd_cols]
  setnames(within_sd_tbl, sd_cols, sub("_within_seed_sd","",sd_cols))
  setorder(within_sd_tbl, qc_label)
  print(within_sd_tbl)
  cat("\n")
}

# =============================================================================
#  B. Friedman test — omnibus across all 6 orderings
#     Block = seed (n=5 independent replicates)
#     Unit  = seed-level mean (averaged across datasets within seed)
# =============================================================================

cat(sprintf("=== B. Friedman omnibus test (blocked by seed, n=%d blocks) ===\n\n",
            n_seeds))
friedman_results <- list()

for (m in metrics_primary) {
  dt_m <- seed_level[!is.na(get(m)), .(qc_label, seed, val=as.numeric(get(m)))]
  # Keep only seeds with all 6 orderings present
  n_present_orderings <- length(unique(dt_m$qc_label))
  complete_seeds <- dt_m[, .N, by=seed][N == n_present_orderings]$seed
  dt_m <- dt_m[seed %in% complete_seeds]
  if (length(complete_seeds) < 3) {
    cat(sprintf("%-22s  skipped (only %d complete seeds)\n", m, length(complete_seeds)))
    next
  }

  mat <- dcast(dt_m, seed ~ qc_label, value.var="val")
  mat[, seed := NULL]
  mat_clean <- mat[complete.cases(mat)]

  ft <- tryCatch(friedman.test(as.matrix(mat_clean)), error=function(e) NULL)
  if (!is.null(ft)) {
    friedman_results[[m]] <- data.table(
      metric   = m,
      n_blocks = nrow(mat_clean),
      chi2     = round(ft$statistic, 3),
      df       = ft$parameter,
      p        = ft$p.value,
      sig      = stars(ft$p.value)
    )
    cat(sprintf("%-22s  n_seeds=%d  chi2=%6.3f  df=%d  p=%.4f  %s\n",
                m, nrow(mat_clean), ft$statistic,
                ft$parameter, ft$p.value, stars(ft$p.value)))
  }
}
cat("\n")

if (length(friedman_results) > 0) {
  fwrite(rbindlist(friedman_results), paste0(opt$out,"_friedman.tsv"), sep="\t")
  cat(sprintf("Friedman results -> %s_friedman.tsv\n\n", opt$out))
}

# =============================================================================
#  C. Wilcoxon signed-rank: each ordering vs QC4
#     Paired by seed (n=5 pairs)
#     Multiple testing: Benjamini-Hochberg FDR
# =============================================================================

cat(sprintf("=== C. Wilcoxon signed-rank vs QC4 (paired by seed, n=%d pairs, BH FDR) ===\n\n",
            n_seeds))
if (n_seeds < 6)
  cat(sprintf("NOTE: With %d pairs the minimum achievable p-value is %.4f.\n",
              n_seeds, 2 * (0.5)^n_seeds))
if (n_datasets > 1)
  cat(sprintf("      Seed means are averaged over %d datasets; datasets are not counted as independent observations.\n", n_datasets))
cat("\n")

wilcox_rows <- list()

for (m in metrics_primary) {
  for (ord in orderings) {
    qc4_dt  <- ref_seed[!is.na(get(m)), .(seed, val_ref =as.numeric(get(m)))]
    comp_dt <- comp_seed[qc_label==ord & !is.na(get(m)),
                         .(seed, val_this=as.numeric(get(m)))]
    paired  <- merge(qc4_dt, comp_dt, by="seed")
    if (nrow(paired) < 3) next

    diffs <- paired$val_this - paired$val_ref
    if (all(diffs == 0, na.rm=TRUE)) next

    wt <- tryCatch(
      wilcox.test(paired$val_this, paired$val_ref, paired=TRUE, exact=FALSE),
      error=function(e) NULL)
    if (is.null(wt)) next

    effect_r <- abs(qnorm(wt$p.value/2)) / sqrt(nrow(paired))
    wilcox_rows[[length(wilcox_rows)+1]] <- data.table(
      metric      = m,
      ordering    = ord,
      n_seeds     = nrow(paired),
      n_datasets  = n_datasets,
      median_diff = round(median(diffs, na.rm=TRUE), 5),
      mean_diff   = round(mean(diffs,   na.rm=TRUE), 5),
      sd_diff     = round(sd(diffs,     na.rm=TRUE), 5),
      W           = wt$statistic,
      p_raw       = wt$p.value,
      effect_r    = round(effect_r, 3)
    )
  }
}

if (length(wilcox_rows) > 0) {
  wilcox_dt <- rbindlist(wilcox_rows)
  wilcox_dt[, p_adj := p.adjust(p_raw, method="BH")]
  wilcox_dt[, sig   := stars(p_adj)]
  setorder(wilcox_dt, metric, ordering)

  cat(sprintf("%-22s  %-6s  %5s  %9s  %9s  %8s  %8s  %s\n",
              "Metric","Order","Seeds","Med_diff","Mean_diff","p_raw","p_BH","sig"))
  cat(strrep("-",82),"\n")
  for (i in seq_len(nrow(wilcox_dt))) {
    r <- wilcox_dt[i]
    cat(sprintf("%-22s  %-6s  %5d  %9.4f  %9.4f  %8.4f  %8.4f  %s\n",
                r$metric, r$ordering, r$n_seeds,
                r$median_diff, r$mean_diff, r$p_raw, r$p_adj, r$sig))
  }
  fwrite(wilcox_dt, paste0(opt$out,"_wilcoxon.tsv"), sep="\t")
  cat(sprintf("\nWilcoxon results -> %s_wilcoxon.tsv\n\n", opt$out))
}

# =============================================================================
#  VISUALISATIONS
#  All plots use seed-level means as data points.
#  Error bars/bands show between-seed SD (genuine biological variability).
# =============================================================================

qc_order <- c("QC1","QC2","QC3","QC4","QC5","QC6")
present  <- intersect(qc_order, unique(seed_level$qc_label))
pal      <- c(QC1="#E41A1C",QC2="#377EB8",QC3="#4DAF4A",
              QC4="#984EA3",QC5="#FF7F00",QC6="#A65628")
col_map  <- pal[present]
ref_col  <- pal["QC4"]

# ── Plot 1: lambda_1000 — seed-level means with individual points ─────────────
if ("lambda_1000" %in% names(seed_level)) {
  cat("Generating lambda_1000 comparison plot...\n")
  png(paste0(opt$out,"_lambda1000_boxplot.png"), width=1200, height=700, res=150)
  par(mar=c(5,5,4,3), bg="white", family="sans")

  l1k_list <- lapply(present, function(q) seed_level[qc_label==q]$lambda_1000)
  names(l1k_list) <- present
  bp <- boxplot(l1k_list, col=adjustcolor(col_map,0.4),
                border=col_map, lwd=1.5,
                ylab=expression(lambda[1000] ~ "(seed-level mean)"),
                main=expression("Genomic Inflation (" * lambda[1000] * ") — Seed-Level Means"),
                sub=sprintf("n=%d independent seeds%s",
                  n_seeds,
                  if(n_datasets>1) sprintf(", each averaged over %d datasets (non-independent)",
                                           n_datasets) else ""),
                cex.main=1.0, col.main="#1F3864", las=1, outline=FALSE)

  # Overlay individual seed points
  for (i in seq_along(present)) {
    y <- seed_level[qc_label==present[i]]$lambda_1000
    points(jitter(rep(i, length(y)), amount=0.12), y,
           pch=21, bg=adjustcolor(col_map[i],0.7),
           col=col_map[i], cex=1.1)
  }
  abline(h=1.0, lty=2, col="#D73027", lwd=1.2)
  abline(v=which(present=="QC4"), col=ref_col, lty=3, lwd=1.5)
  mtext("QC4 (reference)", side=3, at=which(present=="QC4"), cex=0.75, col=ref_col)
  dev.off()
  cat(sprintf("  -> %s_lambda1000_boxplot.png\n", opt$out))
}

# ── Plot 2: Seed interaction plot ─────────────────────────────────────────────
# Lines connect seed-level means across orderings.
# Parallel lines = consistent ordering effect across all seeds.
# Crossing lines = ordering effect depends on specific genotype realisation.
if ("lambda_1000" %in% names(seed_level)) {
  cat("Generating seed interaction plot...\n")
  seeds_uniq <- sort(unique(seed_level$seed))
  seed_pal   <- colorRampPalette(c("#1B7837","#762A83"))(length(seeds_uniq))
  names(seed_pal) <- seeds_uniq

  png(paste0(opt$out,"_seed_interaction.png"), width=1400, height=700, res=150)
  par(mar=c(5,5,4,9), bg="white", family="sans", xpd=FALSE)

  x_pos     <- seq_along(present)
  all_vals  <- seed_level[qc_label %in% present]$lambda_1000
  ylim      <- range(all_vals, na.rm=TRUE)
  ylim      <- ylim + c(-0.05, 0.12)*diff(ylim)

  plot(NULL, xlim=range(x_pos)+c(-0.3,0.3), ylim=ylim,
       xaxt="n", xlab="QC Ordering",
       ylab=expression(lambda[1000] ~ "(seed-level mean)"),
       main=expression("Seed Interaction Plot: " * lambda[1000] * " by QC Ordering"),
       sub=if(n_datasets>1)
         sprintf("Each point = mean across %d datasets within seed (non-independent)",
                 n_datasets) else
         "Each point = single dataset per seed",
       cex.main=1.0, col.main="#1F3864")
  axis(1, at=x_pos, labels=present)

  # Between-seed SD band (grey)
  ord_sum <- seed_level[, .(m=mean(lambda_1000,na.rm=TRUE),
                              s=sd(lambda_1000,na.rm=TRUE)), by=qc_label]
  ord_sum <- ord_sum[match(present, qc_label)]
  polygon(c(x_pos, rev(x_pos)),
          c(ord_sum$m + ord_sum$s, rev(ord_sum$m - ord_sum$s)),
          col=adjustcolor("grey70",0.3), border=NA)
  lines(x_pos, ord_sum$m, col="grey40", lwd=2.5)

  # One line per seed
  for (s in seeds_uniq) {
    s_dat <- seed_level[seed==s][match(present, qc_label)]
    lines(x_pos, s_dat$lambda_1000, col=seed_pal[s], lwd=1.5)
    points(x_pos, s_dat$lambda_1000, col=seed_pal[s], pch=19, cex=1.0)
  }

  abline(v=which(present=="QC4"), col=ref_col, lty=3, lwd=1.2)
  abline(h=1.0, col="#D73027", lty=2, lwd=0.9)

  par(xpd=TRUE)
  legend(max(x_pos)+0.7, mean(ylim),
         legend=c(seeds_uniq, "Overall mean", "Mean \u00b1 SD band"),
         col=c(seed_pal,"grey40",adjustcolor("grey70",0.4)),
         lwd=c(rep(1.5,length(seeds_uniq)),2.5,8),
         lty=1, bty="n", cex=0.72, y.intersp=1.1)
  dev.off()
  cat(sprintf("  -> %s_seed_interaction.png\n", opt$out))
}

# ── Plot 3: Dataset heatmap (only when >1 dataset per seed) ──────────────────
if (n_datasets > 1) {
  cat("Generating dataset heatmap (deviation from QC4)...\n")
  datasets_uniq <- sort(unique(all_metrics$dataset))
  n_dat <- length(datasets_uniq)
  n_ord <- length(present)
  n_met <- length(metrics_primary)

  # Mean across seeds, per ordering × dataset
  hmap_dt <- all_metrics[, lapply(.SD, function(x) mean(as.numeric(x),na.rm=TRUE)),
                          by=.(qc_label, dataset), .SDcols=metrics_primary]

  png(paste0(opt$out,"_dataset_heatmap.png"),
      width=300+200*n_ord, height=200+160*n_met, res=150)
  par(mfrow=c(n_met,1), mar=c(3,10,2,6), bg="white", family="sans")

  for (m in metrics_primary) {
    mat_vals <- matrix(NA, nrow=n_dat, ncol=n_ord,
                       dimnames=list(datasets_uniq, present))
    for (d in datasets_uniq)
      for (q in present)
        mat_vals[d,q] <- hmap_dt[dataset==d & qc_label==q, get(m)]

    # Deviation from QC4 column
    mat_plot <- if ("QC4" %in% present) sweep(mat_vals, 1, mat_vals[,"QC4"], "-") else mat_vals

    # Skip this metric if all values are NA (e.g. ratio metrics where QC4=NA)
    if (all(is.na(mat_plot))) next
    zmax <- max(abs(mat_plot), na.rm=TRUE)
    if (is.na(zmax) || !is.finite(zmax) || zmax==0) zmax <- 1
    pal_div <- colorRampPalette(c("#2166AC","white","#D73027"))(51)

    image(t(mat_plot), col=pal_div, zlim=c(-zmax,zmax),
          xaxt="n", yaxt="n",
          main=paste0(m," (deviation from QC4, averaged across seeds)"),
          cex.main=0.8, col.main="#1F3864")
    axis(1, at=seq(0,1,length.out=n_ord), labels=present, las=2, cex.axis=0.8)
    if (n_dat > 1) {
      axis(2, at=seq(0,1,length.out=n_dat), labels=datasets_uniq, las=2, cex.axis=0.8)
    } else {
      axis(2, at=0.5, labels=datasets_uniq, las=2, cex.axis=0.8)
    }
    for (i in seq_len(n_ord))
      for (j in seq_len(n_dat)) {
        v   <- mat_plot[j,i]
        x_c <- if (n_ord>1) (i-1)/(n_ord-1) else 0.5
        y_c <- if (n_dat>1) (j-1)/(n_dat-1) else 0.5
        if (!is.na(v))
          text(x_c, y_c, sprintf("%.3f",v), cex=0.65,
               col=ifelse(abs(v)>0.6*zmax,"white","#1F3864"))
      }
  }
  dev.off()
  cat(sprintf("  -> %s_dataset_heatmap.png\n", opt$out))
  cat("  NOTE: Heatmap uses means across seeds. Deviations reflect\n")
  cat("        sex-chr representation differences, not independent replication.\n")
} else {
  cat("  Skipping dataset heatmap: only one dataset per seed.\n")
}

# ── Plot 4: Recovery boxplot ──────────────────────────────────────────────────
if ("recovery_gw_pct" %in% names(seed_level)) {
  cat("Generating recovery comparison plot...\n")
  png(paste0(opt$out,"_recovery_boxplot.png"), width=1200, height=700, res=150)
  par(mar=c(5,5,4,3), bg="white", family="sans")

  rec_list <- lapply(present, function(q) seed_level[qc_label==q]$recovery_gw_pct)
  names(rec_list) <- present
  bp <- boxplot(rec_list, col=adjustcolor(col_map,0.4),
                border=col_map, lwd=1.5,
                ylab="% causal SNPs recovered (p < 5e-8)",
                main="Genome-Wide Causal SNP Recovery — Seed-Level Means",
                sub=sprintf("n=%d independent seeds", n_seeds),
                cex.main=1.0, col.main="#1F3864",
                ylim=c(0, max(unlist(rec_list),na.rm=TRUE)*1.15),
                las=1, outline=FALSE)
  for (i in seq_along(present)) {
    y <- seed_level[qc_label==present[i]]$recovery_gw_pct
    points(jitter(rep(i,length(y)),amount=0.12), y,
           pch=21, bg=adjustcolor(col_map[i],0.7), col=col_map[i], cex=1.1)
  }
  abline(v=which(present=="QC4"), col=ref_col, lty=3, lwd=1.5)
  mtext("QC4 (reference)", side=3, at=which(present=="QC4"), cex=0.75, col=ref_col)
  dev.off()
  cat(sprintf("  -> %s_recovery_boxplot.png\n", opt$out))
}

# ── Plot 4b: FP / FDP boxplot across orderings ────────────────────────────────
# Parallel to the recovery boxplot — shows false positive burden per ordering.
# Uses seed-level means as data points (consistent with all other plots).
if ("fdp_gw" %in% names(seed_level) && any(!is.na(seed_level$fdp_gw))) {
  cat("Generating FP/FDP comparison plot...\n")

  png(paste0(opt$out,"_fp_comparison.png"), width=1200, height=700, res=150)
  par(mar=c(5,5,4,3), bg="white", family="sans")

  fdp_list <- lapply(present, function(q) seed_level[qc_label==q]$fdp_gw * 100)
  names(fdp_list) <- present

  bp_fp <- boxplot(fdp_list,
                   col=adjustcolor(col_map, 0.4),
                   border=col_map, lwd=1.5,
                   ylab="False Discovery Proportion at GW threshold (%)",
                   main="False Discovery Proportion (GW) Across QC Orderings",
                   sub=sprintf("n=%d independent seeds  |  FDP = FP / (FP + TP) at p < 5\u00d710\u207b\u2078",
                               n_seeds),
                   cex.main=1.0, col.main="#1F3864",
                   ylim=c(0, max(unlist(fdp_list)*1.2, 5, na.rm=TRUE)),
                   las=1, outline=FALSE)

  # Overlay individual seed points
  for (i in seq_along(present)) {
    y <- seed_level[qc_label==present[i]]$fdp_gw * 100
    points(jitter(rep(i, length(y)), amount=0.12), y,
           pch=21, bg=adjustcolor(col_map[i], 0.7),
           col=col_map[i], cex=1.1)
  }

  abline(h=0, lty=1, col="grey80", lwd=0.8)
  abline(v=which(present=="QC4"), col=ref_col, lty=3, lwd=1.5)
  mtext("QC4 (reference)", side=3, at=which(present=="QC4"), cex=0.75, col=ref_col)
  dev.off()
  cat(sprintf("  -> %s_fp_comparison.png\n", opt$out))
}

# ── Plot 4c: Combined Power vs FDP tradeoff plot ──────────────────────────────
# The most informative single cross-ordering plot.
# X-axis = GW recovery rate (power, %), Y-axis = FDP at GW threshold (%).
# Each point is one ordering's seed-level mean. Error bars = between-seed SD.
# Ideal ordering: top-left (high power, low FDP).
if (all(c("recovery_gw_pct","fdp_gw") %in% names(seed_level))) {
  cat("Generating power vs FDP tradeoff plot...\n")

  tradeoff <- seed_level[, .(
    mean_power = mean(recovery_gw_pct, na.rm=TRUE),
    sd_power   = sd(recovery_gw_pct,   na.rm=TRUE),
    mean_fdp   = mean(fdp_gw * 100,    na.rm=TRUE),
    sd_fdp     = sd(fdp_gw * 100,      na.rm=TRUE),
    n          = sum(!is.na(recovery_gw_pct))
  ), by=qc_label]
  setorder(tradeoff, qc_label)

  # Only keep orderings with data
  tradeoff <- tradeoff[qc_label %in% present]
  tradeoff[, col := col_map[qc_label]]

  x_range <- range(c(0, tradeoff$mean_power + tradeoff$sd_power), na.rm=TRUE)
  y_range <- range(c(0, tradeoff$mean_fdp   + tradeoff$sd_fdp),   na.rm=TRUE)
  x_range[2] <- x_range[2] * 1.15
  y_range[2] <- max(y_range[2] * 1.15, 5)

  png(paste0(opt$out,"_power_vs_fdp.png"), width=900, height=900, res=150)
  par(mar=c(5,5,4,2), bg="white", family="sans")

  plot(NULL,
       xlim=x_range, ylim=y_range,
       xlab="GW Causal SNP Recovery (%, power)",
       ylab="False Discovery Proportion at GW threshold (%)",
       main="Power vs False Discovery Proportion\nper QC Ordering (seed-level means \u00b1 SD)",
       cex.main=0.95, col.main="#1F3864")

  # Ideal quadrant shading: top-right bad, bottom-left ideal
  rect(par("usr")[1], par("usr")[3],
       mean(x_range), mean(y_range),
       col=adjustcolor("#4DAC26", 0.06), border=NA)
  text(mean(c(par("usr")[1], mean(x_range))),
       mean(c(par("usr")[3], mean(y_range))),
       "Ideal\n(high power,\nlow FDP)",
       cex=0.65, col=adjustcolor("#4DAC26", 0.6), font=3)

  # Grid
  abline(h=0, v=0, col="grey80", lwd=0.8)

  # Error bars and points
  for (i in seq_len(nrow(tradeoff))) {
    r <- tradeoff[i]
    if (is.na(r$mean_power) || is.na(r$mean_fdp)) next

    # X error bar (power SD)
    arrows(r$mean_power - r$sd_power, r$mean_fdp,
           r$mean_power + r$sd_power, r$mean_fdp,
           angle=90, code=3, length=0.05, lwd=1.2,
           col=adjustcolor(r$col, 0.6))
    # Y error bar (FDP SD)
    arrows(r$mean_power, r$mean_fdp - r$sd_fdp,
           r$mean_power, r$mean_fdp + r$sd_fdp,
           angle=90, code=3, length=0.05, lwd=1.2,
           col=adjustcolor(r$col, 0.6))

    # Point
    points(r$mean_power, r$mean_fdp,
           pch=ifelse(r$qc_label=="QC4", 23, 21),
           bg=r$col, col=adjustcolor(r$col, 0.8),
           cex=ifelse(r$qc_label=="QC4", 1.8, 1.5))

    # Label
    text(r$mean_power, r$mean_fdp,
         labels=r$qc_label,
         pos=3, cex=0.78, col=r$col, font=2)
  }

  legend("topright", bty="n", cex=0.78,
         legend=c(present, "QC4 reference"),
         col=c(col_map[present], ref_col),
         pch=c(ifelse(present=="QC4", 23, rep(21, length(present))), 23),
         pt.bg=c(col_map[present], ref_col),
         pt.cex=1.3, text.col="#1F3864")

  dev.off()
  cat(sprintf("  -> %s_power_vs_fdp.png\n", opt$out))
}


# ── Plot 4d: TP/FP grouped bar chart across QC orderings ─────────────────────
# Mean recovered causal SNPs (TP) and false positives (FP) per ordering,
# averaged across seeds. Stacked: GW-sig (green) | tested not sig (blue) |
# removed by QC (grey) | false positives (red, right bar).
if (all(c("recovery_gw_sig","fp_gw") %in% names(seed_level))) {
  cat("Generating TP/FP cross-ordering comparison bar chart...\n")

  has_bins <- "bin_removed_by_qc" %in% names(seed_level)

  tp_fp_summary <- seed_level[, {
    lst <- list(
      mean_tp_gw  = mean(recovery_gw_sig, na.rm=TRUE),
      mean_fp_gw  = mean(fp_gw,           na.rm=TRUE),
      sd_tp       = sd(recovery_gw_sig,   na.rm=TRUE),
      sd_fp       = sd(fp_gw,             na.rm=TRUE)
    )
    if (has_bins) {
      lst$mean_missed  <- mean(bin_tested_not_sig, na.rm=TRUE)
      lst$mean_removed <- mean(bin_removed_by_qc,  na.rm=TRUE)
    } else {
      lst$mean_missed  <- mean(pmax(0, causal_tested - recovery_gw_sig), na.rm=TRUE)
      lst$mean_removed <- mean(pmax(0, 50 - causal_tested),              na.rm=TRUE)
    }
    lst
  }, by=qc_label]
  tp_fp_summary <- tp_fp_summary[qc_label %in% present]
  setorder(tp_fp_summary, qc_label)

  n_ord  <- nrow(tp_fp_summary)
  bar_w  <- 0.35
  x_pos  <- seq_len(n_ord)
  y_max  <- max(50, max(tp_fp_summary$mean_fp_gw + tp_fp_summary$sd_fp, na.rm=TRUE)) * 1.3

  png(paste0(opt$out, "_tp_fp_bar.png"), width=1100, height=800, res=150)
  par(mar=c(6, 5, 4, 8), bg="white", family="sans")

  plot(NULL, xlim=c(0.5, n_ord + 0.5), ylim=c(0, y_max),
       xaxt="n", las=1,
       xlab="", ylab="Mean SNP count (averaged across seeds)",
       main=paste0("Causal SNP Recovery and False Positives by QC Ordering\n",
                   "GW threshold (p < 5\u00d710\u207b\u2078), n=", n_seeds, " seeds"),
       cex.main=0.95, col.main="#1F3864")

  for (i in seq_len(n_ord)) {
    r  <- tp_fp_summary[i]

    # Left bar (causal SNPs): stacked segments
    xl <- x_pos[i] - bar_w/2 - 0.02
    base <- 0
    # Removed by QC
    rect(xl-bar_w/2, base, xl+bar_w/2, base+r$mean_removed, col="#CCCCCC", border="white")
    base <- base + r$mean_removed
    # Tested, not significant
    rect(xl-bar_w/2, base, xl+bar_w/2, base+r$mean_missed,  col="#92C5DE", border="white")
    base <- base + r$mean_missed
    # GW significant (TP)
    rect(xl-bar_w/2, base, xl+bar_w/2, base+r$mean_tp_gw,   col="#4DAC26", border="white")
    # TP error bar
    arrows(xl, base + r$mean_tp_gw - r$sd_tp,
           xl, base + r$mean_tp_gw + r$sd_tp,
           angle=90, code=3, length=0.04, lwd=1.2, col="#1F3864")

    # Right bar: false positives
    xr <- x_pos[i] + bar_w/2 + 0.02
    rect(xr-bar_w/2, 0, xr+bar_w/2, r$mean_fp_gw, col="#D73027", border="white")
    if (!is.na(r$sd_fp) && r$sd_fp > 0)
      arrows(xr, max(0, r$mean_fp_gw-r$sd_fp), xr, r$mean_fp_gw+r$sd_fp,
             angle=90, code=3, length=0.04, lwd=1.2, col="#A50026")

    if (r$qc_label=="QC4")
      mtext("ref", side=3, at=x_pos[i], cex=0.65, col=ref_col, line=0.1)
  }

  axis(1, at=x_pos, labels=present, cex.axis=0.9, tick=FALSE, line=0.5)
  mtext("QC Ordering", side=1, line=3.5, cex=0.9, col="#1F3864")
  abline(h=50, lty=2, col="grey50", lwd=0.8)
  mtext("n=50 causal", side=4, at=50, las=2, cex=0.6, col="grey50", line=0.3)

  par(xpd=TRUE)
  legend(n_ord+0.6, y_max*0.95,
         legend=c("GW sig causal (TP)","Tested, not sig","Removed by QC","FP (non-causal GW sig)"),
         fill=c("#4DAC26","#92C5DE","#CCCCCC","#D73027"),
         border="white", bty="n", cex=0.75, y.intersp=1.2)
  dev.off()
  cat(sprintf("  -> %s_tp_fp_bar.png\n", opt$out))
}

if (exists("wilcox_dt") && nrow(wilcox_dt) > 0) {
  cat("Generating significance heatmap...\n")
  all_mets <- unique(wilcox_dt$metric)
  all_ords <- unique(wilcox_dt$ordering)
  heat_mat <- matrix(NA, nrow=length(all_mets), ncol=length(all_ords),
                     dimnames=list(all_mets, all_ords))
  for (i in seq_len(nrow(wilcox_dt)))
    heat_mat[wilcox_dt$metric[i], wilcox_dt$ordering[i]] <- wilcox_dt$p_adj[i]

  png(paste0(opt$out,"_pvalue_heatmap.png"), width=900, height=700, res=150)
  par(mar=c(6,12,4,4), bg="white", family="sans")
  log_heat <- -log10(pmax(heat_mat, 1e-4))
  image(t(log_heat),
        col=colorRampPalette(c("white","#FEE0D2","#FC9272","#DE2D26"))(50),
        xaxt="n", yaxt="n",
        main=sprintf("Wilcoxon BH-adjusted p-values vs QC4\n(-log10; paired by seed, n=%d)",
                     n_seeds),
        col.main="#1F3864", cex.main=0.95)
  axis(1, at=seq(0,1,length.out=ncol(log_heat)),
       labels=colnames(log_heat), las=2, cex.axis=0.9)
  axis(2, at=seq(0,1,length.out=nrow(log_heat)),
       labels=rownames(log_heat), las=2, cex.axis=0.85)
  abline(h=seq(0,1,length.out=nrow(log_heat)), col="grey90", lwd=0.5)
  abline(v=seq(0,1,length.out=ncol(log_heat)), col="grey90", lwd=0.5)
  for (i in seq_len(nrow(log_heat)))
    for (j in seq_len(ncol(log_heat))) {
      pv <- heat_mat[i,j]
      if (!is.na(pv))
        x_pos <- if (ncol(log_heat) > 1) (j-1)/(ncol(log_heat)-1) else 0.5
        y_pos <- if (nrow(log_heat) > 1) (i-1)/(nrow(log_heat)-1) else 0.5
        text(x_pos, y_pos, stars(pv), cex=0.9,
             col=ifelse(pv<0.05,"white","#595959"))
    }
  dev.off()
  cat(sprintf("  -> %s_pvalue_heatmap.png\n", opt$out))
}

# ── Plot 6: Lin's CCC bar chart ───────────────────────────────────────────────
if ("lins_ccc" %in% names(comp_seed) && any(!is.na(comp_seed$lins_ccc))) {
  cat("Generating Lin's CCC bar chart...\n")
  ccc_sum <- comp_seed[!is.na(lins_ccc),
    .(mean_ccc=mean(lins_ccc,na.rm=TRUE),
      sd_ccc  =sd(lins_ccc,  na.rm=TRUE),
      n       =sum(!is.na(lins_ccc))), by=qc_label]
  setorder(ccc_sum, qc_label)

  png(paste0(opt$out,"_lins_ccc.png"), width=900, height=650, res=150)
  par(mar=c(5,5,4,2), bg="white", family="sans")
  bp_c <- barplot(ccc_sum$mean_ccc,
                  names.arg=ccc_sum$qc_label,
                  col=pal[ccc_sum$qc_label],
                  ylim=c(0,1.15),
                  ylab="Lin's CCC vs QC4 (seed-level mean)",
                  main=sprintf("Effect Size Concordance with QC4 (Lin's CCC)\nMean \u00b1 SD across %d seeds", n_seeds),
                  cex.names=0.9, cex.main=0.95, col.main="#1F3864",
                  border="white", las=1)
  arrows(bp_c,
         pmax(0, ccc_sum$mean_ccc - ccc_sum$sd_ccc),
         bp_c,
         pmin(1.1, ccc_sum$mean_ccc + ccc_sum$sd_ccc),
         angle=90, code=3, length=0.06, lwd=1.5, col="#595959")
  abline(h=1.0, lty=2, col="#D73027", lwd=1.2)
  text(bp_c, 0.02, paste0("n=",ccc_sum$n), cex=0.72, col="white", font=2)
  dev.off()
  cat(sprintf("  -> %s_lins_ccc.png\n", opt$out))
}


# =============================================================================
#  PAIRWISE HEATMAPS — All QC ordering pairs
#  For each metric: a 6x6 matrix of mean differences (seed-level means)
#  with Wilcoxon paired tests (BH-corrected across all 15 pairs per metric).
#
#  Fill    = mean difference (row ordering minus column ordering)
#  Cell text = mean diff + significance stars
#  Diagonal  = 0 by definition
#  Matrix is antisymmetric: diff[i,j] = -diff[j,i]
# =============================================================================

cat("\n=== D. Pairwise comparison heatmaps (all ordering pairs) ===\n\n")

# ── Helper: run all pairwise Wilcoxon tests for one metric ───────────────────
pairwise_wilcox <- function(seed_level, metric, orderings) {
  n_ord <- length(orderings)
  diff_mat  <- matrix(0,   nrow=n_ord, ncol=n_ord, dimnames=list(orderings, orderings))
  p_mat     <- matrix(NA,  nrow=n_ord, ncol=n_ord, dimnames=list(orderings, orderings))
  raw_p_vec <- c()
  pair_ids  <- list()

  k <- 0
  for (i in seq_len(n_ord-1)) {
    for (j in (i+1):n_ord) {
      oi <- orderings[i]; oj <- orderings[j]
      di <- seed_level[qc_label==oi & !is.na(get(metric)),
                       .(seed, vi=as.numeric(get(metric)))]
      dj <- seed_level[qc_label==oj & !is.na(get(metric)),
                       .(seed, vj=as.numeric(get(metric)))]
      paired <- merge(di, dj, by="seed")
      if (nrow(paired) < 3) next
      diffs <- paired$vi - paired$vj
      if (all(diffs == 0, na.rm=TRUE)) {
        raw_p_vec <- c(raw_p_vec, 1.0)
      } else {
        wt <- tryCatch(
          wilcox.test(paired$vi, paired$vj, paired=TRUE, exact=FALSE),
          error=function(e) list(p.value=NA))
        raw_p_vec <- c(raw_p_vec, wt$p.value)
      }
      k <- k + 1
      pair_ids[[k]] <- c(oi, oj)
      diff_mat[oi, oj] <-  mean(diffs, na.rm=TRUE)
      diff_mat[oj, oi] <- -mean(diffs, na.rm=TRUE)
    }
  }

  # BH correction across all 15 pairs for this metric
  if (length(raw_p_vec) > 0) {
    adj_p_vec <- p.adjust(raw_p_vec, method="BH")
    for (k2 in seq_along(pair_ids)) {
      oi <- pair_ids[[k2]][1]; oj <- pair_ids[[k2]][2]
      p_mat[oi, oj] <- adj_p_vec[k2]
      p_mat[oj, oi] <- adj_p_vec[k2]   # symmetric
    }
  }
  diag(p_mat) <- NA
  list(diff=diff_mat, p=p_mat, n_pairs=length(raw_p_vec), raw_p=raw_p_vec)
}

# ── Helper: draw one heatmap panel ───────────────────────────────────────────
draw_pairwise_heatmap <- function(diff_mat, p_mat, metric_label,
                                  alpha=0.1, n_seeds=5) {
  n_ord <- nrow(diff_mat)
  ord_labels <- rownames(diff_mat)

  # Colour scale: diverging, centred at 0
  zmax <- max(abs(diff_mat), na.rm=TRUE)
  if (is.na(zmax) || zmax == 0) zmax <- 1
  n_col  <- 101
  pal    <- colorRampPalette(c("#2166AC","#92C5DE","white","#F4A582","#D73027"))(n_col)

  # Significance stars
  stars_fn <- function(p) {
    ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < alpha, "*", ""))))
  }

  par(mar=c(4,4,3.5,5), bg="white", family="sans")
  image(1:n_ord, 1:n_ord, t(diff_mat[n_ord:1, ]),
        col=pal, zlim=c(-zmax, zmax),
        xaxt="n", yaxt="n",
        xlab="", ylab="")

  # Axes — x on bottom, y on left, both show ordering labels
  axis(1, at=1:n_ord, labels=ord_labels,
       las=2, cex.axis=0.9, tick=FALSE, line=-0.5)
  axis(2, at=1:n_ord, labels=rev(ord_labels),
       las=2, cex.axis=0.9, tick=FALSE, line=-0.5)

  # Cell text: mean diff + stars
  for (i in seq_len(n_ord)) {
    for (j in seq_len(n_ord)) {
      row_ord <- ord_labels[n_ord + 1 - j]   # flipped because image() flips rows
      col_ord <- ord_labels[i]
      d <- diff_mat[row_ord, col_ord]
      p <- p_mat[row_ord, col_ord]
      if (i == j) {
        # Diagonal: shade and label "—"
        rect(i-0.5, j-0.5, i+0.5, j+0.5, col="grey88", border=NA)
        text(i, j, "—", cex=0.85, col="grey50")
      } else {
        txt <- paste0(sprintf("%+.3f", d), stars_fn(p))
        text(i, j, txt, cex=0.72,
             col=ifelse(abs(d) > 0.65*zmax, "white", "#1F3864"), font=1)
      }
    }
  }

  # Grid lines
  abline(h=seq(0.5, n_ord+0.5, 1), col="white", lwd=0.8)
  abline(v=seq(0.5, n_ord+0.5, 1), col="white", lwd=0.8)

  # Title
  title(main=sprintf("%s\n(row \u2212 col mean diff across %d seeds)",
                     metric_label, n_seeds),
        cex.main=0.9, col.main="#1F3864", line=1.5)

  # Colour bar on right side
  par(new=TRUE)
  bar_x <- par("usr")[2] + 0.15*diff(par("usr")[1:2])
  ys <- seq(par("usr")[3], par("usr")[4], length.out=n_col+1)
  for (k in seq_len(n_col))
    rect(bar_x, ys[k], bar_x+0.08*diff(par("usr")[1:2]),
         ys[k+1], col=pal[k], border=NA, xpd=TRUE)
  text(bar_x+0.09*diff(par("usr")[1:2]),
       c(par("usr")[3], mean(par("usr")[3:4]), par("usr")[4]),
       labels=sprintf("%+.3f", c(-zmax, 0, zmax)),
       cex=0.65, adj=0, xpd=TRUE, col="#1F3864")

  # Legend for stars
  mtext(sprintf("* BH-adj p<%.2g  (min p=%.4f with n=%d seeds)",
                alpha, 2*(0.5)^n_seeds, n_seeds),
        side=1, line=2.8, cex=0.62, col="#595959")
}

# ── Generate one PNG per metric ───────────────────────────────────────────────
metric_labels <- c(
  lambda_1000       = "Genomic Inflation (\u03bb\u2081\u2030\u2030\u2030)",
  fp_gw             = "False Positives (GW, count)",
  fp_rate_gw        = "FP Rate (GW, per null SNP)",
  fdp_gw            = "False Discovery Proportion (GW)",
  lambda_null       = "Lambda (null variants only)",
  lambda_1000_ratio = "\u03bb\u2081\u2030\u2030\u2030 Ratio vs QC4",
  recovery_gw_pct   = "GW Recovery % (p<5e-8)",
  recovery_sug_pct  = "Suggestive Recovery % (p<1e-5)",
  beta_spearman     = "Beta Spearman r (vs truth)",
  jaccard_snps      = "Jaccard SNP Overlap",
  jaccard_samples   = "Jaccard Sample Overlap"
)

pairwise_results <- list()

for (m in metrics_primary) {
  cat(sprintf("Computing pairwise comparisons: %s\n", m))
  pw <- pairwise_wilcox(seed_level, m, present)
  pairwise_results[[m]] <- pw

  lbl <- if (m %in% names(metric_labels)) metric_labels[[m]] else m

  png(paste0(opt$out, "_pairwise_", m, ".png"),
      width=900, height=870, res=150)
  draw_pairwise_heatmap(pw$diff, pw$p,
                        metric_label=lbl,
                        alpha=0.1, n_seeds=n_seeds)
  dev.off()
  cat(sprintf("  -> %s_pairwise_%s.png\n", opt$out, m))
}

# ── Combined multi-panel PNG (all metrics in one figure) ─────────────────────
n_metrics <- length(metrics_primary)
if (n_metrics > 0) {
  cat("Generating combined pairwise heatmap figure...\n")

  n_cols_fig <- min(3, n_metrics)
  n_rows_fig <- ceiling(n_metrics / n_cols_fig)

  png(paste0(opt$out, "_pairwise_all.png"),
      width=n_cols_fig*800, height=n_rows_fig*820, res=150)
  par(mfrow=c(n_rows_fig, n_cols_fig), oma=c(1,1,3,1))

  for (m in metrics_primary) {
    pw  <- pairwise_results[[m]]
    lbl <- if (m %in% names(metric_labels)) metric_labels[[m]] else m
    draw_pairwise_heatmap(pw$diff, pw$p,
                          metric_label=lbl,
                          alpha=0.1, n_seeds=n_seeds)
  }

  # Blank panels if odd number of metrics
  if (n_metrics %% n_cols_fig != 0) {
    for (k in seq_len(n_cols_fig - n_metrics %% n_cols_fig))
      plot.new()
  }

  mtext("Pairwise QC Ordering Comparison — Mean Difference across Seeds",
        outer=TRUE, cex=1.1, font=2, col="#1F3864", line=1)
  dev.off()
  cat(sprintf("  -> %s_pairwise_all.png\n", opt$out))
}

# ── Write pairwise summary TSV ────────────────────────────────────────────────
pw_rows <- list()
for (m in names(pairwise_results)) {
  pw <- pairwise_results[[m]]
  n_ord <- length(present)
  for (i in seq_len(n_ord-1)) {
    for (j in (i+1):n_ord) {
      oi <- present[i]; oj <- present[j]
      pw_rows[[length(pw_rows)+1]] <- data.table(
        metric     = m,
        ordering_A = oi,
        ordering_B = oj,
        mean_diff_AminusB = round(pw$diff[oi, oj], 5),
        p_adj_BH   = round(pw$p[oi, oj], 4),
        sig_a0.1   = ifelse(!is.na(pw$p[oi,oj]) & pw$p[oi,oj] < 0.1, "*", "")
      )
    }
  }
}
if (length(pw_rows) > 0) {
  pw_dt <- rbindlist(pw_rows)
  fwrite(pw_dt, paste0(opt$out, "_pairwise_tests.tsv"), sep="\t")
  cat(sprintf("Pairwise test results -> %s_pairwise_tests.tsv\n", opt$out))
}

cat("\nAll outputs written with prefix:", opt$out, "\n")
