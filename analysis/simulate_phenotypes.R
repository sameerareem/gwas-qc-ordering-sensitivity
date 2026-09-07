#!/usr/bin/env Rscript
# =============================================================================
#  simulate_phenotypes.R
#  Lung Function Z-Score Phenotype Simulation for GWAS QC Order Benchmarking
#
#  Simulates FEV1, FVC, and FEV1/FVC z-scores as quantitative phenotypes
#  using a biologically realistic model based on GLI-2012 reference equations.
#
#  Phenotype model:
#    z_lung = G*beta + b_age*age_std + b_sex*sex + b_height*height_std
#             + strat_noise + residual
#
#  Covariates output:
#    SEX, AGE (used in GWAS model, not in phenotype generation to avoid
#    conditioning on known covariates during simulation)
#
#  Called from gwas_pipeline_master.py via:
#    Rscript simulate_phenotypes.R --config pheno_config.json
#
#  Required R packages:
#    data.table, optparse, jsonlite
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(jsonlite)
})

# ─────────────────────────────────────────────────────────────────────────────
#  Parse Arguments
# ─────────────────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--config", type = "character", help = "Path to pheno_config.json")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- fromJSON(opt$config)

set.seed(cfg$seed)
cat("[INFO] Lung function z-score phenotype simulation starting\n")
cat("[INFO] Config:\n"); print(cfg)

# ─────────────────────────────────────────────────────────────────────────────
#  1. Load Genotype Data (.bim for SNP list, .fam for individuals)
# ─────────────────────────────────────────────────────────────────────────────

bim_file  <- paste0(cfg$plink_prefix, ".bim")
fam_file  <- paste0(cfg$plink_prefix, ".fam")
pvar_file <- paste0(cfg$plink_prefix, ".pvar")
psam_file <- paste0(cfg$plink_prefix, ".psam")

if (file.exists(pvar_file)) {
  snp_info <- fread(pvar_file, skip = "#CHROM")
  colnames(snp_info)[1:5] <- c("CHR","POS","ID","REF","ALT")
  ind_info <- fread(psam_file)
  plink_format <- "pgen"
  cat(sprintf("[INFO] Loaded %d SNPs from .pvar\n", nrow(snp_info)))
  cat(sprintf("[INFO] Loaded %d individuals from .psam\n", nrow(ind_info)))
} else if (file.exists(bim_file)) {
  snp_info <- fread(bim_file, col.names = c("CHR","ID","CM","POS","A1","A2"))
  ind_info <- fread(fam_file, col.names = c("FID","IID","PAT","MAT","SEX","PHENO"))
  plink_format <- "bed"
  cat(sprintf("[INFO] Loaded %d SNPs from .bim\n", nrow(snp_info)))
  cat(sprintf("[INFO] Loaded %d individuals from .fam\n", nrow(ind_info)))
} else {
  stop("[ERROR] Cannot find PLINK .pvar or .bim file at: ", cfg$plink_prefix)
}

cat(sprintf("[INFO] PLINK format detected: %s\n", plink_format))

N <- nrow(ind_info)
M <- nrow(snp_info)

# Autosomal SNPs only for causal variant selection
# Exclude X, Y, MT, PAR regions
autosomal_snps <- snp_info[
  !(snp_info$CHR %in% c("X","Y","MT","XY","PAR1","PAR2",23,24,25,26)),
]
cat(sprintf("[INFO] N=%d individuals, M=%d total SNPs, %d autosomal\n",
            N, M, nrow(autosomal_snps)))

# ─────────────────────────────────────────────────────────────────────────────
#  2. Simulate Realistic Covariates
# ─────────────────────────────────────────────────────────────────────────────

# ── Sex ──
if ("SEX" %in% names(ind_info)) {
  sex_vec <- ind_info$SEX
} else if ("#SEX" %in% names(ind_info)) {
  sex_vec <- ind_info$`#SEX`
} else {
  cat("[WARN] Sex not found; assigning randomly.\n")
  sex_vec <- sample(c(1L, 2L), size = N, replace = TRUE)
}
# Recode: 1=male -> 0, 2=female -> 1 (regression coding)
sex_coded <- ifelse(sex_vec == 2, 1L, 0L)

# ── Age ──
# Typical GWAS cohort: 40-70 years, roughly normal centred at 55
age_raw <- round(rnorm(N, mean = 55, sd = 8))
age_raw <- pmax(40L, pmin(75L, as.integer(age_raw)))  # clamp to [40, 75]
age_std <- as.vector(scale(age_raw))  # standardised for regression

# ── Height (cm) — correlated with sex ──
# Males: mean 175cm sd 7, Females: mean 162cm sd 6 (approximate EUR values)
height_raw <- ifelse(
  sex_vec == 1,
  rnorm(N, mean = 175, sd = 7),   # male
  rnorm(N, mean = 162, sd = 6)    # female
)
height_std <- as.vector(scale(height_raw))

cat(sprintf("[INFO] Age: mean=%.1f, sd=%.1f, range=[%d,%d]\n",
            mean(age_raw), sd(age_raw), min(age_raw), max(age_raw)))
cat(sprintf("[INFO] Height: mean=%.1f cm, sd=%.1f\n",
            mean(height_raw), sd(height_raw)))

# ─────────────────────────────────────────────────────────────────────────────
#  3. Select Causal SNPs
# ─────────────────────────────────────────────────────────────────────────────

n_causal  <- cfg$n_causal_snps
h2        <- cfg$heritability

causal_idx   <- sample(seq_len(nrow(autosomal_snps)), size = n_causal, replace = FALSE)
causal_snps  <- autosomal_snps[causal_idx, ]

# Effect sizes from Normal(0, sqrt(h2/n_causal))
# This ensures Var(G*beta) ≈ h2 when genotypes are standardised
effect_sd    <- sqrt(h2 / n_causal)
causal_betas <- rnorm(n_causal, mean = 0, sd = effect_sd)

causal_df <- data.table(
  SNP  = causal_snps$ID,
  CHR  = causal_snps$CHR,
  POS  = causal_snps$POS,
  BETA = causal_betas,
  SE   = abs(causal_betas) * 0.1
)

cat(sprintf("[INFO] Selected %d causal SNPs; h2=%.2f; effect_sd=%.4f\n",
            n_causal, h2, effect_sd))

# ─────────────────────────────────────────────────────────────────────────────
#  4. Extract Causal SNP Dosages
# ─────────────────────────────────────────────────────────────────────────────

tmp_dir <- file.path(cfg$out_dir, "tmp_plink_extract")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

causal_snp_list <- file.path(tmp_dir, "causal_snps_list.txt")
fwrite(causal_df[, .(SNP)], causal_snp_list, col.names = FALSE)

dosage_prefix <- file.path(tmp_dir, "causal_dosages")
plink1 <- "plink"

if (plink_format == "bed") {
  plink_args <- c(
    "--bfile",   cfg$plink_prefix,
    "--extract", causal_snp_list,
    "--export",  "A",
    "--out",     dosage_prefix,
    "--allow-no-sex",
    "--allow-extra-chr"
  )
} else {
  plink_args <- c(
    "--pfile",   cfg$plink_prefix,
    "--extract", causal_snp_list,
    "--export",  "A",
    "--out",     dosage_prefix
  )
  plink1 <- "plink2"
}

cat(sprintf("[INFO] Extracting causal SNP dosages with %s ...\n", plink1))
system2(plink1, args = plink_args, stdout = TRUE, stderr = TRUE)

raw_file <- paste0(dosage_prefix, ".raw")
if (!file.exists(raw_file)) {
  stop("[ERROR] PLINK dosage extraction failed — raw file not found.")
}

raw_data  <- fread(raw_file)
iids      <- raw_data$IID
geno_cols <- grep("_[A-Z]$", names(raw_data), value = TRUE)
G         <- as.matrix(raw_data[, ..geno_cols])

# Standardise genotypes
G_std            <- scale(G)
G_std[is.na(G_std)] <- 0

cat(sprintf("[INFO] Dosage matrix: %d individuals x %d causal SNPs\n",
            nrow(G_std), ncol(G_std)))

# ─────────────────────────────────────────────────────────────────────────────
#  5. Simulate Lung Function Z-Score
#
#  Model (based on GLI-2012 reference equation structure):
#
#    z_lung = G_genetic + b_age * age_std + b_sex * sex_coded
#             + b_height * height_std + strat_noise + residual
#
#  Effect sizes reflect published literature:
#    Age:    stronger negative effect (lung function declines with age)
#    Sex:    males have higher FEV1/FVC in absolute terms
#            but we work in z-scores so sex effect is moderate
#    Height: strong positive effect on lung volumes
#
#  The phenotype_subtype in config selects which measure to simulate:
#    "FEV1"      — forced expiratory volume in 1 second
#    "FVC"       — forced vital capacity
#    "FEV1_FVC"  — FEV1/FVC ratio (less height-dependent)
# ─────────────────────────────────────────────────────────────────────────────

pheno_subtype <- cfg$pheno_subtype

# Published approximate effect sizes for lung function z-scores
# (standardised betas from large GWAS meta-analyses)
if (pheno_subtype == "FEV1_FVC") {
  # FEV1/FVC ratio: age effect stronger, height effect weaker
  b_age    <- -0.35   # strong negative: ratio declines with age
  b_sex    <- -0.15   # males have slightly lower FEV1/FVC ratio
  b_height <- -0.05   # ratio is relatively height-independent
} else if (pheno_subtype == "FVC") {
  # FVC: height effect strongest
  b_age    <- -0.25
  b_sex    <- -0.20   # females have lower FVC (negative = female coded as 1)
  b_height <-  0.45   # strong positive height effect
} else {
  # FEV1 (default)
  b_age    <- -0.30
  b_sex    <- -0.18
  b_height <-  0.40
}

cat(sprintf("[INFO] Simulating %s z-score\n", pheno_subtype))
cat(sprintf("[INFO] Effect sizes: age=%.2f, sex=%.2f, height=%.2f\n",
            b_age, b_sex, b_height))

# ── Genetic component ──
genetic_component <- as.vector(G_std %*% causal_betas)
genetic_var       <- var(genetic_component)

if (genetic_var > 0) {
  G_scaled <- genetic_component * sqrt(h2 / genetic_var)
} else {
  G_scaled <- genetic_component
}

# ── Environmental / covariate components ──
age_component    <- b_age    * age_std
sex_component    <- b_sex    * sex_coded
height_component <- b_height * height_std

# ── Population stratification noise ──
strat_noise <- 0.0
if (cfg$strat_variance > 0) {
  n_pops          <- 3
  pop_memberships <- sample(1:n_pops, size = N, replace = TRUE,
                             prob = c(0.60, 0.25, 0.15))
  pop_effects     <- rnorm(n_pops, mean = 0, sd = sqrt(cfg$strat_variance))
  strat_noise     <- as.vector(pop_effects[pop_memberships])
}

# ── Residual noise ──
# Residual variance = 1 - h2 - variance explained by covariates
# We approximate covariate variance and adjust residual accordingly
cov_var      <- var(age_component + sex_component + height_component)
residual_var <- max(0.05, 1 - h2 - cov_var)
noise        <- rnorm(N, mean = 0, sd = sqrt(residual_var))

# ── Compose phenotype ──
pheno_continuous <- G_scaled + age_component + sex_component +
                    height_component + strat_noise + noise

# ── Standardise to z-score (mean=0, sd=1) ──
pheno_z <- as.vector(scale(pheno_continuous))

cat(sprintf("[INFO] %s z-score: mean=%.4f, sd=%.4f, range=[%.2f, %.2f]\n",
            pheno_subtype, mean(pheno_z), sd(pheno_z),
            min(pheno_z), max(pheno_z)))

# Binary phenotype via liability threshold if requested
if (cfg$pheno_type == "binary") {
  threshold <- qnorm(1 - cfg$prevalence)
  pheno_out <- ifelse(pheno_z > threshold, 2L, 1L)
  cat(sprintf("[INFO] Binary phenotype: prevalence=%.2f, cases=%d, controls=%d\n",
              cfg$prevalence, sum(pheno_out == 2L), sum(pheno_out == 1L)))
} else {
  pheno_out <- pheno_z
  cat(sprintf("[INFO] Quantitative %s z-score: mean=%.4f, sd=%.4f\n",
              pheno_subtype, mean(pheno_out), sd(pheno_out)))
}

# ─────────────────────────────────────────────────────────────────────────────
#  6. Write Output Files
# ─────────────────────────────────────────────────────────────────────────────

out_dir <- cfg$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 6a. Phenotype file ──
pheno_dt <- data.table(
  `#FID` = iids,
  IID    = iids,
  PHENO1 = pheno_out
)
pheno_file <- file.path(out_dir, "phenotypes.tsv")
fwrite(pheno_dt, pheno_file, sep = "\t")
cat(sprintf("[INFO] Phenotype file -> %s\n", pheno_file))

# ── 6b. Causal SNP ground-truth file ──
causal_out_file <- file.path(out_dir, "causal_snps.tsv")
fwrite(causal_df, causal_out_file, sep = "\t")
cat(sprintf("[INFO] Causal SNP list -> %s\n", causal_out_file))

# ── 6c. Covariate file — SEX and AGE ──
# PCA not included here — compute post-QC using plink2 --pca approx 10
covar_dt <- data.table(
  `#FID` = iids,
  IID    = iids,
  SEX    = sex_coded,   # 0=male, 1=female (regression coding)
  AGE    = age_raw      # raw age in years (standardise in GWAS model)
)
covar_file <- file.path(out_dir, "covariates.tsv")
fwrite(covar_dt, covar_file, sep = "\t")
cat(sprintf("[INFO] Covariate file (SEX + AGE) -> %s\n", covar_file))

# ── 6d. Simulation summary ──
summary_dt <- data.table(
  parameter = c("N_individuals","N_total_SNPs","N_causal_SNPs",
                "heritability","phenotype_type","phenotype_subtype",
                "sex_effect_size","age_effect_size","height_effect_size",
                "strat_variance","seed"),
  value     = c(N, M, n_causal, h2, cfg$pheno_type, pheno_subtype,
                b_sex, b_age, b_height,
                cfg$strat_variance, cfg$seed)
)
fwrite(summary_dt, file.path(out_dir, "simulation_summary.tsv"), sep = "\t")

cat("\n[INFO] Lung function phenotype simulation complete.\n")
cat(sprintf("[INFO]   Phenotype (%s z-score): %s\n", pheno_subtype, pheno_file))
cat(sprintf("[INFO]   Causal SNPs:            %s\n", causal_out_file))
cat(sprintf("[INFO]   Covariates (SEX, AGE):  %s\n", covar_file))
