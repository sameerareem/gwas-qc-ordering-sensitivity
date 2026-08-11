# Investigating the Sensitivity of GWAS Outcomes to the Order of Data Pre-Processing Steps

**Bachelor of Technology in Biotechnology**  
Sameera Reem Imthiaz (22BBT0122)  
Department of Biotechnology, School of Bio Sciences and Technology  
Vellore Institute of Technology — April 2026

---

## Overview

Genome-wide association studies (GWAS) rely on multi-step quality control (QC) pipelines before association testing. While the individual steps are well established, whether their *sequential order* affects downstream results has received little systematic attention.

This project investigated whether applying the same QC filters in six different sequences produces meaningfully different GWAS outcomes when all other variables are held constant.

**Key finding:** For large genotype datasets with standard QC thresholds, the order in which QC steps are applied is unlikely to be a meaningful source of variability in GWAS results. The dominant source of variation was the simulated dataset, not the pipeline sequence.

---

## Study Design

### Dataset
- 75,000 individuals, ~110,000 SNPs
- Simulated using **msprime** under a European demographic model (`OutOfAfrica_3G09`)
- Quantitative lung function phenotype (FEV1 z-score), SNP-based heritability h² = 0.5, 50 causal variants
- 5 independent replicates (seeds: 42, 123, 456, 789, 999), 3 datasets per seed

### QC Filters Applied (identical thresholds across all pipelines)
| Filter | Threshold |
|--------|-----------|
| Sample missingness (`--mind`) | 5% |
| SNP missingness (`--geno`) | 5% |
| Sex discrepancy check | PLINK `--check-sex` |
| Minor allele frequency (`--maf`) | 1% |
| Hardy-Weinberg equilibrium (`--hwe`) | 1×10⁻⁶ |
| Heterozygosity outliers | ±3 SD from mean F |
| Relatedness (`--king-cutoff`) | 0.125 (KING) |

### The Six QC Orderings
| Pipeline | Step Order |
|----------|-----------|
| QC1 | Mind → Sex → Het → KING → Geno → MAF → HWE |
| QC2 | Geno → MAF → HWE → Mind → Sex → Het → KING |
| QC3 | Geno → Mind → Sex → Het → KING → MAF → HWE |
| QC4 | Geno → Mind → Sex → MAF → HWE → Het → KING |
| QC5 | Sex → Geno → Mind → MAF → HWE → Het → KING |
| QC6 | Sex → Geno → Het → Mind → KING → MAF → HWE |

### Association Testing
- **SAIGE** v1.5.1 (via Docker: `wzhou88/saige:1.5.1`)
- Quantitative trait, inverse-normal transformation, covariates: sex, age

### Outcome Metrics
- Genomic inflation factor (λ_GC, λ₁₀₀₀)
- Causal SNP recovery at genome-wide significance (p < 5×10⁻⁸)
- Effect size concordance (Lin's CCC)
- Jaccard similarity for retained variants and samples
- False discovery proportion

---

## Results Summary

| Metric | Outcome |
|--------|---------|
| Retained samples across pipelines | Differed by fewer than 10 individuals |
| Retained SNPs across pipelines | Differed by fewer than 3 SNPs |
| λ₁₀₀₀ | 0.997 ± 0.01 across all orderings |
| Causal SNP recovery difference | At most 1 SNP between pipelines |
| Lin's CCC (effect sizes) | 0.996–1.000 |
| Jaccard similarity | >0.997 in all but one comparison |

---

## Repository Structure

```
gwas-qc-ordering-sensitivity/
│
├── README.md
│
├── config.yaml                      # Pipeline configuration (sample size, SNP counts, paths)
│
├── data_generation/
│   ├── gwas_pipeline_master.py      # Master Python script: simulates genotypes using msprime
│   └── verify_replicates.py         # Post-generation QC: verifies all replicates are complete
│
├── qc_pipelines/
│   ├── qc1.sh                       # QC ordering 1: Mind → Sex → Het → KING → Geno → MAF → HWE
│   ├── qc2.sh                       # QC ordering 2: Geno → MAF → HWE → Mind → Sex → Het → KING
│   ├── qc3.sh                       # QC ordering 3: Geno → Mind → Sex → Het → KING → MAF → HWE
│   ├── qc4.sh                       # QC ordering 4: Geno → Mind → Sex → MAF → HWE → Het → KING
│   ├── qc5.sh                       # QC ordering 5: Sex → Geno → Mind → MAF → HWE → Het → KING
│   └── qc6.sh                       # QC ordering 6: Sex → Geno → Het → Mind → KING → MAF → HWE
│
├── orchestration/
│   ├── 00_setup_environment.sh      # Dependency checker and directory initialiser (run first)
│   ├── final_run.sh                 # Master script: runs all 6 QC pipelines across all datasets
│   ├── run_compute_metrics.sh       # Runs compute_qc_metrics.R for all seed × dataset combinations
│   └── run_compare_to_qc4.sh       # Runs compare_to_qc4.R aggregating all results
│
└── analysis/
    ├── qc_hetfilter.R               # Identifies heterozygosity outliers (±3 SD)
    ├── simulate_phenotypes.R        # Simulates quantitative lung function phenotype
    ├── compute_lambda.R             # Computes genomic inflation factors (λ_GC, λ₁₀₀₀)
    ├── compute_qc_metrics.R         # Per-ordering metrics: inflation, recovery, CCC, Jaccard
    ├── compute_qc_metrics_1.R       # Variant of metrics computation
    ├── compare_lambda_across_qc.R   # Compares λ across all 6 QC orderings
    ├── compare_to_qc4.R             # Statistical tests comparing all orderings to QC4 baseline
    └── plot_qq.R                    # Generates QQ plots from SAIGE results
```
---

## Software Versions Used

| Tool | Version |
|------|---------|
| PLINK 1.9 | v1.90b6.26 |
| PLINK2 | v2.00a5 |
| SAIGE | v1.5.1 |
| msprime | v1.x (via stdpopsim `OutOfAfrica_3G09`) |
| R | 4.x |
| Python | 3.x |

---

## Acknowledgements

The msprime data simulation pipeline, SAIGE association testing scripts, and R analysis/visualisation scripts were developed with the assistance of AI (Claude, Anthropic). The six PLINK QC pipeline scripts (`qc1.sh`–`qc6.sh`) were independently developed by the author. All code was verified against expected outputs.

---

## Academic Use

This repository accompanies an undergraduate thesis submitted to the Department of Biotechnology, VIT Vellore (April 2026). 