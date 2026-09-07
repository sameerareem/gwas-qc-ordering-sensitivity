#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <plink_prefix> <phenotype_file> [saige_lock]"
    exit 1
fi

DATA=$1
PHENOTYPE=$2
THREADS=13
export OMP_NUM_THREADS=$THREADS

# Lock file path for SAIGE step1 serialisation (passed by final_run.sh)
SAIGE_LOCK="${3:-}"

# SAIGE Docker image — using absolute host paths mounted at /mnt/hdd
# This eliminates /data/ path translation bugs entirely
SAIGE_IMG="wzhou88/saige:1.5.1"
SAIGE_RSCRIPT="/app/.pixi/envs/default/bin/Rscript"
SAIGE_STEP1="/usr/local/bin/step1_fitNULLGLMM.R"
SAIGE_STEP2="/usr/local/bin/step2_SPAtests.R"


QCLOG=${DATA}_qc3_log.txt

echo "GWAS QC tracking log (QC3)" > $QCLOG
echo "Dataset: $DATA" >> $QCLOG
echo "Date: $(date)" >> $QCLOG
echo "Threads used: $THREADS" >> $QCLOG
echo "========================================" >> $QCLOG


log_samples () {
    local STEP=$1
    local FAM=$2

    echo "" >> $QCLOG
    echo "[SAMPLES] STEP: $STEP" >> $QCLOG

    cut -d' ' -f1,2 $FAM | sort --parallel=$THREADS > current_samples_qc3.txt

    if [ -f prev_samples_qc3.txt ]; then
        comm -23 prev_samples_qc3.txt current_samples_qc3.txt > removed_samples_qc3.txt
        echo "Removed samples:" >> $QCLOG
        if [ -s removed_samples_qc3.txt ]; then
            cat removed_samples_qc3.txt >> $QCLOG
        else
            echo "None" >> $QCLOG
        fi
    else
        echo "Initial sample list" >> $QCLOG
    fi

    echo "Samples remaining: $(wc -l < current_samples_qc3.txt)" >> $QCLOG
    mv current_samples_qc3.txt prev_samples_qc3.txt
}

log_snps () {
    local STEP=$1
    local BIM=$2

    echo "" >> $QCLOG
    echo "[SNPS] STEP: $STEP" >> $QCLOG

    cut -f2 $BIM | sort --parallel=$THREADS > current_snps_qc3.txt

    if [ -f prev_snps_qc3.txt ]; then
        comm -23 prev_snps_qc3.txt current_snps_qc3.txt > removed_snps_qc3.txt
        echo "Removed SNPs:" >> $QCLOG
        if [ -s removed_snps_qc3.txt ]; then
            cat removed_snps_qc3.txt >> $QCLOG
        else
            echo "None" >> $QCLOG
        fi
    else
        echo "Initial SNP list" >> $QCLOG
    fi

    echo "SNPs remaining: $(wc -l < current_snps_qc3.txt)" >> $QCLOG
    mv current_snps_qc3.txt prev_snps_qc3.txt
}

echo "Running GWAS QC3 pipeline on dataset: $DATA"
DATASET_DIR="$(pwd)"   # capture before any subshell changes context

# Copy genotypes to /tmp (RAM-backed tmpfs) for fast I/O
# Eliminates HDD read latency during SAIGE step2, saving ~35 min
TMPDIR_QC=$(mktemp -d /tmp/qc3_XXXXXX)
TMP_DATA=$TMPDIR_QC/$(basename $DATA)
cp ${DATA}.bed ${TMP_DATA}.bed
cp ${DATA}.bim ${TMP_DATA}.bim
cp ${DATA}.fam ${TMP_DATA}.fam
echo "[INFO] Genotypes copied to $TMPDIR_QC for RAM-backed I/O"
trap "rm -rf $TMPDIR_QC; echo [INFO] Cleaned up $TMPDIR_QC" EXIT


# INITIAL COUNTS
echo "Initial counts:"
echo "Samples: $(wc -l < ${DATA}.fam)"
echo "SNPs: $(wc -l < ${DATA}.bim)"

log_samples "Raw dataset" ${DATA}.fam
log_snps "Raw dataset" ${DATA}.bim


# STEP 1: Split
plink --bfile $TMP_DATA --autosome --allow-extra-chr --threads $THREADS --make-bed --out ${DATA}_qc3_split_autosomal
plink --bfile $TMP_DATA --chr 23-26 --allow-extra-chr --threads $THREADS --make-bed --out ${DATA}_qc3_split_sexchr


# STEP 2: SNP missingness
plink --bfile $TMP_DATA --geno 0.05 --allow-extra-chr --threads $THREADS --make-bed --out ${DATA}_qc3_step1_genoclean
log_snps "After SNP missingness filtering" ${DATA}_qc3_step1_genoclean.bim

echo "Post SNP-missingness:"
echo "Samples: $(wc -l < ${DATA}_qc3_step1_genoclean.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_step1_genoclean.bim)"


# STEP 3: Sample missingness
plink --bfile ${DATA}_qc3_split_autosomal \
      --mind 0.05 \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_autosomal_mindcheck

cut -d' ' -f1,2 ${DATA}_qc3_split_autosomal.fam | sort --parallel=$THREADS > mind_input_samples_qc3.txt
cut -d' ' -f1,2 ${DATA}_qc3_autosomal_mindcheck.fam | sort --parallel=$THREADS > mind_pass_samples_qc3.txt
comm -23 mind_input_samples_qc3.txt mind_pass_samples_qc3.txt > mind_fail_qc3.txt

plink --bfile ${DATA}_qc3_step1_genoclean \
      --remove mind_fail_qc3.txt \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_step2_mindclean

log_samples "After sample missingness removal" ${DATA}_qc3_step2_mindclean.fam

echo "Post Sample-missingness:"
echo "Samples: $(wc -l < ${DATA}_qc3_step2_mindclean.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_step2_mindclean.bim)"


# STEP 4: Sex check
plink --bfile ${DATA}_qc3_split_sexchr \
      --keep ${DATA}_qc3_step2_mindclean.fam \
      --check-sex \
      --set-hh-missing \
      --allow-extra-chr \
      --threads $THREADS \
      --out ${DATA}_qc3_sexcheck

awk 'NR>1 && $5=="PROBLEM" && $4!=0 {print $1"\t"$2}' ${DATA}_qc3_sexcheck.sexcheck > sex_fail_qc3.txt

plink --bfile ${DATA}_qc3_step2_mindclean \
      --remove sex_fail_qc3.txt \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_step3_sexclean

log_samples "After sex check removal" ${DATA}_qc3_step3_sexclean.fam

echo "Post Sex-check:"
echo "Samples: $(wc -l < ${DATA}_qc3_step3_sexclean.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_step3_sexclean.bim)"


# STEP 5: LD pruning
plink --bfile ${DATA}_qc3_step3_sexclean \
      --indep-pairwise 200 50 0.1 \
      --allow-extra-chr \
      --threads $THREADS \
      --out ${DATA}_qc3_pruned_het


# STEP 6: Heterozygosity
plink --bfile ${DATA}_qc3_step3_sexclean \
      --extract ${DATA}_qc3_pruned_het.prune.in \
      --het \
      --allow-extra-chr \
      --threads $THREADS \
      --out ${DATA}_qc3_het

Rscript qc_hetfilter.R ${DATA}_qc3_het.het $THREADS qc_test_het_fail_qc3.txt


# STEP 7: Remove het outliers
plink --bfile ${DATA}_qc3_step3_sexclean \
      --remove qc_test_het_fail_qc3.txt \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_step4_hetclean

log_samples "After heterozygosity filtering" ${DATA}_qc3_step4_hetclean.fam

echo "Post Heterozygosity:"
echo "Samples: $(wc -l < ${DATA}_qc3_step4_hetclean.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_step4_hetclean.bim)"


# STEP 8: KING
plink2 --bfile ${DATA}_qc3_step4_hetclean \
       --autosome \
       --extract ${DATA}_qc3_pruned_het.prune.in \
       --king-cutoff 0.125 \
       --allow-extra-chr \
       --threads $THREADS \
       --out ${DATA}_qc3_king

plink --bfile ${DATA}_qc3_step4_hetclean \
      --keep ${DATA}_qc3_king.king.cutoff.in.id \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_step5_kingclean

log_samples "After relatedness filtering" ${DATA}_qc3_step5_kingclean.fam

echo "Post Relatedness:"
echo "Samples: $(wc -l < ${DATA}_qc3_step5_kingclean.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_step5_kingclean.bim)"


# STEP 9: MAF
plink --bfile ${DATA}_qc3_step5_kingclean --maf 0.01 --allow-extra-chr --threads $THREADS --make-bed --out ${DATA}_qc3_step6_maf
log_snps "After MAF filtering" ${DATA}_qc3_step6_maf.bim

echo "Post MAF:"
echo "Samples: $(wc -l < ${DATA}_qc3_step6_maf.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_step6_maf.bim)"


# STEP 10: HWE
plink --bfile ${DATA}_qc3_step6_maf --hwe 1e-6 --allow-extra-chr --threads $THREADS --make-bed --out ${DATA}_qc3_final
log_snps "After HWE filtering" ${DATA}_qc3_final.bim

echo "Post HWE:"
echo "Samples: $(wc -l < ${DATA}_qc3_final.fam)"
echo "SNPs: $(wc -l < ${DATA}_qc3_final.bim)"

plink --bfile ${DATA}_qc3_final \
      --autosome \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_final_autosomes

log_snps "Autosomal SNPs for GWAS" ${DATA}_qc3_final_autosomes.bim

# Create LD-pruned plink file for SAIGE step1 (GRM/null model)
# Must be DIFFERENT from step2 bedFile to prevent SAIGE re-running step1
plink --bfile ${DATA}_qc3_final_autosomes \
      --extract ${DATA}_qc3_pruned_het.prune.in \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc3_step1_plink



# STEP 11: SAIGE - Null model
# Serialised via flock: only one null model fits at a time.
# The GLMM PCG solver is memory-bandwidth-bound; 3 simultaneous
# fits at N=66k cause cache thrashing (each 3-5x slower).
# While holding the lock this script uses all 39 cores.
STEP1_THREADS=$(( THREADS * 3 ))
run_step1() {
    # Mount /mnt/hdd directly — no path translation, absolute paths throughout
    docker run --rm \
        -v /mnt/hdd:/mnt/hdd \
        --cpus=$STEP1_THREADS \
        "$SAIGE_IMG" \
        "$SAIGE_RSCRIPT" "$SAIGE_STEP1" \
            --plinkFile="$(pwd)/${DATA}_qc3_final_autosomes" \
            --phenoFile="$PHENOTYPE" \
            --phenoCol=PHENO1 \
            --covarColList=sex,age \
            --sampleIDColinphenoFile=IID \
            --traitType=quantitative \
            --invNormalize=TRUE \
            --ratioCVcutoff=0.01 \
            --LOCO=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE \
            --outputPrefix="$(pwd)/${DATA}_qc3_nullmodel" \
            --nThreads=$STEP1_THREADS
}

# Verify final_autosomes exists before attempting step1
if [ ! -f "${DATA}_qc3_final_autosomes.bed" ]; then
    echo "[ERROR] ${DATA}_qc3_final_autosomes.bed not found — aborting step1"
    exit 1
fi

if [ -n "$SAIGE_LOCK" ]; then
    echo "[$(date +%T)] Waiting for SAIGE step1 lock (qc3)..."
    (
        flock -x 200
        echo "[$(date +%T)] Lock acquired — running SAIGE step1 (qc3)"
        run_step1
    ) 200>"$SAIGE_LOCK"
else
    run_step1
fi



# STEP 12: SAIGE - Association test
# Mount /mnt/hdd directly — absolute paths, no translation
docker run --rm \
    -v /mnt/hdd:/mnt/hdd \
    --cpus=$THREADS \
    "$SAIGE_IMG" \
    "$SAIGE_RSCRIPT" "$SAIGE_STEP2" \
        --bedFile="$(pwd)/${DATA}_qc3_final_autosomes.bed" \
        --bimFile="$(pwd)/${DATA}_qc3_final_autosomes.bim" \
        --famFile="$(pwd)/${DATA}_qc3_final_autosomes.fam" \
        --SAIGEOutputFile="$(pwd)/${DATA}_qc3_saige_result.txt" \
        --minMAF=0.01 \
        --minMAC=20 \
        --markers_per_chunk=2000 \
        --LOCO=FALSE \
        --is_noadjCov=FALSE \
        --GMMATmodelFile="$(pwd)/${DATA}_qc3_nullmodel.rda" \
        --varianceRatioFile="$(pwd)/${DATA}_qc3_nullmodel.varianceRatio.txt" \
        --nThreads=1

# STEP 13: QQ Plot
Rscript plot_qq.R ${DATA}_qc3_saige_result.txt ${DATA}_qc3


# ==========================
# FINAL CLEANUP
# ==========================
rm -f ${DATA}_qc3_split_autosomal.*
rm -f ${DATA}_qc3_split_sexchr.*
rm -f ${DATA}_qc3_step1_genoclean.*
rm -f ${DATA}_qc3_autosomal_mindcheck.*
rm -f ${DATA}_qc3_step2_mindclean.*
rm -f ${DATA}_qc3_step3_sexclean.*
rm -f ${DATA}_qc3_step4_hetclean.*
rm -f ${DATA}_qc3_step5_kingclean.*
rm -f ${DATA}_qc3_step6_maf.*
rm -f ${DATA}_qc3_pruned_het.*
rm -f ${DATA}_qc3_step1_plink.*
rm -f ${DATA}_qc3_het.*
rm -f ${DATA}_qc3_king.*
rm -f ${DATA}_qc3_final_autosomes.*

rm -f mind_fail_qc3.txt sex_fail_qc3.txt qc_test_het_fail_qc3.txt
rm -f prev_samples_qc3.txt prev_snps_qc3.txt current_samples_qc3.txt current_snps_qc3.txt removed_samples_qc3.txt removed_snps_qc3.txt
rm -f mind_input_samples_qc3.txt mind_pass_samples_qc3.txt

echo "========================================"
echo "QC3 pipeline completed successfully."
