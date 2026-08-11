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


QCLOG=${DATA}_qc6_log.txt

echo "GWAS QC tracking log (QC6)" > $QCLOG
echo "Dataset: $DATA" >> $QCLOG
echo "Date: $(date)" >> $QCLOG
echo "Threads used: $THREADS" >> $QCLOG
echo "========================================" >> $QCLOG


log_samples () {
    local STEP=$1
    local FAM=$2

    echo "" >> $QCLOG
    echo "[SAMPLES] STEP: $STEP" >> $QCLOG

    cut -d' ' -f1,2 $FAM | sort --parallel=$THREADS > current_samples_qc6.txt

    if [ -f prev_samples_qc6.txt ]; then
        comm -23 prev_samples_qc6.txt current_samples_qc6.txt > removed_samples_qc6.txt
        echo "Removed samples:" >> $QCLOG
        if [ -s removed_samples_qc6.txt ]; then
            cat removed_samples_qc6.txt >> $QCLOG
        else
            echo "None" >> $QCLOG
        fi
    else
        echo "Initial sample list" >> $QCLOG
    fi

    echo "Samples remaining: $(wc -l < current_samples_qc6.txt)" >> $QCLOG
    mv current_samples_qc6.txt prev_samples_qc6.txt
}

log_snps () {
    local STEP=$1
    local BIM=$2

    echo "" >> $QCLOG
    echo "[SNPS] STEP: $STEP" >> $QCLOG

    cut -f2 $BIM | sort --parallel=$THREADS > current_snps_qc6.txt

    if [ -f prev_snps_qc6.txt ]; then
        comm -23 prev_snps_qc6.txt current_snps_qc6.txt > removed_snps_qc6.txt
        echo "Removed SNPs:" >> $QCLOG
        if [ -s removed_snps_qc6.txt ]; then
            cat removed_snps_qc6.txt >> $QCLOG
        else
            echo "None" >> $QCLOG
        fi
    else
        echo "Initial SNP list" >> $QCLOG
    fi

    echo "SNPs remaining: $(wc -l < current_snps_qc6.txt)" >> $QCLOG
    mv current_snps_qc6.txt prev_snps_qc6.txt
}

cleanup () {
    local PREFIX=$1
    rm -f ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam ${PREFIX}.log
}

echo "Running GWAS QC6 pipeline on dataset: $DATA"
DATASET_DIR="$(pwd)"   # capture before any subshell changes context

# Copy genotypes to /tmp (RAM-backed tmpfs) for fast I/O
# Eliminates HDD read latency during SAIGE step2, saving ~35 min
TMPDIR_QC=$(mktemp -d /tmp/qc6_XXXXXX)
TMP_DATA=$TMPDIR_QC/$(basename $DATA)
cp ${DATA}.bed ${TMP_DATA}.bed
cp ${DATA}.bim ${TMP_DATA}.bim
cp ${DATA}.fam ${TMP_DATA}.fam
echo "[INFO] Genotypes copied to $TMPDIR_QC for RAM-backed I/O"
trap "rm -rf $TMPDIR_QC; echo [INFO] Cleaned up $TMPDIR_QC" EXIT


log_samples "Raw dataset" ${DATA}.fam
log_snps "Raw dataset" ${DATA}.bim
echo "[RAW] Samples: $(wc -l < ${DATA}.fam) | SNPs: $(wc -l < ${DATA}.bim)"


# STEP 1: Split
plink --bfile $TMP_DATA \
      --autosome \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_split_autosomal

plink --bfile $TMP_DATA \
      --chr 23-26 \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_split_sexchr

echo "Post split | Samples: $(wc -l < ${DATA}.fam) | SNPs: $(wc -l < ${DATA}.bim)"


# STEP 2: Sex check
plink --bfile ${DATA}_qc6_split_sexchr \
      --check-sex \
      --set-hh-missing \
      --allow-extra-chr \
      --threads $THREADS \
      --out ${DATA}_qc6_sexcheck

awk 'NR>1 && $5=="PROBLEM" && $4!=0 {print $1"\t"$2}' ${DATA}_qc6_sexcheck.sexcheck > sex_fail_qc6.txt

plink --bfile $TMP_DATA \
      --remove sex_fail_qc6.txt \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step1_sexclean

log_samples "After sex check removal (applied to combined)" ${DATA}_qc6_step1_sexclean.fam
echo "Post sex check | Samples: $(wc -l < ${DATA}_qc6_step1_sexclean.fam) | SNPs: $(wc -l < ${DATA}_qc6_step1_sexclean.bim)"


# STEP 3: SNP missingness
plink --bfile ${DATA}_qc6_step1_sexclean \
      --geno 0.05 \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step2

log_snps "After SNP missingness filtering" ${DATA}_qc6_step2.bim
echo "Post SNP missingness | Samples: $(wc -l < ${DATA}_qc6_step2.fam) | SNPs: $(wc -l < ${DATA}_qc6_step2.bim)"


# STEP 4: LD pruning
plink --bfile ${DATA}_qc6_step2 \
      --indep-pairwise 200 50 0.1 \
      --allow-extra-chr \
      --threads $THREADS \
      --out ${DATA}_qc6_pruned_het


# STEP 5: Heterozygosity
plink --bfile ${DATA}_qc6_step2 \
      --extract ${DATA}_qc6_pruned_het.prune.in \
      --het \
      --allow-extra-chr \
      --threads $THREADS \
      --out ${DATA}_qc6_het

Rscript qc_hetfilter.R ${DATA}_qc6_het.het $THREADS qc_test_het_fail_qc6.txt


# STEP 6: Remove heterozygosity outliers
plink --bfile ${DATA}_qc6_step2 \
      --remove qc_test_het_fail_qc6.txt \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step3_hetclean

log_samples "After heterozygosity filtering" ${DATA}_qc6_step3_hetclean.fam
echo "Post heterozygosity | Samples: $(wc -l < ${DATA}_qc6_step3_hetclean.fam) | SNPs: $(wc -l < ${DATA}_qc6_step3_hetclean.bim)"


# STEP 7: Sample missingness
plink --bfile ${DATA}_qc6_split_autosomal \
      --keep ${DATA}_qc6_step3_hetclean.fam \
      --mind 0.05 \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_autosomal_mindcheck

cut -d' ' -f1,2 ${DATA}_qc6_step3_hetclean.fam | sort --parallel=$THREADS > mind_input_samples_qc6.txt
cut -d' ' -f1,2 ${DATA}_qc6_autosomal_mindcheck.fam | sort --parallel=$THREADS > mind_pass_samples_qc6.txt
comm -23 mind_input_samples_qc6.txt mind_pass_samples_qc6.txt > mind_fail_qc6.txt

plink --bfile ${DATA}_qc6_step3_hetclean \
      --remove mind_fail_qc6.txt \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step4_mindclean

log_samples "After sample missingness removal (applied to combined)" ${DATA}_qc6_step4_mindclean.fam
echo "Post sample missingness | Samples: $(wc -l < ${DATA}_qc6_step4_mindclean.fam) | SNPs: $(wc -l < ${DATA}_qc6_step4_mindclean.bim)"


# STEP 8: KING
plink2 --bfile ${DATA}_qc6_step4_mindclean \
       --autosome \
       --extract ${DATA}_qc6_pruned_het.prune.in \
       --allow-extra-chr \
       --king-cutoff 0.125 \
       --threads $THREADS \
       --out ${DATA}_qc6_king

plink --bfile ${DATA}_qc6_step4_mindclean \
      --keep ${DATA}_qc6_king.king.cutoff.in.id \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step5

log_samples "After relatedness filtering" ${DATA}_qc6_step5.fam
echo "Post KING | Samples: $(wc -l < ${DATA}_qc6_step5.fam) | SNPs: $(wc -l < ${DATA}_qc6_step5.bim)"


# STEP 9: MAF
plink --bfile ${DATA}_qc6_step5 \
      --maf 0.01 \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step6

log_snps "After MAF filtering" ${DATA}_qc6_step6.bim
echo "Post MAF | Samples: $(wc -l < ${DATA}_qc6_step6.fam) | SNPs: $(wc -l < ${DATA}_qc6_step6.bim)"


# STEP 10: HWE
plink --bfile ${DATA}_qc6_step6 \
      --hwe 1e-6 \
      --allow-extra-chr \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_final

log_snps "After HWE filtering" ${DATA}_qc6_final.bim
echo "Post HWE | Samples: $(wc -l < ${DATA}_qc6_final.fam) | SNPs: $(wc -l < ${DATA}_qc6_final.bim)"

plink --bfile ${DATA}_qc6_final \
      --autosome \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_final_autosomes

log_snps "Autosomal SNPs for GWAS" ${DATA}_qc6_final_autosomes.bim

# Create LD-pruned plink file for SAIGE step1 (GRM/null model)
# Must be DIFFERENT from step2 bedFile to prevent SAIGE re-running step1
plink --bfile ${DATA}_qc6_final_autosomes \
      --extract ${DATA}_qc6_pruned_het.prune.in \
      --threads $THREADS \
      --make-bed \
      --out ${DATA}_qc6_step1_plink


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
            --plinkFile="$(pwd)/${DATA}_qc6_final_autosomes" \
            --phenoFile="$PHENOTYPE" \
            --phenoCol=PHENO1 \
            --covarColList=sex,age \
            --sampleIDColinphenoFile=IID \
            --traitType=quantitative \
            --invNormalize=TRUE \
            --ratioCVcutoff=0.01 \
            --LOCO=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE \
            --outputPrefix="$(pwd)/${DATA}_qc6_nullmodel" \
            --nThreads=$STEP1_THREADS
}

# Verify final_autosomes exists before attempting step1
if [ ! -f "${DATA}_qc6_final_autosomes.bed" ]; then
    echo "[ERROR] ${DATA}_qc6_final_autosomes.bed not found — aborting step1"
    exit 1
fi

if [ -n "$SAIGE_LOCK" ]; then
    echo "[$(date +%T)] Waiting for SAIGE step1 lock (qc6)..."
    (
        flock -x 200
        echo "[$(date +%T)] Lock acquired — running SAIGE step1 (qc6)"
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
        --bedFile="$(pwd)/${DATA}_qc6_final_autosomes.bed" \
        --bimFile="$(pwd)/${DATA}_qc6_final_autosomes.bim" \
        --famFile="$(pwd)/${DATA}_qc6_final_autosomes.fam" \
        --SAIGEOutputFile="$(pwd)/${DATA}_qc6_saige_result.txt" \
        --minMAF=0.01 \
        --minMAC=20 \
        --markers_per_chunk=2000 \
        --LOCO=FALSE \
        --is_noadjCov=FALSE \
        --GMMATmodelFile="$(pwd)/${DATA}_qc6_nullmodel.rda" \
        --varianceRatioFile="$(pwd)/${DATA}_qc6_nullmodel.varianceRatio.txt" \
        --nThreads=1

# STEP 13: QQ Plot
Rscript plot_qq.R ${DATA}_qc6_saige_result.txt ${DATA}_qc6



# ==========================
# FINAL CLEANUP 
# ==========================
cleanup ${DATA}_qc6_split_autosomal
cleanup ${DATA}_qc6_split_sexchr
cleanup ${DATA}_qc6_step1_sexclean
cleanup ${DATA}_qc6_step2
cleanup ${DATA}_qc6_autosomal_mindcheck
cleanup ${DATA}_qc6_step3_hetclean
cleanup ${DATA}_qc6_step4_mindclean
cleanup ${DATA}_qc6_step5
cleanup ${DATA}_qc6_step6

rm -f ${DATA}_qc6_sexcheck.sexcheck ${DATA}_qc6_sexcheck.log
rm -f ${DATA}_qc6_pruned_het.prune.in ${DATA}_qc6_pruned_het.prune.out ${DATA}_qc6_pruned_het.log
rm -f ${DATA}_qc6_step1_plink.*
rm -f ${DATA}_qc6_het.het ${DATA}_qc6_het.log
rm -f ${DATA}_qc6_king.king.cutoff.in.id ${DATA}_qc6_king.king.cutoff.out.id ${DATA}_qc6_king.log
rm -f ${DATA}_qc6_final_autosomes.*
rm -f mind_fail_qc6.txt mind_input_samples_qc6.txt mind_pass_samples_qc6.txt
rm -f sex_fail_qc6.txt qc_test_het_fail_qc6.txt
rm -f prev_samples_qc6.txt prev_snps_qc6.txt current_samples_qc6.txt current_snps_qc6.txt removed_samples_qc6.txt removed_snps_qc6.txt

echo "========================================"
echo "QC6 pipeline completed successfully."
