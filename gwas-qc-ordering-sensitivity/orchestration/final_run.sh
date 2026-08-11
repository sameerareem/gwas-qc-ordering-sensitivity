#!/bin/bash
# =============================================================================
#  final_run.sh  —  Master QC orchestration script
#  Runs all 6 QC orderings across all seeds and datasets.
#
#  Key design decisions:
#    - 3 QC pipelines run in PARALLEL per dataset (qc1+qc2+qc3, then qc4+qc5+qc6)
#    - Each pipeline gets its own log file to avoid interleaved output
#    - cd is done inside a subshell to avoid polluting the master working dir
#    - set -e is intentionally NOT used at master level to allow failure capture
#    - Timing recorded per pipeline for bottleneck analysis
# =============================================================================

BASE_DIR="/mnt/hdd/GWAS_synth_data/gwas_synthetic_output"
PIPELINE_DIR="$BASE_DIR"   # qc scripts live in the same directory as output

# One lock file per dataset — serialises SAIGE step1 across the 3 parallel
# QC scripts. The GLMM PCG solver is memory-bandwidth-bound; 3 simultaneous
# fits at N=66k cause cache thrashing and each takes 3-5x longer.
# Running them one at a time with all 39 cores each takes ~10-15 min total
# (3 x ~4 min) vs ~35 min each if run simultaneously.
# Lock file serialises SAIGE step1 within each parallel group.
# Only one null model fits at a time — prevents memory bandwidth
# contention that makes each step1 take 34 min instead of ~5 min.
SAIGE_LOCK=""   # set per-dataset inside the main loop

SEEDS=(seed_42 seed_123 seed_456 seed_789 seed_999)
DATASETS=(dataset_001 dataset_002 dataset_003)

MASTER_LOG="$BASE_DIR/master_run_log.txt"
echo "Master QC run started: $(date)" > "$MASTER_LOG"
echo "========================================" >> "$MASTER_LOG"

# ---------------------------------------------------------------------------
# run_pipeline <pipeline_path> <dataset_dir> <pheno_file> <log_file>
#
# Runs a single QC pipeline inside a subshell so cd does not affect the
# parent script's working directory. Returns exit code of the pipeline.
# ---------------------------------------------------------------------------
run_pipeline() {
    local PIPELINE_PATH="$1"
    local DATASET_DIR="$2"
    local PHENO="$3"
    local LOG="$4"

    # Validate script exists and is readable before launching
    if [ ! -f "$PIPELINE_PATH" ]; then
        echo "[ERROR] Pipeline script not found: $PIPELINE_PATH" >> "$LOG"
        echo "[ERROR] PIPELINE_DIR=$PIPELINE_DIR — check this path is correct" >> "$LOG"
        return 127
    fi
    if [ ! -f "${DATASET_DIR}/genotypes.bed" ]; then
        echo "[ERROR] genotypes.bed not found in $DATASET_DIR" >> "$LOG"
        return 1
    fi

    local LOCK="${5:-}"

    # Subshell: cd is scoped here, does not affect the outer script
    (
        cd "$DATASET_DIR"
        bash "$PIPELINE_PATH" "genotypes" "$PHENO" "$LOCK"
    ) >> "$LOG" 2>&1
    return $?
}

# ---------------------------------------------------------------------------
# run_group <dataset_dir> <pheno_file> <log_dir> <qc_list...>
#
# Runs up to 3 QC pipelines in parallel and waits for all to finish.
# Captures exit codes individually so failures are reported correctly.
# ---------------------------------------------------------------------------
run_group() {
    local DATASET_DIR="$1"
    local PHENO="$2"
    local LOG_DIR="$3"
    local LOCK="$4"
    shift 4
    local QCS=("$@")

    local PIDS=()
    local NAMES=()

    for QC in "${QCS[@]}"; do
        local PIPELINE_PATH="$PIPELINE_DIR/$QC"
        local QC_LOG="$LOG_DIR/${QC%.sh}.log"
        echo "  [LAUNCH] $QC → $QC_LOG" | tee -a "$MASTER_LOG"
        # Write header to individual log for traceability
        {
            echo "=== $QC | $(date) ==="
            echo "Script:   $PIPELINE_PATH"
            echo "Data dir: $DATASET_DIR"
            echo "Pheno:    $PHENO"
            echo "=========================================="
        } > "$QC_LOG"
        run_pipeline "$PIPELINE_PATH" "$DATASET_DIR" "$PHENO" "$QC_LOG" "$LOCK" &
        PIDS+=($!)
        NAMES+=("$QC")
    done

    # Wait for each background job and collect exit codes
    local ALL_OK=true
    for i in "${!PIDS[@]}"; do
        local PID="${PIDS[$i]}"
        local NAME="${NAMES[$i]}"
        if wait "$PID"; then
            echo "  [OK]     $NAME finished" | tee -a "$MASTER_LOG"
        else
            echo "  [FAIL]   $NAME failed — see $LOG_DIR/${NAME%.sh}.log" | tee -a "$MASTER_LOG"
            ALL_OK=false
        fi
    done

    $ALL_OK
}

# ---------------------------------------------------------------------------
# Main loop: seeds × datasets
# ---------------------------------------------------------------------------
for SEED in "${SEEDS[@]}"; do
    for DATASET in "${DATASETS[@]}"; do

        DATASET_DIR="$BASE_DIR/$SEED/replicates/$DATASET"
        PHENO="$DATASET_DIR/phenotypes_updated.tsv"
        LOG_DIR="$DATASET_DIR/qc_logs"
        mkdir -p "$LOG_DIR"

        echo "" | tee -a "$MASTER_LOG"
        echo "========================================"  | tee -a "$MASTER_LOG"
        echo "Processing: $SEED / $DATASET"             | tee -a "$MASTER_LOG"
        echo "Started:    $(date)"                      | tee -a "$MASTER_LOG"
        echo "========================================"  | tee -a "$MASTER_LOG"

        # Validate required files
        if [ ! -f "${DATASET_DIR}/genotypes.bed" ] || \
           [ ! -f "${DATASET_DIR}/genotypes.bim" ] || \
           [ ! -f "${DATASET_DIR}/genotypes.fam" ]; then
            echo "WARNING: Genotype files missing for $SEED/$DATASET — skipping." | tee -a "$MASTER_LOG"
            continue
        fi

        if [ ! -f "$PHENO" ]; then
            echo "WARNING: Phenotype file missing for $SEED/$DATASET — skipping." | tee -a "$MASTER_LOG"
            continue
        fi

        # Symlink R helper scripts into dataset dir (qc scripts expect them in cwd)
        ln -sf "$PIPELINE_DIR/qc_hetfilter.R" "$DATASET_DIR/qc_hetfilter.R"
        ln -sf "$PIPELINE_DIR/plot_qq.R"      "$DATASET_DIR/plot_qq.R"
        ln -sf "$PIPELINE_DIR/compare_lambda_across_qc.R" "$DATASET_DIR/compare_lambda_across_qc.R"

        DATASET_START=$(date +%s)
        DATASET_OK=true

        # Fresh lock file — serialises SAIGE step1 within this dataset group
        SAIGE_LOCK="$LOG_DIR/.saige_step1.lock"
        rm -f "$SAIGE_LOCK"

        # ── Group 1: qc1, qc2, qc3 in parallel ──────────────────────────
        echo "" | tee -a "$MASTER_LOG"
        echo "  [GROUP 1] qc1 qc2 qc3 — launching in parallel" | tee -a "$MASTER_LOG"
        G1_START=$(date +%s)

        if ! run_group "$DATASET_DIR" "$PHENO" "$LOG_DIR" "$SAIGE_LOCK" qc1.sh qc2.sh qc3.sh; then
            echo "  [WARN] One or more of qc1/qc2/qc3 failed — continuing to group 2" | tee -a "$MASTER_LOG"
            DATASET_OK=false
        fi

        G1_END=$(date +%s)
        echo "  [GROUP 1] done in $(( G1_END - G1_START ))s" | tee -a "$MASTER_LOG"

        # ── Group 2: qc4, qc5, qc6 in parallel ──────────────────────────
        echo "" | tee -a "$MASTER_LOG"
        echo "  [GROUP 2] qc4 qc5 qc6 — launching in parallel" | tee -a "$MASTER_LOG"
        G2_START=$(date +%s)

        if ! run_group "$DATASET_DIR" "$PHENO" "$LOG_DIR" "$SAIGE_LOCK" qc4.sh qc5.sh qc6.sh; then
            echo "  [WARN] One or more of qc4/qc5/qc6 failed" | tee -a "$MASTER_LOG"
            DATASET_OK=false
        fi

        G2_END=$(date +%s)
        echo "  [GROUP 2] done in $(( G2_END - G2_START ))s" | tee -a "$MASTER_LOG"

        # ── Compute benchmarking metrics after all 6 QC orderings done ───
        echo "" | tee -a "$MASTER_LOG"
        echo "  [METRICS] Computing QC comparison metrics..." | tee -a "$MASTER_LOG"
        (
            cd "$DATASET_DIR"
            Rscript "$PIPELINE_DIR/compare_lambda_across_qc.R" \
                --results-dir . \
                --baseline-qc qc1 \
                --out "qc_comparison_${SEED}_${DATASET}"
        ) >> "$LOG_DIR/metrics.log" 2>&1 \
        && echo "  [METRICS] Done → qc_comparison_${SEED}_${DATASET}_lambda_comparison.tsv" \
                | tee -a "$MASTER_LOG" \
        || echo "  [METRICS] Warning: metrics script failed — see $LOG_DIR/metrics.log" \
                | tee -a "$MASTER_LOG"

        # ── Cleanup lock file and symlinks ───────────────────────────────
        rm -f "$SAIGE_LOCK"
        rm -f "$DATASET_DIR/qc_hetfilter.R" \
              "$DATASET_DIR/plot_qq.R" \
              "$DATASET_DIR/compare_lambda_across_qc.R"

        DATASET_END=$(date +%s)
        DATASET_ELAPSED=$(( DATASET_END - DATASET_START ))

        if $DATASET_OK; then
            echo "" | tee -a "$MASTER_LOG"
            echo "  [DONE] $SEED/$DATASET completed in ${DATASET_ELAPSED}s ($(( DATASET_ELAPSED/60 ))m $(( DATASET_ELAPSED%60 ))s)" | tee -a "$MASTER_LOG"
        else
            echo "" | tee -a "$MASTER_LOG"
            echo "  [PARTIAL] $SEED/$DATASET finished with errors in ${DATASET_ELAPSED}s — check $LOG_DIR" | tee -a "$MASTER_LOG"
        fi

    done
done

echo "" | tee -a "$MASTER_LOG"
echo "========================================" >> "$MASTER_LOG"
echo "Master QC run finished: $(date)" | tee -a "$MASTER_LOG"
