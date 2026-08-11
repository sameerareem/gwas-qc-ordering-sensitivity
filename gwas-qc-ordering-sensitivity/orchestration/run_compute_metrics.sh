#!/bin/bash
# =============================================================================
#  run_compute_metrics.sh
#  Run compute_qc_metrics.R for every QC ordering across all seed × dataset
#  combinations that have completed SAIGE results.
#
#  Execution order per dataset:
#    1. QC4 first  — no baseline arguments (QC4 is the reference)
#    2. QC1,2,3,5,6 — QC4 metrics passed as baseline values for ratios
#                     and paired statistical tests
#
#  The script is safe to re-run: it skips any ordering whose metrics TSV
#  already exists and is non-empty (use --force to overwrite).
#
#  Usage:
#    bash run_compute_metrics.sh [--force]
#
#  Outputs per ordering per dataset:
#    genotypes_qcN_metrics.tsv
#    genotypes_qcN_metrics_manhattan.png
#    genotypes_qcN_metrics_qqplot.png
#    genotypes_qcN_metrics_recovery.png
#    genotypes_qcN_metrics_effectsize*.png
#    genotypes_qcN_metrics_chr_lambda.png
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR="/mnt/hdd/GWAS_synth_data/gwas_synthetic_output"
PIPELINE_DIR="$BASE_DIR"                          # location of R scripts
SEEDS=(seed_42 seed_123 seed_456 seed_789 seed_999)
DATASETS=(dataset_001 dataset_002 dataset_003)
QC_ORDERS=(qc1 qc2 qc3 qc4 qc5 qc6)
DATA_PREFIX="genotypes"                           # plink prefix inside each dataset dir
CAUSAL_SNP_SUBPATH="phenotypes/causal_snps.tsv"  # relative to seed dir

METRICS_SCRIPT="$PIPELINE_DIR/compute_qc_metrics.R"
LOG_FILE="$BASE_DIR/run_compute_metrics.log"
FORCE=0

# ── Parse arguments ───────────────────────────────────────────────────────────
for arg in "$@"; do
    case $arg in
        --force) FORCE=1 ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# ── Helpers ───────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

# Extract a single numeric value from column $2 of a TSV file $1
# Returns empty string if file missing or column absent
extract_tsv_value() {
    local tsv="$1"
    local col="$2"
    if [ ! -f "$tsv" ]; then echo ""; return; fi
    awk -F'\t' -v col="$col" '
        NR==1 { for(i=1;i<=NF;i++) if($i==col) idx=i }
        NR==2 { if(idx) print $idx }
    ' "$tsv"
}

# ── Preflight checks ──────────────────────────────────────────────────────────
if [ ! -f "$METRICS_SCRIPT" ]; then
    echo "ERROR: compute_qc_metrics.R not found at $METRICS_SCRIPT"
    exit 1
fi

log "=== run_compute_metrics.sh starting ==="
log "Base dir  : $BASE_DIR"
log "Seeds     : ${SEEDS[*]}"
log "Datasets  : ${DATASETS[*]}"
log "QC orders : ${QC_ORDERS[*]}"
log "Force mode: $FORCE"
log ""

n_run=0; n_skip=0; n_fail=0; n_missing=0

# ── Main loop ─────────────────────────────────────────────────────────────────
for SEED in "${SEEDS[@]}"; do
    SEED_DIR="$BASE_DIR/$SEED"

    # causal_snps.tsv lives in seed_*/phenotypes/ (shared across all datasets within a seed)
    CAUSAL_SNPS="$SEED_DIR/$CAUSAL_SNP_SUBPATH"
    if [ ! -f "$CAUSAL_SNPS" ]; then
        log "  [WARN] causal_snps.tsv not found for $SEED — skipping seed"
        log "         Expected: $CAUSAL_SNPS"
        (( n_missing += 6 * ${#DATASETS[@]} )) || true
        continue
    fi

    for DATASET in "${DATASETS[@]}"; do
        DATASET_DIR="$SEED_DIR/replicates/$DATASET"

        if [ ! -d "$DATASET_DIR" ]; then
            log "  [SKIP] Directory not found: $DATASET_DIR"
            (( n_missing += 6 )) || true
            continue
        fi

        log "--- $SEED / $DATASET ---"

        # ── Step 1: QC4 first (reference — no baseline args) ─────────────────
        QC4_SAIGE="$DATASET_DIR/${DATA_PREFIX}_qc4_saige_result.txt"
        QC4_BIM="$DATASET_DIR/${DATA_PREFIX}_qc4_final_autosomes.bim"
        QC4_FAM="$DATASET_DIR/${DATA_PREFIX}_qc4_final_autosomes.fam"
        QC4_METRICS="$DATASET_DIR/${DATA_PREFIX}_qc4_metrics.tsv"
        QC4_SNPS="$DATASET_DIR/${DATA_PREFIX}_qc4_snps.txt"
        QC4_SAMPS="$DATASET_DIR/${DATA_PREFIX}_qc4_samps.txt"

        if [ ! -f "$QC4_SAIGE" ]; then
            log "  [SKIP] QC4 SAIGE results not found — cannot run any ordering"
            log "         Missing: $QC4_SAIGE"
            (( n_missing += 6 )) || true
            continue
        fi

        # Run QC4 metrics (skip if already done and not forcing)
        if [ -s "$QC4_METRICS" ] && [ "$FORCE" -eq 0 ]; then
            log "  [SKIP] QC4 metrics already exist"
        else
            log "  [RUN ] QC4 (reference)"
            if Rscript "$METRICS_SCRIPT" \
                --saige_results "$QC4_SAIGE" \
                --causal_snps   "$CAUSAL_SNPS" \
                --final_bim     "$QC4_BIM" \
                --final_fam     "$QC4_FAM" \
                --qc_label      "QC4" \
                --seed          "$SEED" \
                --dataset       "$DATASET" \
                --out           "$DATASET_DIR/${DATA_PREFIX}_qc4" \
                >> "$LOG_FILE" 2>&1
            then
                log "  [DONE] QC4"
                (( n_run++ )) || true
            else
                log "  [FAIL] QC4"
                (( n_fail++ )) || true
                continue   # without QC4 baseline we cannot run other orderings
            fi
        fi

        # ── Step 2: Extract QC4 baseline values from its metrics TSV ─────────
        QC4_LAMBDA=$(    extract_tsv_value "$QC4_METRICS" "lambda_gc")
        QC4_LAMBDA_1000=$(extract_tsv_value "$QC4_METRICS" "lambda_1000")

        if [ -z "$QC4_LAMBDA" ] || [ -z "$QC4_LAMBDA_1000" ]; then
            log "  [WARN] Could not extract QC4 lambda values from $QC4_METRICS"
            log "         QC4 baselines will be omitted for other orderings"
        else
            log "  [INFO] QC4 baselines: lambda_gc=$QC4_LAMBDA  lambda_1000=$QC4_LAMBDA_1000"
        fi

        # ── Step 3: Write QC4 SNP and sample ID lists for Jaccard ────────────
        if [ ! -s "$QC4_SNPS" ] || [ "$FORCE" -eq 1 ]; then
            if [ -f "$QC4_BIM" ]; then
                awk '{print $2}' "$QC4_BIM" > "$QC4_SNPS"
                log "  [INFO] QC4 SNP list written ($(wc -l < "$QC4_SNPS") SNPs)"
            fi
        fi
        if [ ! -s "$QC4_SAMPS" ] || [ "$FORCE" -eq 1 ]; then
            if [ -f "$QC4_FAM" ]; then
                awk '{print $2}' "$QC4_FAM" > "$QC4_SAMPS"
                log "  [INFO] QC4 sample list written ($(wc -l < "$QC4_SAMPS") samples)"
            fi
        fi

        # ── Step 4: Run all other orderings with QC4 as baseline ─────────────
        for QCN in "${QC_ORDERS[@]}"; do
            [ "$QCN" = "qc4" ] && continue    # already done above

            N="${QCN#qc}"                      # strip "qc" → 1,2,3,5,6
            LABEL="QC${N^}"                    # → QC1, QC2, ... (bash 4+)
            LABEL="QC${N}"                     # simpler: QC1, QC2, QC3, QC5, QC6

            SAIGE_FILE="$DATASET_DIR/${DATA_PREFIX}_qc${N}_saige_result.txt"
            BIM_FILE="$DATASET_DIR/${DATA_PREFIX}_qc${N}_final_autosomes.bim"
            FAM_FILE="$DATASET_DIR/${DATA_PREFIX}_qc${N}_final_autosomes.fam"
            OUT_PREFIX="$DATASET_DIR/${DATA_PREFIX}_qc${N}"
            METRICS_TSV="${OUT_PREFIX}_metrics.tsv"

            # Check SAIGE results exist
            if [ ! -f "$SAIGE_FILE" ]; then
                log "  [MISS] qc${N}: SAIGE results not found — skipping"
                log "         Missing: $SAIGE_FILE"
                (( n_missing++ )) || true
                continue
            fi

            # Skip if already done (unless --force)
            if [ -s "$METRICS_TSV" ] && [ "$FORCE" -eq 0 ]; then
                log "  [SKIP] qc${N}: metrics already exist"
                (( n_skip++ )) || true
                continue
            fi

            log "  [RUN ] qc${N}"

            # Build baseline arguments conditionally
            BASELINE_ARGS=""
            [ -n "$QC4_LAMBDA"      ] && BASELINE_ARGS="$BASELINE_ARGS --baseline_lambda $QC4_LAMBDA"
            [ -n "$QC4_LAMBDA_1000" ] && BASELINE_ARGS="$BASELINE_ARGS --baseline_lambda1000 $QC4_LAMBDA_1000"
            [ -s "$QC4_SNPS"        ] && BASELINE_ARGS="$BASELINE_ARGS --baseline_snps $QC4_SNPS"
            [ -s "$QC4_SAMPS"       ] && BASELINE_ARGS="$BASELINE_ARGS --baseline_samps $QC4_SAMPS"
            [ -f "$QC4_SAIGE"       ] && BASELINE_ARGS="$BASELINE_ARGS --ref_results $QC4_SAIGE"

            # shellcheck disable=SC2086
            if Rscript "$METRICS_SCRIPT" \
                --saige_results "$SAIGE_FILE" \
                --causal_snps   "$CAUSAL_SNPS" \
                --final_bim     "$BIM_FILE" \
                --final_fam     "$FAM_FILE" \
                --qc_label      "$LABEL" \
                --seed          "$SEED" \
                --dataset       "$DATASET" \
                --out           "$OUT_PREFIX" \
                $BASELINE_ARGS \
                >> "$LOG_FILE" 2>&1
            then
                log "  [DONE] qc${N}"
                (( n_run++ )) || true
            else
                log "  [FAIL] qc${N} — check $LOG_FILE for R errors"
                (( n_fail++ )) || true
            fi
        done

        log ""
    done
done

# ── Summary ───────────────────────────────────────────────────────────────────
log "=== run_compute_metrics.sh complete ==="
log "  Completed : $n_run"
log "  Skipped   : $n_skip  (already exist; use --force to rerun)"
log "  Missing   : $n_missing  (SAIGE results not yet available)"
log "  Failed    : $n_fail"
log ""
if [ "$n_fail" -gt 0 ]; then
    log "  Failures logged to: $LOG_FILE"
    log "  Search for '[FAIL]' to locate them"
fi
