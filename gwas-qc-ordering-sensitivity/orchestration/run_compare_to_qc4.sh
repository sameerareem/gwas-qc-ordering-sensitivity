#!/bin/bash
# =============================================================================
#  run_compare_to_qc4.sh
#  Run compare_to_qc4.R once across all seeds and available datasets.
#
#  compare_to_qc4.R searches BASE_DIR recursively for *_metrics.tsv files,
#  aggregates them, and runs all cross-ordering statistical tests and plots.
#  It only needs to be run once — seed-level aggregation is handled internally.
#
#  Prerequisites:
#    run_compute_metrics.sh must have completed first so that
#    *_metrics.tsv files exist under BASE_DIR.
#
#  Usage:
#    bash run_compare_to_qc4.sh [--force]
#
#  --force : rerun even if output files already exist
#
#  Outputs (written to OUTPUT_DIR):
#    qc4_comparison_seed_level.tsv      per-seed aggregated metrics
#    qc4_comparison_all_metrics.tsv     full raw metrics (all runs)
#    qc4_comparison_friedman.tsv        Friedman omnibus test results
#    qc4_comparison_wilcoxon.tsv        Wilcoxon pairwise vs QC4 results
#    qc4_comparison_pairwise_tests.tsv  all-pairs Wilcoxon results
#    qc4_comparison_lambda1000_boxplot.png
#    qc4_comparison_seed_interaction.png
#    qc4_comparison_recovery_boxplot.png
#    qc4_comparison_pvalue_heatmap.png
#    qc4_comparison_lins_ccc.png
#    qc4_comparison_pairwise_<metric>.png  (one per metric)
#    qc4_comparison_pairwise_all.png       (combined figure)
#    qc4_comparison_dataset_heatmap.png    (only if >1 dataset per seed)
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR="/mnt/hdd/GWAS_synth_data/gwas_synthetic_output"
PIPELINE_DIR="$BASE_DIR"
COMPARE_SCRIPT="$PIPELINE_DIR/compare_to_qc4.R"
OUTPUT_DIR="$BASE_DIR/comparison"
OUT_PREFIX="$OUTPUT_DIR/qc4_comparison"
LOG_FILE="$OUTPUT_DIR/run_compare_to_qc4.log"
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

# ── Setup ─────────────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"

log "=== run_compare_to_qc4.sh starting ==="
log "Base dir    : $BASE_DIR"
log "Output dir  : $OUTPUT_DIR"
log "Script      : $COMPARE_SCRIPT"
log "Force mode  : $FORCE"
log ""

# ── Preflight: script exists ──────────────────────────────────────────────────
if [ ! -f "$COMPARE_SCRIPT" ]; then
    log "[ERROR] compare_to_qc4.R not found at: $COMPARE_SCRIPT"
    exit 1
fi

# ── Preflight: count available metrics TSV files ──────────────────────────────
N_METRICS=$(find "$BASE_DIR" -name "*_metrics.tsv" \
            ! -path "$OUTPUT_DIR/*" 2>/dev/null | wc -l)

if [ "$N_METRICS" -eq 0 ]; then
    log "[ERROR] No *_metrics.tsv files found under $BASE_DIR"
    log "        Run run_compute_metrics.sh first to generate per-ordering metrics."
    exit 1
fi

log "[INFO] Found $N_METRICS metrics TSV file(s) across seeds and datasets"

# Show breakdown by seed and QC ordering
log ""
log "Metrics files by seed:"
for SEED_DIR in "$BASE_DIR"/seed_*/; do
    SEED=$(basename "$SEED_DIR")
    N_SEED=$(find "$SEED_DIR" -name "*_metrics.tsv" 2>/dev/null | wc -l)
    if [ "$N_SEED" -gt 0 ]; then
        # Count distinct QC orderings and datasets
        QC_FOUND=$(find "$SEED_DIR" -name "*_metrics.tsv" 2>/dev/null \
                   | xargs -I{} basename {} \
                   | grep -oP 'qc\d+' | sort -u | tr '\n' ' ')
        DS_FOUND=$(find "$SEED_DIR" -name "*_metrics.tsv" 2>/dev/null \
                   | xargs -I{} dirname {} \
                   | xargs -I{} basename {} | sort -u | tr '\n' ' ')
        log "  $SEED : $N_SEED file(s) | QC orderings: ${QC_FOUND}| Datasets: ${DS_FOUND}"
    fi
done
log ""

# ── Check if output already exists ───────────────────────────────────────────
SENTINEL="${OUT_PREFIX}_all_metrics.tsv"
if [ -f "$SENTINEL" ] && [ "$FORCE" -eq 0 ]; then
    log "[SKIP] Output already exists: $SENTINEL"
    log "       Use --force to rerun."
    log ""
    log "Existing outputs in $OUTPUT_DIR:"
    ls -lh "$OUTPUT_DIR"/ 2>/dev/null | grep "qc4_comparison" | \
        awk '{print "  " $NF "  (" $5 ")"}' | tee -a "$LOG_FILE"
    exit 0
fi

# ── Run compare_to_qc4.R ─────────────────────────────────────────────────────
log "[RUN] Starting compare_to_qc4.R"
log "      metrics_dir : $BASE_DIR"
log "      out         : $OUT_PREFIX"
log ""

if Rscript "$COMPARE_SCRIPT" \
    --metrics_dir "$BASE_DIR" \
    --out         "$OUT_PREFIX" \
    >> "$LOG_FILE" 2>&1; then

    log ""
    log "[DONE] compare_to_qc4.R completed successfully"
    log ""

    # ── Report output files ───────────────────────────────────────────────────
    log "Output files:"
    for f in "$OUTPUT_DIR"/qc4_comparison_*; do
        [ -f "$f" ] && log "  $(basename $f)  ($(du -h "$f" | cut -f1))"
    done

    # ── Quick result summary from TSVs ────────────────────────────────────────
    FRIEDMAN_TSV="${OUT_PREFIX}_friedman.tsv"
    WILCOXON_TSV="${OUT_PREFIX}_wilcoxon.tsv"

    if [ -f "$FRIEDMAN_TSV" ]; then
        log ""
        log "Friedman test summary (significant results at p<0.05):"
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next}
                    $col["p"]+0 < 0.05 {
                        printf "  %-24s chi2=%-7s p=%s %s\n",
                        $col["metric"], $col["chi2"], $col["p"], $col["sig"]
                    }' "$FRIEDMAN_TSV" | tee -a "$LOG_FILE" || true
        SIG_FRIEDMAN=$(awk -F'\t' 'NR>1 && $5+0<0.05' "$FRIEDMAN_TSV" | wc -l)
        [ "$SIG_FRIEDMAN" -eq 0 ] && log "  (none at p<0.05)"
    fi

    if [ -f "$WILCOXON_TSV" ]; then
        log ""
        log "Wilcoxon significant comparisons vs QC4 (BH-adj p<0.1):"
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next}
                    $col["p_adj"]+0 < 0.1 {
                        printf "  %-24s %-6s med_diff=%-9s p_adj=%s %s\n",
                        $col["metric"], $col["ordering"],
                        $col["median_diff"], $col["p_adj"], $col["sig"]
                    }' "$WILCOXON_TSV" | tee -a "$LOG_FILE" || true
        SIG_WILCOX=$(awk -F'\t' 'NR>1 && $10+0<0.1' "$WILCOXON_TSV" | wc -l)
        [ "$SIG_WILCOX" -eq 0 ] && log "  (none at BH-adj p<0.1)"
    fi

else
    log ""
    log "[FAIL] compare_to_qc4.R exited with an error"
    log "       Check the log for R error messages:"
    log "       $LOG_FILE"
    log ""
    log "Last 20 lines of log:"
    tail -20 "$LOG_FILE" | tee -a /dev/stderr
    exit 1
fi

log ""
log "=== run_compare_to_qc4.sh complete ==="
log "All outputs written to: $OUTPUT_DIR"
