#!/bin/bash
# =============================================================================
# run_pipeline.sh — CpG methylation data preparation pipeline
# =============================================================================
#
# Usage:
#   bash run_pipeline.sh --preset PRESET [--config FILE]
#
# Atlas presets (Loyfer WGBS):
#   01_byDevLayer              Group by germ/developmental layer
#   02_rmMultSamples           One sample per individual (dedup by PatientID)
#   03_correspMariaTissues     Only tissues also present in Maria's arrays
#   04_maleOnly                Groups where all samples are male
#   05_femaleOnly6gp           Female-only groups, sample 6
#   06_bothsexes6gp            Mixed-sex groups, sample 6
#   08_byTissue                Group by simplified tissue label
#   09_immuneOnly              Immune cells only
#   10_noImmune                No immune cells
#   11_noImmune_sample11gp     No immune, sample 11 groups
#   12_endo                    Endoderm only
#   12_2_endo6gp               Endoderm, sample 6 groups
#   13_meso                    Mesoderm only
#   13_2_meso6gp               Mesoderm, sample 6 groups
#   14_ecto                    Ectoderm only
#   15_pairs_MM                2 males per group
#   16_pairs_FF                2 females per group
#   17_pairs_MF                1 male + 1 female per group
#
# Override any path or threshold via environment variables:
#   DATA_DIR=/my/data bash run_pipeline.sh --preset 04_maleOnly
#   MIN_DATASETS_ATLAS=30 bash run_pipeline.sh --preset 10_noImmune
#
# Or pass a config file:
#   bash run_pipeline.sh --config my_config.sh
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"

# ──────────────────────────────────────────────
#  Default paths — override via env vars or --config
# ──────────────────────────────────────────────

DATA_DIR="${DATA_DIR:-/SAN/ghlab/epigen/Alice/hvCpG_project/data}"
CODE_DIR="${CODE_DIR:-${SCRIPT_DIR}}"

WGBS_DIR="${WGBS_DIR:-${DATA_DIR}/WGBS_human/AtlasLoyfer}"
BETA_PATTERN="${BETA_PATTERN:-${WGBS_DIR}/betaFiles/GSM*.hg38.beta}"
CPG_BED="${CPG_BED:-${WGBS_DIR}/wgbs_tools/references/hg38/CpG.bed.gz}"
META_ATLAS="${META_ATLAS:-${WGBS_DIR}/SupTab1_Loyfer2023.csv}"

GEO_RDS_DIR="${GEO_RDS_DIR:-/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED}"
TCGA_RDS_DIR="${TCGA_RDS_DIR:-/home/alice/tempRDS}"
ARRAY_OUTPUT_BASE="${ARRAY_OUTPUT_BASE:-${DATA_DIR}/Arrays/Maria}"

# Thresholds
MIN_COV="${MIN_COV:-10}"
MIN_SAMPLES="${MIN_SAMPLES:-3}"
MIN_DATASETS_ATLAS="${MIN_DATASETS_ATLAS:-46}"
MIN_DATASETS_ARRAYS="${MIN_DATASETS_ARRAYS:-15}"
MAX_NA="${MAX_NA:-0.2}"
CHUNK_SIZE="${CHUNK_SIZE:-100000}"
LAMBDA_PERCENTILE="${LAMBDA_PERCENTILE:-95}"
RANDOM_SEED="${RANDOM_SEED:-42}"
# Optional: path to a file listing CpG sites to exclude (one chr_pos per line).
# Leave empty to skip exclusion.
EXCLUDE_SITES="${EXCLUDE_SITES:-}"

# ──────────────────────────────────────────────
#  Argument parsing
# ──────────────────────────────────────────────

PRESET=""
CONFIG_FILE=""

usage() {
  grep "^#" "$0" | grep -v "^#!/" | sed "s/^# \{0,2\}//" | head -60
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --preset) PRESET="$2";      shift 2 ;;
    --config) CONFIG_FILE="$2"; shift 2 ;;
    --help|-h) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

[[ -n "$CONFIG_FILE" ]] && { echo "Loading config: $CONFIG_FILE"; source "$CONFIG_FILE"; }

# ──────────────────────────────────────────────
#  Helper: run a named step with logging
# ──────────────────────────────────────────────

_run() {
  local step="$1"; shift
  local log_dir="${OUTPUT_FOLDER}/logs"
  mkdir -p "$log_dir"
  echo ""
  echo "================================================================"
  echo "  $step"
  echo "================================================================"
  "$@" 2>&1 | tee "${log_dir}/${step}.log"
  echo "  [done] $step"
}

# ──────────────────────────────────────────────
#  Shared helper: prepare metadata then run WGBS matrix builder
#
#  Arguments:
#    $1  meta_flags   — extra flags passed to prepare_metadata.py
#    $2  out_prefix   — output prefix (used for filenames AND output folder name)
# ──────────────────────────────────────────────

_run_atlas() {
  local meta_flags="$1"
  local out_prefix="$2"
  local META_OUT="${WGBS_DIR}/${out_prefix}.csv"

  OUTPUT_FOLDER="${WGBS_DIR}/output_${out_prefix}"
  mkdir -p "$OUTPUT_FOLDER"

  _run "01_prepare_metadata" \
    $PYTHON "${CODE_DIR}/prepare_metadata.py" \
      --meta    "$META_ATLAS" \
      --output  "$META_OUT" \
      --seed    "$RANDOM_SEED" \
      $meta_flags

  _run "02_prepare_matrix" \
    $PYTHON "${CODE_DIR}/prepare_beta_matrices.py" \
      --beta_files         "$BETA_PATTERN" \
      --cpg_bed            "$CPG_BED" \
      --output_folder      "$OUTPUT_FOLDER" \
      --meta               "$META_OUT" \
      --minCov             "$MIN_COV" \
      --min_samples        "$MIN_SAMPLES" \
      --min_datasets       "$MIN_DATASETS_ATLAS" \
      --chunk_size         "$CHUNK_SIZE" \
      --lambda_percentile  "$LAMBDA_PERCENTILE" \
      --output_prefix      "$out_prefix" \
      ${EXCLUDE_SITES:+--exclude_sites "$EXCLUDE_SITES"}
}

# ══════════════════════════════════════════════════════════════════════════════
#  Atlas presets — one per original script
# ══════════════════════════════════════════════════════════════════════════════

# Script 1 — Group by developmental/germ layer
run_01_byDevLayer() {
  _run_atlas '--group_col "Germ layer"' "01_byDevLayer"
}

# Script 2 — Remove multiple samples per individual
run_02_rmMultSamples() {
  _run_atlas '--dedup_col "PatientID"' "02_rmMultSamples"
}

# Script 3 — Keep only tissues present in Maria's datasets
run_03_correspMariaTissues() {
  _run_atlas '--keep_col "Found in Maria DS?" --keep_val True' "03_correspMariaTissues"
}

# Script 4 — Male-only groups (all samples in group are male)
run_04_maleOnly() {
  _run_atlas "--sex_filter M" "04_maleOnly"
}

# Script 5 — Female-only groups, randomly sample 6
run_05_femaleOnly6gp() {
  _run_atlas "--sex_filter F --sample_n_groups 6 --min_per_group 3" "05_femaleOnly6gp"
}

# Script 6 — Mixed-sex groups (both M and F present), sample 6
run_06_bothsexes6gp() {
  _run_atlas "--mixed_sex --sample_n_groups 6 --min_per_group 3" "06_bothsexes6gp"
}

# Script 8 — Group by simplified tissue label
run_08_byTissue() {
  _run_atlas '--group_col "Group simplified"' "08_byTissue"
}

# Script 9 — Immune cells only
run_09_immuneOnly() {
  _run_atlas '--keep_col "Immune?" --keep_val True' "09_immuneOnly"
}

# Script 10 — No immune cells
run_10_noImmune() {
  _run_atlas "--exclude_immune" "10_noImmune"
}

# Script 11 — No immune cells, sample 11 groups
run_11_noImmune_sample11gp() {
  _run_atlas "--exclude_immune --sample_n_groups 11 --min_per_group 3" "11_noImmune_sample11gp"
}

# Script 12 — Endoderm only
run_12_endo() {
  _run_atlas "--germ_layer Endo" "12_endo"
}

# Script 12.2 — Endoderm only, sample 6 groups
run_12_2_endo6gp() {
  _run_atlas "--germ_layer Endo --sample_n_groups 6 --min_per_group 3" "12_2_endo6gp"
}

# Script 13 — Mesoderm only
run_13_meso() {
  _run_atlas "--germ_layer Meso" "13_meso"
}

# Script 13.2 — Mesoderm only, sample 6 groups
run_13_2_meso6gp() {
  _run_atlas "--germ_layer Meso --sample_n_groups 6 --min_per_group 3" "13_2_meso6gp"
}

# Script 14 — Ectoderm only
run_14_ecto() {
  _run_atlas "--germ_layer Ecto" "14_ecto"
}

# Script 15 — Paired samples: 2 males per group
run_15_pairs_MM() {
  _run_atlas "--pairs MM" "15_pairs_MM"
}

# Script 16 — Paired samples: 2 females per group
run_16_pairs_FF() {
  _run_atlas "--pairs FF" "16_pairs_FF"
}

# Script 17 — Paired samples: 1 male + 1 female per group
run_17_pairs_MF() {
  _run_atlas "--pairs MF" "17_pairs_MF"
}

# ══════════════════════════════════════════════════════════════════════════════
#  Array presets (Maria RDS datasets)
# ══════════════════════════════════════════════════════════════════════════════

_run_arrays() {
  local extra_flags="$1"
  local out_prefix="$2"
  OUTPUT_FOLDER="${ARRAY_OUTPUT_BASE}/output_${out_prefix}"
  mkdir -p "$OUTPUT_FOLDER"

  _run "02_prepare_matrix" \
    $PYTHON "${CODE_DIR}/prepare_arrays.py" \
      --rds_folders       "$GEO_RDS_DIR" "$TCGA_RDS_DIR" \
      --output_folder     "$OUTPUT_FOLDER" \
      --maxNA             "$MAX_NA" \
      --min_datasets      "$MIN_DATASETS_ARRAYS" \
      --lambda_percentile "$LAMBDA_PERCENTILE" \
      --output_prefix     "$out_prefix" \
      ${EXCLUDE_SITES:+--exclude_sites "$EXCLUDE_SITES"} \
      $extra_flags
}

run_arrays_all() {
  _run_arrays "" "arrays_all"
}

# ══════════════════════════════════════════════════════════════════════════════
#  Dispatch
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "================================================================"
echo "  CpG methylation data preparation pipeline"
echo "  Preset  : ${PRESET:-<none>}"
echo "  Data    : ${DATA_DIR}"
echo "================================================================"

case "$PRESET" in
  # Atlas presets
  01_byDevLayer)           run_01_byDevLayer ;;
  02_rmMultSamples)        run_02_rmMultSamples ;;
  03_correspMariaTissues)  run_03_correspMariaTissues ;;
  04_maleOnly)             run_04_maleOnly ;;
  05_femaleOnly6gp)        run_05_femaleOnly6gp ;;
  06_bothsexes6gp)         run_06_bothsexes6gp ;;
  08_byTissue)             run_08_byTissue ;;
  09_immuneOnly)           run_09_immuneOnly ;;
  10_noImmune)             run_10_noImmune ;;
  11_noImmune_sample11gp)  run_11_noImmune_sample11gp ;;
  12_endo)                 run_12_endo ;;
  12_2_endo6gp)            run_12_2_endo6gp ;;
  13_meso)                 run_13_meso ;;
  13_2_meso6gp)            run_13_2_meso6gp ;;
  14_ecto)                 run_14_ecto ;;
  15_pairs_MM)             run_15_pairs_MM ;;
  16_pairs_FF)             run_16_pairs_FF ;;
  17_pairs_MF)             run_17_pairs_MF ;;
  # Array presets
  arrays_all)              run_arrays_all ;;
  "")
    echo ""
    echo "No preset specified. Available atlas presets:"
    echo "  01_byDevLayer  02_rmMultSamples  03_correspMariaTissues"
    echo "  04_maleOnly    05_femaleOnly6gp   06_bothsexes6gp"
    echo "  08_byTissue    09_immuneOnly      10_noImmune"
    echo "  11_noImmune_sample11gp"
    echo "  12_endo  12_2_endo6gp  13_meso  13_2_meso6gp  14_ecto"
    echo "  15_pairs_MM  16_pairs_FF  17_pairs_MF"
    echo ""
    echo "Array presets:  arrays_all"
    echo ""
    echo "Usage: bash run_pipeline.sh --preset 04_maleOnly"
    exit 1
    ;;
  *)
    echo "Unknown preset: ${PRESET}"
    exit 1
    ;;
esac

echo ""
echo "================================================================"
echo "  Pipeline complete!"
echo "================================================================"
