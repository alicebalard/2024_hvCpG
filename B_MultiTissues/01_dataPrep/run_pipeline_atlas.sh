#!/bin/bash
#$ -N prepAtlas
#$ -S /bin/bash
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=3:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

source /share/apps/source_files/python/python-3.13.0a6.source

# =============================================================================
# run_pipeline.sh — CpG methylation data preparation pipeline
# =============================================================================
#
# Usage:
#   bash run_pipeline.sh --preset PRESET [--config FILE]
#
# Atlas presets (Loyfer WGBS):
#   atlas_general              All data, Source Tissue - Cell type grouping (baseline)
#   02_rmMultSamples           One sample per individual (dedup by PatientID)
#   04_maleOnly                Groups where all samples are male
#   05_femaleOnly6gp           Female-only groups, sample 6
#   06_bothsexes6gp            Mixed-sex groups, sample 6
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
# BASH_SOURCE[0] points to the SGE spool copy of the script, not the
# original — so SCRIPT_DIR would resolve to the spool directory.
# We derive it from the real script path if possible, but CODE_DIR
# can always be overridden via environment variable (recommended on SGE).
_REAL_SCRIPT="$(readlink -f "${BASH_SOURCE[0]}" 2>/dev/null || echo "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$_REAL_SCRIPT")" && pwd)"
PYTHON="${PYTHON:-python3}"

# ──────────────────────────────────────────────
#  Default paths — override via env vars or --config
# ──────────────────────────────────────────────

DATA_DIR="${DATA_DIR:-/SAN/ghlab/epigen/Alice/hvCpG_project/data}"
CODE_DIR="${CODE_DIR:-/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep}"

WGBS_DIR="${WGBS_DIR:-${DATA_DIR}/WGBS_human/AtlasLoyfer}"
BETA_PATTERN="${BETA_PATTERN:-${WGBS_DIR}/betaFiles/GSM*.hg38.beta}"
CPG_BED="${CPG_BED:-${WGBS_DIR}/wgbs_tools/references/hg38/CpG.bed.gz}"
META_ATLAS="${META_ATLAS:-${SCRIPT_DIR}/SupTab1_Loyfer2023_amended.csv}"

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
EXCLUDE_SITES="${EXCLUDE_SITES:-/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/snps_maf_0.01_chr_pos.txt}"

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
  # meta_flags is passed as individual arguments using "$@" after out_prefix.
  # This avoids the classic bash word-splitting bug where flags containing
  # spaces (e.g. --group_col "Germ layer") would be split incorrectly when
  # stored in a plain string variable and expanded unquoted.
  local out_prefix="$1"; shift
  # Remaining positional args are the metadata flags (may be zero args)
  local META_OUT="${WGBS_DIR}/${out_prefix}.csv"

  OUTPUT_FOLDER="${WGBS_DIR}/output_${out_prefix}"
  mkdir -p "$OUTPUT_FOLDER"

  _run "01_prepare_metadata" \
    "$PYTHON" "${CODE_DIR}/prepare_metadata.py" \
      --meta    "$META_ATLAS" \
      --output  "$META_OUT" \
      --seed    "$RANDOM_SEED" \
      "$@"

  _run "02_prepare_matrix" \
    "$PYTHON" "${CODE_DIR}/prepare_beta_matrices.py" \
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

# Script 2 — Remove multiple samples per individual
run_02_rmMultSamples() {
  _run_atlas "02_rmMultSamples" --dedup_col "PatientID"
}

# Script 4 — Male-only groups (all samples in group are male)
run_04_maleOnly() {
  _run_atlas "04_maleOnly" --sex_filter M
}

# Script 5 — Female-only groups, randomly sample 6
run_05_femaleOnly6gp() {
  _run_atlas "05_femaleOnly6gp" --sex_filter F --sample_n_groups 6 --min_per_group 3
}

# Script 6 — Mixed-sex groups (both M and F present), sample 6
run_06_bothsexes6gp() {
  _run_atlas "06_bothsexes6gp" --mixed_sex --sample_n_groups 6 --min_per_group 3
}

# Script 9 — Immune cells only
run_09_immuneOnly() {
  _run_atlas "09_immuneOnly" --keep_col "Immune?" --keep_val True
}

# Script 10 — No immune cells
run_10_noImmune() {
  _run_atlas "10_noImmune" --exclude_immune
}

# Script 11 — No immune cells, sample 11 groups
run_11_noImmune_sample11gp() {
  _run_atlas "11_noImmune_sample11gp" --exclude_immune --sample_n_groups 11 --min_per_group 3
}

# Script 12 — Endoderm only
run_12_endo() {
  _run_atlas "12_endo" --germ_layer Endo
}

# Script 12.2 — Endoderm only, sample 6 groups
run_12_2_endo6gp() {
  _run_atlas "12_2_endo6gp" --germ_layer Endo --sample_n_groups 6 --min_per_group 3
}

# Script 13 — Mesoderm only
run_13_meso() {
  _run_atlas "13_meso" --germ_layer Meso
}

# Script 13.2 — Mesoderm only, sample 6 groups
run_13_2_meso6gp() {
  _run_atlas "13_2_meso6gp" --germ_layer Meso --sample_n_groups 6 --min_per_group 3
}

# Script 14 — Ectoderm only
run_14_ecto() {
  _run_atlas "14_ecto" --germ_layer Ecto
}

# Script 15 — Paired samples: 2 males per group
run_15_pairs_MM() {
  _run_atlas "15_pairs_MM" --pairs MM
}

# Script 16 — Paired samples: 2 females per group
run_16_pairs_FF() {
  _run_atlas "16_pairs_FF" --pairs FF
}

# Script 17 — Paired samples: 1 male + 1 female per group
run_17_pairs_MF() {
  _run_atlas "17_pairs_MF" --pairs MF
}

# General atlas — all data, Source Tissue - Cell type grouping (baseline, no filtering)
run_atlas_general() {
  _run_atlas "atlas_general"
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
  atlas_general)           run_atlas_general ;;
  02_rmMultSamples)        run_02_rmMultSamples ;;
  04_maleOnly)             run_04_maleOnly ;;
  05_femaleOnly6gp)        run_05_femaleOnly6gp ;;
  06_bothsexes6gp)         run_06_bothsexes6gp ;;
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
  "")
    echo ""
    echo "No preset specified. Available presets:"
    echo ""
    echo "  atlas_general          All data, Source Tissue - Cell type grouping (baseline)"
    echo "  02_rmMultSamples       One sample per individual (dedup by PatientID)"
    echo "  04_maleOnly            Groups where all samples are male"
    echo "  05_femaleOnly6gp       Female-only groups, sample 6"
    echo "  06_bothsexes6gp        Mixed-sex groups, sample 6"
    echo "  09_immuneOnly          Immune cells only"
    echo "  10_noImmune            No immune cells"
    echo "  11_noImmune_sample11gp No immune, sample 11 groups"
    echo "  12_endo                Endoderm only"
    echo "  12_2_endo6gp           Endoderm, sample 6 groups"
    echo "  13_meso                Mesoderm only"
    echo "  13_2_meso6gp           Mesoderm, sample 6 groups"
    echo "  14_ecto                Ectoderm only"
    echo "  15_pairs_MM            2 males per group"
    echo "  16_pairs_FF            2 females per group"
    echo "  17_pairs_MF            1 male + 1 female per group"
    echo ""
    echo "Usage: bash run_pipeline_atlas.sh --preset atlas_general"
    echo "       EXCLUDE_SITES=/data/snp.txt bash run_pipeline_atlas.sh --preset 10_noImmune"
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
