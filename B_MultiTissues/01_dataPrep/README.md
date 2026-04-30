# CpG Methylation Data Preparation Pipeline

Modular pipeline for preparing CpG methylation data for hvCpG analysis.
Handles WGBS atlas data (`.beta` files) and Illumina array data (`.RDS` / `.h5`).

---

## File structure

```
cpg_pipeline/
├── cpg_utils.py                       ← Shared core library
├── prepare_metadata.py                ← Step 1: annotate & filter metadata
├── prepare_beta_matrices.py           ← Step 2: WGBS .beta files → HDF5
├── run_pipeline_atlas.sh              ← One-command runner (16 presets)
├── SupTab1_Loyfer2023_amended.csv     ← Single metadata source file
└── README.md
```

---

## Key design principles

**Median SD and lambda are computed AFTER the final CpG selection step.**
`compute_stats_on_filtered_matrix()` in `cpg_utils.py` is the single
point of truth — called at the end of every entry-point script, after
the `min_datasets` coverage filter has been applied.

**All 18 `prepare_metadata_*.py` scripts are replaced by one.**
`prepare_metadata.py` accepts composable flags; every original script
maps to a specific flag combination documented below.

**CpG site exclusion happens at the earliest possible point.**
Pass `--exclude_sites` a plain-text file (one `chr_pos` per line) to
any matrix-building script. For WGBS data the mask is applied before
Pass 1, so excluded sites never accumulate coverage counts. For array
data it is applied immediately after the common CpG set is determined,
before alignment and NA filtering. Either way, excluded sites cannot
appear in the final matrix, statistics, or CpG name list.

---

## Quick start

```bash
# Run any of the analysis cases:
bash run_pipeline_atlas.sh --preset atlas_general
bash run_pipeline_atlas.sh --preset 10_noImmune
bash run_pipeline_atlas.sh --preset 17_pairs_MF

# Override data directory or thresholds:
DATA_DIR=/my/data bash run_pipeline_atlas.sh --preset 12_endo
MIN_DATASETS_ATLAS=30 bash run_pipeline_atlas.sh --preset 06_bothsexes6gp

# Exclude a blacklist of CpG sites (applies to all presets):
EXCLUDE_SITES=/data/snp_cpgs.txt bash run_pipeline_atlas.sh --preset 10_noImmune

# Use a config file:
bash run_pipeline_atlas.sh --config my_config.sh
```

---

## All presets

| Preset | What it does |
|---|---|
| `atlas_general` | All data, Source Tissue - Cell type grouping (baseline) |
| `02_rmMultSamples` | One sample per PatientID |
| `04_maleOnly` | Groups where all samples are male |
| `05_femaleOnly6gp` | Female-only groups, sample 6 |
| `06_bothsexes6gp` | Mixed-sex groups, sample 6 |
| `09_immuneOnly` | Immune cells only |
| `10_noImmune` | No immune cells |
| `11_noImmune_sample11gp` | No immune, sample 11 groups |
| `12_endo` | Endoderm only |
| `12_2_endo6gp` | Endoderm, sample 6 groups |
| `13_meso` | Mesoderm only |
| `13_2_meso6gp` | Mesoderm, sample 6 groups |
| `14_ecto` | Ectoderm only |
| `15_pairs_MM` | Sample 2 males per group |
| `16_pairs_FF` | Sample 2 females per group |
| `17_pairs_MF` | Sample 1M + 1F per group |

---

## `prepare_metadata.py` — all options

Analysis group is always built as `Source Tissue + " - " + Cell type`.

### Row-level filters

| Argument | Description |
|---|---|
| `--exclude_immune` | Drop rows where `Immune? == True` |
| `--germ_layer LAYER` | Keep rows where `Germ layer == LAYER` (Endo / Meso / Ecto) |
| `--dedup_col COL` | Keep first row per unique value of COL (e.g. `PatientID`) |
| `--keep_col COL --keep_val VAL` | Keep rows where COL == VAL. Repeatable. |
| `--drop_col COL --drop_val VAL` | Drop rows where COL == VAL. Repeatable. |

### Group-level filters

| Argument | Description |
|---|---|
| `--sex_filter M\|F` | Keep groups where ALL samples have this sex |
| `--mixed_sex` | Keep groups with BOTH at least one M and one F |
| `--keep_groups_col COL --keep_groups_val VAL` | Keep groups where all rows match |

### Sampling

| Argument | Description |
|---|---|
| `--sample_n_groups N` | Randomly sample up to N eligible groups |
| `--min_per_group N` | Minimum samples per group for eligibility (default: 3) |
| `--pairs MM\|FF\|MF` | Sample exactly 2 individuals per group (MM/FF/MF) |
| `--n_groups N` | Cap on groups when using `--pairs` |
| `--seed N` | Random seed (default: 42) |

### Flag combinations for each preset

```bash
# General (baseline) — no extra flags
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv

# 02 — rmMultSamples
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --dedup_col "PatientID"

# 04 — maleOnly
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --sex_filter M

# 05 — femaleOnly6gp
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --sex_filter F --sample_n_groups 6

# 06 — bothsexes6gp
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --mixed_sex --sample_n_groups 6

# 09 — immuneOnly
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --keep_col "Immune?" --keep_val True

# 10 — noImmune
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --exclude_immune

# 11 — noImmune_sample11groups
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --exclude_immune --sample_n_groups 11

# 12 — endo
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --germ_layer Endo

# 12.2 — endo6gp
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --germ_layer Endo --sample_n_groups 6

# 13 — meso
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --germ_layer Meso

# 13.2 — meso6gp
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --germ_layer Meso --sample_n_groups 6

# 14 — ecto
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --germ_layer Ecto

# 15 — pairs_MM
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --pairs MM

# 16 — pairs_FF
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --pairs FF

# 17 — pairs_MF
python prepare_metadata.py --meta SupTab1_Loyfer2023_amended.csv --pairs MF
```

---

## Step-by-step usage

### Step 1 — Prepare metadata

```bash
python prepare_metadata.py \
  --meta SupTab1_Loyfer2023_amended.csv \
  [filter flags]
```

### Step 2 — Build WGBS beta matrix

```bash
python prepare_beta_matrices.py \
  --beta_files "/data/betaFiles/GSM*.hg38.beta" \
  --cpg_bed   /data/hg38/CpG.bed.gz \
  --output_folder /data/output/atlas_general \
  --meta      atlas_general.csv \
  --minCov    10 \
  --min_samples 3 \
  --min_datasets 46 \
  --chunk_size 100000 \
  --output_prefix all \
  --exclude_sites /data/snp_blacklist.txt   # optional
```

| Argument | Default | Description |
|---|---|---|
| `--beta_files` | required | Glob pattern for `.beta` files |
| `--cpg_bed` | required | CpG BED reference (`.bed` or `.bed.gz`) |
| `--meta` | required | Filtered metadata CSV (output of Step 1) |
| `--minCov` | 10 | Min read coverage per site |
| `--min_samples` | 3 | Min samples per group |
| `--min_datasets` | 46 | Min datasets per CpG |
| `--chunk_size` | 100000 | Rows per processing chunk |
| `--id_pattern` | `r'-([A-Za-z0-9_]+)\.hg38\.beta$'` | Regex to extract sample ID from filename |
| `--lambda_percentile` | 95 | Percentile for lambda = perc/median |
| `--exclude_sites` | — | Text file of `chr_pos` sites to exclude (applied before Pass 1) |

---

## CpG site exclusion

Pass `--exclude_sites` a plain-text file listing sites to remove:

```
# snp_blacklist.txt
# One chr_pos per line. Blank lines and # comments are ignored.
# Gzipped files (.gz) are also accepted.
chr1_12345
chr3_778901
chrX_55102
```

Excluded sites are applied before Pass 1 — they never accumulate coverage counts and cannot appear in the final matrix or statistics.

```bash
EXCLUDE_SITES=/data/snp_blacklist.txt bash run_pipeline_atlas.sh --preset 10_noImmune
```

Or permanently via a config file:

```bash
# my_config.sh
EXCLUDE_SITES=/data/snp_blacklist.txt
MIN_DATASETS_ATLAS=30
```

```bash
bash run_pipeline_atlas.sh --config my_config.sh --preset 10_noImmune
```

---

## `cpg_utils.py` function reference

| Function | Purpose |
|---|---|
| `load_cpg_names_from_bed()` | Parse CpG BED file (gzipped or not) |
| `load_cpg_names_from_h5()` | Read CpG names from an HDF5 file |
| `load_exclude_sites()` | Load exclusion list → `frozenset` of `chr_pos` strings |
| `build_exclude_mask()` | Convert exclusion frozenset → boolean array aligned to reference |
| `load_beta_file()` | Load `.beta` binary → (beta, coverage) arrays |
| `build_sample_to_path_map()` | Map sample IDs → file paths |
| `filter_valid_groups()` | Filter metadata to groups with ≥ N samples |
| `shorten_sample_ids()` | Strip prefix from sample name |
| `build_beta_matrix_chunked()` | Two-pass chunked WGBS matrix builder |
| `compute_median_and_lambda()` | Compute median SD and lambda = perc/median |
| `compute_stats_on_filtered_matrix()` | Stats on final CpG set — canonical call, always after filtration |
| `write_h5_matrix()` | Save matrix + metadata to HDF5 |
| `write_medsd_lambda_tsv()` | Save per-dataset stats to TSV |
| `write_metadata_tsv()` | Save sample metadata to TSV |
| `add_analysis_group()` | Add `Analysis group` = tissue + cell_type |
| `filter_by_column_value()` | Keep/drop rows by column value |
| `deduplicate_by_column()` | Keep first row per unique value of a column |
| `filter_groups_by_column()` | Keep groups where all rows match a value |
| `filter_groups_mixed_sex()` | Keep groups with both M and F present |
| `sample_n_groups()` | Randomly sample N eligible groups |
| `sample_pairs()` | Sample MM / FF / MF pairs per group |
| `normalize_sex()` | Normalise sex strings to 'M' or 'F' |

---

## Adding a new analysis case

1. Add a function in `run_pipeline_atlas.sh`:
```bash
run_18_noImmune_maleOnly() {
  _run_atlas "18_noImmune_maleOnly" --exclude_immune --sex_filter M
}
```

2. Add to the `case "$PRESET"` dispatch block:
```bash
18_noImmune_maleOnly) run_18_noImmune_maleOnly ;;
```

No changes to any Python file required.

---

## Outputs (same for every preset)

| File | Description |
|---|---|
| `{prefix}_matrix_noscale.h5` | HDF5: `matrix`, `cpg_names`, `samples`, `sample_groups` |
| `sample_metadata.tsv` | Sample name and group per sample |
| `{prefix}_medsd_lambda.tsv` | Per-dataset `median_sd` and `lambda` (on final CpGs) |
| `{prefix}_cpg_names.txt` | Retained CpG site IDs |
| `logs/` | Per-step logs |

---

## Dependencies

```bash
pip install numpy pandas h5py bottleneck pyreadr
```
