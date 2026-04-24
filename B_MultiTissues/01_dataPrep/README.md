# CpG Methylation Data Preparation Pipeline

Modular pipeline for preparing CpG methylation data for hvCpG analysis.
Handles WGBS atlas data (`.beta` files) and Illumina array data (`.RDS` / `.h5`).

---

## File structure

```
cpg_pipeline/
├── cpg_utils.py               ← Shared core library
├── prepare_metadata.py        ← Step 1: annotate & filter metadata
├── prepare_beta_matrices.py   ← Step 2a: WGBS .beta files → HDF5
├── prepare_arrays.py          ← Step 2b: RDS / H5 arrays → HDF5
├── run_pipeline.sh            ← One-command runner (18 presets)
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
# Run any of the 18 analysis cases:
bash run_pipeline.sh --preset 04_maleOnly
bash run_pipeline.sh --preset 10_noImmune
bash run_pipeline.sh --preset 17_pairs_MF

# Override data directory or thresholds:
DATA_DIR=/my/data bash run_pipeline.sh --preset 12_endo
MIN_DATASETS_ATLAS=30 bash run_pipeline.sh --preset 06_bothsexes6gp

# Exclude a blacklist of CpG sites (applies to all presets):
EXCLUDE_SITES=/data/snp_cpgs.txt bash run_pipeline.sh --preset 10_noImmune

# Use a config file:
bash run_pipeline.sh --config my_config.sh
```

---

## All presets — original script mapping

| Preset | Original script | What it does |
|---|---|---|
| `01_byDevLayer` | `prepare_metadata_1_byDevLayer.py` | Group by Germ layer column |
| `02_rmMultSamples` | `prepare_metadata_2_rmMultSamples.py` | One sample per PatientID |
| `03_correspMariaTissues` | `prepare_metadata_3_correspMariaTissues.py` | Keep only tissues in Maria's datasets |
| `04_maleOnly` | `prepare_metadata_4_maleOnly6gp.py` | Groups where all samples are male |
| `05_femaleOnly6gp` | `prepare_metadata_5_femaleOnly6gp.py` | Female-only groups, sample 6 |
| `06_bothsexes6gp` | `prepare_metadata_6_bothsexes6gp.py` | Mixed-sex groups, sample 6 |
| `08_byTissue` | `prepare_metadata_8_byTissue.py` | Group by "Group simplified" column |
| `09_immuneOnly` | `prepare_metadata_9_immuneOnly.py` | Immune cells only |
| `10_noImmune` | `prepare_metadata_10_noImmune.py` | No immune cells |
| `11_noImmune_sample11gp` | `prepare_metadata_11_noImmune_sample11groups.py` | No immune, sample 11 groups |
| `12_endo` | `prepare_metadata_12_endo.py` | Endoderm only |
| `12_2_endo6gp` | `prepare_metadata_12_2_endo6gp.py` | Endoderm, sample 6 groups |
| `13_meso` | `prepare_metadata_13_meso.py` | Mesoderm only |
| `13_2_meso6gp` | `prepare_metadata_13_2_meso6gp.py` | Mesoderm, sample 6 groups |
| `14_ecto` | `prepare_metadata_14_ecto.py` | Ectoderm only |
| `15_pairs_MM` | `prepare_metadata_15_pairs_MM.py` | Sample 2 males per group |
| `16_pairs_FF` | `prepare_metadata_16_pairs_FF.py` | Sample 2 females per group |
| `17_pairs_MF` | `prepare_metadata_17_pairs_MF.py` | Sample 1M + 1F per group |
| `arrays_all` | (Maria RDS pipeline) | All Maria array datasets |

---

## `prepare_metadata.py` — all options

### Analysis group construction

| Argument | Description |
|---|---|
| `--group_col COL` | Use a single column directly as group label (e.g. `"Germ layer"` for script 1, `"Group simplified"` for script 8). Mutually exclusive with `--tissue_col`. |
| `--tissue_col COL` | Tissue column for `tissue + " - " + cell_type` formula (default: `"Source Tissue"`). |
| `--cell_col COL` | Cell-type column (default: `"Cell type"`). |

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

### Equivalent flag combinations for each original script

```bash
# Script 1 — byDevLayer
python prepare_metadata.py --meta SupTab1.csv --group_col "Germ layer"

# Script 2 — rmMultSamples
python prepare_metadata.py --meta SupTab1.csv --dedup_col "PatientID"

# Script 3 — correspMariaTissues
python prepare_metadata.py --meta SupTab1.csv \
    --keep_col "Found in Maria DS?" --keep_val True

# Script 4 — maleOnly
python prepare_metadata.py --meta SupTab1.csv --sex_filter M

# Script 5 — femaleOnly6gp
python prepare_metadata.py --meta SupTab1.csv --sex_filter F --sample_n_groups 6

# Script 6 — bothsexes6gp
python prepare_metadata.py --meta SupTab1.csv --mixed_sex --sample_n_groups 6

# Script 8 — byTissue
python prepare_metadata.py --meta SupTab1.csv --group_col "Group simplified"

# Script 9 — immuneOnly
python prepare_metadata.py --meta SupTab1.csv --keep_col "Immune?" --keep_val True

# Script 10 — noImmune
python prepare_metadata.py --meta SupTab1.csv --exclude_immune

# Script 11 — noImmune_sample11groups
python prepare_metadata.py --meta SupTab1.csv \
    --exclude_immune --sample_n_groups 11

# Script 12 — endo
python prepare_metadata.py --meta SupTab1.csv --germ_layer Endo

# Script 12.2 — endo6gp
python prepare_metadata.py --meta SupTab1.csv \
    --germ_layer Endo --sample_n_groups 6

# Script 13 — meso
python prepare_metadata.py --meta SupTab1.csv --germ_layer Meso

# Script 13.2 — meso6gp
python prepare_metadata.py --meta SupTab1.csv \
    --germ_layer Meso --sample_n_groups 6

# Script 14 — ecto
python prepare_metadata.py --meta SupTab1.csv --germ_layer Ecto

# Script 15 — pairs_MM
python prepare_metadata.py --meta SupTab1.csv --pairs MM

# Script 16 — pairs_FF
python prepare_metadata.py --meta SupTab1.csv --pairs FF

# Script 17 — pairs_MF
python prepare_metadata.py --meta SupTab1.csv --pairs MF
```

---

## Step-by-step usage

### Step 2a — WGBS atlas (`.beta` files)

```bash
python prepare_beta_matrices.py \
  --beta_files "/data/betaFiles/GSM*.hg38.beta" \
  --cpg_bed   /data/hg38/CpG.bed.gz \
  --output_folder /data/output/10X \
  --meta      meta_all.csv \
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
| `--meta` | required | Filtered metadata CSV |
| `--minCov` | 10 | Min read coverage per site |
| `--min_samples` | 3 | Min samples per group |
| `--min_datasets` | 46 | Min datasets per CpG |
| `--chunk_size` | 100000 | Rows per processing chunk |
| `--id_pattern` | `r'-([A-Za-z0-9_]+)\.hg38\.beta$'` | Regex to extract sample ID from filename |
| `--lambda_percentile` | 95 | Percentile for lambda = perc/median |
| `--exclude_sites` | — | Text file of `chr_pos` sites to exclude (applied before Pass 1) |

### Step 2b — Illumina arrays (`.RDS` or `.h5` files)

```bash
# From RDS files (GEO + TCGA batches):
python prepare_arrays.py \
  --rds_folders /data/GEO /data/TCGA \
  --output_folder /data/output \
  --min_datasets 15 \
  --output_prefix all \
  --exclude_sites /data/snp_blacklist.txt   # optional

# From H5 files, with a BED reference:
python prepare_arrays.py \
  --h5_files "/data/arrays/*.h5" \
  --cpg_bed /data/hg19/CpG.bed.gz \
  --output_folder /data/output
```

| Argument | Default | Description |
|---|---|---|
| `--rds_folders` | — | One or more dirs containing `.RDS` files |
| `--h5_files` | — | Glob for `.h5` files (can combine with RDS) |
| `--cpg_bed` | optional | Reference BED as master site list |
| `--maxNA` | 0.2 | Max NaN fraction per CpG per dataset |
| `--min_datasets` | 15 | Min datasets per CpG |
| `--lambda_percentile` | 95 | Percentile for lambda |
| `--exclude_sites` | — | Text file of `chr_pos` sites to exclude (applied after CpG intersection, before alignment) |

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

The exclusion is applied **as early as possible** in each pipeline:

| Script | When exclusion is applied |
|---|---|
| `prepare_beta_matrices.py` | Before Pass 1 — excluded sites never accumulate coverage counts |
| `prepare_arrays.py` | After CpG intersection, before alignment and NA filtering |

To apply the same exclusion list to all presets via `run_pipeline.sh`:

```bash
EXCLUDE_SITES=/data/snp_blacklist.txt bash run_pipeline.sh --preset 10_noImmune
```

Or add it permanently to a config file:

```bash
# my_config.sh
EXCLUDE_SITES=/data/snp_blacklist.txt
MIN_DATASETS_ATLAS=30
```

```bash
bash run_pipeline.sh --config my_config.sh --preset 10_noImmune
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
| `load_rds_datasets()` | Load list of `.RDS` files via pyreadr |
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
| `set_analysis_group_from_column()` | Set `Analysis group` from a single column |
| `filter_by_column_value()` | Keep/drop rows by column value |
| `deduplicate_by_column()` | Keep first row per unique value of a column |
| `filter_groups_by_column()` | Keep groups where all rows match a value |
| `filter_groups_mixed_sex()` | Keep groups with both M and F present |
| `sample_n_groups()` | Randomly sample N eligible groups |
| `sample_pairs()` | Sample MM / FF / MF pairs per group |
| `normalize_sex()` | Normalise sex strings to 'M' or 'F' |

---

## Adding a new analysis case

1. Add a function in `run_pipeline.sh`:
```bash
run_15_2_pairs_MM_noImmune() {
  _run_atlas "--exclude_immune --pairs MM" "15_2_pairs_MM_noImmune"
}
```

2. Add to the `case "$PRESET"` dispatch block:
```bash
15_2_pairs_MM_noImmune) run_15_2_pairs_MM_noImmune ;;
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
