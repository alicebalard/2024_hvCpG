# CpG Methylation Data Preparation Pipeline

Prepares WGBS atlas methylation data (`.beta` files) into HDF5 matrices for hvCpG analysis.
One metadata source, one runner, many presets.

---

## What it does (in order)

For each preset, the runner performs three steps:

1. **Filter metadata** (`prepare_metadata.py`) — takes the master table
   `SupTab1_Loyfer2023_amended.csv`, applies the preset's row/group/sampling
   filters, and writes a filtered sample list. Analysis groups are always
   `Source Tissue + " - " + Cell type`.

2. **Build the beta matrix** (`prepare_beta_matrices.py`) — reads the `.beta`
   files for the selected samples, applies coverage and CpG-exclusion masks,
   keeps CpGs covered in enough datasets, and writes an HDF5 matrix.

3. **Compute per-dataset stats** (`cpg_utils.py`) — median SD and lambda are
   computed **once, at the end**, on the final filtered CpG set only.

Everything runs with one command:

```bash
bash run_pipeline_atlas.sh --preset 13_meso
```

Outputs land in `${WGBS_DIR}/output_<preset>/` (see **Outputs** below).

---

## Prerequisites

```bash
pip install numpy pandas h5py bottleneck pyreadr
```

Before any SD-ASM-excluded run, build the SD-ASM blacklist once (see
**Exclusion → both**).

---

## Quick start

```bash
# default run (SNP-excluded, see below)
bash run_pipeline_atlas.sh --preset atlas_general

# override data dir / thresholds
DATA_DIR=/my/data       bash run_pipeline_atlas.sh --preset 12_endo
MIN_DATASETS_ATLAS=30   bash run_pipeline_atlas.sh --preset 06_bothsexes6gp

# submit to the cluster (pass env with -v)
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 13_meso
```

Batch submission lives in `toQsub.txt`.

---

# Details

## Presets

| Preset | What it does |
|---|---|
| `atlas_general` | All data, tissue–celltype grouping (baseline) |
| `02_rmMultSamples` | One sample per PatientID |
| `04_maleOnly` | Groups where all samples are male |
| `05_femaleOnly6gp` | Female-only groups, sampled to 6 |
| `06_bothsexes6gp` | Mixed-sex groups, sampled to 6 |
| `09_immuneOnly` | Immune cells only |
| `10_noImmune` | No immune cells |
| `11_noImmune_sample11gp` | No immune, sampled to 11 groups |
| `12_endo` / `12_2_endo6gp` | Endoderm (all / sampled to 6) |
| `13_meso` / `13_2_meso6gp` | Mesoderm (all / sampled to 6) |
| `14_ecto` | Ectoderm only |
| `15_pairs_MM` / `16_pairs_FF` / `17_pairs_MF` | 2 individuals per group (MM / FF / MF) |

---

## CpG exclusion (`EXCLUDE_MODE`)

CpGs can be masked before matrix building. Controlled by one variable:

| `EXCLUDE_MODE` | Excludes | Output folder suffix |
|---|---|---|
| `snp` *(default)* | SNP-overlapping CpGs (MAF > 0.01) | *(none)* → `output_<preset>` |
| `both` | SNPs **+** all CpGs in SD-ASM regions | `_noSDASM` → `output_<preset>_noSDASM` |
| `none` | nothing | *(none)* |

The suffix means a `both` run never overwrites the corresponding `snp` run —
the two sit side by side, so you can compare with vs without SD-ASM.

Source lists (override paths via `SNP_SITES` / `SDASM_SITES` if needed):

```bash
bash run_pipeline_atlas.sh --preset 13_meso                     # snp (default)
EXCLUDE_MODE=both bash run_pipeline_atlas.sh --preset 13_meso   # snp + SD-ASM
EXCLUDE_MODE=none bash run_pipeline_atlas.sh --preset 13_meso   # nothing
```

### Building the SD-ASM blacklist (once, before `both`)

SD-ASM data are *regions*; the pipeline excludes by *point* (`chr_pos`).
`SD-ASMprep.R` expands regions to the CpGs inside them, using the **same
`CpG.bed`** the pipeline uses (coordinates must match, or exclusion silently
does nothing):

```r
source("B_MultiTissues/01_dataPrep/SD-ASMprep.R")   # builds SDASM_GR
writeSDASM_blacklist(
  SDASM_GR,
  cpg_bed  = "/…/hg38/CpG.bed.gz",              # SAME as pipeline --cpg_bed
  out_file = "gitignore/sdasm_all_chr_pos.txt"  # default SDASM_SITES path
)
```

> **Downstream note:** for `both` runs, read results from
> `output_<preset>_noSDASM` — e.g. in R set `analysis = "13_meso_noSDASM"`.

---

## `prepare_metadata.py` — filters

Groups are always `Source Tissue + " - " + Cell type`.

**Row filters:** `--exclude_immune` · `--germ_layer Endo|Meso|Ecto` ·
`--dedup_col COL` · `--keep_col COL --keep_val VAL` · `--drop_col COL --drop_val VAL`

**Group filters:** `--sex_filter M|F` (all samples one sex) · `--mixed_sex`
(both present) · `--keep_groups_col COL --keep_groups_val VAL`

**Sampling:** `--sample_n_groups N` · `--min_per_group N` (default 3) ·
`--pairs MM|FF|MF` · `--n_groups N` · `--seed N` (default 42)

Preset → flags (representative examples; the rest follow the same pattern):

```bash
02_rmMultSamples   --dedup_col PatientID
04_maleOnly        --sex_filter M
06_bothsexes6gp    --mixed_sex --sample_n_groups 6
10_noImmune        --exclude_immune
13_meso            --germ_layer Meso
13_2_meso6gp       --germ_layer Meso --sample_n_groups 6
17_pairs_MF        --pairs MF
```

---

## `prepare_beta_matrices.py` — key arguments

| Argument | Default | Description |
|---|---|---|
| `--beta_files` | required | Glob for `.beta` files |
| `--cpg_bed` | required | CpG BED reference (`.bed`/`.bed.gz`) |
| `--meta` | required | Filtered metadata (from Step 1) |
| `--minCov` | 10 | Min read coverage per site |
| `--min_samples` | 3 | Min samples per group |
| `--min_datasets` | 46 | Min datasets per CpG |
| `--chunk_size` | 100000 | Rows per chunk |
| `--lambda_percentile` | 95 | lambda = percentile / median SD |
| `--exclude_sites` | — | `chr_pos` blacklist (applied before Pass 1) |

---

## Outputs (per preset, in `output_<preset>[_noSDASM]/`)

| File | Contents |
|---|---|
| `{prefix}_matrix_noscale.h5` | `matrix`, `cpg_names`, `samples`, `sample_groups` |
| `{prefix}_cpg_names.txt` | Retained CpG IDs |
| `{prefix}_medsd_lambda.tsv` | Per-dataset median SD and lambda (final CpGs) |
| `sample_metadata.tsv` | Sample → group |
| `logs/` | Per-step logs |

---

## Adding a preset

In `run_pipeline_atlas.sh`, add a function and a dispatch case — no Python
changes needed:

```bash
run_18_noImmune_maleOnly() { _run_atlas "18_noImmune_maleOnly" --exclude_immune --sex_filter M; }// then:
18_noImmune_maleOnly) run_18_noImmune_maleOnly ;;
```