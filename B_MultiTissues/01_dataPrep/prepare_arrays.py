#!/usr/bin/env python3
"""
prepare_arrays.py — Build filtered CpG beta matrix from Illumina array data
=============================================================================
Entry point for array-based methylation datasets (03_prepDatasetsMaria and
variants in 05_hvCpGalgorithm).

Supports two input formats (can be mixed in one run):
  --rds_folders   One or more directories containing .RDS files (R matrices
                  as produced by Maria's BMIQ + PCs + outlier-removal pipeline)
  --h5_files      Glob pattern for pre-existing .h5 files containing a beta
                  matrix stored under standard keys

Both formats are loaded, CpG sites are intersected (or restricted to a
reference BED if --cpg_bed is provided), and the combined matrix is filtered.

Pipeline
--------
  1. Load all RDS and/or H5 datasets.
  2. Determine the common CpG set (covered in ≥ --min_datasets datasets OR
     using a BED reference as the master list).
  3. Align every dataset to the common CpG order; missing sites → NaN.
  4. Apply per-site NA-fraction filter (--maxNA).
  5. Apply per-group minimum-samples filter.
  6. Build the merged matrix.
  7. Compute per-dataset median SD and lambda on the FINAL filtered CpGs.
  8. Write outputs.

No logit transform is applied — consistent with Maria's approach.

Outputs
-------
  <output_folder>/<prefix>_matrix_noscale.h5
  <output_folder>/sample_metadata.tsv
  <output_folder>/<prefix>_medsd_lambda.tsv
  <output_folder>/<prefix>_cpg_names.txt

Usage examples
--------------
# Load from two RDS folders (GEO + TCGA batches):
python prepare_arrays.py \\
  --rds_folders /data/GEO /data/TCGA \\
  --output_folder /data/output \\
  --min_datasets 15 \\
  --output_prefix all

# Load from H5 files:
python prepare_arrays.py \\
  --h5_files "/data/arrays/*.h5" \\
  --output_folder /data/output \\
  --output_prefix arrays_all

# Mix RDS and H5, restrict to a BED reference:
python prepare_arrays.py \\
  --rds_folders /data/GEO \\
  --h5_files "/data/extra/*.h5" \\
  --cpg_bed /data/hg19/CpG.bed.gz \\
  --output_folder /data/output

Author: Alice Balard
"""

import os
import glob
import argparse
import sys

import numpy as np
import pandas as pd
import h5py
import bottleneck as bn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpg_utils import (
    load_cpg_names_from_bed,
    load_exclude_sites,
    build_exclude_mask,
    load_rds_datasets,
    compute_stats_on_filtered_matrix,
    write_h5_matrix,
    write_medsd_lambda_tsv,
    write_metadata_tsv,
)

# ──────────────────────────────────────────────
#  Arguments
# ──────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Build filtered CpG beta matrix from Illumina array data (RDS or H5).",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)

# Input — at least one of --rds_folders, --rds_files, or --h5_files is required
parser.add_argument("--rds_folders", nargs="+", default=[],
                    help="One or more directories; ALL .RDS files inside are loaded.")
parser.add_argument("--rds_files",   nargs="+", default=[],
                    help="Explicit .RDS file paths. Use when you need a specific subset "
                         "(e.g. --rds_files /data/CD4+_Estonian.RDS /data/CD8+_Estonian.RDS).")
parser.add_argument("--h5_files",    default=None,
                    help="Glob pattern for .h5 files (e.g. '/data/*.h5').")
parser.add_argument("--cpg_bed",     default=None,
                    help="Optional CpG BED reference to define the master site list.")

# HDF5 key names (if input H5 uses non-standard keys)
parser.add_argument("--h5_matrix_key",  default="matrix",    help="H5 key for beta matrix.")
parser.add_argument("--h5_cpg_key",     default="cpg_names", help="H5 key for CpG names.")
parser.add_argument("--h5_sample_key",  default="samples",   help="H5 key for sample names.")

# Output
parser.add_argument("--output_folder",  required=True,
                    help="Directory for output files (created if absent).")
parser.add_argument("--output_prefix",  default="arrays",
                    help="Prefix for output filenames (default: 'arrays').")

# Filtering
parser.add_argument("--maxNA",          type=float, default=0.2,
                    help="Maximum fraction of NaN values per CpG per dataset (default: 0.2).")
parser.add_argument("--min_datasets",   type=int,   default=15,
                    help="Minimum number of datasets a CpG must appear in (default: 15).")
parser.add_argument("--lambda_percentile", type=float, default=95.0,
                    help="Upper percentile for lambda = perc/median (default: 95).")
parser.add_argument("--exclude_sites", default=None, metavar="FILE",
                    help="Optional text file listing CpG sites to exclude "
                         "(one 'chr_pos' identifier per line). Sites are removed "
                         "immediately after the common CpG set is determined — "
                         "before alignment and NA filtering — so they never appear "
                         "in the matrix or statistics.")
parser.add_argument("--max_samples_per_dataset", type=int, default=None, metavar="N",
                    help="Randomly subsample each dataset to at most N individuals "
                         "after loading and before CpG filtering. Datasets with <=N "
                         "samples are unchanged. Use --seed for reproducibility. "
                         "E.g. --max_samples_per_dataset 3 reproduces the 3ind scripts.")
parser.add_argument("--seed", type=int, default=42,
                    help="Random seed for --max_samples_per_dataset (default: 42).")
parser.add_argument("--no_extension",  action="store_true",
                    help="Load files from --rds_folders regardless of extension "
                         "(for raw-cleaned folders where files have no .RDS extension). "
                         "Subdirectories are always skipped.")
parser.add_argument("--exclude_subdir", nargs="+", default=[],
                    help="Subdirectory/file names to skip when using --no_extension "
                         "(e.g. --exclude_subdir CHAMP_Normalization).")

args = parser.parse_args()

if not args.rds_folders and not args.rds_files and not args.h5_files:
    parser.error("Provide at least one of --rds_folders or --h5_files.")

# ──────────────────────────────────────────────
#  Setup
# ──────────────────────────────────────────────

os.makedirs(args.output_folder, exist_ok=True)
p = args.output_prefix

output_h5    = os.path.join(args.output_folder, f"{p}_matrix_noscale.h5")
output_meta  = os.path.join(args.output_folder, "sample_metadata.tsv")
output_stats = os.path.join(args.output_folder, f"{p}_medsd_lambda.tsv")
output_cpgs  = os.path.join(args.output_folder, f"{p}_cpg_names.txt")

print("\n" + "="*62)
print("  prepare_arrays.py")
print("="*62)
print(f"  rds_folders   : {args.rds_folders}")
print(f"  rds_files     : {args.rds_files}")
print(f"  h5_files      : {args.h5_files}")
print(f"  max_samples   : {args.max_samples_per_dataset or '(all)'}")
print(f"  cpg_bed       : {args.cpg_bed}")
print(f"  output_folder : {args.output_folder}")
print(f"  maxNA         : {args.maxNA}")
print(f"  min_datasets  : {args.min_datasets}")
print(f"  exclude_sites : {args.exclude_sites or '(none)'}")
print("="*62 + "\n")

# ══════════════════════════════════════════════════════════════════════════════
#  1. Load all datasets
#     Each entry: (name: str, matrix: np.ndarray float32, cpg_names: list)
# ══════════════════════════════════════════════════════════════════════════════

all_datasets = []   # list of (name, matrix, cpg_names)

# — RDS files ——————————————————————————————————————————————————————————————
for folder in args.rds_folders:
    if args.no_extension:
        # Load all files regardless of extension; skip subdirs and excluded names.
        excl_names = set(args.exclude_subdir)
        rds_files = sorted([
            os.path.join(folder, fn)
            for fn in os.listdir(folder)
            if os.path.isfile(os.path.join(folder, fn)) and fn not in excl_names
        ])
    else:
        rds_files = sorted(glob.glob(os.path.join(folder, "*.RDS")) +
                           glob.glob(os.path.join(folder, "*.rds")))
    if not rds_files:
        print(f"  Warning: no RDS files found in: {folder}")
        continue
    print(f"\n  Loading {len(rds_files)} RDS file(s) from: {folder}")
    all_datasets.extend(load_rds_datasets(rds_files))

# — Explicit RDS file list (--rds_files) ————————————————————————————————
if args.rds_files:
    present = [f for f in args.rds_files if os.path.isfile(f)]
    missing = [f for f in args.rds_files if not os.path.isfile(f)]
    for m in missing:
        print(f"  Warning: explicit RDS file not found: {m}")
    if present:
        print(f"\n  Loading {len(present)} explicit RDS file(s).")
        all_datasets.extend(load_rds_datasets(present))

# — H5 files ———————————————————————————————————————————————————————————————
if args.h5_files:
    h5_paths = sorted(glob.glob(args.h5_files))
    if not h5_paths:
        print(f"  ⚠️  No H5 files matched: {args.h5_files}")
    else:
        print(f"\n  Loading {len(h5_paths)} H5 file(s).")
        for path in h5_paths:
            name = os.path.splitext(os.path.basename(path))[0]
            with h5py.File(path, "r") as f:
                mat      = f[args.h5_matrix_key][:].astype(np.float32)
                raw_cpg  = f[args.h5_cpg_key][:]
                raw_smp  = f[args.h5_sample_key][:]
            cpg_list = [s.decode("utf-8") if isinstance(s, bytes) else s for s in raw_cpg]
            smp_list = [s.decode("utf-8") if isinstance(s, bytes) else s for s in raw_smp]
            # Store sample names on the dataset for bookkeeping
            all_datasets.append((name, mat, cpg_list, smp_list))
            print(f"    Loaded {name}: {mat.shape[0]:,} CpGs × {mat.shape[1]:,} samples")

if not all_datasets:
    raise RuntimeError("No datasets were loaded. Check your input paths.")

print(f"\n  Total datasets loaded: {len(all_datasets)}")

# ══════════════════════════════════════════════════════════════════════════════
#  2. Normalise dataset tuples
#     RDS datasets have no pre-existing sample list → auto-generate names.
#     Ensure all tuples are (name, matrix, cpg_names, sample_names).
# ══════════════════════════════════════════════════════════════════════════════

normalised = []
for entry in all_datasets:
    if len(entry) == 3:
        name, mat, cpg_names = entry
        # Auto-generate sample names: datasetname_0, datasetname_1, ...
        smp_names = [f"{name}_{i}" for i in range(mat.shape[1])]
    else:
        name, mat, cpg_names, smp_names = entry
    normalised.append((name, mat, cpg_names, smp_names))

# Per-dataset sample cap (--max_samples_per_dataset)
# Applied after loading but before CpG intersection so coverage counts
# and statistics reflect only the retained samples.
if args.max_samples_per_dataset is not None:
    rng    = np.random.default_rng(args.seed)
    capped = []
    for name, mat, cpg_names, smp_names in normalised:
        n, cap = mat.shape[1], args.max_samples_per_dataset
        if n > cap:
            idx       = rng.choice(n, size=cap, replace=False)
            mat       = mat[:, idx]
            smp_names = [smp_names[i] for i in idx]
            print(f"  {name}: subsampled {n} -> {cap} samples.")
        capped.append((name, mat, cpg_names, smp_names))
    normalised = capped
    print(f"  max_samples_per_dataset={args.max_samples_per_dataset}: "
          f"total samples after cap: {sum(len(s) for _,_,_,s in normalised):,}")

# ══════════════════════════════════════════════════════════════════════════════
#  3. Determine common CpG set
# ══════════════════════════════════════════════════════════════════════════════

if args.cpg_bed:
    # Reference BED defines the master list; intersect to only keep sites
    # actually present in every dataset.
    ref_cpgs   = load_cpg_names_from_bed(args.cpg_bed)
    ref_set    = set(ref_cpgs)
    # Count presence across datasets
    cpg_counts = {}
    for _, _, cpg_names, _ in normalised:
        for c in set(cpg_names):
            if c in ref_set:
                cpg_counts[c] = cpg_counts.get(c, 0) + 1
    common_cpgs = [c for c in ref_cpgs if cpg_counts.get(c, 0) >= args.min_datasets]
    print(f"\n  Reference BED: {len(ref_cpgs):,} sites")
    print(f"  After min_datasets={args.min_datasets} filter: {len(common_cpgs):,} sites")
else:
    # Count CpG presence across all datasets
    cpg_counts = {}
    for _, _, cpg_names, _ in normalised:
        for c in set(cpg_names):
            cpg_counts[c] = cpg_counts.get(c, 0) + 1
    # Use ordering from the first dataset that contains each CpG
    first_cpg_order = []
    seen = set()
    for _, _, cpg_names, _ in normalised:
        for c in cpg_names:
            if c not in seen:
                first_cpg_order.append(c)
                seen.add(c)
    common_cpgs = [c for c in first_cpg_order if cpg_counts[c] >= args.min_datasets]
    print(f"\n  CpGs covered in ≥{args.min_datasets} datasets: {len(common_cpgs):,}")

if not common_cpgs:
    raise RuntimeError(
        f"No CpGs survived the min_datasets={args.min_datasets} filter. "
        "Lower --min_datasets or check your input files."
    )

# Remove excluded sites from the common CpG list at the earliest possible
# point — before alignment and NA filtering — so they never appear in the
# matrix or statistics.
if args.exclude_sites:
    excl_sites  = load_exclude_sites(args.exclude_sites)
    n_before    = len(common_cpgs)
    common_cpgs = [c for c in common_cpgs if c not in excl_sites]
    print(f"  Exclusion list: removed {n_before - len(common_cpgs):,} sites "
          f"-> {len(common_cpgs):,} remain.")
    if not common_cpgs:
        raise RuntimeError("No CpGs remain after applying the exclusion list.")

cpg_index_map = {c: i for i, c in enumerate(common_cpgs)}
n_sites = len(common_cpgs)

# ══════════════════════════════════════════════════════════════════════════════
#  4–6. Align, filter, and merge into final matrix
# ══════════════════════════════════════════════════════════════════════════════

print("\n  ── Aligning datasets and applying NA filter ──")

all_sample_names:  list        = []
all_sample_groups: list        = []
matrix_blocks:     list        = []   # list of (n_sites, n_samples_in_dataset) arrays

for name, mat, cpg_names, smp_names in normalised:

    # Build row-index array to align this dataset to common_cpgs
    ds_cpg_map = {c: i for i, c in enumerate(cpg_names)}
    row_idx    = np.array(
        [ds_cpg_map.get(c, -1) for c in common_cpgs], dtype=np.int64
    )

    # Aligned matrix: common_cpgs × n_samples; missing rows → NaN
    aligned = np.full((n_sites, mat.shape[1]), np.nan, dtype=np.float32)
    present = row_idx >= 0
    aligned[present, :] = mat[row_idx[present], :]

    # Apply per-CpG NA-fraction filter
    na_frac = np.mean(np.isnan(aligned), axis=1)   # fraction of samples that are NaN
    aligned[na_frac > args.maxNA, :] = np.nan

    matrix_blocks.append(aligned)
    all_sample_names.extend(smp_names)
    all_sample_groups.extend([name] * len(smp_names))

    n_valid = int(np.sum(~np.all(np.isnan(aligned), axis=1)))
    print(f"    {name}: {aligned.shape[1]} samples, "
          f"{n_valid:,}/{n_sites:,} CpGs pass maxNA={args.maxNA}")

# Merge horizontally
print("\n  Merging all blocks...")
final_matrix = np.hstack(matrix_blocks).astype(np.float32)   # (n_sites, n_total_samples)

# ══════════════════════════════════════════════════════════════════════════════
#  7. Statistics on the final filtered matrix
#     This is correct — computed AFTER all filtration steps.
# ══════════════════════════════════════════════════════════════════════════════

print("\n  ── Computing median SD and lambda on final filtered CpGs ──")
group_medians, group_lambdas = compute_stats_on_filtered_matrix(
    final_matrix, all_sample_groups, args.lambda_percentile
)

# ══════════════════════════════════════════════════════════════════════════════
#  8. Write outputs
# ══════════════════════════════════════════════════════════════════════════════

write_h5_matrix(output_h5, final_matrix, common_cpgs, all_sample_names, all_sample_groups)
write_metadata_tsv(output_meta, all_sample_names, all_sample_groups)
write_medsd_lambda_tsv(output_stats, group_medians, group_lambdas)

with open(output_cpgs, "w") as f:
    f.write("\n".join(common_cpgs) + "\n")
print(f"  Saved CpG name list → {output_cpgs}")

print(f"\n🎉 Done!")
print(f"   Matrix  : {final_matrix.shape[0]:,} CpGs × {final_matrix.shape[1]:,} samples")
print(f"   Output  : {args.output_folder}")
