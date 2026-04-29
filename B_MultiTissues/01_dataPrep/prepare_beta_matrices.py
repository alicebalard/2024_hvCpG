#!/usr/bin/env python3
"""
prepare_beta_matrices.py — Build filtered CpG beta matrix from WGBS .beta files
=================================================================================
Entry point for Loyfer-style WGBS atlas data (04_prepAtlas).

Pipeline
--------
  1. Load metadata (must already contain 'Analysis group' column — run
     prepare_metadata.py first).
  2. Select groups with ≥ --min_samples samples.
  3. Two-pass chunked matrix builder (see cpg_utils.build_beta_matrix_chunked):
       Pass 1 — count per-site dataset coverage (no I/O to final matrix).
       Site selection — keep CpGs covered in ≥ --min_datasets groups.
       Pass 2 — write filtered beta values into the final matrix.
  4. Compute per-dataset median SD and lambda on the FINAL filtered CpGs.
  5. Write outputs.

Outputs
-------
  <output_folder>/all_matrix_noscale.h5   — HDF5 matrix
  <output_folder>/sample_metadata.tsv
  <output_folder>/all_medsd_lambda.tsv    — statistics
  <output_folder>/cpg_names.txt

Usage
-----
python prepare_beta_matrices.py \\
  --beta_files "/data/betaFiles/GSM*.hg38.beta" \\
  --cpg_bed   /data/hg38/CpG.bed.gz \\
  --output_folder /data/output/10X \\
  --meta      /data/metadata_with_analysis_group.csv \\
  --minCov    10 \\
  --min_samples 3 \\
  --min_datasets 46 \\
  --chunk_size 100000

Author: Alice Balard
"""

import os
import glob
import argparse
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpg_utils import (
    load_cpg_names_from_bed,
    load_exclude_sites,
    build_exclude_mask,
    build_sample_to_path_map,
    filter_valid_groups,
    shorten_sample_ids,
    build_beta_matrix_chunked,
    write_h5_matrix,
    write_medsd_lambda_tsv,
    write_metadata_tsv,
)

# ──────────────────────────────────────────────
#  Arguments
# ──────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Build filtered CpG beta matrix from WGBS .beta files.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)

# Input
parser.add_argument("--beta_files",    required=True,
                    help="Glob pattern for .beta files, e.g. '/data/GSM*.hg38.beta'.")
parser.add_argument("--cpg_bed",       required=True,
                    help="CpG BED reference file (.bed or .bed.gz). Columns: chr, pos.")
parser.add_argument("--meta",          required=True,
                    help="Metadata CSV with 'Analysis group' and 'Sample name' columns.")

# Output
parser.add_argument("--output_folder", required=True,
                    help="Directory for output files (created if absent).")
parser.add_argument("--output_prefix", default="all",
                    help="Prefix for output filenames (default: 'all').")

# Filtering thresholds
parser.add_argument("--minCov",        type=int,   default=10,
                    help="Minimum read coverage to include a beta value (default: 10).")
parser.add_argument("--min_samples",   type=int,   default=3,
                    help="Minimum samples per group (default: 3).")
parser.add_argument("--min_datasets",  type=int,   default=46,
                    help="Minimum datasets a CpG must be covered in (default: 46).")
parser.add_argument("--lambda_percentile", type=float, default=95.0,
                    help="Upper percentile for lambda = perc/median (default: 95).")

# Performance
parser.add_argument("--chunk_size",    type=int,   default=100_000,
                    help="CpG rows per processing chunk (default: 100,000).")

# Column name overrides (for non-standard metadata)
parser.add_argument("--group_col",     default="Analysis group",
                    help="Metadata column for group labels.")
parser.add_argument("--sample_col",    default="Sample name",
                    help="Metadata column for sample names.")
parser.add_argument("--id_pattern",    default=r"-([A-Za-z0-9_]+)\.hg38\.beta$",
                    help=r"Regex to extract sample ID from beta filename. "
                         r"One capture group required. "
                         r"Default: r'-([A-Za-z0-9_]+)\.hg38\.beta$'")
parser.add_argument("--id_sep",        default="-",
                    help="Separator used to strip prefix from sample names (default: '-').")
parser.add_argument("--exclude_sites", default=None, metavar="FILE",
                    help="Optional text file listing CpG sites to exclude "
                         "(one 'chr_pos' identifier per line). Sites are removed "
                         "at Pass 1 — before any coverage counting — so they never "
                         "appear in the matrix or statistics.")

args = parser.parse_args()

# ──────────────────────────────────────────────
#  Setup
# ──────────────────────────────────────────────

os.makedirs(args.output_folder, exist_ok=True)

output_h5    = os.path.join(args.output_folder, "all_matrix_noscale.h5")
output_meta  = os.path.join(args.output_folder, "sample_metadata.tsv")
output_stats = os.path.join(args.output_folder, "all_medsd_lambda.tsv")
output_cpgs  = os.path.join(args.output_folder, "cpg_names.txt")

print("\n" + "="*62)
print("  prepare_beta_matrices.py")
print("="*62)
print(f"  beta_files    : {args.beta_files}")
print(f"  cpg_bed       : {args.cpg_bed}")
print(f"  meta          : {args.meta}")
print(f"  output_folder : {args.output_folder}")
print(f"  minCov        : {args.minCov}")
print(f"  min_samples   : {args.min_samples}")
print(f"  min_datasets  : {args.min_datasets}")
print(f"  chunk_size    : {args.chunk_size:,}")
print(f"  exclude_sites : {args.exclude_sites or '(none)'}")
print("="*62 + "\n")

# ──────────────────────────────────────────────
#  Load inputs
# ──────────────────────────────────────────────

beta_files = sorted(glob.glob(args.beta_files))
if not beta_files:
    raise FileNotFoundError(f"No beta files matched pattern: {args.beta_files}")
print(f"Found {len(beta_files):,} beta files.")

meta = pd.read_csv(args.meta)
for required in [args.group_col, args.sample_col]:
    if required not in meta.columns:
        raise ValueError(f"Metadata missing required column: '{required}'")

cpg_names = load_cpg_names_from_bed(args.cpg_bed)
nr_sites   = len(cpg_names)

# Build site-exclusion mask (None if no list provided)
if args.exclude_sites:
    excl_sites   = load_exclude_sites(args.exclude_sites)
    exclude_mask = build_exclude_mask(cpg_names, excl_sites)
else:
    exclude_mask = None

sample_to_path        = build_sample_to_path_map(beta_files, id_pattern=args.id_pattern)
samples_per_group     = filter_valid_groups(meta, args.group_col, args.sample_col, args.min_samples)
samples_per_group_short = shorten_sample_ids(samples_per_group, sep=args.id_sep)

# ──────────────────────────────────────────────
#  Core computation  (stats computed AFTER filtration inside this call)
# ──────────────────────────────────────────────

matrix, group_medians, group_lambdas, all_sample_names, all_sample_groups, final_mask = \
    build_beta_matrix_chunked(
        samples_per_group_short = samples_per_group_short,
        sample_to_path          = sample_to_path,
        nr_sites                = nr_sites,
        min_cov                 = args.minCov,
        min_samples_per_group   = args.min_samples,
        min_datasets            = args.min_datasets,
        chunk_size              = args.chunk_size,
        lambda_percentile       = args.lambda_percentile,
        exclude_mask            = exclude_mask,
    )

selected_cpgs = [cpg_names[i] for i, keep in enumerate(final_mask) if keep]

# ──────────────────────────────────────────────
#  Write outputs
# ──────────────────────────────────────────────

write_h5_matrix(output_h5, matrix, selected_cpgs, all_sample_names, all_sample_groups)
write_metadata_tsv(output_meta, all_sample_names, all_sample_groups)
write_medsd_lambda_tsv(output_stats, group_medians, group_lambdas)

with open(output_cpgs, "w") as f:
    f.write("\n".join(selected_cpgs) + "\n")
print(f"  Saved CpG name list → {output_cpgs}")

print(f"\n🎉 Done!")
print(f"   Matrix  : {matrix.shape[0]:,} CpGs × {matrix.shape[1]:,} samples")
print(f"   Output  : {args.output_folder}")
