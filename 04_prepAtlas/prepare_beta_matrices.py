#!/usr/bin/env python3
"""
Prepare Loyfer WGBS Atlas Data (No Scaling)
---------------------------------
This script:
- Filters cell/tissue groups with at least MIN_SAMPLES_PER_GROUP samples
- Loads beta files (.beta)
- Builds a CpG-by-sample matrix for each group
- Masks coverage < minCov
- Saves: (a) filtered matrix, (b) median of per-CpG row SD per dataset, (c) lambda per dataset, (d) CpG names
- Selects CpGs covered in ≥ MIN_DATASETS datasets
Author: Alice Balard
"""

import os
import glob
import numpy as np
import pandas as pd
import re
import h5py
import bottleneck as bn
import gzip
import argparse

###############
## Arguments ##
###############

parser = argparse.ArgumentParser(description="Prepare Loyfer WGBS Atlas Data")
parser.add_argument("--beta_files", required=True, help="Files pattern for beta files")
parser.add_argument("--cpg_bed", required=True, help="Path to CpG bed file (chr pos extracols)")
parser.add_argument("--output_folder", required=True, help="Output folder")
parser.add_argument("--meta", required=True, help="Metadata CSV file")
parser.add_argument("--minCov", type=int, default=10, help="Minimum coverage")
parser.add_argument("--min_samples", type=int, default=3, help="Minimum samples per group")
parser.add_argument("--min_datasets", type=int, default=46, help="Minimum datasets per CpG")
parser.add_argument("--chunk_size", type=int, default=100000, help="(Unused in optimized write) Row chunk size if needed")
args = parser.parse_args()

# Resolve args
beta_files = glob.glob(args.beta_files)
cpg_bed = args.cpg_bed
output_folder = args.output_folder
meta = pd.read_csv(args.meta)
minCov = args.minCov
MIN_SAMPLES_PER_GROUP = args.min_samples
MIN_DATASETS = args.min_datasets
CHUNK_SIZE = args.chunk_size  # kept for compatibility, not used in column-wise write

#################
## Main script ##
#################

os.chdir(output_folder)

output_file = "all_matrix_noscale.h5"
metadata_file = "sample_metadata.tsv"
output_file_medsd_lambda = "all_medsd_lambda.tsv"

# Filter groups with enough samples
group_counts = meta.groupby("Analysis group").size()
valid_groups = group_counts[group_counts >= MIN_SAMPLES_PER_GROUP].index.tolist()
print(f"Found {len(valid_groups)} groups with >= {MIN_SAMPLES_PER_GROUP} samples.")

samples_per_group = {
    g: meta.loc[meta["Analysis group"] == g, "Sample name"].tolist()
    for g in valid_groups
}
samples_per_group_short = {
    g: [s.split("-")[-1] for s in samples]
    for g, samples in samples_per_group.items()
}

# Build a fast map from short sample id to beta file path (avoid repeated regex scanning)
sample_to_path = {}
dup_count = 0
for fn in beta_files:
    m = re.search(r"-([A-Za-z0-9_]+)\.hg38\.beta$", os.path.basename(fn))
    if m:
        key = m.group(1)
        if key in sample_to_path:
            dup_count += 1
        sample_to_path[key] = fn
if dup_count:
    print(f"⚠️ Warning: {dup_count} duplicate sample IDs in beta files; last one wins.")

# Load CpG names
cpg_names_start = []
with gzip.open(cpg_bed, "rt") as f:  # "rt" = read text mode
    for line in f:
        chrom, pos = line.strip().split("\t")[:2]
        pos = int(pos)
        cpg_names_start.append(f"{chrom}_{pos}")

NR_SITES = len(cpg_names_start)
print(f"✅ Loaded {NR_SITES:,} CpG names.")

# Helper: Load beta + coverage
def load_beta(path):
    arr = np.fromfile(path, dtype=np.uint8).reshape((-1, 2))
    meth = arr[:, 0]
    cov = arr[:, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
    return beta, cov

# Pass 1: Compute coverage counts (no HDF5 writes)
dataset_pass_count = np.zeros(NR_SITES, dtype=np.uint16)
sample_names = []
sample_groups = []
all_sample_names = []
all_sample_groups = []
group_medians = {}
group_lambdas = {}
sample_idx = 0

print("Pass 1: Computing coverage counts...")
for group, samples in samples_per_group_short.items():
    group_betas = []
    group_added = 0
    for s in samples:
        path = sample_to_path.get(s)
        if path is None:
            print(f"⚠️ No beta file found for short sample ID: {s}")
            continue
        beta, cov = load_beta(path)
        if len(beta) != NR_SITES:
            raise ValueError(f"Mismatch: {s} has {len(beta)} CpGs, expected {NR_SITES}")
        beta[cov < minCov] = np.nan
        group_betas.append(beta)

        # Track sample info (strings stored as UTF-8 later)
        sample_names.append(s)
        sample_groups.append(group)
        all_sample_names.append(s)
        all_sample_groups.append(group)
        sample_idx += 1
        group_added += 1

    if group_added == 0:
        print(f"⚠️ Skipping {group}: no valid samples")
        continue

    # Stack CpGs x group_samples
    mat = np.column_stack(group_betas).astype(np.float32)

    # Enforce per-group coverage rule
    valid_counts = np.sum(~np.isnan(mat), axis=1)
    mat[valid_counts < MIN_SAMPLES_PER_GROUP, :] = np.nan

    # Update global dataset pass counts
    coverage_counts = np.sum(~np.isnan(mat), axis=1)
    dataset_pass_count[coverage_counts >= MIN_SAMPLES_PER_GROUP] += 1

    # Group SD stats
    row_sds = bn.nanstd(mat, axis=1)
    if np.all(np.isnan(row_sds)):
        median_sd = np.nan
        lambda_value = np.nan
    else:
        median_sd = float(np.nanmedian(row_sds))
        perc95 = float(np.nanpercentile(row_sds, 95))
        lambda_value = (perc95 / median_sd) if (median_sd and np.isfinite(median_sd)) else np.nan

    group_medians[group] = median_sd
    group_lambdas[group] = lambda_value
    print(f"✅ {group}: median_sd = {median_sd:.4f}, lambda = {lambda_value:.4f}")

# Apply CpG selection
final_mask = dataset_pass_count >= MIN_DATASETS
selected_cpgs = np.array(cpg_names_start)[final_mask]
n_sel = int(np.sum(final_mask))
n_cols = sample_idx
print(f"✅ Selected {n_sel:,} CpGs covered in >= {MIN_SAMPLES_PER_GROUP} samples in >= {MIN_DATASETS} datasets.")

# Pass 2: Write filtered matrix directly to HDF5 (column-wise; load each beta once)
output_path = os.path.join(output_folder, output_file)
with h5py.File(output_path, "w") as h5f:
    matrix_dset = h5f.create_dataset(
        "matrix",
        shape=(n_sel, n_cols),
        dtype="float32",
        compression="gzip"
    )

    print("Pass 2: Writing filtered matrix (column-wise)...")
    col_idx = 0
    for group, samples in samples_per_group_short.items():
        for s in samples:
            path = sample_to_path.get(s)
            if path is None:
                # keep column as NaN (already initialized by HDF5 zeros; set to NaN explicitly)
                matrix_dset[:, col_idx] = np.full(n_sel, np.nan, dtype=np.float32)
                col_idx += 1
                continue
            beta, cov = load_beta(path)
            if len(beta) != NR_SITES:
                raise ValueError(f"Mismatch on write: {s} has {len(beta)} CpGs, expected {NR_SITES}")
            beta[cov < minCov] = np.nan
            # write only selected rows to this column
            matrix_dset[:, col_idx] = beta[final_mask]
            col_idx += 1

    # Save metadata
    h5f.create_dataset("cpg_names", data=selected_cpgs.astype("S"))
    h5f.create_dataset("samples", data=[str(s) for s in sample_names], dtype=h5py.string_dtype(encoding='utf-8'))
    h5f.create_dataset("sample_groups", data=[str(g) for g in sample_groups], dtype=h5py.string_dtype(encoding='utf-8'))
    
print("\n✅ Validation:")
print(f"Matrix shape: {n_sel:,} CpGs x {n_cols:,} samples")
print(f"Saved filtered matrix to: {output_path}")

# Save metadata TSV
meta_df = pd.DataFrame({"sample": all_sample_names, "dataset": all_sample_groups})
meta_df.to_csv(metadata_file, sep="\t", index=False)
print(f"✅ Saved metadata to: {metadata_file}")

# Save medians and lambdas
df = pd.DataFrame({
    "dataset": list(group_medians.keys()),
    "median_sd": list(group_medians.values()),
    "lambda": [group_lambdas[k] for k in group_medians.keys()]
})
df.to_csv(output_file_medsd_lambda, sep="\t", index=False)
print(f"✅ Saved medians and lambdas to TSV: {output_file_medsd_lambda}")

print("\n🎉 All done.")
