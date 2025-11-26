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

###############
## Arguments ##
###############

import argparse

parser = argparse.ArgumentParser(description="Prepare Loyfer WGBS Atlas Data")
parser.add_argument("--beta_files", required=True, help="Files pattern for beta files")
parser.add_argument("--cpg_bed", required=True, help="Path to CpG bed file (chr pos extracols)")
parser.add_argument("--output_folder", required=True, help="Output folder")
parser.add_argument("--meta", required=True, help="Metadata CSV file")
parser.add_argument("--minCov", type=int, default=10, help="Minimum coverage")
parser.add_argument("--min_samples", type=int, default=3, help="Minimum samples per group")
parser.add_argument("--min_datasets", type=int, default=46, help="Minimum datasets per CpG")
parser.add_argument("--chunk_size", type=int, default=100000, help="Chunk size for processing")

args = parser.parse_args()

# Then replace hardcoded variables with args.*
beta_files = glob.glob(args.beta_files)
cpg_bed = args.cpg_bed
output_folder = args.output_folder
meta = pd.read_csv(args.meta)
minCov = args.minCov
MIN_SAMPLES_PER_GROUP = args.min_samples
MIN_DATASETS = args.min_datasets
CHUNK_SIZE = args.chunk_size

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

# Load CpG names
cpg_names = []
with gzip.open(cpg_bed, "rt") as f: # "rt" = read text mode
    for line in f:
        chrom, pos = line.strip().split("\t")[:2]
        pos = int(pos)
        cpg_names.append(f"{chrom}_{pos}")

NR_SITES = len(cpg_names)
print(f"Loaded {NR_SITES:,} CpG names.")

# Helper: Load beta + coverage
def load_beta(path):
    arr = np.fromfile(path, dtype=np.uint8).reshape((-1, 2))
    meth = arr[:, 0]
    cov = arr[:, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
    return beta, cov

# Prepare HDF5
all_valid_samples = [s for group in samples_per_group_short.values() for s in group]
output_path = os.path.join(output_folder, output_file)
h5f = h5py.File(output_path, "w")

matrix_dset = h5f.create_dataset(
    "matrix",
    shape=(NR_SITES, len(all_valid_samples)),
    dtype="float32",
    compression="gzip"
)

# Track CpG coverage across datasets
dataset_pass_count = np.zeros(NR_SITES, dtype=np.uint16)

sample_names = []
sample_groups = []
sample_idx = 0
all_sample_names = []
all_sample_groups = []
group_medians = {}
group_lambdas = {}

# Process groups
for group, samples in samples_per_group_short.items():
    print(f"🔄 Processing {group} ({len(samples)} samples)")
    group_start_idx = sample_idx
    for s in samples:
        matches = [fn for fn in beta_files if re.search(f"-{s}\\.hg38\\.beta$", fn)]
        if not matches:
            print(f"⚠️ No beta file found for: {s}")
            continue
        beta, cov = load_beta(matches[0])
        if len(beta) != NR_SITES:
            raise ValueError(f"Mismatch: {s} has {len(beta)} CpGs, expected {NR_SITES}")
        beta[cov < minCov] = np.nan
        matrix_dset[:, sample_idx] = beta
        sample_names.append(s.encode())
        sample_groups.append(group.encode())
        all_sample_names.append(s)
        all_sample_groups.append(group)
        sample_idx += 1

    # Group-level stats
    group_sample_indices = list(range(group_start_idx, sample_idx))
    if len(group_sample_indices) == 0:
        print(f"!!Skipping {group}: no valid samples")
        continue
    mat = matrix_dset[:, group_sample_indices]
    valid_counts = np.sum(~np.isnan(mat), axis=1)
    mat[valid_counts < MIN_SAMPLES_PER_GROUP, :] = np.nan
    matrix_dset[:, group_sample_indices] = mat

    # Update CpG coverage tracker
    coverage_counts = np.sum(~np.isnan(mat), axis=1)
    mask = coverage_counts >= MIN_SAMPLES_PER_GROUP
    dataset_pass_count[mask] += 1

    # Compute SD stats
    row_sds = bn.nanstd(mat, axis=1)
    median_sd = np.nanmedian(row_sds)
    percentile_95 = np.nanpercentile(row_sds, 95)
    lambda_value = percentile_95 / median_sd
    group_medians[group] = median_sd
    group_lambdas[group] = lambda_value
    print(f"{group}: median_sd = {median_sd:.4f}, lambda = {lambda_value:.4f}")

# Apply CpG selection
final_mask = dataset_pass_count >= MIN_DATASETS
selected_cpgs = np.array(cpg_names)[final_mask]
print(f"Selected {len(selected_cpgs):,} CpGs covered in >= {MIN_SAMPLES_PER_GROUP} samples in >= {MIN_DATASETS} datasets.")

# Save filtered matrix in chunks
sel_rows = np.where(final_mask)[0]
n_sel = len(sel_rows)
n_cols = sample_idx

filtered_dset = h5f.create_dataset("filtered_matrix", shape=(n_sel, n_cols),
                                   dtype="float32", compression="gzip")

for start in range(0, n_sel, CHUNK_SIZE):
    end = min(start + CHUNK_SIZE, n_sel)
    rows = sel_rows[start:end]
    filtered_dset[start:end, :] = matrix_dset[rows, :n_cols]

# Replace original matrix with filtered one
del h5f["matrix"]
h5f.move("filtered_matrix", "matrix")

# Save metadata in HDF5
h5f.create_dataset("selected_cpg_names", data=selected_cpgs.astype("S"))
max_len = max(len(s) for s in sample_names)
h5f.create_dataset("samples", data=np.array(sample_names, dtype=f"S{max_len}"))
max_len = max(len(s) for s in sample_groups)
h5f.create_dataset("sample_groups", data=np.array(sample_groups, dtype=f"S{max_len}"))
h5f.close()
print(f"Saved all samples to: {output_path}")

print("\nValidation:")
print(f"Matrix shape: {n_sel:,} CpGs x {n_cols:,} samples")
print(f"Number of selected CpGs: {len(selected_cpgs):,}")
print(f"Number of samples: {n_cols:,}")

# Save metadata TSV
meta_df = pd.DataFrame({"sample": all_sample_names, "dataset": all_sample_groups})
meta_df.to_csv(metadata_file, sep="\t", index=False)
print(f"Saved metadata to: {metadata_file}")

# Save medians and lambdas
df = pd.DataFrame({
    "dataset": list(group_medians.keys()),
    "median_sd": list(group_medians.values()),
    "lambda": [group_lambdas[k] for k in group_medians.keys()]
})
df.to_csv(output_file_medsd_lambda, sep="\t", index=False)
print(f"Saved medians and lambdas to TSV: {output_file_medsd_lambda}")

print("\nAll done.")

