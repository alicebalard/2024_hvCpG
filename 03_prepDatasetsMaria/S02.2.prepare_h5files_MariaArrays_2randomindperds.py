#!/usr/bin/env python3
"""
- Load .RDS files from two batches (GEO & TCGA)
- Merge them
- Save everything into a large HDF5 file (all_matrix_noscale.h5) along with sample names, groups, CpG names, and summary statistics (median_sd, lambda)
- NB: no logit transform for these ones to stick to Maria's approach

Author: Alice Balard
"""
# 🧩 Setup

import os
import glob
import pyreadr
import pandas as pd
import numpy as np
import h5py
import bottleneck as bn

# --- PARAMETERS ---
folder1 = "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/"
folder2 = "/home/alice/tempRDS/"
output_folder = "/home/alice/arraysh5_3ind"
os.makedirs(output_folder, exist_ok=True)

output_path = os.path.join(output_folder, "all_matrix_noscale.h5")
metadata_file = os.path.join(output_folder, "sample_metadata.tsv")
output_file_medsd_lambda = os.path.join(output_folder, "all_medsd_lambda.tsv")

# --- STEP 1: Load all RDS files ---
def load_rds_matrix(file_list):
    result = []
    for f in file_list:
        print(f"📖 Reading: {f}")
        read_result = pyreadr.read_r(f)
        df = next(iter(read_result.values()))
        mat = df.values
        rownames = df.index.tolist()
        result.append((os.path.splitext(os.path.basename(f))[0], mat, rownames))
    return result

rds_files1 = glob.glob(os.path.join(folder1, "*.RDS"))
rds_files2 = glob.glob(os.path.join(folder2, "*.RDS"))

rds_list_mat1 = load_rds_matrix(rds_files1)
rds_list_mat2 = load_rds_matrix(rds_files2)

rds_list_mat = rds_list_mat1 + rds_list_mat2
print(f"✅ Total datasets loaded: {len(rds_list_mat)}")

# Keep only 2 random individuals (columns) per dataset ---
np.random.seed(42)  # fix seed for reproducibility
reduced_list_mat = []
for name, mat, rownames in rds_list_mat:
    if mat.shape[1] > 2:  # only reduce if more than 2 individuals
        cols = np.random.choice(mat.shape[1], size=2, replace=False)
        mat = mat[:, cols]
    reduced_list_mat.append((name, mat, rownames))

rds_list_mat = reduced_list_mat
print("✅ Reduced to 2 individuals per dataset")

# --- STEP 2: Identify common CpGs (covered in ≥15 datasets) ---
all_cpgs = []
for _, _, rownames in rds_list_mat:
    all_cpgs.extend(rownames)

cpg_counts = pd.Series(all_cpgs).value_counts()
common_cpgs = set(cpg_counts[cpg_counts >= 15].index.tolist())
sorted_common_cpgs = sorted(common_cpgs)
print(f"✅ {len(sorted_common_cpgs):,} CpGs covered in ≥15 datasets")

# --- STEP 3: Initialize HDF5 ---
NR_SITES = len(sorted_common_cpgs)
all_valid_samples = []
for name, mat, rownames in rds_list_mat:
    all_valid_samples.extend([f"{name}_{i}" for i in range(mat.shape[1])])

h5f = h5py.File(output_path, "w")
noscaled_dset = h5f.create_dataset(
    "matrix",
    shape=(NR_SITES, len(all_valid_samples)),
    dtype="float32",
    compression="gzip"
)

sample_names = []
sample_groups = []
all_sample_names = []
all_sample_groups = []
group_medians = {}
group_lambdas = {}

sample_idx = 0
cpg_index_map = {cpg: i for i, cpg in enumerate(sorted_common_cpgs)}

# --- STEP 4: Write filtered data to HDF5 ---
for group_name, mat, rownames in rds_list_mat:
    df = pd.DataFrame(mat, index=rownames)
    df_filtered = df.loc[df.index.intersection(sorted_common_cpgs)]
    if df_filtered.shape[0] < NR_SITES:
        missing = sorted_common_cpgs.copy()
        df_filtered = df_filtered.reindex(sorted_common_cpgs)
    mat_filtered = df_filtered.values.astype(np.float32)
    # Write to HDF5
    n_samples = mat_filtered.shape[1]
    noscaled_dset[:, sample_idx:sample_idx + n_samples] = mat_filtered
    # Metadata
    for i in range(n_samples):
        sample_name = f"{group_name}_{i}"
        sample_names.append(sample_name)
        sample_groups.append(group_name)
        all_sample_names.append(sample_name)
        all_sample_groups.append(group_name)
    # Stats
    row_sds = bn.nanstd(mat_filtered, axis=1)
    median_sd = np.nanmedian(row_sds)
    percentile_95 = np.nanpercentile(row_sds, 95)
    lambda_value = percentile_95 / median_sd
    group_medians[group_name] = median_sd
    group_lambdas[group_name] = lambda_value
    print(f"✅ {group_name}: median_sd = {median_sd:.4f}, lambda = {lambda_value:.4f}")
    sample_idx += n_samples

# --- STEP 5: Finalize HDF5 ---
max_len = max(len(s) for s in sample_names)
dt = h5py.string_dtype(encoding='utf-8')
h5f.create_dataset("samples", data=np.array(sample_names, dtype=dt))

max_len = max(len(s) for s in sample_groups)
dt = h5py.string_dtype(encoding='utf-8')
h5f.create_dataset("sample_groups", data=np.array(sample_groups, dtype=dt))

h5f.create_dataset("cpg_names", data=np.array(sorted_common_cpgs, dtype="S"))

h5f.close()
print(f"✅ Saved matrix to: {output_path}")

# --- STEP 6: Save metadata ---
meta_df = pd.DataFrame({
    "sample": all_sample_names,
    "dataset": all_sample_groups
})
meta_df.to_csv(metadata_file, sep="\t", index=False)
print(f"✅ Saved metadata to: {metadata_file}")

# --- STEP 7: Save medians/lambdas ---
df = pd.DataFrame({
    "dataset": list(group_medians.keys()),
    "median_sd": list(group_medians.values()),
    "lambda": [group_lambdas[k] for k in group_medians]
})
df.to_csv(output_file_medsd_lambda, sep="\t", index=False)
print(f"✅ Saved stats to: {output_file_medsd_lambda}")

print("\n🎉 All done.")
