#!/usr/bin/env python3
"""
- Load .RDS files from two batches (GEO & TCGA)
- Merge them
- Create a reduced HDF5 matrix where each CpG has exactly `nSamples × nDatasets` non-NaN values
- Save associated metadata and stats.

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
import random

def process_reduced_h5(nSamples, nDatasets, output_root="/home/alice/arraysh5_reducedMimicAtlas"):
    # --- PARAMETERS ---
    folder1 = "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/"
    folder2 = "/home/alice/tempRDS/"

    # Create dynamic output folder
    output_folder = f"{output_root}_{nSamples}samples_{nDatasets}datasets"
    os.makedirs(output_folder, exist_ok=True)

    output_path = os.path.join(output_folder, "all_scaled_matrix.h5")
    metadata_file = os.path.join(output_folder, "all_metadata.tsv")
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
    for name, mat, _ in rds_list_mat:
        all_valid_samples.extend([f"{name}_{i}" for i in range(mat.shape[1])])
    
    h5f = h5py.File(output_path, "w")
    scaled_dset = h5f.create_dataset(
        "scaled_matrix",
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
    
    # --- STEP 4: Map dataset → sample matrix column indices ---
    dataset_to_indices = {}
    sample_idx = 0
    for group_name, mat, _ in rds_list_mat:
        n_samples = mat.shape[1]
        dataset_to_indices[group_name] = list(range(sample_idx, sample_idx + n_samples))
        sample_idx += n_samples
    
    # --- STEP 5: Build reduced matrix (optimized) ---
    reduced_matrix = np.full((NR_SITES, len(all_valid_samples)), np.nan, dtype=np.float32)
    
    # Precompute CpG -> datasets map and dataset name -> rowname lookup
    cpg_to_datasets = {cpg: [] for cpg in sorted_common_cpgs}
    group_to_mat = {}
    group_to_rownames = {}
    
    for group_name, mat, rownames in rds_list_mat:
        group_to_mat[group_name] = mat
        group_to_rownames[group_name] = {c: i for i, c in enumerate(rownames)}
        for cpg in rownames:
            if cpg in cpg_to_datasets:
                cpg_to_datasets[cpg].append(group_name)
    
    # Fill reduced matrix
    for i, cpg in enumerate(sorted_common_cpgs):
        datasets = cpg_to_datasets[cpg]
        if len(datasets) < nDatasets:
            continue
        selected_datasets = random.sample(datasets, nDatasets)
        for group_name in selected_datasets:
            mat = group_to_mat[group_name]
            row_idx = group_to_rownames[group_name].get(cpg, None)
            if row_idx is None:
                continue
            sample_indices = dataset_to_indices[group_name]
            if len(sample_indices) < nSamples:
                continue
            selected_samples = random.sample(sample_indices, nSamples)
            for j in selected_samples:
                local_idx = j - sample_indices[0]
                value = mat[row_idx, local_idx]
                if pd.notna(value):
                    reduced_matrix[i, j] = value
                    
    # --- STEP 6: Write to HDF5 ---
    scaled_dset[:, :] = reduced_matrix
    
    # --- STEP 7: Metadata ---
    sample_idx = 0
    for group_name, mat, _ in rds_list_mat:
        n_samples = mat.shape[1]
        for i in range(n_samples):
            sample_name = f"{group_name}_{i}"
            sample_names.append(sample_name)
            sample_groups.append(group_name)
            all_sample_names.append(sample_name)
            all_sample_groups.append(group_name)
            sample_idx += n_samples
    
    h5f.create_dataset("samples", data=np.array(sample_names, dtype=h5py.string_dtype('utf-8')))
    h5f.create_dataset("sample_groups", data=np.array(sample_groups, dtype=h5py.string_dtype('utf-8')))
    h5f.create_dataset("cpg_names", data=np.array(sorted_common_cpgs, dtype="S"))
    
    h5f.close()
    print(f"✅ Saved reduced matrix to: {output_path}")
    
    # --- STEP 8: Save metadata ---
    meta_df = pd.DataFrame({
        "sample": all_sample_names,
        "dataset": all_sample_groups
    })
    meta_df.to_csv(metadata_file, sep="\t", index=False)
    print(f"✅ Saved metadata to: {metadata_file}")
    
    # --- STEP 9: Save medians/lambdas ---
    cpg_to_index = {cpg: i for i, cpg in enumerate(sorted_common_cpgs)}  # Precompute once
    
    for group_name, mat, rownames in rds_list_mat:
        cpg_idx = [cpg_to_index[cpg] for cpg in rownames if cpg in cpg_to_index]
        if not cpg_idx:
            continue
        col_indices = dataset_to_indices[group_name]
        subset = reduced_matrix[np.ix_(cpg_idx, col_indices)]
        if np.all(np.isnan(subset)):
            median_sd = np.nan
            lambda_value = np.nan
        else:
            row_sds = bn.nanstd(subset, axis=1)
            median_sd = np.nanmedian(row_sds)
            percentile_95 = np.nanpercentile(row_sds, 95)
            lambda_value = percentile_95 / median_sd if median_sd != 0 else np.nan
            group_medians[group_name] = median_sd
            group_lambdas[group_name] = lambda_value
            print(f"📊 {group_name}: median_sd = {median_sd:.4f}, lambda = {lambda_value:.4f}")
    
    df = pd.DataFrame({
        "dataset": list(group_medians.keys()),
        "median_sd": list(group_medians.values()),
        "lambda": [group_lambdas[k] for k in group_medians]
    })
    df.to_csv(output_file_medsd_lambda, sep="\t", index=False)
    print(f"✅ Saved stats to: {output_file_medsd_lambda}")
    
    print("\n🎉 All done.")

# Run for different parameters:

## NB: for memory resons, do it by block
for s in range(3,5):
    for n in range(3, 31):  # last is exclusive
        print(f"🚀 Running for nDatasets = {n}")
        process_reduced_h5(nSamples=s, nDatasets=n)
