#!/usr/bin/env python3
"""
# Prepare Loyfer WGBS Atlas Data
---------------------------------
This script:
 1️⃣ Filters cell/tissue groups with at least 3 samples
 2️⃣ Loads beta files (.beta)
 3️⃣ Builds a CpG-by-sample matrix for each group
 4️⃣ Masks coverage < 20
 5️⃣ Applies logit transform: log2(p / (1 - p)) with clipping
 6️⃣ Saves: (a) matrix, (b) median of per-CpG row SD per data set, (c) lambda per dataset, (d) CpG names

Author: Alice Balard
"""

# 🧩 Setup
import os
import glob
import numpy as np
import pandas as pd
import re
import h5py
import bottleneck as bn

os.chdir('/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer')
output_folder = "datasets_prepared"
os.makedirs(output_folder, exist_ok=True)  # Create folder if it doesn't exist

NR_SITES = 29401795  # Known length of CpGs list
epsilon = 1e-6       # For safe logit transform

# 📂 1️⃣ Read metadata & filter valid groups
meta = pd.read_csv("SupTab1_Loyfer2023.csv")

group_counts = meta.groupby("Source Tissue").size()
valid_groups = group_counts[group_counts >= 3].index.tolist()

print(f"✅ Found {len(valid_groups)} groups with ≥ 3 samples.")

# Build dict: {group: [sample1, sample2, ...]}
samples_per_group = {
    g: meta.loc[meta["Source Tissue"] == g, "Sample name"].tolist()
    for g in valid_groups
}

samples_per_group_short = {
    g: [s.split("-")[-1] for s in samples]
    for g, samples in samples_per_group.items()
}

# 🧬 2️⃣ Load CpG names
with open("hg38CpGpos_Loyfer2023.txt") as f:
    cpg_names = [line.strip() for line in f]

assert len(cpg_names) == NR_SITES, f"Expected {NR_SITES} CpGs, got {len(cpg_names)}"

print(f"✅ Loaded {len(cpg_names):,} CpG names.")

# 💡 3️⃣ Helper: Load beta + coverage
def load_beta(path):
    """
    A .beta file is a binary file: NR_SITES rows × 2 columns:
      - [0]: # methylated reads (uint8)
      - [1]: total coverage (uint8)
    """
    arr = np.fromfile(path, dtype=np.uint8).reshape((-1, 2))
    meth = arr[:, 0]
    cov = arr[:, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        beta = np.where(cov == 0, np.nan, meth / cov)
    return beta, cov

# 📊 4️⃣ For each group: build matrix, logit transform, save
# Collect all beta files once (outside the loop!)
all_files = glob.glob("betaFiles/GSM*.hg38.beta")

for group, samples in samples_per_group_short.items():
    print(f"🔄 Processing {group} ({len(samples)} samples)")
    betas = []
    for s in samples:
        # Match any file that ends with -ZXXXXXXX.hg38.beta
        matches = [fn for fn in all_files if re.search(f"-{s}\\.hg38\\.beta$", fn)]
        if not matches:
            print(f"⚠️  No beta file found for: {s}")
            continue
        fn = matches[0]
        beta, cov = load_beta(fn)
        if len(beta) != NR_SITES:
            raise ValueError(f"Mismatch: {s} has {len(beta)} CpGs, expected {NR_SITES}")
        # Mask low coverage
        beta[cov < 20] = np.nan
        # Logit transform with clipping
        p = np.clip(beta, epsilon, 1 - epsilon)
        scaled = np.log2(p / (1 - p))
        betas.append(scaled)
    if not betas:
        print(f"⚠️  Skipping {group}: no valid samples")
        continue
    # 🧬 build matrix: rows = CpGs, columns = samples
    mat = np.column_stack(betas)
    ## Count valid (non-NaN) values per CpG (row)
    valid_counts = np.sum(~np.isnan(mat), axis=1)
    ## Find rows with fewer than 3 valid values
    rows_to_mask = valid_counts < 3
    ## Mask entire rows with NaN
    mat[rows_to_mask, :] = np.nan
    ## ✅ Check shape
    print("Matrix shape:", mat.shape)
    ## Should be (number_of_CpGs, number_of_samples)
    print(f"Rows (CpGs): {mat.shape[0]:,}")
    print(f"Columns (Samples): {mat.shape[1]}")
    ## 🔍 Peek at the first few rows and columns
    print("First 5 rows, all columns:")
    print(mat[:5, :])    
    #mat = np.vstack(betas).T
    # Save matrix
    with h5py.File(os.path.join(output_folder, f"{group}_scaled_matrix.h5"), "w") as f:
        f.create_dataset("scaled_matrix", data=mat, compression="gzip")
        f.create_dataset("samples", data=np.array(samples, dtype='S'))  # Save as bytes
        f.create_dataset("cpg_names", data=np.array(cpg_names, dtype='S'))
    print(f"✅ Saved: {group}_scaled_matrix.h5")
    # Calculate row SDs
    row_sds = bn.nanstd(mat, axis=1)
    # Calculate median SD
    median_sd = np.nanmedian(row_sds)
    # Save median SD only
    with h5py.File(os.path.join(output_folder, f"{group}_median_sd.h5"), "w") as f:
        f.create_dataset("median_sd", data=np.array(median_sd))
        print(f"✅ Saved: {group}_median_sd.h5 [median_sd = {median_sd:.4f}]")
    # Calculate lambda: (95th percentile / median)
    percentile_95 = np.nanpercentile(row_sds, 95)
    lambda_value = percentile_95 / median_sd
    # Save lambda only
    with h5py.File(os.path.join(output_folder, f"{group}_lambda.h5"), "w") as f:
        f.create_dataset("lambda", data=np.array(lambda_value))
        print(f"✅ Saved: {group}_lambda.h5 [lambda = {lambda_value:.4f}]")

# ---- Save shared CpG names once ----
with h5py.File(os.path.join(output_folder, "CpG_names.h5"), "w") as f:
    f.create_dataset("cpg_names", data=np.array(cpg_names, dtype='S'))

## File names have spaces: correct that
for filename in os.listdir(output_folder):
    if " " in filename:
        # New name: replace spaces with underscores
        new_filename = filename.replace(" ", "_")
        # Full paths
        old_path = os.path.join(output_folder, filename)
        new_path = os.path.join(output_folder, new_filename)
        # Rename
        os.rename(old_path, new_path)
        print(f"Renamed: {filename} -> {new_filename}")
    
print("\n🎉 All done.")
