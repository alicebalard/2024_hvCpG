#!/usr/bin/env python3
"""
Combine:
- R logic: load Maria’s .RDS-like files, merge them, keep CpGs covered in ≥15 datasets
- WGBS steps 3, 5, 6: build CpG-by-sample matrix, logit transform, save: matrix, median SD, lambda, CpG names
"""

import os
import glob
import pandas as pd
import numpy as np
import h5py
import bottleneck as bn
import pyreadr

# ------------------------------
# 📂 1) Load Maria’s CpG list
# ------------------------------
marias_hvcpgs = pd.read_csv(
    "/home/your_user/2024_hvCpG/03_prepDatasetsMaria/Derakhshan2022_ST5_hvCpG.txt"
)
print(f"✅ Loaded Maria’s hvCpG list: {marias_hvcpgs.shape[0]} rows")

# ------------------------------
# 📂 2) Find all .RDS files (batch 1)
# ------------------------------
folder1 = os.path.expanduser(
    "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/"
)
rds_files1 = glob.glob(os.path.join(folder1, "*.RDS"))
print(f"✅ Found {len(rds_files1)} .RDS files")

# ------------------------------
# 📂 3) Load all .RDS into list of matrices
# ------------------------------
rds_list_mat1 = []
for f in rds_files1:
    print(f"Reading: {f}")
    read_result = pyreadr.read_r(f)
    # Take the first object in the .RDS file:
    df = next(iter(read_result.values()))
    mat = df.values
    rownames = df.index.tolist()
    rds_list_mat1.append((os.path.splitext(os.path.basename(f))[0], mat, rownames))

# ------------------------------
# 📂 4) Load TCGA .RDS batch
# ------------------------------
tcga_file = "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/TCGA/TCGA_BMIQ_age_sex_PC_adjusted_OUTLIERS_REMOVED_round2.RDS"
tcga_result = pyreadr.read_r(tcga_file)

rds_list_mat2 = []
for key, df in tcga_result.items():
    mat = df.values
    rownames = df.index.tolist()
    rds_list_mat2.append((key, mat, rownames))

# Merge both
rds_list_mat = rds_list_mat1 + rds_list_mat2

print(f"✅ Total datasets loaded: {len(rds_list_mat)}")

# ------------------------------
# 📌 5) Keep only CpGs covered in ≥15 datasets
# ------------------------------
all_cpgs = []
for _, _, rownames in rds_list_mat:
    all_cpgs.extend(rownames)

cpg_counts = pd.Series(all_cpgs).value_counts()
common_cpgs = set(cpg_counts[cpg_counts >= 15].index.tolist())

print(f"✅ {len(common_cpgs):,} CpGs covered in ≥15 datasets")

# ------------------------------
# 📌 6) Filter each matrix for common CpGs + rebuild as aligned matrix
# ------------------------------
Maria_filtered_list_mat = []
for name, mat, rownames in rds_list_mat:
    df = pd.DataFrame(mat, index=rownames)
    df_filtered = df[df.index.isin(common_cpgs)]
    Maria_filtered_list_mat.append((name, df_filtered))

print(f"✅ Filtered all matrices")

# ------------------------------
# 📌 7) Combine into one CpG-by-sample matrix
# ------------------------------
# The rows = common CpGs (in same order!)
# The columns = samples (one column per dataset)

common_cpgs_sorted = sorted(list(common_cpgs))

all_columns = []
for name, df in Maria_filtered_list_mat:
    # Align rows
    aligned = df.reindex(common_cpgs_sorted)
    col = aligned.values[:, 0] if aligned.shape[1] == 1 else aligned.mean(axis=1).values
    all_columns.append(col)

big_matrix = np.column_stack(all_columns)
print(f"✅ Combined matrix: shape = {big_matrix.shape}")

# ------------------------------
# 📌 8) Logit transform with clipping
# ------------------------------
epsilon = 1e-6
p = np.clip(big_matrix, epsilon, 1 - epsilon)
scaled_matrix = np.log2(p / (1 - p))

# ------------------------------
# 📌 9) Save: matrix, median SD, lambda, CpG names
# ------------------------------
output_folder = "/home/your_user/2024_hvCpG/Maria_prepared"
os.makedirs(output_folder, exist_ok=True)

output_name = "Maria_array_combined"

# Save scaled matrix
with h5py.File(os.path.join(output_folder, f"{output_name}_scaled_matrix.h5"), "w") as f:
    f.create_dataset("scaled_matrix", data=scaled_matrix, compression="gzip")
    f.create_dataset("cpg_names", data=np.array(common_cpgs_sorted, dtype='S'))

print(f"✅ Saved scaled matrix: {output_name}_scaled_matrix.h5")

# Median of per-row SDs
row_sds = bn.nanstd(scaled_matrix, axis=1)
median_sd = np.nanmedian(row_sds)
percentile_95 = np.nanpercentile(row_sds, 95)
lambda_value = percentile_95 / median_sd

# Save median SD
with h5py.File(os.path.join(output_folder, f"{output_name}_median_sd.h5"), "w") as f:
    f.create_dataset("median_sd", data=np.array(median_sd))

# Save lambda
with h5py.File(os.path.join(output_folder, f"{output_name}_lambda.h5"), "w") as f:
    f.create_dataset("lambda", data=np.array(lambda_value))

print(f"✅ Saved median_sd: {median_sd:.4f} | lambda: {lambda_value:.4f}")

print("\n🎉 All done.")
