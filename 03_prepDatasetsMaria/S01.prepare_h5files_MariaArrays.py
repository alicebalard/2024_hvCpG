#!/usr/bin/env python3
"""
Combine:
- R logic: load Maria’s .RDS-like files, merge them, keep CpGs covered in ≥15 datasets
- WGBS steps 3, 5, 6: build CpG-by-sample matrix, logit transform, save: matrix, median SD, lambda, CpG names
"""
import os
import glob
import pyreadr
import pandas as pd
import numpy as np
import h5py
import bottleneck as bn

output_folder = os.path.expanduser("/SAN/ghlab/epigen/Alice/hvCpG_project/data/arrays_human/27dsh5files")

# 1) Paths
folder = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/arrays_human/30datasetsMaria"
os.makedirs(output_folder, exist_ok=True)

# 2) Find .RDS files
rds_files = glob.glob(os.path.join(folder, "*.RDS"))
print(f"✅ Found {len(rds_files)} .RDS files")

# 3) Load all .RDS into list
rds_list_mat = []
for f in rds_files:
    print(f"Reading: {f}")
    read_result = pyreadr.read_r(f)
    df = next(iter(read_result.values()))
    mat = df.values
    rownames = df.index.tolist()
    rds_list_mat.append((os.path.splitext(os.path.basename(f))[0], mat, rownames))

# 4) Keep CpGs covered in ≥15 datasets
all_cpgs = []
for _, _, rownames in rds_list_mat:
    all_cpgs.extend(rownames)

cpg_counts = pd.Series(all_cpgs).value_counts()
common_cpgs = set(cpg_counts[cpg_counts >= 15].index.tolist())
sorted_common_cpgs = sorted(common_cpgs)
print(f"✅ {len(common_cpgs):,} CpGs covered in ≥15 datasets")

# Save the global list to file
with open(os.path.join(output_folder, "sorted_common_cpgs.txt"), "w") as f:
    for cpg in sorted_common_cpgs:
        f.write(f"{cpg}\n")

# 5) Process each dataset separately — keep all samples (no mean!)
for name, mat, rownames in rds_list_mat:
    df = pd.DataFrame(mat, index=rownames)

    # Align CpGs to the full global set — insert NaNs for missing CpGs
    df_aligned = df.reindex(sorted_common_cpgs)

    # Logit transform (clip to avoid infs)
    epsilon = 1e-6
    p = np.clip(df_aligned.values, epsilon, 1 - epsilon)
    scaled_matrix = np.log2(p / (1 - p))  # logit transform

    # Row-wise SDs (ignoring NaNs)
    row_sds = bn.nanstd(scaled_matrix, axis=1)
    median_sd = np.nanmedian(row_sds)
    perc95 = np.nanpercentile(row_sds, 95)
    lambda_value = perc95 / median_sd if median_sd != 0 else np.nan

    # Save scaled matrix (CpGs x Samples)
    with h5py.File(os.path.join(output_folder, f"{name}_scaled_matrix.h5"), "w") as f:
        f.create_dataset("scaled_matrix", data=scaled_matrix, compression="gzip")
        f.create_dataset("cpg_names", data=np.array(sorted_common_cpgs, dtype='S'))
        f.create_dataset("sample_names", data=np.array(df_aligned.columns.tolist(), dtype='S'))

    # Save median SD
    with h5py.File(os.path.join(output_folder, f"{name}_median_sd.h5"), "w") as f:
        f.create_dataset("median_sd", data=np.array(median_sd))

    # Save lambda
    with h5py.File(os.path.join(output_folder, f"{name}_lambda.h5"), "w") as f:
        f.create_dataset("lambda", data=np.array(lambda_value))

    print(f"  Saved for {name} | CpGs: {scaled_matrix.shape[0]} | Samples: {scaled_matrix.shape[1]} | median_sd: {median_sd:.4f} | lambda: {lambda_value:.4f}")

print("\n🎉 All done — 1 H5 set per dataset saved.")

### FOR LSHTM:
### ------------------------------
### 📂 2) Find all .RDS files (batch 1)
### ------------------------------
##folder1 = os.path.expanduser(
##    "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/"
##)
##rds_files1 = glob.glob(os.path.join(folder1, "*.RDS"))
##print(f"✅ Found {len(rds_files1)} .RDS files")

### ------------------------------
### 📂 3) Load all .RDS into list of matrices
### ------------------------------
##rds_list_mat1 = []
##for f in rds_files1:
##    print(f"Reading: {f}")
##    read_result = pyreadr.read_r(f)
##    # Take the first object in the .RDS file:
##    df = next(iter(read_result.values()))
##    mat = df.values
##    rownames = df.index.tolist()
##    rds_list_mat1.append((os.path.splitext(os.path.basename(f))[0], mat, rownames))

### ------------------------------
### 📂 4) Load TCGA .RDS batch
### ------------------------------
##tcga_file = "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/TCGA/TCGA_BMIQ_age_sex_PC_adjusted_OUTLIERS_REMOVED_round2.RDS"
##tcga_result = pyreadr.read_r(tcga_file)
##
##rds_list_mat2 = []
##for key, df in tcga_result.items():
##    mat = df.values
##    rownames = df.index.tolist()
##    rds_list_mat2.append((key, mat, rownames))
##
### Merge both
##rds_list_mat = rds_list_mat1 + rds_list_mat2
##
##print(f"✅ Total datasets loaded: {len(rds_list_mat)}")
