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

output_folder = os.path.expanduser("/home/alice/arraysh5files") # to scp into ing-s1 soon after
os.makedirs(output_folder, exist_ok=True)

# 1) Find .RDS batch 1 GEO
folder1 = "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/"
rds_files1 = glob.glob(os.path.join(folder1, "*.RDS"))
print(f"✅ Found {len(rds_files1)} .RDS files")

# Load all .RDS into list
rds_list_mat1 = []
for f in rds_files1:
    print(f"Reading: {f}")
    read_result = pyreadr.read_r(f)
    df = next(iter(read_result.values()))
    mat = df.values
    rownames = df.index.tolist()
    rds_list_mat1.append((os.path.splitext(os.path.basename(f))[0], mat, rownames))

# 2) Find .RDS batch 2 TCGA files
## !!!! NB: done in R before, rm after : prepTCGAforpy.R
folder2 = "/home/alice/tempRDS/"
rds_files2 = glob.glob(os.path.join(folder2, "*.RDS"))
if len(rds_files2) == 0:
    print("⚠️  WARNING!!! No .RDS files found — in R, run `prepTCGAforpy`")
else:
    print(f"✅ Found {len(rds_files2)} .RDS files")

# Load all .RDS into list
rds_list_mat2 = []
for f in rds_files2:
    print(f"Reading: {f}")
    read_result = pyreadr.read_r(f)
    df = next(iter(read_result.values()))
    mat = df.values
    rownames = df.index.tolist()
    rds_list_mat2.append((os.path.splitext(os.path.basename(f))[0], mat, rownames))

# 3) Merge both
rds_list_mat = rds_list_mat1 + rds_list_mat2

print(f"✅ Total datasets loaded: {len(rds_list_mat)}")

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
<<<<<<< HEAD
=======
## ✅ 406,036 CpGs covered in ≥15 datasets
>>>>>>> 9e5ca0b9b22fe1eb6b6eb68d730b65a13c66fc72

# 5) Process each dataset separately — keep all samples (no mean!)
for name, mat, rownames in rds_list_mat:
    df = pd.DataFrame(mat, index=rownames)

    # Align CpGs to the full global set — insert NaNs for missing CpGs
    df_aligned = df.reindex(sorted_common_cpgs)

<<<<<<< HEAD
    # Logit transform (clip to avoid infs)
    epsilon = 1e-6
    p = np.clip(df_aligned.values, epsilon, 1 - epsilon)
    scaled_matrix = np.log2(p / (1 - p))  # logit transform
=======
    ## NB: DO NOT logit transform, data where already transformed
    scaled_matrix = df_aligned.values
>>>>>>> 9e5ca0b9b22fe1eb6b6eb68d730b65a13c66fc72

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
