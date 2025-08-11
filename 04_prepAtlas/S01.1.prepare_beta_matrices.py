#!/usr/bin/env python3
"""
# Prepare Loyfer WGBS Atlas Data
---------------------------------
This script:
 1️⃣ Filters cell/tissue groups with at least 3 samples
 2️⃣ Loads beta files (.beta)
 3️⃣ Builds a CpG-by-sample matrix for each group
 4️⃣ Masks coverage < 10
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

output_folder = "/home/alice/Documents/Project_hvCpG/10X" ## in cluster "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X"
os.chdir(output_folder)
input_folder = "../betaFiles"

output_file = "all_scaled_matrix.h5"
metadata_file = "sample_metadata.tsv"
output_file_medsd_lambda = "all_medsd_lambda.tsv"

minCov = 10 # We will mask sites for which the coverage is below this

# 📂 1️⃣ Read metadata & filter valid groups

meta = pd.read_csv("../SupTab1_Loyfer2023.csv")

# Create composite group
meta["Composite Group"] = meta["Source Tissue"] + " - " + meta["Cell type"]

# Count samples per composite group
group_counts = meta.groupby("Composite Group").size()

# Keep groups with ≥ 3 samples
valid_groups = group_counts[group_counts >= 3].index.tolist()

print(f"✅ Found {len(valid_groups)} composite groups (Source Tissue + Cell type) with ≥ 3 samples.")

# Build dict: {group: [sample1, sample2, ...]}
samples_per_group = {
    g: meta.loc[meta["Composite Group"] == g, "Sample name"].tolist()
    for g in valid_groups
}

samples_per_group_short = {
    g: [s.split("-")[-1] for s in samples]
    for g, samples in samples_per_group.items()
}

# 🧬 2️⃣ Load CpG names

with open("../hg38CpGpos_Loyfer2023.txt") as f:
    cpg_names = [line.strip() for line in f]

print(f"✅ Loaded {len(cpg_names):,} CpG names.")
NR_SITES = len(cpg_names)

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
        beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
    return beta, cov

# 📊 4️⃣ For each group: build matrix, CLIP THEN logit transform (avoid massive issues for large 0s or 1s transformed).
# Save one file for all groups

# Collect all beta files once (outside the loop!)
all_files = glob.glob("../betaFiles/GSM*.hg38.beta")

## Store everything
# Count total number of valid samples first (flat list of all samples)

all_valid_samples = [s for group in samples_per_group_short.values() for s in group]
output_path = os.path.join(output_folder, "all_scaled_matrix.h5")
h5f = h5py.File(output_path, "w")
scaled_dset = h5f.create_dataset(
    "scaled_matrix",
    shape=(NR_SITES, len(all_valid_samples)),
    dtype="float32",
    compression="gzip"
)
sample_names = []
sample_groups = []
sample_idx = 0  # to track sample index across all groups

# To collect per-group stats
all_sample_names = []
all_sample_groups = []
group_medians = {}
group_lambdas = {}

for group, samples in samples_per_group_short.items():
    print(f"🔄 Processing {group} ({len(samples)} samples)")
    group_start_idx = sample_idx  # mark where this group starts
    for s in samples:
        matches = [fn for fn in all_files if re.search(f"-{s}\\.hg38\\.beta$", fn)]
        if not matches:
            print(f"⚠️  No beta file found for: {s}")
            continue
        fn = matches[0]
        beta, cov = load_beta(fn)
        if len(beta) != NR_SITES:
            raise ValueError(f"Mismatch: {s} has {len(beta)} CpGs, expected {NR_SITES}")
        ## Clip at 0.01 then logit
        epsilon = 0.01
        beta_clipped = np.clip(beta, epsilon, 1 - epsilon)
        scaled = np.log2(beta_clipped / (1 - beta_clipped)).astype(np.float32)
        ## Mask low coverage sites (we mask but keep the order of CpGs)
        scaled[cov < minCov] = np.nan
        scaled_dset[:, sample_idx] = scaled
        sample_names.append(s.encode())
        sample_groups.append(group.encode())
        all_sample_names.append(s)
        all_sample_groups.append(group)
        sample_idx += 1

    # Extract matrix for this group from the HDF5 dataset
    group_sample_indices = list(range(group_start_idx, sample_idx))
    if len(group_sample_indices) == 0:
        print(f"⚠️  Skipping {group}: no valid samples")
        continue
    mat = scaled_dset[:, group_sample_indices]
    # Mask CpGs with <3 non-NaNs
    valid_counts = np.sum(~np.isnan(mat), axis=1)
    mat[valid_counts < 3, :] = np.nan

    # ✅ Save masked matrix back to dataset
    scaled_dset[:, group_sample_indices] = mat

    # Compute stats
    row_sds = bn.nanstd(mat, axis=1)
    median_sd = np.nanmedian(row_sds)
    percentile_95 = np.nanpercentile(row_sds, 95)
    lambda_value = percentile_95 / median_sd
    group_medians[group] = median_sd
    group_lambdas[group] = lambda_value
    print(f"✅ {group}: median_sd = {median_sd:.4f}, lambda = {lambda_value:.4f}")

# Finalize HDF5
max_len = max(len(s) for s in sample_names)
dt = f'S{max_len}'
h5f.create_dataset("samples", data=np.array(sample_names, dtype=dt))

max_len = max(len(s) for s in sample_groups)
dt = f'S{max_len}'
h5f.create_dataset("sample_groups", data=np.array(sample_groups, dtype=dt))

h5f.create_dataset("cpg_names", data=np.array(cpg_names, dtype="S"))

h5f.close()
print(f"✅ Saved all samples to: {output_path}")

# Save metadata file
meta_df = pd.DataFrame({
    "sample": all_sample_names,
    "dataset": all_sample_groups
})
meta_df.to_csv(metadata_file, sep="\t", index=False)
print(f"✅ Saved metadata to: {metadata_file}")

# Save SDs and lambdas

# Combine the dictionaries into a DataFrame
df = pd.DataFrame({
    "dataset": list(group_medians.keys()),
    "median_sd": list(group_medians.values()),
    "lambda": [group_lambdas[k] for k in group_medians.keys()]
})

# Save to TSV
df.to_csv(output_file_medsd_lambda, sep="\t", index=False)

print(f"✅ Saved medians and lambdas to TSV: {output_file_medsd_lambda}")

print("\n🎉 All done.")

## 10X with correct clipping, 7th August

##✅ Found 46 composite groups (Source Tissue + Cell type) with ≥ 3 samples.
##✅ Loaded 29,401,795 CpG names.
##🔄 Processing Abdominal Subcut. - Adipocytes (3 samples)
##✅ Abdominal Subcut. - Adipocytes: median_sd = 0.6063, lambda = 2.6899
##🔄 Processing Bladder - Epithelium (5 samples)
##✅ Bladder - Epithelium: median_sd = 1.1999, lambda = 1.8013
##🔄 Processing Blood - B cells (3 samples)
##✅ Blood - B cells: median_sd = 0.7681, lambda = 2.3146
##🔄 Processing Blood - Granulocytes (3 samples)
##✅ Blood - Granulocytes: median_sd = 0.7269, lambda = 2.3427
##🔄 Processing Blood - Monocytes (3 samples)
##✅ Blood - Monocytes: median_sd = 0.8485, lambda = 2.1827
##🔄 Processing Blood - NK (3 samples)
##✅ Blood - NK: median_sd = 0.8558, lambda = 2.2026
##🔄 Processing Blood - T central memory CD4 (3 samples)
##✅ Blood - T central memory CD4: median_sd = 0.6959, lambda = 2.5560
##🔄 Processing Blood - T cytotoxic (CD8+) cells (3 samples)
##✅ Blood - T cytotoxic (CD8+) cells: median_sd = 0.7793, lambda = 2.2730
##🔄 Processing Blood - T effector cell CD8 (3 samples)
##✅ Blood - T effector cell CD8: median_sd = 0.7139, lambda = 2.5509
##🔄 Processing Blood - T effector memory CD4 (3 samples)
##✅ Blood - T effector memory CD4: median_sd = 0.6629, lambda = 2.6953
##🔄 Processing Blood - T helper(CD4+) cells (3 samples)
##✅ Blood - T helper(CD4+) cells: median_sd = 0.7060, lambda = 2.4016
##🔄 Processing Bone marrow - Erythrocyte progenitors (3 samples)
##✅ Bone marrow - Erythrocyte progenitors: median_sd = 0.4534, lambda = 2.8974
##🔄 Processing Brain - Neuronal (10 samples)
##✅ Brain - Neuronal: median_sd = 1.2014, lambda = 1.8241
##🔄 Processing Brain - Oligodendrocytes (4 samples)
##✅ Brain - Oligodendrocytes: median_sd = 0.8917, lambda = 1.9653
##🔄 Processing Breast - Basal epithelial (4 samples)
##✅ Breast - Basal epithelial: median_sd = 1.0155, lambda = 1.8505
##🔄 Processing Breast - Luminal epithelial (3 samples)
##✅ Breast - Luminal epithelial: median_sd = 0.8777, lambda = 2.1080
##🔄 Processing Colon - Endocrine (3 samples)
##✅ Colon - Endocrine: median_sd = 0.9731, lambda = 2.0978
##🔄 Processing Colon - Epithelium (5 samples)
##✅ Colon - Epithelium: median_sd = 0.9189, lambda = 1.9414
##🔄 Processing Endometrium - Epithelium (3 samples)
##✅ Endometrium - Epithelium: median_sd = 0.7785, lambda = 2.3460
##🔄 Processing Fallopien tubes - Epithelium (3 samples)
##✅ Fallopien tubes - Epithelium: median_sd = 0.7369, lambda = 2.2802
##🔄 Processing Gastric antrum - Epithelium (3 samples)
##✅ Gastric antrum - Epithelium: median_sd = 0.6390, lambda = 2.5203
##🔄 Processing Gastric body - Epithelium (3 samples)
##✅ Gastric body - Epithelium: median_sd = 0.6542, lambda = 2.4930
##🔄 Processing Gastric fundus - Epithelium (3 samples)
##✅ Gastric fundus - Epithelium: median_sd = 0.6133, lambda = 2.5540
##🔄 Processing Heart - Cardiomyocyte (4 samples)
##✅ Heart - Cardiomyocyte: median_sd = 0.8948, lambda = 2.0003
##🔄 Processing Heart - Fibroblast (4 samples)
##✅ Heart - Fibroblast: median_sd = 0.7425, lambda = 2.4008
##🔄 Processing Kidney glomerular - Endothelium (3 samples)
##✅ Kidney glomerular - Endothelium: median_sd = 0.7020, lambda = 2.3966
##🔄 Processing Kidney glomerular - Podocyte (3 samples)
##✅ Kidney glomerular - Podocyte: median_sd = 0.7658, lambda = 2.2720
##🔄 Processing Kidney tubular - Endothelium (3 samples)
##✅ Kidney tubular - Endothelium: median_sd = 0.7441, lambda = 2.4076
##🔄 Processing Kidney tubular - Epithelium (3 samples)
##✅ Kidney tubular - Epithelium: median_sd = 0.8099, lambda = 2.3028
##🔄 Processing Liver - Hepatocyte (6 samples)
##✅ Liver - Hepatocyte: median_sd = 0.9215, lambda = 1.9534
##🔄 Processing Lung alveolar - Endothelium (3 samples)
##✅ Lung alveolar - Endothelium: median_sd = 0.7575, lambda = 2.2733
##🔄 Processing Lung alveolar - Epithelium (3 samples)
##✅ Lung alveolar - Epithelium: median_sd = 0.6348, lambda = 2.4780
##🔄 Processing Lung bronchus - Epithelium (3 samples)
##✅ Lung bronchus - Epithelium: median_sd = 0.6711, lambda = 2.4284
##🔄 Processing Lung interstitial - Macrophages (3 samples)
##✅ Lung interstitial - Macrophages: median_sd = 0.7719, lambda = 2.2833
##🔄 Processing Pancreas - Acinar (4 samples)
##✅ Pancreas - Acinar: median_sd = 0.7151, lambda = 2.3649
##🔄 Processing Pancreas - Alpha (3 samples)
##✅ Pancreas - Alpha: median_sd = 0.6694, lambda = 2.5248
##🔄 Processing Pancreas - Beta (3 samples)
##✅ Pancreas - Beta: median_sd = 0.7447, lambda = 2.4184
##🔄 Processing Pancreas - Delta (3 samples)
##✅ Pancreas - Delta: median_sd = 0.7897, lambda = 2.3703
##🔄 Processing Pancreas - Duct (4 samples)
##✅ Pancreas - Duct: median_sd = 0.9136, lambda = 2.0156
##🔄 Processing Pancreas - Endothelium (4 samples)
##✅ Pancreas - Endothelium: median_sd = 0.8505, lambda = 2.1393
##🔄 Processing Prostate - Epithelium (4 samples)
##✅ Prostate - Epithelium: median_sd = 0.9489, lambda = 1.9747
##🔄 Processing Small intestine - Epithelium (3 samples)
##✅ Small intestine - Epithelium: median_sd = 0.7015, lambda = 2.4121
##🔄 Processing Thyroid - Epithelium (3 samples)
##✅ Thyroid - Epithelium: median_sd = 0.7273, lambda = 2.3525
##🔄 Processing Tongue - Epithelium (4 samples)
##✅ Tongue - Epithelium: median_sd = 0.9290, lambda = 1.9847
##🔄 Processing Tonsil palatine - Epithelium (3 samples)
##✅ Tonsil palatine - Epithelium: median_sd = 0.6162, lambda = 2.6465
##🔄 Processing Vascular saphenous - Endothelium (3 samples)
##✅ Vascular saphenous - Endothelium: median_sd = 0.7071, lambda = 2.3832
##✅ Saved all samples to: /home/alice/Documents/Project_hvCpG/10X/all_scaled_matrix.h5
##✅ Saved metadata to: sample_metadata.tsv
##✅ Saved medians and lambdas to TSV: all_medsd_lambda.tsv
##
##🎉 All done.
