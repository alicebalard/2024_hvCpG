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

output_folder = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer"
os.chdir(output_folder)
input_folder = "betaFiles"
output_file = "all_scaled_matrix.h5"
metadata_file = "sample_metadata.tsv"
output_file_medsd_lambda = "all_medsd_lambda.tsv"

epsilon = 1e-6       # For safe logit transform
minCov = 10 # We will mask sites for which the coverage is below this

# 📂 1️⃣ Read metadata & filter valid groups

meta = pd.read_csv("SupTab1_Loyfer2023.csv")

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

with open("hg38CpGpos_Loyfer2023.txt") as f:
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

# 📊 4️⃣ For each group: build matrix, logit transform. Save one file for all groups

# Collect all beta files once (outside the loop!)
all_files = glob.glob("betaFiles/GSM*.hg38.beta")

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
        p = np.clip(beta, epsilon, 1 - epsilon)
        scaled = np.log2(p / (1 - p)).astype(np.float32)
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

##🔄 Processing Abdominal Subcut. - Adipocytes (3 samples)
##✅ Abdominal Subcut. - Adipocytes: median_sd = 0.6171, lambda = 12.5324
##🔄 Processing Bladder - Epithelium (5 samples)
##✅ Bladder - Epithelium: median_sd = 2.1960, lambda = 3.7061
##🔄 Processing Blood - B cells (3 samples)
##✅ Blood - B cells: median_sd = 0.8454, lambda = 9.3145
##🔄 Processing Blood - Granulocytes (3 samples)
##✅ Blood - Granulocytes: median_sd = 0.7851, lambda = 9.9101
##🔄 Processing Blood - Monocytes (3 samples)
##✅ Blood - Monocytes: median_sd = 0.9880, lambda = 8.0698
##🔄 Processing Blood - NK (3 samples)
##✅ Blood - NK: median_sd = 0.9469, lambda = 8.4149
##🔄 Processing Blood - T central memory CD4 (3 samples)
##✅ Blood - T central memory CD4: median_sd = 0.7114, lambda = 11.1253
##🔄 Processing Blood - T cytotoxic (CD8+) cells (3 samples)
##✅ Blood - T cytotoxic (CD8+) cells: median_sd = 0.8222, lambda = 9.5628
##🔄 Processing Blood - T effector cell CD8 (3 samples)
##✅ Blood - T effector cell CD8: median_sd = 0.7314, lambda = 10.8723
##🔄 Processing Blood - T effector memory CD4 (3 samples)
##✅ Blood - T effector memory CD4: median_sd = 0.6701, lambda = 11.8299
##🔄 Processing Blood - T helper(CD4+) cells (3 samples)
##✅ Blood - T helper(CD4+) cells: median_sd = 0.7559, lambda = 10.3149
##🔄 Processing Bone marrow - Erythrocyte progenitors (3 samples)
##✅ Bone marrow - Erythrocyte progenitors: median_sd = 0.4560, lambda = 15.5837
##🔄 Processing Brain - Neuronal (10 samples)
##✅ Brain - Neuronal: median_sd = 5.2185, lambda = 1.5295
##🔄 Processing Brain - Oligodendrocytes (4 samples)
##✅ Brain - Oligodendrocytes: median_sd = 1.0246, lambda = 7.8146
##🔄 Processing Breast - Basal epithelial (4 samples)
##✅ Breast - Basal epithelial: median_sd = 6.4379, lambda = 1.2633
##🔄 Processing Breast - Luminal epithelial (3 samples)
##✅ Breast - Luminal epithelial: median_sd = 1.0460, lambda = 7.5503
##🔄 Processing Colon - Endocrine (3 samples)
##✅ Colon - Endocrine: median_sd = 1.2429, lambda = 6.5422
##🔄 Processing Colon - Epithelium (5 samples)
##✅ Colon - Epithelium: median_sd = 1.0272, lambda = 7.7236
##🔄 Processing Endometrium - Epithelium (3 samples)
##✅ Endometrium - Epithelium: median_sd = 0.8331, lambda = 9.4733
##🔄 Processing Fallopien tubes - Epithelium (3 samples)
##✅ Fallopien tubes - Epithelium: median_sd = 0.8093, lambda = 9.5993
##🔄 Processing Gastric antrum - Epithelium (3 samples)
##✅ Gastric antrum - Epithelium: median_sd = 0.6647, lambda = 11.5654
##🔄 Processing Gastric body - Epithelium (3 samples)
##✅ Gastric body - Epithelium: median_sd = 0.6811, lambda = 11.3344
##🔄 Processing Gastric fundus - Epithelium (3 samples)
##✅ Gastric fundus - Epithelium: median_sd = 0.6306, lambda = 12.0879
##🔄 Processing Heart - Cardiomyocyte (4 samples)
##✅ Heart - Cardiomyocyte: median_sd = 1.0055, lambda = 7.9813
##🔄 Processing Heart - Fibroblast (4 samples)
##✅ Heart - Fibroblast: median_sd = 0.7509, lambda = 10.3942
##🔄 Processing Kidney glomerular - Endothelium (3 samples)
##✅ Kidney glomerular - Endothelium: median_sd = 0.7310, lambda = 10.6118
##🔄 Processing Kidney glomerular - Podocyte (3 samples)
##✅ Kidney glomerular - Podocyte: median_sd = 0.8150, lambda = 9.6349
##🔄 Processing Kidney tubular - Endothelium (3 samples)
##✅ Kidney tubular - Endothelium: median_sd = 0.7835, lambda = 10.0943
##🔄 Processing Kidney tubular - Epithelium (3 samples)
##✅ Kidney tubular - Epithelium: median_sd = 0.8881, lambda = 8.8115
##🔄 Processing Liver - Hepatocyte (6 samples)
##✅ Liver - Hepatocyte: median_sd = 1.0220, lambda = 7.6088
##🔄 Processing Lung alveolar - Endothelium (3 samples)
##✅ Lung alveolar - Endothelium: median_sd = 0.7868, lambda = 9.8271
##🔄 Processing Lung alveolar - Epithelium (3 samples)
##✅ Lung alveolar - Epithelium: median_sd = 0.6640, lambda = 11.4899
##🔄 Processing Lung bronchus - Epithelium (3 samples)
##✅ Lung bronchus - Epithelium: median_sd = 0.7015, lambda = 10.9161
##🔄 Processing Lung interstitial - Macrophages (3 samples)
##✅ Lung interstitial - Macrophages: median_sd = 0.8251, lambda = 9.4686
##🔄 Processing Pancreas - Acinar (4 samples)
##✅ Pancreas - Acinar: median_sd = 0.7339, lambda = 10.5206
##🔄 Processing Pancreas - Alpha (3 samples)
##✅ Pancreas - Alpha: median_sd = 0.6935, lambda = 11.2484
##🔄 Processing Pancreas - Beta (3 samples)
##✅ Pancreas - Beta: median_sd = 0.7783, lambda = 10.1910
##🔄 Processing Pancreas - Delta (3 samples)
##✅ Pancreas - Delta: median_sd = 0.8431, lambda = 9.5016
##🔄 Processing Pancreas - Duct (4 samples)
##✅ Pancreas - Duct: median_sd = 1.0930, lambda = 7.3992
##🔄 Processing Pancreas - Endothelium (4 samples)
##✅ Pancreas - Endothelium: median_sd = 0.8912, lambda = 9.0293
##🔄 Processing Prostate - Epithelium (4 samples)
##✅ Prostate - Epithelium: median_sd = 1.0966, lambda = 7.2404
##🔄 Processing Small intestine - Epithelium (3 samples)
##✅ Small intestine - Epithelium: median_sd = 0.7487, lambda = 10.3884
##🔄 Processing Thyroid - Epithelium (3 samples)
##✅ Thyroid - Epithelium: median_sd = 0.7770, lambda = 10.0721
##🔄 Processing Tongue - Epithelium (4 samples)
##✅ Tongue - Epithelium: median_sd = 1.0171, lambda = 7.8011
##🔄 Processing Tonsil palatine - Epithelium (3 samples)
##✅ Tonsil palatine - Epithelium: median_sd = 0.6276, lambda = 12.2909
##🔄 Processing Vascular saphenous - Endothelium (3 samples)
##✅ Vascular saphenous - Endothelium: median_sd = 0.7455, lambda = 10.3830
