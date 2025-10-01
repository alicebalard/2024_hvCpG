#!/usr/bin/env python3
"""
# Prepare Loyfer WGBS Atlas Data (No Scaling)
---------------------------------
This script:
 1️⃣ Filters cell/tissue groups with at least 3 samples. NEW: take only CD4+CD8+
 2️⃣ Loads beta files (.beta)
 3️⃣ Builds a CpG-by-sample matrix for each group. NEW: select only the 450k CpGs (to compare with array data)
 4️⃣ Masks coverage < 10
 5️⃣ Saves: (a) matrix, (b) median of per-CpG row SD per data set, (c) lambda per dataset, (d) CpG names
 6️⃣ NEW: Select CpGs covered in the 2 cell types
Author: Alice Balard
"""

import os
import glob
import numpy as np
import pandas as pd
import re
import h5py
import bottleneck as bn

###########################
## Part 1: prepare files ##
###########################

output_folder = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer//10X_CD4+CD8+"
os.chdir(output_folder)
input_folder = "../betaFiles"

output_file = "all_matrix_noscale.h5"
metadata_file = "sample_metadata.tsv"
output_file_medsd_lambda = "all_medsd_lambda.tsv"

minCov = 10  # Mask sites with coverage below this

# 📂 1️⃣ Read metadata & filter valid groups
meta = pd.read_csv("../SupTab1_Loyfer2023.csv")

# Keep first row
first_row = meta.iloc[[0], :]

# NEW: Keep only rows where Cell type contains CD4+ or CD8+
mask = meta["Cell type"].str.contains("CD4\\+|CD8\\+", regex=True)
cd4_cd8_rows = meta[mask]

# Combine first row and CD4/CD8 rows
meta = pd.concat([first_row, cd4_cd8_rows]).drop_duplicates()

## Same as before
meta["Composite Group"] = meta["Source Tissue"] + " - " + meta["Cell type"]
group_counts = meta.groupby("Composite Group").size()
valid_groups = group_counts[group_counts >= 3].index.tolist()

print(f"✅ Found {len(valid_groups)} composite groups (Source Tissue + Cell type) with ≥ 33samples.")

samples_per_group = {
    g: meta.loc[meta["Composite Group"] == g, "Sample name"].tolist()
    for g in valid_groups
}
samples_per_group_short = {
    g: [s.split("-")[-1] for s in samples]
    for g, samples in samples_per_group.items()
}

# 🧬 2️⃣ Load CpG names
cpg_bed = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz"

cpg_names = []
with gzip.open(cpg_bed, "rt") as f:  # "rt" = read text mode
    for line in f:
        chrom, pos, _ = line.strip().split("\t")[:3]
        pos = int(pos)
        cpg_names.append(f"{chrom}_{pos}-{pos+1}")

NR_SITES = len(cpg_names)
print(f"✅ Loaded {NR_SITES:,} CpG names.")

# 💡 3️⃣ Helper: Load beta + coverage
def load_beta(path):
    arr = np.fromfile(path, dtype=np.uint8).reshape((-1, 2))
    meth = arr[:, 0]
    cov = arr[:, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
    return beta, cov

# 📊 4️⃣ Build matrix without scaling
all_files = glob.glob("../betaFiles/GSM*.hg38.beta")
all_valid_samples = [s for group in samples_per_group_short.values() for s in group]

output_path = os.path.join(output_folder, output_file)
h5f = h5py.File(output_path, "w")
matrix_dset = h5f.create_dataset(
    "matrix",
    shape=(NR_SITES, len(all_valid_samples)),
    dtype="float32",
    compression="gzip"
)

sample_names = []
sample_groups = []
sample_idx = 0
all_sample_names = []
all_sample_groups = []
group_medians = {}
group_lambdas = {}

for group, samples in samples_per_group_short.items():
    print(f"🔄 Processing {group} ({len(samples)} samples)")
    group_start_idx = sample_idx
    for s in samples:
        matches = [fn for fn in all_files if re.search(f"-{s}\\.hg38\\.beta$", fn)]
        if not matches:
            print(f"⚠️  No beta file found for: {s}")
            continue
        beta, cov = load_beta(matches[0])
        if len(beta) != NR_SITES:
            raise ValueError(f"Mismatch: {s} has {len(beta)} CpGs, expected {NR_SITES}")
        # Mask low coverage
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
        print(f"⚠️  Skipping {group}: no valid samples")
        continue
    mat = matrix_dset[:, group_sample_indices]
    valid_counts = np.sum(~np.isnan(mat), axis=1)
    mat[valid_counts < 3, :] = np.nan
    matrix_dset[:, group_sample_indices] = mat

    row_sds = bn.nanstd(mat, axis=1)
    median_sd = np.nanmedian(row_sds)
    percentile_95 = np.nanpercentile(row_sds, 95)
    lambda_value = percentile_95 / median_sd
    group_medians[group] = median_sd
    group_lambdas[group] = lambda_value
    print(f"✅ {group}: median_sd = {median_sd:.4f}, lambda = {lambda_value:.4f}")

# Finalize HDF5
max_len = max(len(s) for s in sample_names)
h5f.create_dataset("samples", data=np.array(sample_names, dtype=f"S{max_len}"))
max_len = max(len(s) for s in sample_groups)
h5f.create_dataset("sample_groups", data=np.array(sample_groups, dtype=f"S{max_len}"))
h5f.create_dataset("cpg_names", data=np.array(cpg_names, dtype="S"))
h5f.close()
print(f"✅ Saved all samples to: {output_path}")

# Save metadata
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

##############################################################################
## Part 2: select CpGs covered in the 2 cell types AND only for array (NEW) ##
##############################################################################

h5_path = "all_matrix_noscale.h5"
X = 2  # Min number of datasets in which CpG is covered by ≥3 samples
min_samples = 3
chunk_size = 100_000  # Number of CpGs to process at once

with h5py.File(h5_path, "r") as h5f:
    matrix = h5f["matrix"]
    samples = h5f["samples"][:].astype(str)
    sample_groups = h5f["sample_groups"][:].astype(str)
    cpg_names = h5f["cpg_names"][:]
    NR_SITES = matrix.shape[0]

    # NEW: Load list of CpGs of interest
    filter_path = "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/hg38array.txt"
    with open(filter_path) as f:
        allowed_cpgs = set(line.strip() for line in f if line.strip())

    print(f"✅ Loaded {len(allowed_cpgs)} CpGs from hg38array.txt")

    # Build a boolean mask for CpGs to keep
    allowed_mask = np.isin(cpg_names, list(allowed_cpgs))

    # Get sample indices per group
    df = pd.DataFrame({"sample": samples, "group": sample_groups})
    group_to_indices = df.groupby("group").indices

    # Tracker: how many datasets each CpG is covered in (≥3 samples)
    dataset_pass_count = np.zeros(NR_SITES, dtype=np.uint16)
    for group, indices in group_to_indices.items():
        print(f"⏳ Processing group: {group} with {len(indices)} samples")
        for start in range(0, NR_SITES, chunk_size):
            end = min(start + chunk_size, NR_SITES)
            submatrix_chunk = matrix[start:end, indices]  # shape: (chunk_size, num_samples_in_group)

            # Count non-NaNs per row (CpG)
            coverage_counts = np.sum(~np.isnan(submatrix_chunk), axis=1)

            # Mask: CpGs with ≥ min_samples
            mask = coverage_counts >= min_samples

            # Increment global counter
            dataset_pass_count[start:end][mask] += 1

    # Final CpG selection; NEW: only select the cpg corresponding to array
    final_mask = (dataset_pass_count >= X) & allowed_mask
    selected_cpgs = np.array(cpg_names)[final_mask]

    print(f"✅ Selected {len(selected_cpgs):,} CpGs (≥{min_samples} samples in ≥{X} datasets, and present in allowed list).")


    # Convert from bytes -> str
    selected_cpgs_str = [c.decode() for c in selected_cpgs]

    # Save
    np.savetxt(
        f"selected_cpgs_min{min_samples}_in{X}_datasets.txt",
        selected_cpgs_str,
        fmt="%s"
    )
    
    print(f"📝 Saved CpG list to: selected_cpgs_min{min_samples}_in{X}_datasets.txt")

print("\n🎉 All done.")


