import h5py
import numpy as np
import pandas as pd
from collections import defaultdict
from tqdm import tqdm

########
## 5X ##
########

# === Paths ===
h5_path = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/5X/all_scaled_matrix.h5"

# === Load HDF5 contents ===
with h5py.File(h5_path, "r") as h5f:
    cpg_names = [x.decode() for x in h5f["cpg_names"][:]]
    sample_groups = [x.decode() for x in h5f["sample_groups"][:]]
    mat = h5f["scaled_matrix"][:]  # shape: [CpG, Sample]

# === Build mapping: group → sample indices ===
group_to_indices = defaultdict(list)
for idx, group in enumerate(sample_groups):
    group_to_indices[group].append(idx)

# === Init result arrays ===
covered_5 = np.zeros(mat.shape[0], dtype=np.uint16)

# === Iterate over groups ===
print("🔍 Counting datasets with ≥3 covered samples per CpG...")

for group, col_idxs in tqdm(group_to_indices.items(), desc="Groups"):
    submat = mat[:, col_idxs]  # shape = [CpG, samples in this dataset]
    is_cov5 = ~np.isnan(submat)
    count_5 = np.sum(is_cov5, axis=1)
    # Add 1 if ≥3 samples for this CpG in the group
    covered_5 += count_5 >= 3

#########
## 10X ##
#########

# === Paths ===
h5_path = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/all_scaled_matrix.h5"

# === Load HDF5 contents ===
with h5py.File(h5_path, "r") as h5f:
    cpg_names = [x.decode() for x in h5f["cpg_names"][:]]
    sample_groups = [x.decode() for x in h5f["sample_groups"][:]]
    mat = h5f["scaled_matrix"][:]  # shape: [CpG, Sample]

# === Build mapping: group → sample indices ===
group_to_indices = defaultdict(list)
for idx, group in enumerate(sample_groups):
    group_to_indices[group].append(idx)

# === Init result arrays ===
covered_10 = np.zeros(mat.shape[0], dtype=np.uint16)

# === Iterate over groups ===
print("🔍 Counting datasets with ≥3 covered samples per CpG...")

for group, col_idxs in tqdm(group_to_indices.items(), desc="Groups"):
    submat = mat[:, col_idxs]  # shape = [CpG, samples in this dataset]
    is_cov10 = ~np.isnan(submat)
    count_10 = np.sum(is_cov10, axis=1)
    # Add 1 if ≥3 samples for this CpG in the group
    covered_10 += count_10 >= 3

# === Make frequency tables ===
import numpy as np
import pandas as pd

# Create frequency table (number of CpGs per coverage count)
unique, counts = np.unique(covered_5, return_counts=True)

# Assemble into a dataframe
freq_table = pd.DataFrame({
    "datasets_covered_in": unique,
    "num_CpGs": counts
})

# Save to TSV
freq_table.to_csv("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/CpG_coverage_freqtable5X.tsv", sep="\t", index=False)

# Create frequency table (number of CpGs per coverage count)
unique, counts = np.unique(covered_10, return_counts=True)

# Assemble into a dataframe
freq_table = pd.DataFrame({
    "datasets_covered_in": unique,
    "num_CpGs": counts
})

# Save to TSV
freq_table.to_csv("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/CpG_coverage_freqtable10X.tsv", sep="\t", index=False)
