import numpy as np
import h5py
import pandas as pd
import os

## My directory
os.chdir("/home/alice/Documents/10X")

h5_path = "all_scaled_matrix.h5"
X = 46  # Min number of datasets in which CpG is covered by ≥3 samples
min_samples = 3
chunk_size = 100_000  # Number of CpGs to process at once

with h5py.File(h5_path, "r") as h5f:
    matrix = h5f["scaled_matrix"]
    samples = h5f["samples"][:].astype(str)
    sample_groups = h5f["sample_groups"][:].astype(str)
    cpg_names = h5f["cpg_names"][:].astype(str)
    NR_SITES = matrix.shape[0]

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

    # Final CpG selection
    final_mask = dataset_pass_count >= X
    selected_cpgs = np.array(cpg_names)[final_mask]
    print(f"✅ Selected {len(selected_cpgs):,} CpGs covered in ≥{min_samples} samples in ≥{X} datasets.")

    # Save result
    np.savetxt(f"selected_cpgs_min{min_samples}_in{X}_datasets.txt", selected_cpgs, fmt="%s")
    print(f"📝 Saved CpG list to: selected_cpgs_min{min_samples}_in{X}_datasets.txt")
