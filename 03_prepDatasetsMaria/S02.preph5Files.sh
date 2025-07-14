import h5py
import numpy as np
import os

epsilon = 1e-6
output_folder = "./scaled_matrices"
os.makedirs(output_folder, exist_ok=True)

# Example: all your HDF5 matrix files
matrix_files = glob.glob("path/to/output/*.h5")

# Load them all
matrices = {}
for f in matrix_files:
    with h5py.File(f, "r") as h5:
        mat = h5["beta_matrix"][:]
        cpgs = h5["cpg_names"][:].astype(str)
        samples = h5["sample_names"][:].astype(str)
        matrices[os.path.basename(f)] = (mat, cpgs, samples)

# Now you can loop by group — suppose you define groups by filename pattern or metadata
groups = {
    "Group1": ["Dataset1_matrix.h5", "Dataset2_matrix.h5"],
    "Group2": ["Dataset3_matrix.h5"],
    # etc.
}

for group, files in groups.items():
    print(f"🔄 Processing {group}")
    betas = []

    for f in files:
        mat, cpgs, samples = filtered_matrices[f]
        # Logit transform
        p = np.clip(mat, epsilon, 1 - epsilon)
        logit = np.log2(p / (1 - p))
        betas.append(logit)

    # Combine
    combined = np.column_stack(betas)

    # Mask rows with <3 valid samples
    valid_counts = np.sum(~np.isnan(combined), axis=1)
    combined[valid_counts < 3, :] = np.nan

    # Save scaled matrix
    out_file = os.path.join(output_folder, f"{group}_scaled_matrix.h5")
    with h5py.File(out_file, "w") as f:
        f.create_dataset("scaled_matrix", data=combined, compression="gzip")
        f.create_dataset("cpg_names", data=np.array(common_cpgs, dtype="S"))
    print(f"✅ Saved: {out_file}")

    # SD stats
    row_sds = np.nanstd(combined, axis=1)
    median_sd = np.nanmedian(row_sds)
    percentile_95 = np.nanpercentile(row_sds, 95)
    lambda_value = percentile_95 / median_sd

    with h5py.File(out_file.replace("_scaled_matrix.h5", "_median_sd.h5"), "w") as f:
        f.create_dataset("median_sd", data=np.array(median_sd))

    with h5py.File(out_file.replace("_scaled_matrix.h5", "_lambda.h5"), "w") as f:
        f.create_dataset("lambda", data=np.array(lambda_value))

    print(f"✅ Saved SD: median = {median_sd:.4f}, lambda = {lambda_value:.4f}")
