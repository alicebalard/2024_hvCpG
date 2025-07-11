##############################
## Data prep Maria datasets ##
## Hosted on LSHTM server of Matt Silver

## On Maria
MariasCpGs <- read.csv("~/2024_hvCpG/03_prepDatasetsMaria/Derakhshan2022_ST5_hvCpG.txt")

folder1 <- "/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/"
folder1 <- normalizePath("/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/GEO/BMIQ + 10 PCs + age + sex OUTLIERS REMOVED/")
rds_files1 <- list.files(path = folder1, pattern = "\\.RDS$", full.names = TRUE)

# Read all .rds files into a list and convert each to a matrix
rds_list_mat1 <- lapply(rds_files1, function(file) {
  as.matrix(readRDS(file))
})
## Name the list elements by file names (without extensions)
names(rds_list_mat1) <- gsub("\\.RDS$", "", basename(rds_files1))

rds_list_mat2 <- readRDS("/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/TCGA/TCGA_BMIQ_age_sex_PC_adjusted_OUTLIERS_REMOVED_round2.RDS")
rds_list_mat2 <- lapply(rds_list_mat2, function(x) {
  as.matrix(x)
})

rds_list_mat <- c(rds_list_mat1, rds_list_mat2)

# Keep only the array background (covered in 15 datasets out of 30)
all_cpgs <- unlist(lapply(rds_list_mat, rownames))
cpg_counts <- table(all_cpgs)
common_cpgs <- names(cpg_counts[cpg_counts >= 15])

Maria_filtered_list_mat <- lapply(rds_list_mat, function(mat) {
  mat[rownames(mat) %in% common_cpgs, ]
})

rm(rds_list_mat1,rds_list_mat2, rds_list_mat, common_cpgs,cpg_counts,all_cpgs, folder1, rds_files1)


# 📊 For each group: load, logit transform, save
# Collect all beta files once (outside the loop!)
all_files = glob.glob("XX")

for group, samples in samples_per_group_short.items():
    print(f"🔄 Processing {group} ({len(samples)} samples)")
    betas = []
    for s in samples:
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
