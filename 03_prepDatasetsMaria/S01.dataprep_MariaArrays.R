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

library(rhdf5)

# Save each matrix as an HDF5 file with CpG names
output_folder <- "~/outputHDF5_temp"

dir.create(output_folder, showWarnings = FALSE)

for (name in names(Maria_filtered_list_mat)) {
  mat <- Maria_filtered_list_mat[[name]]
  file <- file.path(output_folder, paste0(name, "_matrix.h5"))
  
  h5createFile(file)
  h5write(mat, file, "beta_matrix")
  h5write(rownames(mat), file, "cpg_names")
  h5write(colnames(mat), file, "sample_names")
  
  cat("✅ Saved:", file, "\n")
}