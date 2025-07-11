##############################
## Data prep Maria datasets ##
## Hosted on LSHTM server of Matt Silver

library(dplyr)

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

my_list_mat_Mariads <- lapply(rds_list_mat, function(mat) {
  mat[rownames(mat) %in% common_cpgs, ]
})

cpgnames <- unique(unlist(sapply(my_list_mat_Mariads, row.names)))
cpgnames <- cpgnames[order(cpgnames)]

## Load mQTL-matched controls
cistrans_GoDMC_hvCpG_matched_control <- 
  read.table("~/2024_hvCpG/03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

table(cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% MariasCpGs$CpG)
cpgnames[cpgnames %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
cpgnames[cpgnames %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length

sub_cistrans_GoDMC_hvCpG_matched_control <- cistrans_GoDMC_hvCpG_matched_control[
  cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% cpgnames &
    cistrans_GoDMC_hvCpG_matched_control$controlCpG_name %in% cpgnames,]

rm(rds_list_mat1,rds_list_mat2, rds_list_mat, common_cpgs, cpg_counts,all_cpgs, folder1, rds_files1, 
   cistrans_GoDMC_hvCpG_matched_control)

message("This script loads: \nMariasCpGs: list of hvCpG from Derakhshan2022 (ST5)\nmy_list_mat_Mariads: list of datasets\ncpgnames: all the cpgs covered in the datasets after filtration\nsub_cistrans_GoDMC_hvCpG_matched_control: cpgnames and matching mQTL controls (if exist)")
