## NB: done in R before, rm after                                                                                                                                                       
# Load the full list from file2                                                                                                                                                         
file2 <- readRDS("/mnt/old_user_accounts/p3/maria/PhD/Data/datasets/TCGA/TCGA_BMIQ_age_sex_PC_adjusted_OUTLIERS_REMOVED_round2.RDS")

# Choose output folder                                                                                                                                                                  
output_dir <- "/home/alice/tempRDS"
dir.create(output_dir, showWarnings = FALSE)

# Save each element (e.g., bladder, breast) into its own RDS file                                                                                                                       
for (name in names(file2)) {
  path <- file.path(output_dir, paste0(name, ".RDS"))
  saveRDS(file2[[name]], path)
  message("Saved: ", path)
}
