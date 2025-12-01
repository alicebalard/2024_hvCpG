library(here)
source(here("05_hvCpGalgorithm/quiet_library.R"))

source(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R"))
hvCpGandControls <- prephvCpGandControls(codeDir = here())

# Load the names of CpGs
# Assuming they are stored in the HDF5 file as a dataset called "cpg_names"
cpg_names <- h5read("~/Documents/Project_hvCpG/10X/all_matrix_noscale.h5", "cpg_names")

# Find the positions/indices of the CpGs in hvCpG and matching control
col_indices <- match(c(hvCpGandControls$DerakhshanhvCpGs_names, hvCpGandControls$mQTLcontrols_names),
                     cpg_names)
# Remove NAs in case some CpGs are not found
col_indices <- col_indices[!is.na(col_indices)]

# Extract only the desired columns from the matrix
# Assuming the matrix is stored as "matrix" in the HDF5 file
subset_matrix <- h5read(file_path, "matrix", index = list(NULL, col_indices))

# Check result
dim(subset_matrix)
head(subset_matrix)

## TBC (trop long)

## same but names(rds_list_mat) 



# For each dataset
sd_diff <- lapply(names(rds_list_mat), function(dataset_name, selectN = NULL) {
  
  if (!is.null(selectN) && selectN < ncol(mat)) {
    mat <- mat[, 1:selectN, drop = FALSE]
  }
  
  mat <- rds_list_mat[[dataset_name]]
  
  hv_cpgs <- sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name
  ctrl_cpgs <- sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name
  
  # Only keep pairs that exist in the dataset
  keep <- hv_cpgs %in% rownames(mat) & ctrl_cpgs %in% rownames(mat)
  hv_cpgs <- hv_cpgs[keep]
  ctrl_cpgs <- ctrl_cpgs[keep]
  
  # Calculate SD
  hv_sds <- rowSds(as.matrix(mat[hv_cpgs, , drop = FALSE]), na.rm = TRUE)
  ctrl_sds <- rowSds(as.matrix(mat[ctrl_cpgs, , drop = FALSE]), na.rm = TRUE)
  
  # Difference
  hv_sds - ctrl_sds
})

names(sd_diff) <- names(rds_list_mat)

inner95_df <- lapply(names(sd_diff), function(dataset_name) {
  diffs <- sd_diff[[dataset_name]]
  range95 <- quantile(diffs, c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    dataset = dataset_name,
    lower_95 = range95[1],
    upper_95 = range95[2],
    median_diff = median(diffs, na.rm = TRUE)
  )
}) %>% bind_rows()

# Compute inner 95% for first 3 individuals
sd_diff_3 <- lapply(names(rds_list_mat), function(dataset_name) {
  mat <- rds_list_mat[[dataset_name]]
  mat <- mat[, 1:min(3, ncol(mat)), drop = FALSE]  # first 3 columns
  
  hv_cpgs <- sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name
  ctrl_cpgs <- sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name
  
  # Only keep pairs that exist in the dataset
  keep <- hv_cpgs %in% rownames(mat) & ctrl_cpgs %in% rownames(mat)
  hv_cpgs <- hv_cpgs[keep]
  ctrl_cpgs <- ctrl_cpgs[keep]
  
  hv_sds <- rowSds(as.matrix(mat[hv_cpgs, , drop = FALSE]), na.rm = TRUE)
  ctrl_sds <- rowSds(as.matrix(mat[ctrl_cpgs, , drop = FALSE]), na.rm = TRUE)
  
  hv_sds - ctrl_sds
})

names(sd_diff_3) <- names(rds_list_mat)

inner95_df_3 <- lapply(names(sd_diff_3), function(dataset_name) {
  diffs <- sd_diff_3[[dataset_name]]
  range95 <- quantile(diffs, c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    dataset = dataset_name,
    lower_95 = range95[1],
    upper_95 = range95[2],
    median_diff = median(diffs, na.rm = TRUE)
  )
}) %>% bind_rows()

# Add a type column to distinguish
inner95_df$type <- "all individuals"
inner95_df_3$type <- "first 3 individuals"

# Combine
plot_df <- rbind(inner95_df, inner95_df_3)

# Set dodge width
dodge_width <- 0.4

ggplot(plot_df, aes(x = dataset, y = median_diff, fill = type)) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), 
                width = 0.2, position = position_dodge(width = dodge_width), color = "black") +
  geom_point(shape = 21, size = 3, stroke = 0.5,
             position = position_dodge(width = dodge_width)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = c(0.85, 1), legend.title = element_blank(),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(y = "SD difference (hvCpG - control) ± inner 95% range",
       x = NULL,
       title = "CpG SD differences across datasets")