## Mars 2025
## Alice Balard
## Identify hvCpGs in arrays

library(dplyr)
library(data.table)
library(matrixStats)

## Files untared from Maria's datasets folder
folder_path <- "/SAN/ghlab/pophistory/Alice/hvCpG_project/data/arrays_human/30datasetsMaria"

# List all .rds files in the folder
rds_files <- list.files(path = folder_path, pattern = "\\.RDS$", full.names = TRUE)

# Read all .rds files into a list and convert each to a matrix
## NB: needs big memory allocation to translate into matrix
rds_list_mat <- lapply(rds_files, function(file) {
    as.matrix(readRDS(file))
})

## Save in RDS for now
# saveRDS(rds_list_mat, "/SAN/ghlab/pophistory/Alice/hvCpG_project/data/arrays_human/rds_list_27.RDS")

rds_list_mat <- readRDS("/SAN/ghlab/pophistory/Alice/hvCpG_project/data/arrays_human/rds_list_27.RDS")

# Optionally, name the list elements by file names (without extensions)
names(rds_list_mat) <- gsub("\\.RDS$", "", basename(rds_files))

##     GSM1870951 GSM1870952 GSM1870953
## cg00000957 0.88650894 0.87378629 0.90568321
## cg00001349 0.88900837 0.84949440 0.89971964
## cg00001583 0.06226518 0.05918149 0.05000655

########################
## Part 1: preprocessing

## tbc see Maria's code

## array background 406 306 CpGs covered in at least 15 of the 30 datasets used in this study.

## Identify the array background
all_cpgs <- unlist(lapply(rds_list_mat, rownames))
cpg_counts <- table(all_cpgs)
common_cpgs <- names(cpg_counts[cpg_counts >= 15])
rm(all_cpgs)

# Keep only the array background                                                                                                  
for (i in 1:length(rds_list_mat)){
    rds_list_mat[[i]] <- rds_list_mat[[i]][rownames(rds_list_mat[[i]]) %in% common_cpgs,]
}

##############################################
## Part 2: identify hvCpGs based on thresholds

## Within each dataset, calculate the CpGs variance, and keep the top 5%
filtered_matrices <- lapply(rds_list_mat, function(mat) {
  # Step 1: Calculate row variances
  row_variances <- rowVars(mat, na.rm = T)
  
  # Step 2: Find the threshold for the top 5% variances
  threshold <- quantile(row_variances, probs = 0.95)
  
  # Step 3: Subset the matrix to keep only rows with variance above the threshold
  mat[row_variances > threshold, , drop = FALSE]
})

## Select hvCpG as top 5 CpG covered in >= 65% of the ds
all_cpgs_top <- unlist(lapply(filtered_matrices, rownames))
cpg_counts_top <- table(all_cpgs_top)

cpg_counts <- data.frame(cpg_counts) %>% rename("cpgs"="all_cpgs", "all_cpgs"="Freq")
cpg_counts_top <- data.frame(cpg_counts_top) %>% rename("cpgs"="all_cpgs_top", "all_cpgs_top"="Freq")

cpg_counts_full <- merge(cpg_counts, cpg_counts_top)
cpg_counts_full <- cpg_counts_full[cpg_counts_full$all_cpgs_top >= 0.65*cpg_counts_full$all_cpgs, , drop = FALSE]

## hvCpGs:
# rownames(cpg_counts_full)

nrow(cpg_counts_full)                                                                                                            
## 4074
## Maria has 4143 with 3 more datasets

##################################################
## Part 3: identify hvCpGs based on Bayesian stats

## We fix some parameters:
lambda = 5 # multiplicative factor of variance to be a hvCpG
p1 = 0.65 ## is a CpG hv in at least 65% of datasets?

## Within each dataset, calculate the hvCpG
hv_k <- lapply(rds_list_mat, function(mat){
    
    ## Scale
    mat = log2(mat/(1-mat))
    
    ## Calculate row variances
    row_variances <- rowVars(mat, na.rm = T)

    ## Calculate the median variance in this dataset
    medVar <- median(row_variances, na.rm = T)

    ## Is the variance higher than lambda medVar?
    names(row_variances[row_variances >= lambda * medVar])
})

## a hv in dataset k is a hv overall if hv in at least threshold
find_common_elements <- function(vector_list, threshold) {
  # Combine all vectors into a single vector
  all_elements <- unlist(vector_list)
  
  # Count occurrences of each unique element
  element_counts <- table(all_elements)
  
  # Calculate the threshold count
  threshold_count <- length(vector_list) * threshold
  
  # Find elements that meet or exceed the threshold
  common_elements <- names(element_counts[element_counts >= threshold_count])
  
  return(common_elements)
}

hv <- find_common_elements(vector_list = hv_k, threshold = p1)


