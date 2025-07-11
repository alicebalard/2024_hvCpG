###################################################################################
## Script to prepare Atlas objects used in our algorithm, once, and for each CpG ##
###################################################################################

library(rhdf5)

## 📁 Folder with all data prepared in S01
folder <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets_prepared"

## 🔍 List all scaled_matrix files
matrix_files <- list.files(folder, pattern = "_scaled_matrix\\.h5$", full.names = TRUE)

## 📦 Storage lists
scaled_rows <- list()
median_sds <- list()
lambdas <- list()

message("NB: for each CpG, we need at least 3 samples covered in a group for this group to be considered, otherwise the remaining one or 2 samples are masked by Nas.")
## 🔄 Loop over each matrix file [to do ONCE]
for (matrix_file in matrix_files) {
    ## Extract group name from file name
    group <- sub("_scaled_matrix\\.h5$", "", basename(matrix_file))
    ## ---- Read median_sd ----
    median_sd_file <- file.path(folder, paste0(group, "_median_sd.h5"))
    median_sd <- h5read(median_sd_file, "median_sd")
    median_sds[[group]] <- median_sd
    ## ---- Read lambda ----
    lambda_file <- file.path(folder, paste0(group, "_lambda.h5"))
    lambda <- h5read(lambda_file, "lambda")
    lambdas[[group]] <- lambda   
    cat(sprintf("✅ %s: median_sd=%.4f, lambda=%.4f\n",
                group, median_sd, lambda))
}

str(median_sds)
str(lambdas)

## To do for each CpG 
## 🧬 YpG of interest


pos <- 1 ## to loop over

## source_scaled_mat_1CpG <- function(pos){}

cpg_index <- 0 + pos ## Python is in base 0

## Create empty list to hold rows
scaled_rows <- list()

## 🔄 Loop over each matrix file
for (matrix_file in matrix_files) {
  ## Extract group name from file name
  group <- sub("_scaled_matrix\\.h5$", "", basename(matrix_file))
  
  ## Read only this row (1 x N)
  row <- h5read(matrix_file, "scaled_matrix", index = list(NULL, cpg_index))
  
  ## Store in list with group name
  scaled_rows[[group]] <- t(row)
}

