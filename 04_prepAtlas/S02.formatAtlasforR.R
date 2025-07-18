###################################################################################
## Script to prepare Atlas objects used in our algorithm, once, and for each CpG ##
###################################################################################

## Should output:
## (1) median sd sigma k (1 sigma per dataset)
## (2) lambdas as a list (95th percentile sd/ median sd) (1 lambda per dataset)
## (3) a function able to retrieve the methylation values for a given CpG

library(rhdf5)

## 📁 Folder with all data prepared in S01
folder <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets_prepared"

## 🔍 List all scaled_matrix files
matrix_files <- list.files(folder, pattern = "_scaled_matrix\\.h5$", full.names = TRUE)

## 📦 Storage lists
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

##str(median_sds)
##str(lambdas)

## Function of 🧬  CpG of interestt
source_scaled_mat_1CpG <- function(pos) {
    cpg_index <- pos
    scaled_rows <- list()

    for (matrix_file in matrix_files) {
        ## 🗝️ SAFE: no explicit open/close of hd5
        row <- rhdf5::h5read(matrix_file, "scaled_matrix", index = list(NULL, cpg_index))
        group <- sub("_scaled_matrix\\.h5$", "", basename(matrix_file))
        scaled_rows[[group]] <- t(row)
    }

    return(scaled_rows)
}
## e.g. source_scaled_mat_1CpG(1)
