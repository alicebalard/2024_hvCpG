############################################################################################
## Script to prepare Maria's arrays objects used in our algorithm, once, and for each CpG ##
############################################################################################                                                               
## This outputs:                                                                  
## (1) median sd sigma k (1 sigma per dataset)                                    
## (2) lambdas as a list (95th percentile sd/ median sd) (1 lambda per dataset)   
## (3) a function able to retrieve the methylation values for a given CpG
library(rhdf5)

## 📁 Folder with all data prepared in S01                                          
folder <- "~/arraysh5files/"
# folder <- "/mnt/ing-s1/alice_data/arraysh5files/" ## doesn't work!

## 🔍 List all scaled_matrix files                                                  
matrix_files <- list.files(folder, pattern = "_scaled_matrix\\.h5$", full.names = TRUE)

## 📦 Storage lists                                                         
median_sds <- list()
lambdas <- list()

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

## Function of 🧬  CpG of interest
source_scaled_mat_1CpG <- function(pos) {
    scaled_rows <- list()
   
    for (matrix_file in matrix_files) {
        # Read the vector of CpG names
        cpg_names <- h5read(matrix_file, "cpg_names")
        cpg_name <- cpg_names[pos]
        
        # Read one row (CpG) from the matrix (as a column vector)
        row <- h5read(matrix_file, "scaled_matrix", index = list(NULL,pos))
        # Transpose and format
        row_df <- as.data.frame(t(row))
        rownames(row_df) <- cpg_name
        
        # Extract group name
        group <- sub("_scaled_matrix\\.h5$", "", basename(matrix_file))
        scaled_rows[[group]] <- row_df
    }
    
    return(scaled_rows)
}

## e.g. 
source_scaled_mat_1CpG(1)
