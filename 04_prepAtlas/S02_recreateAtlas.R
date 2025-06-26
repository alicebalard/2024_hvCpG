## Recreate Atlas data at the format needed for my algorithm

dataset_dir <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets"

## List beta and coverage files, ensure they match by group
beta_files <- list.files(dataset_dir, pattern = "^dataset_.*\\.rds$", full.names = TRUE)
coverage_files <- list.files(dataset_dir, pattern = "^coverage_.*\\.rds$", full.names = TRUE)

## CpG names
cpgnames <- readRDS(paste0(dataset_dir, "/cpg_names.rds"))

## Keep sites with coverage 20+, in at least 3 samples
recreateAtlas <- function(coverageMin = 20, NIND = 3){

    ## Extract group names (e.g., "Adipocytes", "Blood_Cauc", etc.)
    get_group = function(x) sub("^(dataset|coverage)_(.*)\\.rds$", "\\2", basename(x))
    groups = intersect(get_group(beta_files), get_group(coverage_files))

    ## Make named vectors for easy matching
    beta_files_named = setNames(beta_files, get_group(beta_files))
    coverage_files_named = setNames(coverage_files, get_group(coverage_files))

    process_in_chunks <- function(beta_mat, coverage_mat, group, chunk_size = 50000) {
        rownames(beta_mat) = cpgnames
        ## Mask beta values where coverage < coverageMin
        beta_mat[coverage_mat < coverageMin] = NA

        n_rows = nrow(beta_mat)
        keep_rows = logical(n_rows)

        for (i in seq(1, n_rows, by = chunk_size)) {
            idx = i:min(i + chunk_size - 1, n_rows)
            keep_rows[idx] = rowSums(!is.na(beta_mat[idx, , drop = FALSE])) >= NIND
        }

        beta_mat_filtered = beta_mat[keep_rows, , drop = FALSE]
        return(beta_mat_filtered)
    }

    result_list = vector("list", length(groups))
    names(result_list) = groups

    for (i in seq_along(groups)) {
        grp = groups[i]
        beta_mat = readRDS(beta_files_named[[grp]])
        coverage_mat = readRDS(coverage_files_named[[grp]])

        result_list[[i]] = process_in_chunks(beta_mat, coverage_mat, grp)

                                        # Clean up memory
        rm(beta_mat, coverage_mat)
        gc()
    }
    
    return(result_list)
}

## system.time(Atlas_data <- recreateAtlas()) ## how long does it take to recreate the dataset in R?
## 40 min with 100G needed
##     user   system  elapsed 
## 2147.568  179.553 2342.715

## table(sapply(Atlas_data, ncol))
## 3  4  5  6  8 10 11 13 19 22 
## 11  9  4  1  2  1  2  1  1  1 

## sapply(Atlas_data, ncol)
##       Adipocytes        Bladder-Ep           Blood-B      Blood-Granul 
##                3                 5                 5                 3 
## Blood-Mono+Macro          Blood-NK           Blood-T   Breast-Basal-Ep 
##               11                 3                22                 4 
##Breast-Luminal-Ep          Colon-Ep       Endothelium        Eryth-prog 
##                3                 8                19                 3 
##     Fallopian-Ep        Gastric-Ep      Head-Neck-Ep      Heart-Cardio 
##                3                11                13                 4 
##      Heart-Fibro         Kidney-Ep         Liver-Hep     Lung-Ep-Alveo 
##                4                 8                 6                 4 
##     Lung-Ep-Bron            Neuron         Oligodend    Ovary+Endom-Ep 
##                3                10                 4                 4 
##  Pancreas-Acinar    Pancreas-Alpha     Pancreas-Beta    Pancreas-Delta 
##                4                 3                 3                 3 
##    Pancreas-Duct       Prostate-Ep      Small-Int-Ep       Smooth-Musc 
##                4                 4                 5                 5 
##       Thyroid-Ep 
##                3 


 

