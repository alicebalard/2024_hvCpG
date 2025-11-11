## hvCpG algorithm (batched HDF5 loading)
## Alice Balard — Updated Aug 2025

###########
## Setup ##
###########

quiet_library <- function(pkg) {
  # Check if installed
  installed <- requireNamespace(pkg, quietly = TRUE)
  
  # Install if missing
  if (!installed) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      suppressMessages(suppressWarnings(
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
      ))
    }
    
    cran_pkgs <- suppressMessages(available.packages(repos = "https://cloud.r-project.org"))
    if (pkg %in% rownames(cran_pkgs)) {
      suppressMessages(suppressWarnings(
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      ))
    } else {
      suppressMessages(suppressWarnings(
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      ))
    }
  }
  
  # Load silently
  suppressPackageStartupMessages(
    suppressMessages(
      suppressWarnings(
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      )
    )
  )
  
  # Get version after loading
  version <- as.character(utils::packageVersion(pkg))
  
  # Print only our message
  cat(sprintf("Load package %s v%s\n", pkg, version))
}

quiet_library_all <- function(pkgs) {
  invisible(lapply(pkgs, quiet_library))
}

quiet_library_all(c("dplyr", "data.table", "matrixStats", "reshape2", "ggrepel",
                    "parallel", "rhdf5", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "tidyr", 
                    "ggplot2", "GenomicRanges", "rtracklayer"))
## NB: not all libraries are necessary; to clean when packaging

###############
## Data load ##
###############

prepData <- function(analysis, dataDir, subsetMetadata) {   
  if (grepl("MariasarraysREDUCED", analysis)) {
    x <- sub("MariasarraysREDUCED", "", analysis)
    basepath <- file.path("/home/alice/arraysh5_reducedMimicAtlas", x)
    metapath <- file.path(basepath, "all_metadata.tsv")
    if (!file.exists(metapath)) {
      stop("❌ Run 03_prepDatasetsMaria/S02.prepare_REDUCED_arrays_mimicAtlas.py")
    }
    metadata <- read.table(metapath, sep = "\t", header = TRUE)
    medsd_lambdas <- read.table(file.path(basepath, "all_medsd_lambda.tsv"), sep = "\t", header = TRUE)
    cpg_names_all <- h5read(file.path(basepath, "all_scaled_matrix.h5"), "cpg_names")
    h5file <- file.path(basepath, "all_scaled_matrix.h5")        
  } else {
    metadata <- read.table(file.path(dataDir, "sample_metadata.tsv"), sep = "\t", header = TRUE)
    medsd_lambdas <- read.table(file.path(dataDir, "all_medsd_lambda.tsv"), sep = "\t", header = TRUE)
    cpg_names_all <- h5read(file.path(dataDir, "all_matrix_noscale.h5"), "cpg_names")
    h5file <- file.path(dataDir, "all_matrix_noscale.h5")
  }
  ## Possible subset of samples or datasets
  if (!isFALSE(subsetMetadata)) {
    message("The algorithm runs on a subset of metadata")
    metadata = subsetMetadata
  }
  return(list(
    metadata = metadata,
    medsd_lambdas = medsd_lambdas,
    cpg_names_all = cpg_names_all,
    h5file = h5file
  ))
}

###########################################
## Likelihood function for a given CpG j ##
###########################################

getLogLik_oneCpG_optimized_fast <- function(Mdf, metadata, dataset_groups, ds_params, p0, p1, alpha) {
  
  # samples <- metadata$sample # to rm? useless?
  datasets <- unique(metadata$dataset)
  
  # Precompute p0/p1 matrix and mixture probs
  p0p1_mat <- matrix(c(p0, 1 - p1, 1 - p0, p1), nrow = 2, byrow = TRUE)
  proba_hvCpG_vec <- c(1 - alpha, alpha)
  
  log_P_Mj <- 0
  
  for (k in datasets) {
    Mij_vals <- as.numeric(Mdf[dataset_groups[[k]], , drop = FALSE])
    if (length(Mij_vals) < 3 || all(is.na(Mij_vals))) next
    
    # Precompute mean and SDs
    mu_jk <- mean(Mij_vals, na.rm = TRUE)
    
    ## Cal precomputed sds for both cases
    params =  ds_params[k, ]
    
    # Vectorized density
    norm_probs <- matrix(0, nrow = length(Mij_vals), ncol = 2)
    norm_probs[,1] <- dnorm(Mij_vals, mu_jk, params$sd0)
    norm_probs[,2] <- dnorm(Mij_vals, mu_jk, params$sd1)
    
    # Compute zjk_probs safely
    zjk_probs <- array(0, dim = c(length(Mij_vals), 2, 2))
    for (zjk in 0:1) {
      zjk_probs[,, zjk+1] <- cbind(
        norm_probs[, zjk+1] * p0p1_mat[zjk+1,1],
        norm_probs[, zjk+1] * p0p1_mat[zjk+1,2]
      )
    }
    
    # Sum across latent states and mixture
    col_sums <- apply(zjk_probs, c(1,3), sum)
    dataset_loglik <- sum(log(rowSums(col_sums %*% proba_hvCpG_vec)))
    if (!is.finite(dataset_loglik)) dataset_loglik <- 0
    
    log_P_Mj <- log_P_Mj + dataset_loglik
  }
  
  return(log_P_Mj)
}

################################
## Optimisation per CpG ##
################################

runOptim1CpG_gridrefine <- function(Mdf, metadata, dataset_groups, ds_params, p0, p1) {
  # Step 1. Coarse grid search (0, 0.05, 0.1...)
  grid <- seq(0, 1, length.out = 21)
  logliks <- vapply(grid, function(a) {
    getLogLik_oneCpG_optimized_fast(Mdf, metadata, dataset_groups, ds_params, p0, p1, a)
  }, numeric(1))
  
  best_idx <- which.max(logliks)
  alpha_start <- grid[best_idx]
  
  # Step 2. Local refinement with Brent, only in neighborhood
  lower <- ifelse(best_idx == 1, 0, grid[best_idx - 1])
  upper <- ifelse(best_idx == length(grid), 1, grid[best_idx + 1])
  
  resOpt <- optim(
    par = alpha_start,
    fn = function(alpha) {
      getLogLik_oneCpG_optimized_fast(Mdf, metadata, dataset_groups, ds_params, p0, p1, alpha)
    },
    method = "Brent",
    lower = lower, upper = upper,
    control = list(fnscale = -1)
  )
  return(resOpt$par)
}

#########################################
## Batch loading + parallel processing ##
#########################################

getAllOptimAlpha_parallel_batch_fast <- function(cpg_names_vec, NCORES, p0, p1, prep, batch_size = 1000, Nds) {
  metadata       <- prep$metadata
  cpg_names_all  <- prep$cpg_names_all
  h5file         <- prep$h5file
  medsd_lambdas  <- prep$medsd_lambdas
  
  ## Precompute dataset-level parameters
  ds_params = medsd_lambdas %>%
    dplyr::select(dataset, median_sd, lambda) %>%
    dplyr::mutate(sd0 = pmax(median_sd, 1e-4),
                  sd1 = pmax(lambda * median_sd, 1e-4)) %>%
    as.data.frame()
  rownames(ds_params) = ds_params$dataset
  
  ## Build a list of row indices grouped by dataset
  dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)
  
  # Read sample names once
  samples <- h5read(h5file, "samples")
  
  # Map CpG names to indices
  cpg_indices <- match(cpg_names_vec, cpg_names_all)
  if (anyNA(cpg_indices)) {
    stop("Some CpG names not found in HDF5: ", paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
  }
  
  # Container for results
  all_results <- vector("list", length(cpg_indices))
  
  # Split into batches
  batches <- split(cpg_indices, ceiling(seq_along(cpg_indices) / batch_size))
  
  for (b in seq_along(batches)) {
    message(sprintf(
      "📦 Loading batch %d / %d (%d CpGs) at %s",
      b, length(batches), length(batches[[b]]), Sys.time()
    ))
    
    # Load block of matrix (samples × CpGs) with rhdf5::h5read direct slice
    col_batches <- batches[[b]]
    M_batch <- h5read(
      file   = h5file,
      name   = "matrix",
      index  = list(NULL, col_batches)  # NULL = all rows, subset columns
    )
    
    # Assign dimnames
    rownames(M_batch) <- samples
    colnames(M_batch) <- cpg_names_all[col_batches]
    
    # Reorder rows to match metadata
    M_batch <- M_batch[metadata$sample, , drop = FALSE]
    
    sample_to_dataset <- setNames(metadata$dataset, metadata$sample)
    
    # Split CpGs into chunks (not one per worker)
    idx_split <- split(seq_len(ncol(M_batch)), cut(seq_len(ncol(M_batch)), NCORES, labels = FALSE))
    
    # Run in parallel over chunks
    chunk_results <- mclapply(idx_split, function(idx) {
      sapply(idx, function(i) {
        Mdf <- M_batch[, i, drop = FALSE]
        
        ## Require at least Nds datasets with data
        datasets_present <- unique(sample_to_dataset[names(Mdf[!is.na(Mdf), ])])
        if (length(datasets_present) < Nds) return(NA_real_)
        
        res <- tryCatch(
          runOptim1CpG_gridrefine(Mdf = Mdf, metadata = metadata, dataset_groups = dataset_groups,
                                  ds_params = ds_params, p0 = p0, p1 = p1),
          error = function(e) NA_real_
        )
        return(res)
      })
    }, mc.cores = NCORES)
    
    batch_results <- unlist(chunk_results, use.names = FALSE)
    
    # Store results (original positions)
    all_results[match(batches[[b]], cpg_indices)] <- batch_results
  }
  
  # Build final matrix
  my_matrix <- matrix(unlist(all_results), ncol = 1)
  rownames(my_matrix) <- cpg_names_vec
  colnames(my_matrix) <- "alpha"
  
  return(my_matrix)
}

######################
## Top-level runner ##
######################
##   dataDir: if  testLocalPC = "~/Documents/Project_hvCpG/10X/",
##   Atlas10X = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/",
##   Maria = "/home/alice/arraysh5"

runAndSave_fast <- function(analysis, cpg_names_vec, resultDir, NCORES, p0, p1, overwrite = FALSE,
                            batch_size = 1000, dataDir, skipsave=FALSE, Nds=3, subsetMetadata = FALSE) {
  t <- Sys.time()
  prep <- prepData(analysis, dataDir, subsetMetadata)
  message("Preparing the data took ", round(Sys.time() - t), " seconds")
  
  obj_name <- paste0("results_", analysis, "_", length(cpg_names_vec), "CpGs_", p0, "p0_", p1, "p1")
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)
  
  if (!grepl("/$", resultDir)) resultDir <- paste0(resultDir, "/")
  if (!dir.exists(resultDir)) {
    dir.create(resultDir, recursive = TRUE)
    message("New result directory ", resultDir, " created")
  }
  
  file_name <- paste0(resultDir, obj_name, ".RData")
  if (!overwrite && file.exists(file_name)) {
    message("⚠️ File already exists: ", file_name)
    return(invisible(NULL))
  }
  
  # Run batch + parallel processing
  result <- getAllOptimAlpha_parallel_batch_fast(
    cpg_names_vec = cpg_names_vec, NCORES = NCORES,
    p0 = p0, p1 = p1, prep = prep, batch_size = batch_size, Nds = Nds
  )
  
  assign(obj_name, result, envir = .GlobalEnv)
  
  if (!skipsave) {
    message("Saving to file: ", file_name)
    save(list = obj_name, file = file_name, envir = .GlobalEnv)
    message("💾 Result saved successfully.")
  }
}

###################################################
###### Part 2: (log) Probability per tissue #######
###################################################

getLogPhv_oneCpG_byTissue <- function(Mdf, metadata, dataset_groups, ds_params, p1) {
    
    datasets <- unique(metadata$dataset)
    
    ## Run for each dataset, and save in an list
    vector_datasets_logphv <- vector()
    for (k in datasets) {
        ## Safe extraction of the correct individual values
        sample_names <- metadata$sample[dataset_groups[[k]]]
        Mij_vals <- as.numeric(Mdf[sample_names, , drop = FALSE])    

        if (length(Mij_vals) < 3 || all(is.na(Mij_vals))) next 
        
                                        # Precompute mean and SDs
        mu_jk <- mean(Mij_vals, na.rm = TRUE)
        
        ## Call precomputed sds for both cases
        params =  ds_params[k, ]
        
        ## NB: select only 3 individuals per dataset!Makes it comparable across datasets with different N
        Mij_vals_3 <- sample(Mij_vals, 3, replace = F)
        
        ## Vectorized density
        norm_probs <- matrix(0, nrow = length(Mij_vals_3), ncol = 2)
        norm_probs[,1] <- dnorm(Mij_vals_3, mu_jk, params$sd0)
        norm_probs[,2] <- dnorm(Mij_vals_3, mu_jk, params$sd1)
        
        ## Compute log(Phv|Dk) as the sum of logs on all individuals
        dataset_logphv <- sum(log((norm_probs[,1]*p1)/(norm_probs[,1]*p1 + norm_probs[,2]*(1-p1))))
        vector_datasets_logphv[[k]] <- dataset_logphv
    }
    return(vector_datasets_logphv) # named vector of a logPr per dataset for this CpG
}

## Run and cbind all CpG to have a matrix CpG x logPrhv
getALL_LogPhv_oneCpG_byTissue_batch <- function(cpg_names_vec, NCORES, p1, prep, batch_size = 1000, Nds) {
    metadata       <- prep$metadata
    cpg_names_all  <- prep$cpg_names_all
    h5file         <- prep$h5file
    medsd_lambdas  <- prep$medsd_lambdas
    
    ## Precompute dataset-level parameters
    ds_params = medsd_lambdas %>%
        dplyr::select(dataset, median_sd, lambda) %>%
        dplyr::mutate(sd0 = pmax(median_sd, 1e-4),
                      sd1 = pmax(lambda * median_sd, 1e-4)) %>%
        as.data.frame()
    rownames(ds_params) = ds_params$dataset
    
    ## Build a list of row indices grouped by dataset
    dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)
    
    ## Read sample names once
    samples <- h5read(h5file, "samples")
    
    ## Map CpG names to indices
    cpg_indices <- match(cpg_names_vec, cpg_names_all)
    if (anyNA(cpg_indices)) {
        stop("Some CpG names not found in HDF5: ", paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
    }
    
    ## Container for results
    all_results <- vector("list", length(cpg_indices))
    
    ## Split into batches
    batches <- split(cpg_indices, ceiling(seq_along(cpg_indices) / batch_size))
    
    for (b in seq_along(batches)) {
        message(sprintf(
            "📦 Loading batch %d / %d (%d CpGs) at %s",
            b, length(batches), length(batches[[b]]), Sys.time()
        ))
        
     ## Load block of matrix (samples × CpGs) with rhdf5::h5read direct slice
        col_batches <- batches[[b]]
        M_batch <- h5read(
            file   = h5file,
            name   = "matrix",
            index  = list(NULL, col_batches)  # NULL = all rows, subset columns
        )
        
        ## Assign dimnames
        rownames(M_batch) <- samples
        colnames(M_batch) <- cpg_names_all[col_batches]
        
        ## Reorder rows to match metadata
        M_batch <- M_batch[metadata$sample, , drop = FALSE]
        
        sample_to_dataset <- setNames(metadata$dataset, metadata$sample)
        
        ## Split CpGs into chunks (not one per worker)
        idx_split <- split(seq_len(ncol(M_batch)), cut(seq_len(ncol(M_batch)), NCORES, labels = FALSE))
        
        ## Run in parallel over chunks
        chunk_results <- mclapply(idx_split, function(idx) {
            sapply(idx, function(i) {
                Mdf <- M_batch[, i, drop = FALSE]
                
                ## Require at least Nds datasets with data
                datasets_present <- unique(sample_to_dataset[names(Mdf[!is.na(Mdf), ])])
                if (length(datasets_present) < Nds) return(NA_real_)
                
                res <- tryCatch(
                    getLogPhv_oneCpG_byTissue(Mdf = Mdf, metadata = metadata, dataset_groups = dataset_groups,
                                              ds_params = ds_params, p1 = p1),
                    error = function(e) NA_real_
                )
                return(res) # named vector for one CpG
            })
        }, mc.cores = NCORES)
        
        ## Convert per-chunk results into a clean list of named vectors
        batch_results <- do.call(c, chunk_results)
        
        ## Make sure each CpG has a named vector of same length (datasets)
        dataset_names <- rownames(ds_params)
        batch_matrix <- matrix(NA_real_, nrow = length(batch_results), ncol = length(dataset_names),
                               dimnames = list(names(batch_results), dataset_names))
        
        ## Fill in each CpG's log-probabilities
        for (i in seq_along(batch_results)) {
            v <- batch_results[[i]]
            if (is.null(v)) next
            batch_matrix[i, names(v)] <- unlist(v)
        }
        
        ## Append to global results (as list of matrices)
        all_results[[b]] <- batch_matrix
    }
    
    ## Combine batches by rows (CpGs)
    my_matrix <- do.call(rbind, all_results)
    
    return(my_matrix)
}

runAndSave_tissueAnalysis <- function(analysis, cpg_names_vec, resultDir, NCORES, p1, overwrite = FALSE,
                                      batch_size = 1000, dataDir, skipsave=FALSE, Nds=3, subsetMetadata = FALSE) {
    t <- Sys.time()
    prep <- prepData(analysis, dataDir, subsetMetadata)
    message("Preparing the data took ", round(Sys.time() - t), " seconds")
    
    obj_name <- paste0("results_bytissue_", analysis, "_", length(cpg_names_vec), "CpGs_", p1, "p1")
    obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)
    
    if (!grepl("/$", resultDir)) resultDir <- paste0(resultDir, "/")
    if (!dir.exists(resultDir)) {
        dir.create(resultDir, recursive = TRUE)
        message("New result directory ", resultDir, " created")
    }
    
    file_name <- paste0(resultDir, obj_name, ".RData")
    if (!overwrite && file.exists(file_name)) {
        message("⚠️ File already exists: ", file_name)
        return(invisible(NULL))
    }
    
    ## Run batch + parallel processing
    result <- getALL_LogPhv_oneCpG_byTissue_batch(
        cpg_names_vec = cpg_names_vec, NCORES = NCORES,
        p1 = p1, prep = prep, batch_size = batch_size, Nds = Nds
    )
    
    assign(obj_name, result, envir = .GlobalEnv)
    
    if (!skipsave) {
        message("Saving to file: ", file_name)
        save(list = obj_name, file = file_name, envir = .GlobalEnv)
        message("💾 Result saved successfully.")
    }
}
