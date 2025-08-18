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

quiet_library_all(c("dplyr", "data.table", "matrixStats", "ggplot2", "reshape2", "ggrepel",
              "parallel", "rhdf5"))

###############
## Data load ##
###############

prepData <- function(analysis, dataDir) {   
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

getLogLik_oneCpG_optimized <- function(Mdf, metadata, medsd_lambdas, p0, p1, alpha) {
  datasets <- unique(metadata$dataset)
  x <- na.omit(as.data.frame(Mdf))
  x$sample <- rownames(x)
  long_x <- reshape2::melt(x, id.vars = "sample", variable.name = "CpG", value.name = "value")
  long_x <- left_join(long_x, metadata, by = "sample")
  datasets <- datasets[datasets %in% long_x$dataset]
  
  log_P_Mj <- 0
  p0p1_mat <- matrix(c(p0, 1 - p1, 1 - p0, p1), nrow = 2, byrow = TRUE)
  proba_hvCpG_vec <- c(1 - alpha, alpha)
  
  for (k in datasets) {
    samples_in_k <- metadata$sample[metadata$dataset %in% k]
    samples_in_M <- intersect(samples_in_k, rownames(Mdf))
    Mij_vals <- as.numeric(Mdf[samples_in_M, , drop = FALSE])
    if (length(Mij_vals) < 3 || all(is.na(Mij_vals))) next
    
    mu_jk <- mean(Mij_vals, na.rm = TRUE)
    sd_k <- medsd_lambdas[medsd_lambdas$dataset %in% k, "median_sd"]
    lambda_k <- medsd_lambdas[medsd_lambdas$dataset %in% k, "lambda"]
    sd_values <- pmax(c(sd_k, lambda_k * sd_k), 1e-4)
    
    norm_probs <- tryCatch({
      matrix(c(
        dnorm(Mij_vals, mu_jk, sd_values[1]),
        dnorm(Mij_vals, mu_jk, sd_values[2])
      ), ncol = 2)
    }, error = function(e) matrix(1, nrow = length(Mij_vals), ncol = 2))
    
    zjk_probs <- array(dim = c(length(Mij_vals), 2, 2))
    for (zjk in 0:1) {
      zjk_probs[,,zjk+1] <- cbind(
        norm_probs[, zjk+1] * p0p1_mat[zjk+1, 1],
        norm_probs[, zjk+1] * p0p1_mat[zjk+1, 2]
      )
    }
    
    col_sums <- apply(zjk_probs, c(1, 3), sum)
    dataset_loglik <- sum(log(rowSums(col_sums %*% proba_hvCpG_vec)))
    if (!is.finite(dataset_loglik)) dataset_loglik <- 0
    
    log_P_Mj <- log_P_Mj + dataset_loglik
  }
  return(log_P_Mj)
}

################################
## Optimisation per CpG ##
################################

fun2optim <- function(par, ...) {
  alpha <- if (length(par) == 1) par else par["alpha"]
  getLogLik_oneCpG_optimized(alpha = alpha, ...)
}

runOptim1CpG <- function(Mdf, metadata, medsd_lambdas, p0, p1) {
  start_alphas <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  results <- lapply(start_alphas, function(start_alpha) {
    resOpt <- optim(
      par = start_alpha,
      fn = fun2optim,
      Mdf = Mdf, metadata = metadata,
      medsd_lambdas = medsd_lambdas,
      p0 = p0, p1 = p1,
      method = "Brent",
      lower = 0, upper = 1,
      control = list(fnscale = -1)
    )
    list(par = resOpt$par, value = resOpt$value)
  })
  best_idx <- which.max(sapply(results, function(x) x$value))
  return(results[[best_idx]]$par)
}

#########################################
## Batch loading + parallel processing ##
#########################################

getAllOptimAlpha_parallel_batch <- function(cpgPos_vec, NCORES, p0, p1, prep, batch_size = 1000) {
    metadata <- prep$metadata
    medsd_lambdas <- prep$medsd_lambdas
    cpg_names_all <- prep$cpg_names_all
    h5file <- prep$h5file
    
    all_results <- vector("list", length(cpgPos_vec))
    
    ## Split into batches
    batches <- split(cpgPos_vec, ceiling(seq_along(cpgPos_vec) / batch_size))
    
    for (b in seq_along(batches)) {
        message(sprintf(
            "📦 Loading batch %d / %d (%d CpGs) at %s",
            b,
            length(batches),
            length(batches[[b]]),
            Sys.time()
        ))
        ## Read multiple CpGs at once from HDF5
        M_batch <- h5read(h5file, "matrix", index = list(NULL, batches[[b]]))
        rownames(M_batch) <- metadata$sample
        colnames(M_batch) <- cpg_names_all[batches[[b]]]
        
        ## Process each CpG in parallel
        batch_results <- mclapply(seq_along(batches[[b]]), function(i) {
            Mdf <- M_batch[, i, drop = FALSE]
            if (length(table(setNames(metadata$dataset, metadata$sample)[names(Mdf[!is.na(Mdf),])])) < 3) {
                return(NA_real_)
            }
            tryCatch({
                runOptim1CpG(Mdf = Mdf, metadata = metadata, medsd_lambdas = medsd_lambdas, p0 = p0, p1 = p1)
            }, error = function(e) NA_real_)
        }, mc.cores = NCORES)
        
        ## Store results
        all_results[batches[[b]]] <- batch_results
    }
    
    ## Build final matrix
    my_matrix <- matrix(unlist(all_results), ncol = 1)
    rownames(my_matrix) <- cpg_names_all[cpgPos_vec]
    colnames(my_matrix) <- "alpha"
    return(my_matrix)
}

######################
## Top-level runner ##
######################

##   dataDir: if  testLocalPC = "~/Documents/Project_hvCpG/10X/",
##   Atlas10X = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/",
##   Maria = "/home/alice/arraysh5"

runAndSave <- function(analysis, cpgPos_vec, resultDir, NCORES, p0, p1, overwrite = FALSE, batch_size = 1000, dataDir) {
  prep <- prepData(analysis, dataDir)
  
  obj_name <- paste0("results_", analysis, "_", length(cpgPos_vec), "CpGs_", p0, "p0_", p1, "p1")
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)
  if (!grepl("/$", resultDir)) resultDir <- paste0(resultDir, "/")
  if (!dir.exists(resultDir)) stop("Result directory does not exist: ", resultDir)
  file_name <- paste0(resultDir, obj_name, ".RData")
  
  if (!overwrite && file.exists(file_name)) {
    message("⚠️ File already exists: ", file_name)
    return(invisible(NULL))
  }
  
  result <- getAllOptimAlpha_parallel_batch(
    cpgPos_vec = cpgPos_vec, NCORES = NCORES, 
    p0 = p0, p1 = p1, prep = prep, batch_size = batch_size
  )
  
  assign(obj_name, result, envir = .GlobalEnv)
  message("Saving to file: ", file_name)
  save(list = obj_name, file = file_name, envir = .GlobalEnv)
  message("💾 Result saved successfully.")
}
