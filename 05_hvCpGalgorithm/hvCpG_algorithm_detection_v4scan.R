## hvCpG algorithm
## Alice Balard
## July 2025

###########
## Setup ##
###########

packages <- c("dplyr", "data.table", "matrixStats", "ggplot2", "reshape2","ggrepel",
              "parallel", "rhdf5")

## Install any packages that are not already installed:
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if(length(to_install)) install.packages(to_install)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))
rm(packages, to_install)

###############
## Data load ##
###############

## We need:

## (1) metadata:
## sample | dataset

## (2) median sd sigma k (1 sigma per dataset) & lambdas (95th percentile sd/ median sd) (1 lambda per dataset)
## dataset | median_sd | lambda

## (3) a function able to retrieve the methylation values for a given CpG
## sample | M value (scaled or not)

## (4) the full cpg names (should be retrieved in same h5 files as for (3)
library(rhdf5)

prepData <- function(analysis) {
  # Define base directories
  base_dirs <- list(
    testLocalPC = "~/Documents/10X/",
    Atlas10X = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/",
    Atlas5X = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/5X/"
  )
  # Special case for Maria's reduced arrays
  is_marias_reduced <- grepl("MariasarraysREDUCED", analysis)
  
  if (analysis %in% names(base_dirs)) {
    mydir <- base_dirs[[analysis]]
    metadata <- read.table(file.path(mydir, "sample_metadata.tsv"), sep = "\t", header = TRUE)
    medsd_lambdas <- read.table(file.path(mydir, "all_medsd_lambda.tsv"), sep = "\t", header = TRUE)
    cpg_names_all <- h5read(file.path(mydir, "all_scaled_matrix.h5"), "cpg_names")
    
    source_M_1CpG <- function(cpgpos) {
      M <- h5read(file.path(mydir, "all_scaled_matrix.h5"), "scaled_matrix", index = list(NULL, cpgpos))
      rownames(M) <- metadata$sample
      colnames(M) <- cpg_names_all[cpgpos]
      return(M)
    }
    
  } else if (is_marias_reduced) {
    # Extract suffix from analysis string
    x <- sub("MariasarraysREDUCED", "", analysis)
    basepath <- file.path("/home/alice/arraysh5_reducedMimicAtlas", x)
    metapath <- file.path(basepath, "all_metadata.tsv")
    if (!file.exists(metapath)) {
      stop("❌ Run 03_prepDatasetsMaria/S02.prepare_REDUCED_arrays_mimicAtlas.py")
    }
    metadata <- read.table(metapath, sep = "\t", header = TRUE)
    medsd_lambdas <- read.table(file.path(basepath, "all_medsd_lambda.tsv"), sep = "\t", header = TRUE)
    source_M_1CpG <- function(cpgpos) {
      M <- h5read(file.path(basepath, "all_scaled_matrix.h5"), "scaled_matrix", index = list(NULL, cpgpos))
      rownames(M) <- metadata$sample
      return(M)
    }
    cpg_names_all <- h5read(file.path(basepath, "all_scaled_matrix.h5"), "cpg_names")
  } else {
    stop("❌ Unknown analysis type.")
  }
  
  return(list(
    metadata = metadata,
    medsd_lambdas = medsd_lambdas,
    cpg_names_all = cpg_names_all,
    source_M_1CpG = source_M_1CpG
  ))
}

##############
## Run algo ##

# Maximum likelihood analysis
## The equation is: $$ \log\left(P(M_j)\right) = \sum_{i=1}^{n} \log\left( \sum_{Z_j=0}^{1} \left( \sum_{Z_{j,k}=0}^{1} P(M_{i,j} \mid Z_{j,k}) \times P(Z_{j,k} \mid Z_j) \right) \times P(Z_j) \right) $$             

###########################################
## Likelihood function for a given CpG j ##

getLogLik_oneCpG_optimized <- function(Mdf, metadata, medsd_lambdas, p0, p1, alpha){
  datasets <- unique(metadata$dataset)
  
  # ## TEST if unlog?
  # unlog <- function(x) {
  #   odds <- 2^x
  #   beta <- odds / (1 + odds)
  #   beta
  # }
  # Mdf = unlog(Mdf)
  
  ## Remove empty datasets
  x = na.omit(as.data.frame(Mdf))
  x$sample = rownames(x)
  long_x = reshape2::melt(x, id.vars = "sample", variable.name = "CpG", value.name = "value")
  long_x = left_join(long_x, metadata, by = "sample")
  datasets = datasets[datasets %in% long_x$dataset]
  
  ## initialise loglik
  log_P_Mj <- 0
  
  ## Precompute probability matrix
  p0p1_mat <- matrix(c(p0, 1 - p1, 1 - p0, p1), nrow = 2, byrow = TRUE)
  proba_hvCpG_vec <- c(1 - alpha, alpha)
  
  for (k in datasets) {
    ## Extract M values for the good samples
    Mij_vals <- Mdf[metadata$sample[metadata$dataset %in% k], ]
    if (all(is.na(Mij_vals))) next
    Mij_vals <- na.omit(as.numeric(Mij_vals))
    if (length(Mij_vals) == 0) next
    
    ## Mean of methylation value in this dataset for this CpG:
    mu_jk <- mean(Mij_vals, na.rm = TRUE)
    
    ## Sd_k and lamda_k for this k
    sd_k <- medsd_lambdas[medsd_lambdas$dataset %in% k,"median_sd"]
    lambda_k <- medsd_lambdas[medsd_lambdas$dataset %in% k,"lambda"]
    
    ## 2 possible sd, if CpG or hvCpG
    #sd_values <- c(sd_k, lambda_k * sd_k)
    ## New: clamp sd to avoid division by 0
    sd_values <- pmax(c(sd_k, lambda_k * sd_k), 1e-4)
    
    
    ## Calculate normal probabilities
    norm_probs <- tryCatch({
      matrix(c(
        dnorm(Mij_vals, mu_jk, sd_values[1]),
        dnorm(Mij_vals, mu_jk, sd_values[2])
      ), ncol = 2)
    }, error = function(e) {
      warning(paste("DNORM error in dataset", k))
      matrix(1, nrow = length(Mij_vals), ncol = 2)
    })
    
    ## Zjk probs: 2 hidden states, 2 configurations each
    zjk_probs <- array(dim = c(length(Mij_vals), 2, 2))
    
    for (zjk in 0:1) {
      zjk_probs[,,zjk+1] <- cbind(
        norm_probs[, zjk+1] * p0p1_mat[zjk+1, 1],
        norm_probs[, zjk+1] * p0p1_mat[zjk+1, 2]
      )
    }
    
    ## Combine, sum, log-sum-exp
    col_sums <- apply(zjk_probs, c(1, 3), sum)
    
    dataset_loglik <- sum(log(rowSums(col_sums %*% proba_hvCpG_vec)))
    
    if (!is.finite(dataset_loglik)) {
      warning(paste("Non-finite loglik in dataset", k))
      dataset_loglik <- 0
    }
    
    log_P_Mj <- log_P_Mj + dataset_loglik
  }
  return(log_P_Mj)
}

## We optimise over the parameter alpha, and specify everything else
################################
## Helper functions for optim ##
# fun2optim <- function(par, ...) {
#   alpha <- if (length(par) == 1) par else par["alpha"]
#   getLogLik_oneCpG_optimized(alpha = alpha, ...)
# }

## to rm
fun2optim <- function(alpha, Mdf, metadata, medsd_lambdas, p0, p1) {
  loglik <- getLogLik_oneCpG_optimized(Mdf, metadata, medsd_lambdas, p0, p1, alpha)
  
  # Save alpha and loglik to global trace
  alpha_trace[[length(alpha_trace) + 1]] <<- list(alpha = alpha, loglik = loglik)
  
  return(loglik)
}

####################################
## Run optimizer for a single CpG ##
runOptim1CpG <- function(Mdf, metadata, medsd_lambdas, p0, p1) {
  # Define multiple starting points for alpha
  start_alphas <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  # Run optim from each starting point and store results
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
  # Pick the result with the best (maximum) log-likelihood
  best_idx <- which.max(sapply(results, function(x) x$value))
  best_alpha <- results[[best_idx]]$par
  return(best_alpha)
}

#####################################################
## Parallelise optimizer for all CpG with mclapply ##

getAllOptimAlpha_parallel <- function(cpgPos_vec, NCORES, p0, p1, prep) {
  
  metadata <- prep$metadata
  medsd_lambdas <- prep$medsd_lambdas
  cpg_names_all <- prep$cpg_names_all
  source_M_1CpG <- prep$source_M_1CpG
  
  safe_run <- function(cpgPos) {
    Mdf = source_M_1CpG(cpgPos)
    ## to rm, test
    print(Mdf)
    ## ✅ “Proceed only if there are at least 3 matrices that contain at least one non-NA value.”
    if (length(table(setNames(metadata$dataset, metadata$sample)[names(Mdf[!is.na(Mdf),])])) < 3) {
      message(sprintf("CpG at index %s not covered in enough (>= 3) datasets", cpgPos))
      return(NA_real_)
    } else {
      tryCatch({
        runOptim1CpG(Mdf=Mdf, metadata = metadata, medsd_lambdas = medsd_lambdas, p0 = p0, p1 = p1
        )
      },
      error = function(e) {
        message(sprintf("CpG %s failed: %s", cpgPos, e$message))
        NA_real_
      })
    }
  }
  
  ## Do the parallel run OUTSIDE safe_run!
  res <- parallel::mclapply(cpgPos_vec, safe_run, mc.cores = NCORES, mc.preschedule = FALSE)
  
  ## Return matrix
  my_matrix <- matrix(unlist(res), ncol = 1)
  
  ## Name with CpG names
  rownames(my_matrix) <- cpg_names_all[cpgPos_vec]
  colnames(my_matrix) <- "alpha"
  return(my_matrix)
}

###################################
## Top-level: run & save results ##

runAndSave <- function(analysis, cpgPos_vec, resultDir, NCORES, p0, p1, overwrite = FALSE) {
  prep = prepData(analysis)
  
  ## Generate a safe object name
  obj_name <- paste0("results_", analysis, "_", length(cpgPos_vec), "CpGs_", p0, "p0_", p1, "p1")
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)
  
  ## Ensure directory ends with /
  if (!grepl("/$", resultDir)) resultDir <- paste0(resultDir, "/")
  if (!dir.exists(resultDir)) stop("Result directory does not exist: ", resultDir)
  
  file_name <- paste0(resultDir, obj_name, ".RData")
  
  if (overwrite == FALSE){
    ## Run only if file does not exist (by default)
    if (file.exists(file_name)) {
      message("⚠️  File already exists: ", file_name)
      message("⏭️  Skipping run to avoid overwriting.")
      return(invisible(NULL))
    }
  }
  result <- getAllOptimAlpha_parallel(cpgPos_vec = cpgPos_vec, NCORES = NCORES, 
                                      p0 = p0, p1 = p1, prep = prep)
  
  ## Assign to global env so user has it
  assign(obj_name, result, envir = .GlobalEnv)
  
  message("Saving to file: ", file_name)
  save(list = obj_name, file = file_name, envir = .GlobalEnv)
  message("💾 Result saved successfully.")
}
