## hvCpG algorithm
## Alice Balard
## June 2025

## Input data for our algorithm needs to be in the format list of dataframes:
# filtered_list_mat$Blood_Cauc %>% head
#            GSM1870951  GSM1870952  GSM1870953 GSM1870954 GSM1870955  GSM1870956 GSM1870957
# cg00000957 0.88650894 0.873786294 0.905683210 0.89756560 0.93243437 0.934731380 0.90938323
# cg00001349 0.88900837 0.849494403 0.899719641 0.83146357 0.85704137 0.872804810 0.91467193
# cg00001583 0.06226518 0.059181490 0.050006552 0.07168150 0.04283758 0.044865908 0.03723361

###########
## Setup ##

## Specific libPath for R --vanilla
.libPaths(c(
  "/opt/R/packages/lib_4.4.1",
  "/home/alice/R/x86_64-pc-linux-gnu-library/4.4",
  "/usr/lib/R/library"
))

packages <- c("dplyr", "data.table", "matrixStats", "ggplot2", "reshape2","ggrepel",
              "parallel")

# Install any packages that are not already installed
to_install <- setdiff(packages, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))
rm(packages, to_install)

##############
## Run algo ##

# Maximum likelihood analysis
## The equation is: $$ \log\left(P(M_j)\right) = \sum_{i=1}^{n} \log\left( \sum_{Z_j=0}^{1} \left( \sum_{Z_{j,k}=0}^{1} P(M_{i,j} \mid Z_{j,k}) \times P(Z_{j,k} \mid Z_j) \right) \times P(Z_j) \right) $$

## Functions to calculate different parameters of the datasets:
## These functions take as imput a list of matrices containing our datasets
# It outputs
# * scaled_list_mat: a list of scaled matrices
# * mu_jk: named list for each dataset of vectors containing the mean methylation for each CpG j
# * sigma_k: named list for each dataset of elements = the median sd for all CpGs

scale_my_list <- function(our_list_datasets){
  lapply(our_list_datasets, function(k){ ## Scale all the datasets in list
    k = log2(k/(1-k))
  })
}

calc_mu_jk <- function(scaled_list_mat){
  lapply(scaled_list_mat, function(k){rowMeans(k, na.rm = T)}) ## list of mean methylation for each CpG j in all datasets
}

calc_sigma_k <- function(scaled_list_mat){
  lapply(scaled_list_mat, function(k){
    sd_j <- rowSds(k, na.rm = T)  ## Calculate the row (=per CpG j) sd
    return(median(sd_j, na.rm = T))  ## Calculate the median sd in dataset k 
  })
}

calc_lambdas <- function(scaled_list_mat) {
    ## Calculate row SDs for each dataset matrix
    all_sd_jk <- sapply(scaled_list_mat, matrixStats::rowSds, na.rm = TRUE)
    ## Compute lambda values (95th percentile / median) for each dataset
    lambdas <- sapply(all_sd_jk, \(x) quantile(x, 0.95, na.rm = TRUE) / median(x, na.rm = TRUE))
    ## Clean names and return
    names(lambdas) <- gsub(".95%", "", names(lambdas))
    return(lambdas)
}

## Likelihood function for a given CpG j
## The arguments are:
##   * my_list_mat: a list of matrices, i.e. our datasets
## 
## For the function to optimise itself:
## * j: a CpG (default=first CpG)
## * lambdas: a vector of lambda, one per dataset (multiplicative factor of the sd_k for hvCpG)
## * p0: proba of true negative
## * p1: proba of true positive
## * alpha: probability of a CpG site to be a hvCpG

getLogLik_oneCpG_optimized <- function(scaled_list_mat, mu_jk_list, sigma_k_list,
                                       j, p0, p1, alpha, lambdas) {
    
    ## Ensure the number of lambda parameters matches scaled_list_mat's length
    if (length(lambdas) != length(scaled_list_mat)) {
        stop("Number of lambda arguments must match length(scaled_list_mat)")
    }
    
    ## Precompute all necessary values upfront
    datasets <- names(scaled_list_mat)
    log_P_Mj <- 0 ## initialise
    
    ## Precompute probability matrices
    p0p1_mat <- matrix(c(p0, 1-p1, 1-p0, p1), nrow=2, byrow=TRUE)
    
    proba_hvCpG_vec <- c(1-alpha, alpha)
    
    for (k in datasets) {
        ## 1. SAFE ROW INDEXING
        row_idx <- which(rownames(scaled_list_mat[[k]]) == j)
        if (length(row_idx) == 0) next  ## Skip missing CpGs
        
        ## 2. SAFE VALUE EXTRACTION
        Mij_vals <- tryCatch({
            scaled_list_mat[[k]][row_idx, ]
        }, error = function(e) numeric(0))
        
        ## Remove NAs and check length
        Mij_vals <- na.omit(Mij_vals)
        if (length(Mij_vals) == 0) next
        
        ## 3. MISSING VALUE HANDLING FOR MU
        mu_jk <- mu_jk_list[[k]][j]
        if (is.na(mu_jk)) {
            warning(paste("Missing mu_jk for", j, "in dataset", k))
            next
        }
        
        ## 4. SAFE STANDARD DEVIATION
        sd_k <- sigma_k_list[[k]]
        if (is.na(sd_k) || sd_k <= 0) {
            warning(paste("Invalid sd_k for dataset", k))
            next
        }
        
        ## 5. LAMBDA VALIDATION
        lambda_k <- lambdas[k]
        if (is.na(lambda_k) || lambda_k <= 0) {
            warning(paste("Invalid lambda_k for dataset", k))
            next
        }
        
        ## 6. VECTORIZED CALCULATIONS WITH CHECKS
        sd_values <- c(sd_k, lambda_k * sd_k)
        if (any(is.na(sd_values)) || any(sd_values <= 0)) {
            warning(paste("Invalid sd_values for", j, "in dataset", k))
            next
        }
        
        ## 7. PROBABILITY CALCULATIONS: proba(M(i,j) knowing Z(j,k) 0 or 1)
        ## a matrix of 2 columns (sd or lambda sd), one row per individual
        norm_probs <- tryCatch({
            matrix(c(
                dnorm(Mij_vals, mu_jk, sd_values[1]),
                dnorm(Mij_vals, mu_jk, sd_values[2])
            ), ncol=2)
        }, error = function(e) {
            warning(paste("DNORM error for", j, "in dataset", k))
            matrix(1, nrow=length(Mij_vals), ncol=2)  ## Fallback to neutral values
        })
        
        ## 8. ARRAY OPERATIONS: double the previous matrix to create two, for Zj=0 or 1
        zjk_probs <- array(dim = c(length(Mij_vals), 2, 2))
        for (zjk in 0:1) {
            zjk_probs[,,zjk+1] <- cbind(norm_probs[, zjk+1] * p0p1_mat[zjk+1, 1],
                                        norm_probs[, zjk+1] * p0p1_mat[zjk+1, 2])
        }
        ## zjk_probs[,,1] <- cbind(norm_probs[, 1] * p0p1_mat[1, 1], ### M with sd_k * p0
        ##                             norm_probs[, 1] * p0p1_mat[1, 2])## M with sd_k * 1-p1
        ## zjk_probs[,,2] <- cbind(norm_probs[, 2] * p0p1_mat[2, 1],## M with lambda_k * sd_k * 1-p0
        ##                             norm_probs[, 2] * p0p1_mat[2, 2])## M with lambda_k * sd_k * p1
        
        ## 9. LOGSUMEXP IMPLEMENTATION
        ## Sum columns across depth dimension
        col_sums <- apply(zjk_probs, c(1,3), sum)  ## c(1,3) applies the function sum to both rows (1) and the third dimension (3)
        ## apply(X, c(1,3), sum) = X[,1,] + X[,2,] + X[,3,] = sum of the columns within each row and third dim level
        ## for us, does row sum of
        ## [(M with sd_k * p0) + (M with sd_k * 1-p1)] +
        ## [(M with lambda_k * sd_k * 1-p0) + (M with lambda_k * sd_k * p1)]
        ## Which corresponds to:
        ## P(Mj|Zjk=0) * P(Zjk=0|Zj=0) + P(Mj|Zjk=0) * P(Zjk=0|Zj=1) +
        ## P(Mj|Zjk=1) * P(Zjk=1|Zj=0) + P(Mj|Zjk=1) * P(Zjk=1|Zj=1)
        ## Multiply columns proba_hvCpG_vec and sum, then log, then sum for all individuals
        
        ## How to handle log(0)? We replace -Inf by the next lowest value
        
        ## Debugging corner START
        ## print(sum(log(rowSums(col_sums %*% proba_hvCpG_vec))))
        ## Debugging corner END
        
        dataset_loglik <- sum(log(rowSums(col_sums %*% proba_hvCpG_vec)))
        
        if (!is.finite(dataset_loglik)) {
            warning(paste("Non-finite loglik for", j, "in dataset", k))
            dataset_loglik <- 0  ## Neutral value for problematic calculations
        }
        log_P_Mj <- log_P_Mj + dataset_loglik
    }
    return(log_P_Mj)
}

## Optim to find alpha
## https://www.magesblog.com/post/2013-03-12-how-to-use-optim-in-r/

## We optimise over the parameter alpha, and specify everything else

## Core helper: logit version of fun2optim
fun2optim_logit <- function(theta, scaled_list_mat, mu_jk_list, sigma_k_list,
                            j, p0, p1, lambdas) {
  alpha <- 1 / (1 + exp(-theta))  # logit transform
  names(alpha) <- "alpha"
  fun2optim(par = alpha,
            scaled_list_mat = scaled_list_mat,
            mu_jk_list = mu_jk_list,
            sigma_k_list = sigma_k_list,
            j = j, p0 = p0, p1 = p1, lambdas = lambdas)
}

## Core helper: direct version
fun2optim <- function(par, scaled_list_mat, mu_jk_list, sigma_k_list,
                      j, p0, p1, lambdas) {
  alpha <- par["alpha"]
  getLogLik_oneCpG_optimized(
    scaled_list_mat = scaled_list_mat,
    mu_jk_list = mu_jk_list,
    sigma_k_list = sigma_k_list,
    j = j, p0 = p0, p1 = p1, alpha = alpha, lambdas = lambdas
  )
}

## Run optimizer for a single CpG
runOptim1CpG <- function(CpG, scaled_list_mat, mu_jk_list, sigma_k_list,
                         lambdas, optimMeth, p0, p1) {
  if (optimMeth == "Nelder-Mead") {
##      par_init <- 0  # theta
      par_init <- -9.2  # theta for alpha = 0.0001
      resOpt <- optim(
      par = par_init,
      fn = fun2optim_logit,
      scaled_list_mat = scaled_list_mat,
      mu_jk_list = mu_jk_list,
      sigma_k_list = sigma_k_list,
      j = CpG, p0 = p0, p1 = p1, lambdas = lambdas,
      method = "Nelder-Mead",
      control = list(fnscale = -1) ## maximize instead of minimize
    )
    alpha_hat <- 1 / (1 + exp(-resOpt$par))  # back-transform
    return(alpha_hat)

  } else if (optimMeth == "L-BFGS-B") {
    par_init <- c(alpha = 0)
    par_lower <- c(alpha = 0)
    par_upper <- c(alpha = 1)
    resOpt <- optim(
      par = par_init,
      fn = fun2optim,
      scaled_list_mat = scaled_list_mat,
      mu_jk_list = mu_jk_list,
      sigma_k_list = sigma_k_list,
      j = CpG, p0 = p0, p1 = p1, lambdas = lambdas,
      method = "L-BFGS-B",
      lower = par_lower,
      upper = par_upper,
      control = list(fnscale = -1)
    )
    return(resOpt$par)
  } else {
      stop("Unknown optimization method: ", optimMeth)
  }
}

getAllOptimAlpha_parallel <- function(my_list_mat, cpgvec, optimMeth, NCORES, p0, p1) {
<<<<<<< HEAD
    ## Prepare data inside the function
=======
    message("🚀 Scaling datasets...")
>>>>>>> 9f54ab5809dc895c7c9c309976e28718b6851eea
    scaled_list_mat <- scale_my_list(my_list_mat)
    message("📊 Computing mu_jk, sigma_k, lambdas...")
    mu_jk_list <- calc_mu_jk(scaled_list_mat)
    sigma_k_list <- calc_sigma_k(scaled_list_mat)
    lambdas <- calc_lambdas(scaled_list_mat)
<<<<<<< HEAD
    message("Inputs prepared")

=======
    message("✅ Inputs ready. Running parallel optimization...")
>>>>>>> 9f54ab5809dc895c7c9c309976e28718b6851eea
    ## Use safe wrapper with tryCatch
    safe_run <- function(CpG) {
        tryCatch(
            runOptim1CpG(
                CpG,
                scaled_list_mat = scaled_list_mat,
                mu_jk_list = mu_jk_list,
                sigma_k_list = sigma_k_list,
                lambdas = lambdas,
                optimMeth = optimMeth,
                p0 = p0,
                p1 = p1
            ),
            error = function(e) {
            message(sprintf("CpG %s failed: %s", CpG, e$message))
            NA_real_
        }
        )
    }

    ## Run in parallel
    res <- mclapply(cpgvec, safe_run, mc.cores = NCORES, mc.preschedule = FALSE)

    ## Make sure length is same
    names(res) <- cpgvec

    my_matrix <- matrix(unlist(res), ncol = 1)
    rownames(my_matrix) <- names(res)
    colnames(my_matrix) <- "alpha_hat"
    return(my_matrix)    
}

## Top-level: run & save results
runAndSave <- function(my_list_mat, cpgvec, resultDir, optimMeth, NCORES, p0, p1) {
  # Generate a safe object name
  var_mylist <- deparse(substitute(my_list_mat))
  var_cpgvec <- deparse(substitute(cpgvec))
  obj_name <- paste0(
    "results_", optimMeth, "_",
    var_mylist, "_",
    var_cpgvec, "_",
    p0, "p0_", p1, "p1"
  )
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)

  # Ensure directory ends with /
  if (!grepl("/$", resultDir)) resultDir <- paste0(resultDir, "/")
  if (!dir.exists(resultDir)) stop("Result directory does not exist: ", resultDir)

  file_name <- paste0(resultDir, obj_name, ".RData")

  # Check if file exists first
  if (file.exists(file_name)) {
    message("File already exists: ", file_name)
    message("Skipping run to avoid overwriting.")
    return(invisible(NULL))
  }

  # Run only if file does not exist
  result <- getAllOptimAlpha_parallel(
    my_list_mat, cpgvec, optimMeth = optimMeth, NCORES = NCORES, p0 = p0, p1 = p1
  )

  # Assign to global env so user has it
  assign(obj_name, result, envir = .GlobalEnv)

  message("Saving to file: ", file_name)
  save(list = obj_name, file = file_name, envir = .GlobalEnv)
  message("💾 Result saved successfully.")
}
