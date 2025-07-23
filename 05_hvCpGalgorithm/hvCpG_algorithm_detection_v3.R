## hvCpG algorithm
## Alice Balard
## July 2025

## Input data for our algorithm needs to be in the R format list of dataframes.
## To run, this script needs:
## scaled_list_mat --> function of a CpG that will retrieve methylation in all datasets for this CpG
#            GSM1870951  GSM1870952  GSM1870953 GSM1870954 GSM1870955  GSM1870956 GSM1870957
# cg00000957 0.88650894 0.873786294 0.905683210 0.89756560 0.93243437 0.934731380 0.90938323
## sigma k --> median sd per dataset for all CpGs
## lambda k --> lambda per dataset considering all CpGs

## For Atlas data on UCL server:
## source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/S02.formatAtlasforR.R")
## Outputs:
## median_sds (sigma_k)
## lambdas
## source_scaled_mat_1CpG(1 to length(cpgvec)

## For Maria's array on LSHTM server:
## source("~/2024_hvCpG/03_prepDatasetsMaria/S02.formatArraysforR.R")
## Outputs:
## median_sds (sigma_k)
## lambdas
## source_scaled_mat_1CpG(1 to length(cpgvec)

###########
## Check ##
# 🔍 List of required variables
required_vars <- c("median_sds", "lambdas")

# ❌ Check for missing variables
missing_vars <- required_vars[!sapply(required_vars, exists)]

# 🛑 If any are missing, stop and print informative error
if (length(missing_vars) > 0) {
  stop(sprintf("❌ Missing required variable(s): %s", paste(missing_vars, collapse = ", ")))
}

###########
## Setup ##
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

###########################################
## Likelihood function for a given CpG j ##
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
                                       p0, p1, alpha, lambdas) {
    ## lambdas must match length
    if (length(lambdas) != length(scaled_list_mat)) {
        stop("Number of lambda arguments must match length(scaled_list_mat)")
    }

    datasets <- names(scaled_list_mat)
    log_P_Mj <- 0

    ## Precompute probability matrix
    p0p1_mat <- matrix(c(p0, 1 - p1, 1 - p0, p1), nrow = 2, byrow = TRUE)
    proba_hvCpG_vec <- c(1 - alpha, alpha)

    for (k in datasets) {

        ## Get the only row: the whole matrix is 1-row already
        Mij_vals <- scaled_list_mat[[k]]
        if (is.null(Mij_vals) || length(Mij_vals) == 0) next
        Mij_vals <- as.numeric(Mij_vals)
        Mij_vals <- na.omit(Mij_vals)
        if (length(Mij_vals) == 0) next
        mu_jk <- mu_jk_list[[k]]
        if (is.na(mu_jk)) {
            warning(paste("Missing mu_jk for dataset", k))
            next
        }
        sd_k <- sigma_k_list[[k]]
        if (is.na(sd_k) || sd_k <= 0) {
            warning(paste("Invalid sd_k for dataset", k))
            next
        }
        lambda_k <- lambdas[[k]]
        if (is.na(lambda_k) || lambda_k <= 0) {
            warning(paste("Invalid lambda_k for dataset", k))
            next
        }
        sd_values <- c(sd_k, lambda_k * sd_k)
        if (any(is.na(sd_values)) || any(sd_values <= 0)) {
            warning(paste("Invalid sd_values in dataset", k))
            next
        }

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

## https://www.magesblog.com/post/2013-03-12-how-to-use-optim-in-r/

## We optimise over the parameter alpha, and specify everything else

################################
## Helper functions for optim ##
fun2optim <- function(par, ...) {
  alpha <- if (length(par) == 1) par else par["alpha"]
  getLogLik_oneCpG_optimized(alpha = alpha, ...)
}

####################################
## Run optimizer for a single CpG ##

runOptim1CpG <- function(scaled_list_mat, mu_jk_list, sigma_k_list,
                         lambdas, p0, p1) {
    resOpt <- optim(
        par = 0.001,  # initial alpha
        fn = fun2optim,
        scaled_list_mat = scaled_list_mat,
        mu_jk_list = mu_jk_list,
        sigma_k_list = sigma_k_list,
        p0 = p0, p1 = p1, lambdas = lambdas,
        method = "Brent",
        lower = 0, upper = 1,
        control = list(fnscale = -1)
    )
    return(resOpt$par)
}

#####################################################
## Parallelise optimizer for all CpG with mclapply ##

getAllOptimAlpha_parallel <- function(cpgvec, NCORES, p0, p1) {

    safe_run <- function(cpgPos) {
        scaled_list_mat <- source_scaled_mat_1CpG(cpgPos)
        
        ## ✅ “Proceed only if there are at least 3 matrices that contain at least one non-NA value.”
        ## Logical vector: TRUE if matrix has at least one non-NA value
        has_data <- sapply(scaled_list_mat, function(m) any(!is.na(m)))

        ## Count how many matrices have real data
        num_with_data <- sum(has_data)

        ## Proceed only if >= 3 have real data
        if (num_with_data < 3) {
            message(sprintf("CpG at index %s not covered in enough (>= 3) datasets", cpgPos))
            return(NA_real_)
        } else {
            mu_jk_list <- lapply(scaled_list_mat, function(mat) rowMeans(mat, na.rm = TRUE))
            tryCatch({
                runOptim1CpG(
                    scaled_list_mat = scaled_list_mat,
                    mu_jk_list = mu_jk_list,
                    sigma_k_list = median_sds,
                    lambdas = lambdas,
                    p0 = p0,
                    p1 = p1
                )
            },
            error = function(e) {
                message(sprintf("CpG %s failed: %s", cpgPos, e$message))
                NA_real_
            })
        }
    }

    ## Do the parallel run OUTSIDE safe_run!
    res <- parallel::mclapply(cpgvec,
                         safe_run,
                         mc.cores = NCORES,
                         mc.preschedule = FALSE
                     )

    ## Name & return matrix
    names(res) <- cpgvec
    my_matrix <- matrix(unlist(res), ncol = 1)
    rownames(my_matrix) <- names(res)
    colnames(my_matrix) <- "alpha"
    return(my_matrix)
}

###################################
## Top-level: run & save results ##

runAndSave <- function(datasetName, cpgvec, resultDir, NCORES, p0, p1) {
    ## Generate a safe object name
    var_cpgvec <- deparse(substitute(cpgvec))
    obj_name <- paste0(
        "results_", datasetName, "_",  var_cpgvec, "_",
        p0, "p0_", p1, "p1")
    obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)

    ## Ensure directory ends with /
    if (!grepl("/$", resultDir)) resultDir <- paste0(resultDir, "/")
    if (!dir.exists(resultDir)) stop("Result directory does not exist: ", resultDir)

    file_name <- paste0(resultDir, obj_name, ".RData")

    ## Check if file exists first
    if (file.exists(file_name)) {
        message("File already exists: ", file_name)
        message("Skipping run to avoid overwriting.")
        return(invisible(NULL))
    }

    ## Run only if file does not exist
    result <- getAllOptimAlpha_parallel(
        cpgvec, NCORES = NCORES, p0 = p0, p1 = p1
    )

    ## Assign to global env so user has it
    assign(obj_name, result, envir = .GlobalEnv)

    message("Saving to file: ", file_name)
    save(list = obj_name, file = file_name, envir = .GlobalEnv)
    message("💾 Result saved successfully.")
}
