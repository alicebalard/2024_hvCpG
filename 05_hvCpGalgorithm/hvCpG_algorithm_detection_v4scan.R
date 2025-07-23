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

## <<- assigns the variable to the global environment (or parent environment) rather than the local function scope. Be careful with it!
prepData <- function(analysis){
    if (analysis == "Atlas"){
        metadata <<- read.table("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/sample_metadata.tsv", sep ="\t", header = TRUE)

        medsd_lambdas <<- read.table("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/all_medsd_lambda.tsv", sep = "\t", header = TRUE)

        ## Function of one CpG of interest:
        source_M_1CpG <<- function(cpgpos) {
            M = rhdf5::h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/all_scaled_matrix.h5", "scaled_matrix", index = list(NULL, cpgpos))
            rownames(M) = metadata$sample
            return(M)
        }

        ## Names of all CpG in the dataset:
        cpg_names_all <<- rhdf5::h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/all_scaled_matrix.h5", "cpg_names")
    }
}

##############
## Run algo ##

# Maximum likelihood analysis
## The equation is: $$ \log\left(P(M_j)\right) = \sum_{i=1}^{n} \log\left( \sum_{Z_j=0}^{1} \left( \sum_{Z_{j,k}=0}^{1} P(M_{i,j} \mid Z_{j,k}) \times P(Z_{j,k} \mid Z_j) \right) \times P(Z_j) \right) $$             

###########################################
## Likelihood function for a given CpG j ##

getLogLik_oneCpG_optimized <- function(Mdf, metadata, medsd_lambdas, p0, p1, alpha){
    datasets <- unique(metadata$dataset)
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
        sd_values <- c(sd_k, lambda_k * sd_k)

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
fun2optim <- function(par, ...) {
  alpha <- if (length(par) == 1) par else par["alpha"]
  getLogLik_oneCpG_optimized(alpha = alpha, ...)
}

####################################
## Run optimizer for a single CpG ##

runOptim1CpG <- function(Mdf, metadata, medsd_lambdas, p0, p1) {
    resOpt <- optim(
        par = 0.001,  # initial alpha
        fn = fun2optim,
        Mdf = Mdf, metadata = metadata,
        medsd_lambdas = medsd_lambdas,
        p0 = p0, p1 = p1,
        method = "Brent",
        lower = 0, upper = 1,
        control = list(fnscale = -1)
    )
    return(resOpt$par)
}

#####################################################
## Parallelise optimizer for all CpG with mclapply ##

getAllOptimAlpha_parallel <- function(cpgPos_vec, NCORES, p0, p1) {

    safe_run <- function(cpgPos) {
        Mdf = source_M_1CpG(cpgPos)

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

runAndSave <- function(analysis, cpgPos_vec, resultDir, NCORES, p0, p1) {
    ## Prepare data in the environment:
    prepData(analysis)

    ## Generate a safe object name
    obj_name <- paste0("results_", analysis, "_", length(cpgPos_vec), "CpGs_", p0, "p0_", p1, "p1")
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
    result <- getAllOptimAlpha_parallel(cpgPos_vec = cpgPos_vec, NCORES = NCORES, p0 = p0, p1 = p1)

    ## Assign to global env so user has it
    assign(obj_name, result, envir = .GlobalEnv)

    message("Saving to file: ", file_name)
    save(list = obj_name, file = file_name, envir = .GlobalEnv)
    message("💾 Result saved successfully.")
}
