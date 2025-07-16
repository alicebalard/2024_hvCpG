#!/bin/bash
#$ -N runhvCpGMaria
#$ -S /bin/bash
#$ -pe smp 20
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=10:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

## If needed, to prepare the files
## source /share/apps/source_files/python/python-3.13.0a6.source
## python3 /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/03_prepDatasetsMaria/S01.prepare_h5files_MariaArrays.py

R --vanilla <<EOF

myNthreads=20 ## specify here

## Load algorithm (30sec)
system.time(source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/hvCpG_algorithm_detection_v3.R"))

## Load data & functions specific to Atlas (2 sec)
system.time(source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/S02.formatAtlasforR.R"))




## Load cpg list (~1 min 30) 29,401,795 CpG names
system.time(cpg_names <- h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets_prepared/CpG_names.h5", "cpg_names"))

length(cpg_names); head(cpg_names)
##[1] 29401795
##[1] "chr1_10469-10470" "chr1_10471-10472" "chr1_10484-10485" "chr1_10489-10490"
##[5] "chr1_10493-10494" "chr1_10497-10498"

system.time(runAndSave("Atlas", cpgvec = head(cpg_names, 10000),p0=0.80, p1=0.65, NCORES=myNthreads,
		       resultDir="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas/"))

## 100, 5G, 1C = 87 sec
## 1000, 5G, 1C = 471.074 sec
## 1000 CpGs, 4G, 16C = 57 sec
## 10000 CpGs, 4G, 25C = 759 sec
## 10000 CpGs, 5G, 20C = 3806 sec (but coverage is now 10 min, not 20 min, so much more to process!)

message("Done!")

q()
EOF



############ STOP!
                                                                                          
###########
## NB: this script has to run on R outside of Rstudio, for compatibility issues with parallelisation

## Specific libPath for R --vanilla
.libPaths(c(
  "/opt/R/packages/lib_4.4.1",
  "/home/alice/R/x86_64-pc-linux-gnu-library/4.4",
  "/usr/lib/R/library"
))

setwd("/home/alice/2024_hvCpG/")

## Maria's data (Hosted on LSHTM server of Matt Silver)
source("03_prepDatasetsMaria/dataprep_MariaArrays.R")

cat(paste0("Prepare Maria's dataset and source functions (v1) for optimisation:\n"))
source("05_hvCpGalgorithm/hvCpG_algorithm_detection_v1.R")
cat("Functions sourced.\n")

## Set up my number of threads
myNthreads=50

##########################################################################
## Part 1. Run hvCpG on 3453 hvCpG from Maria and matched mQTL controls ##
##########################################################################
hvCpGandControls <- c(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, 
                     sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name)

## Nelder-Mead or L-BFGS-B?
runAndSave(my_list_mat = my_list_mat_Mariads, cpgvec = hvCpGandControls,
           optimMeth="Nelder-Mead", NCORES=myNthreads, p0=0.95, p1=0.65, 
           resultDir="/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/")

runAndSave(my_list_mat = my_list_mat_Mariads, cpgvec = hvCpGandControls,
           optimMeth="L-BFGS-B", NCORES=myNthreads, p0=0.95, p1=0.65, 
           resultDir="/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/")
## We choose Nelder-Mead, as L-BFGS-B produces weird zero values (see Rmd plots).

library(foreach)
library(doParallel)

# Number of *parallel tasks* (each task uses 5 cores)
NWORKERS <- 10

cl <- makeCluster(NWORKERS)
registerDoParallel(cl)

# Your grid
params <- seq(0.4, 1, 0.05)
grid <- expand.grid(p0 = params, p1 = params)

# No need to assign to `results` — just run for side effects
foreach(i = 1:nrow(grid), .packages = packages) %dopar% {
  row <- grid[i, ]
  runAndSave(
    my_list_mat = my_list_mat_Mariads,
    cpgvec = hvCpGandControls,
    optimMeth = "Nelder-Mead",
    NCORES = 5,
    p0 = row$p0,
    p1 = row$p1,
    resultDir = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/"
  )
}

stopCluster(cl)

## On Rstudio
## myres <- readRDS("/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/results_Nelder-Mead_minitest_0.95.RDS")
## head(myres)
## nrow(myres)

######################################################################
## Part 2. Subset Maria's data following Atlas data, to check power ##
######################################################################

## Extract column counts from Atlas data
## (done in S02 as it needs loads of ram and time!)
## n_matrices <- length(my_list_mat_Mariads)
## ncols_list1 <- sapply(my_list_mat_Mariads, ncol)
n_matrices <- 33
ncols_list1 <- c(3, 5, 5, 3, 11, 3, 22, 4, 3, 8, 19, 3, 3, 11, 13, 4, 4, 8, 6, 4, 3, 10, 4, 4, 4, 3, 3, 3, 4, 4, 5, 5, 3)

## Track which matrices in my_list_mat_Mariads have enough columns
eligible_indices <- lapply(ncols_list1, function(nc) {
  which(sapply(my_list_mat_Mariads, function(m) ncol(m) >= nc))
})

## Create the output list
## Create the output list with names preserved
my_list_mat_Mariads_mimicAtlas <- vector("list", n_matrices)
used_indices <- integer()

for (i in seq_len(n_matrices)) {
  options <- setdiff(eligible_indices[[i]], used_indices)
  if (length(options) == 0) {
    stop(sprintf("No more available matrices for required %d columns at position %d", ncols_list1[i], i))
  }
  chosen_idx <- sample(options, 1)
  used_indices <- c(used_indices, chosen_idx)
  mat <- my_list_mat_Mariads[[chosen_idx]]
  selected_cols <- sample(ncol(mat), ncols_list1[i])
  my_list_mat_Mariads_mimicAtlas[[i]] <- mat[, selected_cols, drop = FALSE]
  ## Set the name in the new list to match the chosen matrix name
  names(my_list_mat_Mariads_mimicAtlas)[i] <- names(my_list_mat_Mariads)[chosen_idx]
}

## Check
sapply(my_list_mat_Mariads_mimicAtlas, ncol)

## Run the algorithm on the reduced matrices: will we manage to find them?

## Define your fixed list of (p0, p1) pairs
my_p_combinations <- list(c(0.45, 0.6),## best AUC, closest to Maria's data; low specificity and mid sensitivity
                          c(0.95, 0.65), ## my original idea, high specificity but mid sensitivity
                          c(0.9, 0.9)) ## high specificity and sensitivity
# Loop over each pair
for (pair in my_p_combinations) {
  p0_val <- pair[1]
  p1_val <- pair[2]

  message(sprintf("Running for p0 = %.2f | p1 = %.2f", p0_val, p1_val))

  runAndSave(
    my_list_mat = my_list_mat_Mariads_mimicAtlas,
    cpgvec = test3000CpGsvec,
    optimMeth = "Nelder-Mead",
    NCORES = myNthreads,
    p0 = p0_val,
    p1 = p1_val,
    resultDir = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/"
  )
}
