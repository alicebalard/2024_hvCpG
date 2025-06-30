###########
## NB: this script has to run on R outside of Rstudio, for compatibility issues with parallelisation

setwd("/home/alice/2024_hvCpG/")

## Maria's data (Hosted on LSHTM server of Matt Silver)
source("03_prepDatasetsMaria/dataprep_MariaArrays.R")
my_list_mat_Mariads <- Maria_filtered_list_mat; rm(Maria_filtered_list_mat)
cpgnames <- unique(unlist(sapply(my_list_mat_Mariads, row.names)))
cpgnames <- cpgnames[order(cpgnames)]

cat(paste0("Prepare Maria's dataset and source functions for optimisation:\n"))
source("05_hvCpGalgorithm/hvCpG_algorithm_detection_v1.R")
cat("Functions sourced.\n")

## Load mQTL-matched controls
cistrans_GoDMC_hvCpG_matched_control <- read.table("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

table(cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% MariasCpGs$CpG)
cpgnames[cpgnames %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
cpgnames[cpgnames %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length

##########################################################################
## Part 1. Run hvCpG on 1500 hvCpG from Maria and matched mQTL controls ##
##########################################################################

set.seed(123)
sub1500 <- cistrans_GoDMC_hvCpG_matched_control[sample(1:nrow(cistrans_GoDMC_hvCpG_matched_control), 1500), ]
test3000CpGsvec <- c(sub1500$hvCpG_name, sub1500$controlCpG_name)

## Nelder-Mead or L-BFGS-B?
runAndSave(my_list_mat = my_list_mat_Mariads, cpgvec = test3000CpGsvec,
           optimMeth="Nelder-Mead", NCORES=8, p0=0.95, p1=0.65, resultDir="/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/")

runAndSave(my_list_mat = my_list_mat_Mariads, cpgvec = test3000CpGsvec,
           optimMeth="L-BFGS-B", NCORES=8, p0=0.95, p1=0.65, resultDir="/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/")
## We choose Nelder-Mead, as L-BFGS-B produces weird zero values (see Rmd plots).

## Vary p0 true negative and p1 true positive:
params <- seq(0.4, 1, 0.05)
grid <- expand.grid(p0 = params, p1 = params)

apply(grid, 1, function(row) {
    runAndSave(
        my_list_mat = my_list_mat_Mariads,
        cpgvec = test3000CpGsvec,
        optimMeth = "Nelder-Mead",
        NCORES = 8,
        p0 = as.numeric(row["p0"]),
        p1 = as.numeric(row["p1"]),
        resultDir = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/"
    )
})

## On Rstudio
## myres <- readRDS("/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/results_Nelder-Mead_minitest_0.95.RDS")
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
runAndSave(
  my_list_mat = my_list_mat_Mariads_mimicAtlas,
  cpgvec = test3000CpGsvec,
  optimMeth = "Nelder-Mead",
  NCORES = 8,
  p0 = 0.45,
  p1 = 0.6,
  resultDir = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/"
)
