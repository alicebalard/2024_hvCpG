setwd("~/2024_hvCpG/")

## Load data & functions specific to Maria's arrays (to do before)
system.time(source("03_prepDatasetsMaria/S02.formatArraysforR.R"))

## Load algorithm (30sec)
system.time(source("05_hvCpGalgorithm/hvCpG_algorithm_detection_v3.R"))

length(cpg_names); head(cpg_names)
## [1] 406036
# [1] "cg00000029" "cg00000108" "cg00000109" "cg00000165" "cg00000236" "cg00000289"

## Load hvCpG of Maria and matching mQTL
cistrans_GoDMC_hvCpG_matched_control <- read.table("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

## Load full cpg names
cpg_names <- readLines("~/arraysh5files/sorted_common_cpgs.txt")

## In which positions are they?
cpg_names[cpg_names %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
cpg_names[cpg_names %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length

sub_cistrans_GoDMC_hvCpG_matched_control <- cistrans_GoDMC_hvCpG_matched_control[
  cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% cpg_names &
    cistrans_GoDMC_hvCpG_matched_control$controlCpG_name %in% cpg_names,]

## Positions of targets in cpg_names:
match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names)
match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names)

source_scaled_mat_1CpG(pos = match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names)[2])

pos2check <- c(match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names),
               match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names))

system.time(runAndSave("Maria", cpgvec = pos2check, p0=0.80, p1=0.65, NCORES=30, resultDir="05_hvCpGalgorithm/resultsDir/Mariads/"))
## 774 sec = 3 minutes

# Your grid
params <- seq(0.4, 1, 0.05)
grid <- expand.grid(p0 = params, p1 = params)

apply(grid, 1, function(row) {
  runAndSave("Maria", cpgvec = pos2check, p0=as.numeric(row["p0"]), p1=as.numeric(row["p1"]), NCORES=30,
             resultDir="05_hvCpGalgorithm/resultsDir/Mariads/")
})





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
