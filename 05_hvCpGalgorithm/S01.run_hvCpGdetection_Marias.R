setwd("~/2024_hvCpG/")

## Load algorithm
system.time(source("05_hvCpGalgorithm/hvCpG_algorithm_detection_v4scan.R"))

## Load hvCpG of Maria and matching mQTL
cistrans_GoDMC_hvCpG_matched_control <- read.table("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

## Names of all CpG in the dataset:                                                                                                                                                   
cpg_names_all <<- rhdf5::h5read(paste0("/home/alice/arraysh5/all_scaled_matrix.h5"), "cpg_names")
length(cpg_names_all); head(cpg_names_all)
## [1] 406036
# [1] "cg00000029" "cg00000108" "cg00000109" "cg00000165" "cg00000236" "cg00000289"

## In which positions are they?
cpg_names_all[cpg_names_all %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
cpg_names_all[cpg_names_all %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length

sub_cistrans_GoDMC_hvCpG_matched_control <- cistrans_GoDMC_hvCpG_matched_control[
  cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% cpg_names_all &
    cistrans_GoDMC_hvCpG_matched_control$controlCpG_name %in% cpg_names_all,]

## Positions of targets in cpg_names:
match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all)
match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names_all)

## test
source_M_1CpG(cpgpos = match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all)[2])

pos2check <- c(match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all),
               match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names_all))

system.time(runAndSave(analysis = "Maria", cpgPos_vec = pos2check, resultDir="05_hvCpGalgorithm/resultsDir/Mariads/",
                       NCORES=30, p0=0.80, p1=0.65))
## algo v3: 774 sec = 13 minutes
## algo v4: 120.994 = 2 minutes. 6 times faster.

##################################
## Comparisons WORKS (rm since) ##
##################################
##load("05_hvCpGalgorithm/resultsDir/Mariads/test/results_Maria_6906CpGs_0_8p0_0_65p1.RData")
##resV4 = results_Maria_6906CpGs_0_8p0_0_65p1
##load("05_hvCpGalgorithm/resultsDir/Mariads/results_Maria_pos2check_0_8p0_0_65p1.RData")
##resV3 = results_Maria_pos2check_0_8p0_0_65p1
##table(round(resV3,3) == round(resV4,3))
#### TRUE 
## 6906 

## Your grid
params <- seq(0.4, 1, 0.05)
grid <- expand.grid(p0 = params, p1 = params)

apply(grid, 1, function(row) {
  runAndSave(analysis = "Maria", cpgPos_vec = pos2check, resultDir="05_hvCpGalgorithm/resultsDir/Mariads/",
                       NCORES=30, p0=as.numeric(row["p0"]), p1=as.numeric(row["p1"]))
})

## On Rstudio
## myres <- readRDS("/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/results_Nelder-Mead_minitest_0.95.RDS")
## head(myres)
## nrow(myres)

######################################################################
## Part 2. Subset Maria's data following Atlas data, to check power ##
######################################################################

## In Atlas data, a lot of CpG have only 3 datasets with 3 samples

###############
## Done before: python3 S02.prepare_REDUCED_arrays_mimicAtlas.py 
###############

for (s in 3:5){
    for (d in 3:30){
        if (file.exists(paste0("/home/alice/arraysh5_reducedMimicAtlas_", s, "samples_", d, "datasets/all_scaled_matrix.h5")))
        system.time(runAndSave(analysis = paste0("MariasarraysREDUCED_", s, "samples_", d, "datasets"), cpgPos_vec = pos2check,
                               resultDir="05_hvCpGalgorithm/resultsDir/Mariads/", NCORES=30, p0=0.80, p1=0.65))
    }
}

for (d in 3:30){
    if (file.exists(paste0("/home/alice/arraysh5_reducedMimicAtlas_allsamples_", d, "datasets/all_scaled_matrix.h5")))
        system.time(runAndSave(analysis = paste0("MariasarraysREDUCED_allsamples_", d, "datasets"), cpgPos_vec = pos2check,
                               resultDir="05_hvCpGalgorithm/resultsDir/Mariads/", NCORES=30, p0=0.80, p1=0.65))
}


###############
## To do after:
###############

## rm -rf /home/alice/arraysh5_reducedMimicAtlas_*samples*

