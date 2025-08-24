setwd("~/2024_hvCpG/")

## Load algorithm
source("05_hvCpGalgorithm/hvCpG_algorithm_detection_v5batches.R")

## Names of all CpG in the dataset:                                                                                                                                                   
cpg_names_all <<- rhdf5::h5read(paste0("/home/alice/arraysh5/all_matrix_noscale.h5"), "cpg_names")
length(cpg_names_all); head(cpg_names_all)
## [1] 394240
## [1] "cg00000029" "cg00000108" "cg00000109" "cg00000165" "cg00000236"
## [6] "cg00000289"

system.time(runAndSave(analysis = "arrayAll_algov5",
                       cpg_names_vec = cpg_names_all,
                       dataDir = "/home/alice/arraysh5",
                       resultDir="05_hvCpGalgorithm/resultsDir/Arrays/",
                       p0=0.80, p1=0.65, NCORES=30, batch_size = 10000))

#Saving to file: 05_hvCpGalgorithm/resultsDir/Arrays/results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1.RData
#💾 Result saved successfully.
#      user     system    elapsed 
#797097.912   2733.189  28188.744 
