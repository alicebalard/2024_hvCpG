#!/usr/bin/env Rscript
# 27th August
# Optimise algo v5 --> v6

###############################
message("Define parameters...")
dataDir <- "~/Documents/Project_hvCpG/10X/"
nslots <- 5

message(paste0("We work in ", dataDir, " with ", nslots, " cores."))
codeDir = "~/Documents/GIT/2024_hvCpG"
setwd(codeDir)

##################################
message("Source the algorithm...")
source(file.path(codeDir, "05_hvCpGalgorithm/hvCpG_algorithm_detection_v5batches.R"))

###########################
message("Load the cpgs which are covered in 26 cells...")
cpg_46 <- read.table(file.path(dataDir, "selected_cpgs_min3_in46_datasets.txt"))$V1

########################################################
message(paste0("Run algorithm on ", length(head(cpg_46)), " CpGs..."))

resDir = file.path(codeDir, "05_hvCpGalgorithm/resultsDir/Atlas10X_test/")

system.time(runAndSave_fast(
  analysis = "test",
  cpg_names_vec = head(cpg_46,100),
  resultDir = resDir,
  NCORES = 5,
  p0 = 0.80,
  p1 = 0.65, 
  batch_size = 100,
  dataDir = dataDir, overwrite = T)
)

## 100 CpGs
## 75.046 
## now 68.505

load("05_hvCpGalgorithm/resultsDir/Atlas10X_test/results_test_100CpGs_0_8p0_0_65p1.RData")
head(results_test_100CpGs_0_8p0_0_65p1,10)
# alpha
# chr1_17452-17453   1.230424e-01
# chr1_17478-17479   6.474096e-09
# chr1_17483-17484   6.474096e-09
# chr1_17492-17493   6.474096e-09
# chr1_17562-17563   1.507041e-02
# chr1_17571-17572   6.474096e-09
# chr1_99166-99167   5.367075e-01
# chr1_101831-101832 1.319087e-01
# chr1_125414-125415 5.887608e-01
# chr1_127492-127493 3.011410e-01

# elapsed
# 33.187 
# 17.305 
# 17.592
# 16.945
## 100 CpGs
# 59.712 
#55.982
# 54.156

## 10CpGs
# 24.967 
# 24.161 
#19.087 