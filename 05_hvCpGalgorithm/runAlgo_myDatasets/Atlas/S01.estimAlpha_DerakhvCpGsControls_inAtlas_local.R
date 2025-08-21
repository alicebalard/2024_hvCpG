#!/usr/bin/env Rscript
# 12th August
# Run algo v5 for hvCpGs and controls, without transformation

###############################
message("Define parameters...")
dataDir <- "~/Documents/Project_hvCpG/10X/"
nslots <- 8

message(paste0("We work in ", dataDir, " with ", nslots, " cores."))

codeDir = "~/Documents/GIT/2024_hvCpG"
resDir = file.path(codeDir, "05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/")

##################################
message("Source the algorithm...")
source(file.path(codeDir, "05_hvCpGalgorithm/hvCpG_algorithm_detection_v5batches.R"))

###########################
message("Load the cpgs...")
cpg_names_all <- h5read(file.path(dataDir, "all_matrix_noscale.h5"), "cpg_names")

###################################
## Which are covered in 26 cells ##
cpg_46 <- read.table(file.path(dataDir, "selected_cpgs_min3_in46_datasets.txt"))$V1

#############################
## Prep hvCpG and controls ##
source("runAlgo_myDatasets/Atlas/prephvCpGandControls.R")
hvCpGandControls <- prephvCpGandControls(codeDir, cpg_names_all, cpg_46)
DerakhshanhvCpGs_names_filtered = hvCpGandControls$DerakhshanhvCpGs_names_filtered
mQTLcontrols_names_filtered = hvCpGandControls$mQTLcontrols_names_filtered
rm(hvCpGandControls)

########################################################
message(paste0("Run algorithm on the ", length(DerakhshanhvCpGs_names_filtered), " hvCpG covered in 46 cells..."))

system.time(runAndSave(
  analysis = "hvCpGsMariav5",
  cpg_names_vec = DerakhshanhvCpGs_names_filtered,
  resultDir = resDir,
  NCORES = nslots,
  p0 = 0.80,
  p1 = 0.65, 
  batch_size = 1000,
  dataDir = dataDir, overwrite = T)
)
# New result directory ~/Documents/GIT/2024_hvCpG//05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/ created
# 📦 Loading batch 1 / 2 (1000 CpGs) at 2025-08-21 09:40:05.751693
# 📦 Loading batch 2 / 2 (663 CpGs) at 2025-08-21 09:53:24.189292
# Saving to file: ~/Documents/GIT/2024_hvCpG//05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1.RData
# 💾 Result saved successfully.
###########################################################
message(paste0("Run algorithm on the ", length(mQTLcontrols_names_filtered), " controls covered in 46 cells..."))

system.time(runAndSave(
  analysis = "mQTLcontrolsv5",
  cpg_names_vec = mQTLcontrols_names_filtered,
  resultDir = resDir,
  NCORES = nslots,
  p0 = 0.80,
  p1 = 0.65, 
  batch_size = 1000,
  dataDir = dataDir, overwrite = T)
)
# 📦 Loading batch 1 / 2 (1000 CpGs) at 2025-08-21 10:04:19.439537
# 📦 Loading batch 2 / 2 (474 CpGs) at 2025-08-21 10:17:43.350317
# Saving to file: ~/Documents/GIT/2024_hvCpG//05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1.RData
# 💾 Result saved successfully.

######################
### Explore results ##

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1.RData")
reshvCpG = results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1; rm(results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1)

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1.RData")
resControls = results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1; rm(results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1)

df = rbind(data.frame(alpha = reshvCpG, type = "hvCpG Derakhshan"),
           data.frame(alpha = resControls, type = "mQTL controls"))

ggplot(df, aes(x = type, y = alpha)) +
  geom_jitter(aes(fill=type), pch=21, size = 3, alpha = .1)+ 
  geom_violin(aes(col=type))+
  geom_boxplot(aes(col=type), width = .1) + 
  theme_minimal(base_size = 14) +                                                    
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Probability of being a hvCpG") 

reshvCpG[rownames(reshvCpG) %in% c("chr22_27438478-27438479", "chr1_1944364-1944365"),]
# chr22_27438478-27438479   6.474096e-09 hvCpG Derakhshan

## Test 1 CpG at a time
system.time(runAndSave(
  analysis = "mQTLcontrolsv5",
  cpgPos_vec = which(cpg_names_all %in% "chr22_27438478-27438479"),
  resultDir = resDir,
  NCORES = nslots,
  p0 = 0.80,
  p1 = 0.65, 
  batch_size = 1000,
  dataDir = dataDir)
)

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_mQTLcontrolsv5_1CpGs_0_8p0_0_65p1.RData")
results_mQTLcontrolsv5_1CpGs_0_8p0_0_65p1
# alpha
# chr22_27438478-27438479 0.4418279

#### ISSUE WHEN DONE IN BATCHES!!
head(DerakhshanhvCpGs_positions, 50)

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_50CpGs_0_8p0_0_65p1.RData")
results_hvCpGsMariav5_50CpGs_0_8p0_0_65p1

head(cbind(head(reshvCpG, 50), results_hvCpGsMariav5_50CpGs_0_8p0_0_65p1),10)
# alpha      alpha
# chr1_1944341-1944342     0.4806861 0.48068614
# chr1_1944364-1944365     0.7561706 0.75617061
# chr1_41918976-41918977   0.3396265 0.33962655
# chr1_43007301-43007302   0.5823137 0.58231367
# chr1_44213119-44213120   0.6402617 0.64026169
# chr1_58248633-58248634   0.6443377 0.17959387
# chr1_58577536-58577537   0.5371961 0.55077547
# chr1_62783528-62783529   0.5018367 0.12677161
# chr1_62783538-62783539   0.4238029 0.55125701

system.time(runAndSave(
  analysis = "hvCpGsMariav5",
  cpg_names_vec = c("chr1_41918976-41918977","chr1_1944341-1944342", "chr1_62783538-62783539"),
  resultDir = resDir,
  NCORES = nslots,
  p0 = 0.80,
  p1 = 0.65, 
  batch_size = 1000,
  dataDir = dataDir, overwrite = T)
)

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_3CpGs_0_8p0_0_65p1.RData")
results_hvCpGsMariav5_3CpGs_0_8p0_0_65p1

# alpha
# chr1_1944341-1944342   0.4806861
# chr1_62783538-62783539 0.4985959
# chr1_41918976-41918977 0.1795939

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_100CpGs_0_8p0_0_65p1.RData")
results_hvCpGsMariav5_100CpGs_0_8p0_0_65p1

cbind(head(results_hvCpGsMariav5_100CpGs_0_8p0_0_65p1, 50), results_hvCpGsMariav5_50CpGs_0_8p0_0_65p1)
