#!/usr/bin/env Rscript
# 27th August
# Optimise algo v5

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
message("Load the cpgs...")
cpg_names_all <- h5read(file.path(dataDir, "all_matrix_noscale.h5"), "cpg_names")

###################################
## Which are covered in 26 cells ##
cpg_46 <- read.table(file.path(dataDir, "selected_cpgs_min3_in46_datasets.txt"))$V1

#############################
## Prep hvCpG and controls ##
source("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R")
hvCpGandControls <- prephvCpGandControls(codeDir, cpg_46)
DerakhshanhvCpGs_names_filtered = hvCpGandControls$DerakhshanhvCpGs_names_filtered
mQTLcontrols_names_filtered = hvCpGandControls$mQTLcontrols_names_filtered
dictionary <- hvCpGandControls$dictionary
rm(hvCpGandControls)

########################################################
message(paste0("Run algorithm on the ", length(DerakhshanhvCpGs_names_filtered), " hvCpG covered in 46 cells..."))

resDir = file.path(codeDir, "05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/")

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
# 📦 Loading batch 1 / 2 (1000 CpGs) at 2025-08-21 15:56:00.907993
# 📦 Loading batch 2 / 2 (450 CpGs) at 2025-08-21 17:05:44.917912
# Saving to file: ~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_1450CpGs_0_8p0_0_65p1.RData
# 💾 Result saved successfully.
# user   system  elapsed 
# 9245.875   79.822 5052.667 
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
  dataDir = dataDir)
)

# 📦 Loading batch 1 / 2 (1000 CpGs) at 2025-08-21 14:02:06.679142
# 📦 Loading batch 2 / 2 (474 CpGs) at 2025-08-21 14:20:51.707893
# Saving to file: ~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1.RData
# 💾 Result saved successfully.
## 5 threads, 7G/thread, batch 1000CpGs --> 15min for 1000CpGs

######################
### Explore results ##

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_hvCpGsMariav5_1450CpGs_0_8p0_0_65p1.RData")
reshvCpG = results_hvCpGsMariav5_1450CpGs_0_8p0_0_65p1; rm(results_hvCpGsMariav5_1450CpGs_0_8p0_0_65p1)

load("~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X_hvCpGandControls/results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1.RData")
resControls = results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1; rm(results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1)

summary(reshvCpG)
summary(resControls)

df = rbind(data.frame(alpha = reshvCpG, type = "hvCpG Derakhshan"),
           data.frame(alpha = resControls, type = "mQTL controls"))

ggplot(df, aes(x = type, y = alpha)) +
  geom_jitter(aes(fill=type), pch=21, size = 3, alpha = .1)+ 
  geom_violin(aes(col=type))+
  geom_boxplot(aes(col=type), width = .1) + 
  theme_minimal(base_size = 14) +                                                    
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Probability of being a hvCpG") 

summary(lm(alpha ~ type, data = df))
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.188740   0.005650  33.406  < 2e-16 ***
#   typemQTL controls -0.020765   0.007958  -2.609  0.00912 **  

testCpGs <- c(rownames(head(df[df$type %in% "hvCpG Derakhshan" & df$alpha < 0.1,],1)),
              rownames(head(df[df$type %in% "mQTL controls" & df$alpha < 0.1,],1)),
              rownames(head(df[df$type %in% "hvCpG Derakhshan" & df$alpha > 0.9,],1)),
              rownames(head(df[df$type %in% "mQTL controls" & df$alpha > 0.9,],1)))

###################################################
## Check if results make sense based on raw data ##
source("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/checkRawDataSD.R")
checkRawDataSD(testCpGs)
df <- df[rownames(df) %in% testCpGs, ]
df$alpha <- round(df$alpha, 2)
df

##############################################################
## Let's have a look at the surprising values in array data ##

## hypo variable in atlas/ hyper in array
df$probe = dictionary[match(rownames(df), dictionary$hg38),"illu450k"]

df

##############################################################################
## Median SD (& bootstrapped 95% CI) of hvCpGs and controls in all datasets ##
library(matrixStats)
library(boot)
library(dplyr)

sdselected <- function(cpg_indices, B = 1000) {
  
  ## Extract sub matrix
  mat = h5read("~/Documents/Project_hvCpG/10X/all_matrix_noscale.h5", "matrix", 
               index = list(NULL, cpg_indices))
  rownames(mat) <- h5read("~/Documents/Project_hvCpG/10X/all_matrix_noscale.h5", "samples")
  colnames(mat) <- cpg_names_all[cpg_indices]
  
  metadata = read.table(file.path(dataDir, "sample_metadata.tsv"), sep = "\t", header = TRUE)
  
  results <- metadata %>%
    group_by(dataset) %>%
    group_modify(~{
      samples <- .x$sample
      submat <- mat[rownames(mat) %in% samples, , drop = FALSE]
      
      if (nrow(submat) == 0) {
        return(tibble(median_sd = NA_real_, lower = NA_real_, upper = NA_real_))
      }
      
      row_sds <- rowSds(as.matrix(submat), na.rm = TRUE)
      
      # If all SDs are identical, CI cannot be computed
      if (length(unique(row_sds)) == 1) {
        return(tibble(
          median_sd = unique(row_sds),
          lower = NA_real_,
          upper = NA_real_
        ))
      }
      
      boot_median <- function(data, indices) median(data[indices], na.rm = TRUE)
      boot_res <- boot(data = row_sds, statistic = boot_median, R = B)
      
      ci <- tryCatch(boot.ci(boot_res, type = "perc"), error = function(e) NULL)
      
      if (!is.null(ci) && !is.null(ci$percent)) {
        lower_ci <- ci$percent[4]
        upper_ci <- ci$percent[5]
      } else {
        lower_ci <- NA_real_
        upper_ci <- NA_real_
      }
      
      tibble(
        median_sd = median(row_sds, na.rm = TRUE),
        lower = lower_ci,
        upper = upper_ci
      )
    }) %>%
    ungroup()
  results
}

df <- rbind(data.frame(
  medianSD = sdselected(cpg_indices = match(DerakhshanhvCpGs_names_filtered, cpg_names_all)),
  type="hvCpG"),
  data.frame(
    medianSD = sdselected(cpg_indices = match(mQTLcontrols_names_filtered, cpg_names_all)),
    type="mQTLcontrol"))

ggplot(df, aes(x = medianSD.dataset, y = medianSD.median_sd, fill = type)) +
  geom_pointrange(
    aes(ymin = medianSD.lower, ymax = medianSD.upper),
    shape = 21, size = 0.4, stroke = 0.5
  ) +
  geom_point(pch = 21, size = 2) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = c(0.2, 0.8),
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(title = NULL)) +
  xlab("") +
  ylab("Median SD ± 95% CI")
