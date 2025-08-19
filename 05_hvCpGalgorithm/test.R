source("loadMyLibs.R")

##################################
message("Source the algorithm...")
source("hvCpG_algorithm_detection_v5batches.R")

###########################
message("Load the cpgs...")
cpg_names_all <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                        "cpg_names")

resDir = "/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/"

indexVec = which(cpg_names_all %in% c("chr1_1944364-1944365", "chr1_778669-778670", 
                                      "chr11_86301834-86301835", "chr22_27438478-27438479"))

system.time(runAndSave(
  analysis = "test",
  cpgPos_vec = indexVec,
  resultDir = resDir,
  NCORES = 1,
  p0 = 0.80,
  p1 = 0.65, 
  batch_size = 100,
  dataDir = "/home/alice/Documents/Project_hvCpG/10X/", overwrite = T)
)

load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_test_4CpGs_0_8p0_0_65p1.RData")
res = results_test_4CpGs_0_8p0_0_65p1

round(res,3)

samples <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                        "samples")

sample_groups <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                           "sample_groups")

subRawData <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
       "matrix", index = list(NULL, indexVec))

colnames(subRawData) = cpg_names_all[indexVec]

subRawData <- as.data.frame(subRawData)

subRawData$samples = samples
subRawData$sample_groups = sample_groups

# Exclude metadata columns
num_cols <- setdiff(names(subRawData), c("samples", "sample_groups"))

# Step 1: compute SD per group per CpG site
df_long <- subRawData %>%
  select(all_of(num_cols), sample_groups) %>%
  pivot_longer(cols = -sample_groups, names_to = "CpG", values_to = "value") %>%
  group_by(sample_groups, CpG) %>%
  summarise(sd = sd(value, na.rm = TRUE), .groups = "drop")

# Step 2: median of SDs per group
med <- df_long %>%
  group_by(CpG) %>%
  summarise(median_sd = median(sd, na.rm = TRUE))

# Step 3: histogram of SD distributions per group
ggplot(df_long, aes(x = sd)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ CpG, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of SD per CpG",
       x = "Standard deviation across samples",
       y = "Count")+
  geom_vline(data = med, aes(xintercept = median_sd))

round(res,3)

## Check if my full run is correct
load("resultsDir/Atlas10X/Atlas_batch001/results_Atlas10X_100000CpGs_0_8p0_0_65p1.RData")

rownames(results_Atlas_100000CpGs_0_8p0_0_65p1)[indexVec]
results_Atlas_100000CpGs_0_8p0_0_65p1[indexVec]


load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1.RData")
reshvCpG = results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1; rm(results_hvCpGsMariav5_1663CpGs_0_8p0_0_65p1)

load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1.RData")
resControls = results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1; rm(results_mQTLcontrolsv5_1474CpGs_0_8p0_0_65p1)

reshvCpG[rownames(reshvCpG) %in% cpg_names_all[indexVec],]

resControls[rownames(resControls) %in% cpg_names_all[indexVec],]
