#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

data_dir <- args[1] ## where the data is
task_id <- as.integer(args[2]) ## which task on the array
chunk_size <- as.integer(args[3]) ## what is the size of the chunk
batch_size <- as.integer(args[4]) ## how many CpGs are loaded at once

## Run on CpGs which are covered in all 46 cells 
cpg_46 <- read.table(file.path(data_dir, "selected_cpgs_min3_in46_datasets.txt"))$V1

## Batch
start_idx <- (task_id - 1) * chunk_size + 1
end_idx <- min(task_id * chunk_size, length(cpg_46))

if (start_idx > length(cpg_46)) {
  stop("Start index beyond end of CpG list.")
}

subset_cpgs <- cpg_46[start_idx:end_idx]

#############******Prepare 4 subsets of data*********#############

# Load metadata to use for the algorithm
metadata <- read.table(file.path(data_dir, "sample_metadata.tsv"), sep = "\t", header = TRUE)

# Load data with sex
meta2 <- read.csv("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/SupTab1_Loyfer2023.csv")

# Create the dataset identifier
meta2$CompositeGroup <- paste(meta2$Source.Tissue, meta2$Cell.type, sep = " - ")

##***## Case 1: only males datasets

## Identify which groups contain only males
only_male_groups <- meta2 %>%
  dplyr::group_by(CompositeGroup) %>%
  dplyr::summarise(all_male = all(sex == "M")) %>%
  dplyr::filter(all_male) %>%
  dplyr::pull(CompositeGroup)

# Subset metadata to those groups
meta_male_only <- meta2 %>%
  dplyr::filter(CompositeGroup %in% only_male_groups)

# Extract short sample IDs (the Z000000T7 part)
meta_male_only$ShortSampleID <- sub(".*-(Z[0-9A-Za-z]+)$", "\\1", meta_male_only$Sample.name)

mySubsetMetadata_1 <- metadata[metadata$sample %in% meta_male_only$ShortSampleID,]

##***## Case 2: all but only males datasets

# Subset metadata to those groups
meta_no_male_only <- meta2 %>%
  dplyr::filter(!CompositeGroup %in% only_male_groups)

# Extract short sample IDs (the Z000000T7 part)
meta_no_male_only$ShortSampleID <- sub(".*-(Z[0-9A-Za-z]+)$", "\\1", meta_no_male_only$Sample.name)

mySubsetMetadata_2 <- metadata[metadata$sample %in% meta_no_male_only$ShortSampleID,]

##***## Case 3: only females datasets

## Identify which groups contain only females
only_female_groups <- meta2 %>%
  dplyr::group_by(CompositeGroup) %>%
  dplyr::summarise(all_female = all(sex == "F")) %>%
  dplyr::filter(all_female) %>%
  dplyr::pull(CompositeGroup)

# Subset metadata to those groups
meta_female_only <- meta2 %>%
  dplyr::filter(CompositeGroup %in% only_female_groups)

# Extract short sample IDs (the Z000000T7 part)
meta_female_only$ShortSampleID <- sub(".*-(Z[0-9A-Za-z]+)$", "\\1", meta_female_only$Sample.name)

mySubsetMetadata_3 <- metadata[metadata$sample %in% meta_female_only$ShortSampleID,]

##***## Case 4: all but only females datasets

# Subset metadata to those groups
meta_no_female_only <- meta2 %>%
  dplyr::filter(!CompositeGroup %in% only_female_groups)

# Extract short sample IDs (the Z000000T7 part)
meta_no_female_only$ShortSampleID <- sub(".*-(Z[0-9A-Za-z]+)$", "\\1", meta_no_female_only$Sample.name)

mySubsetMetadata_4 <- metadata[metadata$sample %in% meta_no_female_only$ShortSampleID,]

#############***************#############

######### Run for 4 subsets of data #########

runByGroup <- function(mySubsetMetadata, subdir){

    make_result_dir <- function(subdir, task_id) {
        file.path(
            "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas",
            subdir,
            sprintf("Atlas_batch%03d", task_id)
        )
    }
    
    ## Output directory
    result_dir <- make_result_dir("AtlasMalesOnly", task_id)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

    message(paste0("If new, results will be saved in dir: ", result_dir))

    ## Load algorithm
    source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/hvCpG_algorithm_detection_v6.R")

    ## Run
    myNthreads <- as.numeric(Sys.getenv("NSLOTS", unset = "1"))  # Use all cores

    message("Run algo:")

    system.time(runAndSave_fast(
        analysis = "Atlas10X",
        cpg_names_vec = subset_cpgs,
        dataDir = data_dir,
        resultDir = result_dir,
        NCORES = myNthreads,
        p0 = 0.80,
        p1 = 0.65,
        batch_size = batch_size,
        subsetMetadata = mySubsetMetadata)
        )
}

runByGroup(mySubsetMetadata = mySubsetMetadata_1, subdir = "Atlas10X_onlymales")

runByGroup(mySubsetMetadata = mySubsetMetadata_2, subdir = "Atlas10X_allbutonlymales")

runByGroup(mySubsetMetadata = mySubsetMetadata_3, subdir = "Atlas10X_onlyfemales")

runByGroup(mySubsetMetadata = mySubsetMetadata_4, subdir = "Atlas10X_allbutonlyfemales")
