#!/bin/bash
#$ -N prepAtlasData
#$ -S /bin/bash
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -pe smp 48
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

## 17 June 2025

## Alice Balard
## Preparation of data from the article: https://www.nature.com/articles/s41586-022-05580-6#Sec35

## Input data for our algorithm needs to be in the format:
## filtered_list_mat$Blood_Cauc %>% head
##            GSM1870951  GSM1870952  GSM1870953 GSM1870954 GSM1870955  GSM1870956 GSM1870957
## cg00000957 0.88650894 0.873786294 0.905683210 0.89756560 0.93243437 0.934731380 0.90938323
## cg00001349 0.88900837 0.849494403 0.899719641 0.83146357 0.85704137 0.872804810 0.91467193
## cg00001583 0.06226518 0.059181490 0.050006552 0.07168150 0.04283758 0.044865908 0.03723361

## In beta files, "chr1 10468 10470 3 18" means a base in 10468, then C in 10469 then G in 10470

## read beta files in different languages: https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md
## beta file is the simplest representation of methylation data. It is a binary file with a fixed size of 2 * 8 * NR_SITES bytes (~54MB for hg19), holding a matrix of uint8 values with dimensions of (NR_SITES x 2). For each of the NR_SITES CpG sites in the genome, it holds 2 values: the #meth and #covered. Meaning, the i'th row in the matrix corresponds to the i'th CpG site:
## * meth: the number of times the i'th site (CpGi) site is observed in a methylated state.
## * coverage: the total number of times i'th site (CpGi) is observed. #coverage==0 is equivalent to a missing value (NaN).
## CpGi's beta value is obtained by dividing #meth / #coverage.

## ❗Note: this uint8 format limits values to range [0,255]. In case a CpG site appears over 255 times in the raw data (e.g. bam), its representation is normalized to this range. For example, if site CpG100 appeared 510 times, from which 100 times it was methylated, the 100'th row will be (50, 255). We find it negligible for most cases, but if you find it important, use the lbeta format (wgbstools pat2beta --lbeta ARGS...)

echo "Starting R inline script..."

Rscript --vanilla - << 'EOF'

library(data.table)

DIR <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/"
OUTDIR <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets"

## Read CpG names (only once)
cpg_names <- fread(paste0(DIR, "hg38CpGpos_Loyfer2023.txt"), header = FALSE)[[1]]

## Read metadata as data.table for speed and efficiency
metadata <- fread(paste0(DIR, "SupTab1_Loyfer2023.csv"))
names(metadata) <- gsub(" ", ".", names(metadata))

## Select organs with at least 3 people
metadata <- metadata[, .SD[.N >= 3], by = Group]

## List all beta files
files <- list.files(paste0(DIR, "betaFiles"), full.names = TRUE)
metadata[, file := files[match(sub(".*-", "", metadata$Sample.name), sub(".hg38.beta", "", sub(".*-", "", files)))]]
metadata <- metadata[!is.na(file)]

## Functions to read beta and coverage
read_beta <- function(file) {
    file_size <- file.info(file)$size
    raw <- readBin(file, "integer", n = file_size, size = 1, signed = FALSE)
    mat <- matrix(raw, ncol = 2, byrow = TRUE)
    beta <- mat[,1] / mat[,2]
    beta[is.nan(beta) | is.infinite(beta)] <- NA_real_
    return(beta)
}

read_coverage <- function(file) {
    file_size <- file.info(file)$size
    raw <- readBin(file, "integer", n = file_size, size = 1, signed = FALSE)
    mat <- matrix(raw, ncol = 2, byrow = TRUE)
    return(mat[,2])
}

## Process and save both beta and coverage matrices for each group
for (grp in unique(metadata$Group)) {
    coverage_file <- file.path(OUTDIR, paste0("coverage_", grp, ".rds"))
    if (!file.exists(coverage_file)) {    
        group_samples <- metadata[Group == grp]
        group_matrix <- do.call(
            cbind,
            lapply(group_samples$file, read_beta)
        )
        colnames(group_matrix) <- group_samples$Sample.name

        coverage_matrix <- do.call(
            cbind,
            lapply(group_samples$file, read_coverage)
        )
        colnames(coverage_matrix) <- group_samples$Sample.name

        ## Save each group as its own RDS file (beta values)
        saveRDS(
            group_matrix,
            file = file.path(OUTDIR, paste0("dataset_", grp, ".rds")),
            compress = "xz"
        )
        ## Save each group as its own RDS file (coverage values)
        saveRDS(
            coverage_matrix,
            file = file.path(OUTDIR, paste0("coverage_", grp, ".rds")),
            compress = "xz"
        )
    }
}

## Save cpg_names and metadata separately
saveRDS(cpg_names, file = file.path(OUTDIR, "cpg_names.rds"), compress = "xz")
saveRDS(metadata, file = file.path(OUTDIR, "metadata.rds"), compress = "xz")

EOF

echo "Done!"
