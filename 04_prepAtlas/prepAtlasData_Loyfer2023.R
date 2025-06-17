## 4 June 2025
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

# Optimized version of prepAtlasData_Loyfer2023.R (no dataframes, too heavy)
DIR <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/"

## Read the CpG names (store once, use references)
cpg_names_path <- paste0(DIR, "hg38CpGpos_Loyfer2023.txt")
if(!exists("cpg_names")) {  # Only load if not already in memory
  cpg_names <- readLines(cpg_names_path)
}

## Read in metadata
metadata <- read.csv(paste0(DIR, "SupTab1_Loyfer2023.csv"))

## Filter groups with ≥3 samples
valid_groups <- names(which(table(metadata$Group) >= 3))
metadata <- metadata[metadata$Group %in% valid_groups, ]

## File handling
files <- list.files(paste0(DIR, "betaFiles"), full.names = TRUE)
sample_names <- as.character(sub(".*-", "", metadata$Sample.name))

# Create sample-file mapping
sample_file_map <- setNames(
  lapply(sample_names, function(sname) {
    f <- grep(sname, files, value = TRUE, fixed = TRUE)
    if(length(f) == 1) f else NA
  }),
  sample_names
)
sample_file_map <- sample_file_map[!is.na(sample_file_map)]

## Process files using matrix storage
process_file <- function(file) {
  file_size <- file.info(file)$size
  content <- readBin(file, "integer", n = file_size, size = 1, signed = FALSE)
  matrix(content, ncol = 2, byrow = TRUE)
}

## Process files in batches per group
final_output <- list()
sample_groups <- metadata$Group[match(names(sample_file_map), sample_names)]

df <- data.frame(sample = names(sample_file_map), group = sample_groups)

for(current_group in unique(df$group)) {
    group_samples <- df$sample[df$group %in% current_group]

    ## Process samples in current group
    group_data <- lapply(sample_file_map[group_samples], function(f) {
        mat = process_file(f)
        beta = mat[,1] / mat[,2]  # Calculate beta values
        beta[beta == Inf] = NA
        return(beta)
    })

    ## Convert to matrix (more memory efficient than data.frame)
    group_matrix <- do.call(cbind, group_data)
    dimnames(group_matrix) <- list(NULL, group_samples)  # No row names stored
    final_output[[current_group]] <- group_matrix
}

# Save optimized structure
saveRDS(list(
  data = final_output,
  cpg_names = cpg_names,
  metadata = metadata
), file = paste0(DIR, "listDatasetsLoyfer2023.RDS"))
