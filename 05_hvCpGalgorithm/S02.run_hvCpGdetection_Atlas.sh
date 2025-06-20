#!/bin/bash
#$ -N runhvCpGAtlas
#$ -S /bin/bash
#$ -l tmem=100G
#$ -l h_vmem=100G
#$ -l h_rt=240:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

R --vanilla <<EOF

setwd("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/")

###############
## Data load ##
###############

## Atlas data (Hosted on UCL cs server; NB need high ram (100Gb) to read!)
## Outputs cpgnames and my_list_mat:
source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/S02_recreateAtlas.R")
my_list_mat <- recreateAtlas() ## 40 min with 100G needed
rm(coverage_files, beta_files)

######################
## Launch algorithm ##
######################

cat(paste0("Prepare ", which, " dataset and source functions for optimisation:\n"))
source("hvCpG_algorithm_detection.R")
cat("Functions sourced.\n")

## Run algorithm on Atlas data:
result <- getAllOptimAlpha_parallel(cpgvec = cpgnames, optimMeth="Nelder-Mead", NCORES=8, p0=0.95)
save(result, test3000CpGsvec, MariasCpGs, file = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/results_NM_Atlas_0.95.RDA")

q()
EOF

##################################################
## Are Maria's hvCpG covered by the Atlas data? ##
##################################################
#### Liftover hg38-hg19
### 1. Required packages
##if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
##BiocManager::install(c("rtracklayer", "GenomicRanges", "minfi"))
##
##library(GenomicRanges)
##library(rtracklayer)
##library(minfi)  # for Illumina manifest
##
### 2. Parse to GRanges
##parse_coords_to_granges <- function(coord_strings) {
##  parts <- regmatches(coord_strings, regexec("^(chr[^_]+)_([0-9]+)-([0-9]+)$", coord_strings))
##  parsed <- do.call(rbind, lapply(parts, function(x) x[2:4]))
##  parsed <- as.data.frame(parsed, stringsAsFactors = FALSE)
##  parsed$start <- as.integer(parsed$start)
##  parsed$end <- as.integer(parsed$end)
##  
##  GRanges(seqnames = parsed[,1],
##          ranges = IRanges(start = parsed[,2], end = parsed[,3]),
##          strand = "*")
##}
##gr_hg38 <- parse_coords_to_granges(cpgnames)
##
### 3. LiftOver to hg19
##chain_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz"
##chain_file <- "hg38ToHg19.over.chain"
##if (!file.exists(chain_file)) {
##  download.file(chain_url, destfile = paste0(chain_file, ".gz"))
##  R.utils::gunzip(paste0(chain_file, ".gz"))
##}
##chain <- import.chain(chain_file)
##gr_hg19 <- unlist(liftOver(gr_hg38, chain))
##
### 4. Load Illumina 450K annotation
##data(Locations)  # from `minfiData`/`minfi` package
### If not loaded, try:
### annotation <- getAnnotation(getPlatform("IlluminaHumanMethylation450k"))
##annotation <- getAnnotation(IlluminaHumanMethylation450k)
##
### 5. Match hg19 positions to Illumina CpGs
### Match by chromosome and start position
##gr_anno <- GRanges(seqnames = annotation$chr,
##                   ranges = IRanges(start = annotation$pos, width = 1),
##                   names = rownames(annotation))
##
### Use `findOverlaps`
##hits <- findOverlaps(gr_hg19, gr_anno, type = "equal")
##
### 6. Extract matched CpG IDs
##matched_cpg_ids <- names(gr_anno)[subjectHits(hits)]
##matched_cpg_ids
##
