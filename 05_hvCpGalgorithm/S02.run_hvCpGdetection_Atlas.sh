#!/bin/bash
#$ -N runhvCpGAtlas
#$ -S /bin/bash
#$ -pe smp 10
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=100:00:00
#$ -m abe ## Send an email if abort, begins or end
#$ -M alice.cam.balard@gmail.com
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

source /share/apps/source_files/python/python-3.13.0a6.source

python3 /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/S01.prepare_beta_matrices.py

R --vanilla <<EOF

myNthreads=16 ## specify here

## Load algorithm (30sec)
system.time(source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/hvCpG_algorithm_detection_v3.R"))

## Load data & functions specific to Atlas (2 sec)
system.time(source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/S02.formatAtlasforR.R"))

## Load cpg list (1 min) 29,401,795 CpG names
system.time(cpg_names <- h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets_prepared/Abdominal_Subcut._scaled_matrix.h5", "cpg_names"))

length(cpg_names); head(cpg_names)
##[1] 29401795
##[1] "chr1_10469-10470" "chr1_10471-10472" "chr1_10484-10485" "chr1_10489-10490"
##[5] "chr1_10493-10494" "chr1_10497-10498"

system.time(runAndSave("Atlas", cpgvec = head(cpg_names, 1000),p0=0.80, p1=0.65, NCORES=myNthreads,
		       resultDir="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas/"))

##100, 5G, 1C = 87 sec
##1000, 5G, 1C = 
## 825 CpGs, 5G, 1C in 1403 sec

message("Done!")

q()
EOF

#########################
## end test ## **********

#####################################
## Run the algorithm in each chunk ##
#####################################
#### List all .rds files
##chunk_files <- list.files(output_dir, pattern = "^filtered_chunk_\\d+\\.rds$", full.names = TRUE)
##chunk_files <- sort(chunk_files)  # Ensure correct order
##
### Process each chunk one at a time
##for (f in chunk_files) {
##  message("Reading: ", f)
##  chunk_data <- readRDS(f)
##
##  # `chunk_data` is a list:
##  #   chunk_data[["Adipocytes"]]  -> matrix for Adipocytes
##  #   chunk_data[["Blood_Cauc"]]  -> matrix for Blood Caucasian, etc.
##
##  ## Do your downstream analysis here!
##
##  rm(chunk_data); gc()  # Free RAM
##}
##
##cat("Run finished!")
##
##q()
##EOF
##
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
