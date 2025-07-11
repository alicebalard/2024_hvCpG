NSLOTS=10

echo "Number of slots allocated: $NSLOTS"

R --vanilla <<EOF

myNthreads <- as.integer(Sys.getenv("NSLOTS"))
cat("Number of threads inside R:", myNthreads, "\n")

## Atlas data was preprocessed in 04/S02 and cut in chunks of 500k CpG
## Filter: CpGs with at least 20 coverage in at least 3 individuals and at least half of the datasets
input_dir <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/filtered_chunks/"

## Source the function:
source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/hvCpG_algorithm_detection_v3.R")

#########################
## test ## **************
chunk_idx <- 1

## Build filename (must match your sprintf format!)
file_path <- file.path(input_dir, sprintf("filtered_chunk_%04d.rds", chunk_idx))

## Load it
chunk_data <- readRDS(file_path)

cpgnames <- unique(unlist(sapply(chunk_data, row.names)))
cpgnames <- cpgnames[order(cpgnames)]

#### test subset 50k: 8 cores, 10Gb/core 2h30 for 50k pos
#### test subset 50k: 12 cores, 10Gb/core: 30min (must have bugged) ou 2h20?
#### test all 500k: 12 cores, 10Gb/core

### Take first 50 CpG names
cpg_subset <- cpgnames[1:50]
chunk_data_subset50 <- lapply(
  chunk_data,
 function(mat) {
    rows_to_keep <- intersect(rownames(mat), cpg_subset)
    mat[rows_to_keep, , drop = FALSE]
  }
)

system.time(runAndSave(my_list_mat = chunk_data_subset50, cpgvec = cpg_subset,
           optimMeth="Nelder-Mead", NCORES=myNthreads, p0=0.95, p1=0.65, 
resultDir="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas/"))

#system.time(runAndSave(my_list_mat = chunk_data, cpgvec = cpgnames,                                                                                                        
#           optimMeth="Nelder-Mead", NCORES=myNthreads, p0=0.95, p1=0.65,                                                                                                               
#resultDir="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas/"))     

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
