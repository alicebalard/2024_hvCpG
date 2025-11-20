#!/usr/bin/env Rscript

## Load hyperVarMeth (installed in my CS HPC home R.4.2 directory)
library(hyperVarMeth)
#####################

args <- commandArgs(trailingOnly = TRUE)

data_dir <- args[1] ## where the data is
task_id <- as.integer(args[2]) ## which task on the array
chunk_size <- as.integer(args[3]) ## what is the size of the chunk
batch_size <- as.integer(args[4]) ## how many CpGs are loaded at once
res_dir <- as.integer(args[5]) ## where the results should be stored

## Run on CpGs which are covered in all 46 cells 
cpg_46 <- read.table(file.path(data_dir, "selected_cpgs_min3_in46_datasets.txt"))$V1

## Batch
start_idx <- (task_id - 1) * chunk_size + 1
end_idx <- min(task_id * chunk_size, length(cpg_46))

if (start_idx > length(cpg_46)) {
  stop("Start index beyond end of CpG list.")
}

subset_cpgs <- cpg_46[start_idx:end_idx]

# Output directory
result_dir <- sprintf(paste0(res_dir, "Atlas_batch%03d"), task_id)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

message(paste0("If new, results will be saved in dir: ", result_dir))

## Run
myNthreads <- as.numeric(Sys.getenv("NSLOTS", unset = "1"))  # Use all cores

message("Run algo:")

system.time(hyperVarMeth::runAndSave_fast(
    analysis = "Atlas10X",
    cpg_names_vec = subset_cpgs,
    dataDir = data_dir,
    resultDir = result_dir,
    NCORES = myNthreads,
    p0 = 0.80,
    p1 = 0.65,
    batch_size = batch_size)
)
