#!/usr/bin/env Rscript

## Load hyperVarMeth (installed in my CS HPC home R.4.2 directory)
library(hyperVarMeth)
#####################

args <- commandArgs(trailingOnly = TRUE)

data_dir <- args[1] ## where the data is
task_id <- as.integer(args[2]) ## which task on the array
chunk_size <- as.integer(args[3]) ## what is the size of the chunk
batch_size <- as.integer(args[4]) ## how many CpGs are loaded at once
res_dir <- args[5] ## where the results should be stored

## Run on CpGs which are covered
cpg_names <- rhdf5::h5read(file.path(data_dir, "all_matrix_noscale.h5"), "cpg_names")

## Batch
start_idx <- (task_id - 1) * chunk_size + 1
end_idx <- min(task_id * chunk_size, length(cpg_names))

if (start_idx > length(cpg_names)) {
  stop("Start index beyond end of CpG list.")
}

subset_cpgs <- cpg_names[start_idx:end_idx]

# Output directory
result_dir <- sprintf(paste0(res_dir, "Atlas_batch%03d"), task_id)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

message(paste0("If new, results will be saved in dir: ", result_dir))

## Run
myNthreads <- as.numeric(Sys.getenv("NSLOTS", unset = "1"))  # Use all cores

message("Run algo:")

system.time(hyperVarMeth::runAndSave_tissueAnalysis(
  analysis = "Atlas10X_tissue0.5",
  cpg_names_vec = subset_cpgs,
  resultDir = result_dir,
  NCORES = myNthreads,
  p1 = 0.5,
  batch_size = batch_size,
  dataDir = data_dir)
)
