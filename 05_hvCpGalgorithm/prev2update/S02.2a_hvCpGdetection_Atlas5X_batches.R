#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
chunk_size <- as.integer(args[2])

# Load required objects
cpg_names_all <- rhdf5::h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/5X/all_scaled_matrix.h5", "cpg_names")

start_idx <- (task_id - 1) * chunk_size + 1
end_idx <- min(task_id * chunk_size, length(cpg_names_all))

if (start_idx > length(cpg_names_all)) {
  stop("Start index beyond end of CpG list.")
}

subset_cpgs_pos <- start_idx:end_idx

# Output directory
result_dir <- sprintf("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas5X/Atlas_batch%03d", task_id)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

## Load algorithm
system.time(source("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/hvCpG_algorithm_detection_v4scan.R"))                                                            

## Run
myNthreads <- as.numeric(Sys.getenv("NSLOTS", unset = "1"))  # Use all cores

runAndSave(
  analysis = "Atlas5X",
  cpgPos_vec = subset_cpgs_pos,
  resultDir = result_dir,
  NCORES = myNthreads,
  p0 = 0.80,
  p1 = 0.65
)
