#!/usr/bin/env Rscript

## Load hyperVarMeth (Installed in ing-p5 in ‘/home/alice/R/x86_64-pc-linux-gnu-library/4.4’)
library(hyperVarMeth)
#####################

args <- commandArgs(trailingOnly = TRUE)

data_dir <- args[1] ## where the data is
batch_size <- as.integer(args[2]) ## how many CpGs are loaded at once
res_dir <- args[3] ## where the results should be stored
nthreads <- as.integer(args[4]) ## number of threads

## Run on CpGs which are covered
cpg_names <- rhdf5::h5read(file.path(data_dir, "all_matrix_noscale.h5"), "cpg_names")

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
message(paste0("If new, results will be saved in dir: ", res_dir))

message("Run algo:")
system.time(hyperVarMeth::runAndSave_tissueAnalysis(
  analysis = "Array_tissue0.5",
  cpg_names_vec = cpg_names,
  resultDir = res_dir,
  NCORES = nthreads,
  p1 = 0.5,
  batch_size = batch_size,
  dataDir = data_dir)
)
