#!/bin/bash
#$ -N recreateAtlasDataChunks
#$ -S /bin/bash
#$ -l tmem=50G
#$ -l h_vmem=50G
#$ -l h_rt=48:00:00
#$ -pe smp 8 # Request 8 cores
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

echo "Starting chunked R inline script with half-group filter..."

Rscript --vanilla - << 'EOF'

library(future)
library(future.apply)

## Allow larger globals
options(future.globals.maxSize = 2 * 1024^3)  # 2GB

plan(multisession, workers = 8)

dataset_dir <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/datasets"
output_dir <- "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/filtered_chunks/"
chunk_size <- 500000
coverageMin <- 20
NIND <- 3

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

beta_files <- list.files(dataset_dir, pattern = "^dataset_.*\\.rds$", full.names = TRUE)
coverage_files <- list.files(dataset_dir, pattern = "^coverage_.*\\.rds$", full.names = TRUE)
cpgnames <- readRDS(file.path(dataset_dir, "cpg_names.rds"))

get_group <- function(x) sub("^(dataset|coverage)_(.*)\\.rds$", "\\2", basename(x))
groups <- intersect(get_group(beta_files), get_group(coverage_files))
beta_files_named <- setNames(beta_files, get_group(beta_files))
coverage_files_named <- setNames(coverage_files, get_group(coverage_files))

total_cpgs <- length(cpgnames)
n_chunks <- ceiling(total_cpgs / chunk_size)

# Pre-index to avoid sending cpgnames
rownames_list <- lapply(groups, function(grp) {
  list(
    beta = beta_files_named[[grp]],
    coverage = coverage_files_named[[grp]]
  )
})
names(rownames_list) <- groups

for (chunk_idx in seq_len(n_chunks)) {
  idx_start <- (chunk_idx - 1) * chunk_size + 1
  idx_end <- min(chunk_idx * chunk_size, total_cpgs)
  rows_chunk <- idx_start:idx_end

  out_file <- file.path(output_dir, sprintf("filtered_chunk_%04d.rds", chunk_idx))
  if (file.exists(out_file)) {
    message(sprintf("⏭️  Chunk %d/%d already exists. Skipping.", chunk_idx, n_chunks))
    next
  }

  message(sprintf("Processing chunk %d/%d: rows %d-%d",
                  chunk_idx, n_chunks, idx_start, idx_end))

  chunk_pass_matrix <- future_lapply(groups, function(grp) {
    beta_mat <- readRDS(rownames_list[[grp]]$beta)
    coverage_mat <- readRDS(rownames_list[[grp]]$coverage)

    rownames(beta_mat) <- cpgnames
    rownames(coverage_mat) <- cpgnames

    beta_chunk <- beta_mat[rows_chunk, , drop = FALSE]
    coverage_chunk <- coverage_mat[rows_chunk, , drop = FALSE]

    beta_chunk[coverage_chunk < coverageMin] <- NA
    keep_rows <- rowSums(!is.na(beta_chunk)) >= NIND

    rm(beta_mat, coverage_mat, beta_chunk, coverage_chunk)
    gc()
    keep_rows
  })

  pass_counts <- Reduce(`+`, chunk_pass_matrix)
  final_keep_rows <- pass_counts >= length(groups) / 2
  rows_final <- rows_chunk[final_keep_rows]

  message(sprintf("➡️ %d CpGs pass in ≥ half groups for this chunk.",
                  length(rows_final)))

  if (length(rows_final) == 0) {
    message("⚠️ Chunk ", chunk_idx, " has no rows passing filter. Skipping save.")
    next
  }

  chunk_group_list <- future_lapply(groups, function(grp) {
    beta_mat <- readRDS(rownames_list[[grp]]$beta)
    coverage_mat <- readRDS(rownames_list[[grp]]$coverage)

    rownames(beta_mat) <- cpgnames
    rownames(coverage_mat) <- cpgnames

    beta_chunk <- beta_mat[rows_final, , drop = FALSE]
    coverage_chunk <- coverage_mat[rows_final, , drop = FALSE]

    beta_chunk[coverage_chunk < coverageMin] <- NA
    keep_rows_group <- rowSums(!is.na(beta_chunk)) >= NIND

    beta_chunk_filtered <- beta_chunk[keep_rows_group, , drop = FALSE]

    rm(beta_mat, coverage_mat, beta_chunk, coverage_chunk)
    gc()
    beta_chunk_filtered
  })

  names(chunk_group_list) <- groups

  out_file <- file.path(output_dir, sprintf("filtered_chunk_%04d.rds", chunk_idx))
  saveRDS(chunk_group_list, out_file, compress = "xz")
  message("✅ Saved: ", out_file)
}

message("✅✅✅ All chunked processing done.")

EOF

echo "Done!"

#### sapply(Atlas_data, ncol)
####       Adipocytes        Bladder-Ep           Blood-B      Blood-Granul 
####                3                 5                 5                 3 
#### Blood-Mono+Macro          Blood-NK           Blood-T   Breast-Basal-Ep 
####               11                 3                22                 4 
####Breast-Luminal-Ep          Colon-Ep       Endothelium        Eryth-prog 
####                3                 8                19                 3 
####     Fallopian-Ep        Gastric-Ep      Head-Neck-Ep      Heart-Cardio 
####                3                11                13                 4 
####      Heart-Fibro         Kidney-Ep         Liver-Hep     Lung-Ep-Alveo 
####                4                 8                 6                 4 
####     Lung-Ep-Bron            Neuron         Oligodend    Ovary+Endom-Ep 
####                3                10                 4                 4 
####  Pancreas-Acinar    Pancreas-Alpha     Pancreas-Beta    Pancreas-Delta 
####                4                 3                 3                 3 
####    Pancreas-Duct       Prostate-Ep      Small-Int-Ep       Smooth-Musc 
####                4                 4                 5                 5 
####       Thyroid-Ep 
####                3 
##
