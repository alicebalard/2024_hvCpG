#!/bin/bash

## Prepare the data for each analysis

## Excluding SNPs + SD-ASM (core analyses)
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset atlas_general
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 02_rmMultSamples
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 12_endo
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 12_2_endo6gp
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 13_meso
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 13_2_meso6gp
qsub -v EXCLUDE_MODE=both run_pipeline_atlas.sh --preset 14_ecto

## Run with only SNPs excluded, to test downstream the enrichment in SNPs in top 1% CpGs
qsub -v EXCLUDE_MODE=snp run_pipeline_atlas.sh --preset atlas_general
qsub -v EXCLUDE_MODE=snp run_pipeline_atlas.sh --preset 12_endo
qsub -v EXCLUDE_MODE=snp run_pipeline_atlas.sh --preset 13_meso
qsub -v EXCLUDE_MODE=snp run_pipeline_atlas.sh --preset 14_ecto

## Run with no site excluded, to test downstream the enrichment in SNPs and SD-ASM in top 1% CpGs
qsub -v EXCLUDE_MODE=none run_pipeline_atlas.sh --preset atlas_general
qsub -v EXCLUDE_MODE=none run_pipeline_atlas.sh --preset 12_endo
qsub -v EXCLUDE_MODE=none run_pipeline_atlas.sh --preset 13_meso
qsub -v EXCLUDE_MODE=none run_pipeline_atlas.sh --preset 14_ecto
