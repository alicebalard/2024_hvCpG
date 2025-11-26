#!/bin/bash
#$ -N prepAtlas
#$ -S /bin/bash
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

source /share/apps/source_files/python/python-3.13.0a6.source

######################
## Prepare metadata ##
######################

DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/"

META_IN="${DIR}SupTab1_Loyfer2023.csv"
META_OUT="${DIR}SupTab1_Loyfer2023.with_analysis_group.csv"

python3 /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/prepare_metadata.py --meta "$META_IN" --output "$META_OUT"

###########################
# Run preparation script ##
###########################

python3 /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/prepare_beta_matrices.py \
  --beta_files "${DIR}betaFiles/GSM*.hg38.beta" \
  --cpg_bed "${DIR}wgbs_tools/references/hg38/CpG.bed.gz" \
  --output_folder "${DIR}10X" \
  --meta "$META_OUT" \
  --minCov 10 \
  --min_samples 3 \
  --min_datasets 46 \
  --chunk_size 100000
