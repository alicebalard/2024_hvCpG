#!/bin/bash
#$ -N runTissuealgo_atlas
#$ -S /bin/bash
#$ -pe smp 1
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-100
#$ -tc 30

## arg 1: where the starting data is
DATA_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/output_atlas_general/"

## arg 2: $SGE_TASK_ID

## arg 3: How big are chunks sent to arrays?
CHUNK_SIZE=250000

## arg 4: How many CpGs are loaded at the same time?
BATCH_SIZE=10000

## arg 5: where the results should be stored
RES_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/resultsDir_gitIgnored/Atlas/Atlas10X_tissueAnalysis/"

echo "**** Job $JOB_NAME.$SGE_TASK_ID started at $(date) ****"

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/02_runAlgo/Atlas/tissueSpec/S01.1_runTissueSpecificAlgo.R "$DATA_DIR" "$SGE_TASK_ID" "$CHUNK_SIZE" "$BATCH_SIZE" "$RES_DIR"

cp /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/output_atlas_general/all_medsd_lambda.tsv $RES_DIR/.

echo "**** Job $JOB_NAME.$SGE_TASK_ID finished at $(date) ****"
