#!/bin/bash
#$ -N runalgo5_atlas_100k.10T.5G
#$ -S /bin/bash
#$ -pe smp 10
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-240

CHUNK_SIZE=100000 ## How big are chunks sent to arrays?
BATCH_SIZE=10000 ## How many CpGs are loaded at the same time?

echo "**** Job $JOB_NAME.$SGE_TASK_ID started at $(date) ****"

DATA_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/"

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/S02.1_runalgov5_atlas_cscluster.R "$DATA_DIR" "$SGE_TASK_ID" "$CHUNK_SIZE" "$BATCH_SIZE"

echo "**** Job $JOB_NAME.$SGE_TASK_ID finished at $(date) ****"











