#!/bin/bash
#$ -N runalgo6_atlas_100k.5T.10G
#$ -S /bin/bash
#$ -pe smp 5
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-240
<<<<<<< HEAD
=======
#$ -tc 30
>>>>>>> 499f91a94f459808c6acfc15b74437747cb6f4bc

CHUNK_SIZE=100000 ## How big are chunks sent to arrays?
BATCH_SIZE=10000 ## How many CpGs are loaded at the same time?

echo "**** Job $JOB_NAME.$SGE_TASK_ID started at $(date) ****"

DATA_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/"

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/S02.1_runalgov6_atlas_cscluster.R "$DATA_DIR" "$SGE_TASK_ID" "$CHUNK_SIZE" "$BATCH_SIZE"

echo "**** Job $JOB_NAME.$SGE_TASK_ID finished at $(date) ****"
