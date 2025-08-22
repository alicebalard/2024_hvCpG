#!/bin/bash
#$ -N runalgo5_atlas_bychunk
#$ -S /bin/bash
#$ -pe smp 15
#$ -l tmem=6G
#$ -l h_vmem=6G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -l tscratch=50G ## use of scratch
#$ -t 1-231 ## For 23,036,026 CpGs in chunks of 100k
#$ -tc 100

CHUNK_SIZE=100000

echo "**** Job $JOB_NAME started at $(date) ****"

## Create scratch directory for in/out issues
mkdir -p /scratch0/abalard/$JOB_ID

TEMPDIR="/scratch0/abalard/$JOB_ID"

echo "Copy data in scratch..."

orig_dataDir="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/"
cp $orig_dataDir/* $TEMPDIR/. ## Copy in scratch

## Will be sent to different tasks in an array
Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/R02.estimAlpha_DerakhvCpGsControls_inAtlas.R $TEMPDIR $SGE_TASK_ID $CHUNK_SIZE

## Rm scratch data
function finish {
    rm -rf $TEMPDIR
}

trap finish EXIT ERR INT TERM

echo "**** Job $JOB_NAME finished at $(date) ****"










