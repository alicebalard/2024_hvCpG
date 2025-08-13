#!/bin/bash
#$ -N runalgo5_atlas_hvcpgControls
#$ -S /bin/bash
#$ -pe smp 5
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -l tscratch=50G ## use of scratch

echo "**** Job $JOB_NAME started at $(date) ****"

## Create scratch directory for in/out issues
mkdir -p /scratch0/abalard/$JOB_ID

TEMPDIR="/scratch0/abalard/$JOB_ID"

echo "Copy data in scratch..."

orig_dataDir="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/"
cp $orig_dataDir/* $TEMPDIR/. ## Copy in scratch

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/R02.estimAlpha_DerakhvCpGsControls_inAtlas.R $TEMPDIR $NSLOTS

## Rm scratch data
function finish {
    rm -rf $TEMPDIR
}

trap finish EXIT ERR INT TERM

echo "**** Job $JOB_NAME finished at $(date) ****"








