#!/bin/bash
#$ -N runhvCpGAtlas_batches
#$ -S /bin/bash
#$ -pe smp 5
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-295 ## For 29,401,795 CpGs in batches of 100k
#$ -tc 100

CHUNK_SIZE=100000

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/S02.a_hvCpGdetection_Atlas_batches.R $SGE_TASK_ID $CHUNK_SIZE
