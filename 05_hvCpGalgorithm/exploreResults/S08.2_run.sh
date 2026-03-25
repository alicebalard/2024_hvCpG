#!/bin/bash
#$ -N runrGREAT
#$ -S /bin/bash
#$ -pe smp 10
#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/monitor_resources.sh &

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/exploreResults/S08_rGREAT_pchuckle.R
