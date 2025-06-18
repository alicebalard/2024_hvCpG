#!/bin/bash
#$ -N prepAtlasData
#$ -S /bin/bash
#$ -l tmem=50G
#$ -l h_vmem=50G
#$ -l h_rt=10:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/prepAtlasData_Loyfer2023.R
