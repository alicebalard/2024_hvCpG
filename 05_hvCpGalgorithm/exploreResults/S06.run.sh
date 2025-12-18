#!/bin/bash
#$ -N runannotrGREAT
#$ -S /bin/bash
#$ -pe smp 10
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/exploreResults/S06_analysesExtraAtlas.R
