#!/bin/bash
#$ -N prepAtlas
#$ -S /bin/bash
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

source /share/apps/source_files/python/python-3.13.0a6.source

python3 /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/S01.1.prepare_beta_matrices.py

