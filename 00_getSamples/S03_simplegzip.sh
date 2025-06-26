#!/bin/bash
#$ -N simpleGzip
#$ -S /bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -pe smp 8
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

################################################
## DL files from a SRA list: previous step (S1)
DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/00RawFastq/dataset2"
cd $DIR

gzip *fastq ## skip already gzipped ones by default

