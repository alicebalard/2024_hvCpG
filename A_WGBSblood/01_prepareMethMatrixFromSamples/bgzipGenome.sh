#!/bin/bash
#$ -N bgzip
#$ -S /bin/bash
#$ -pe smp 1
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=01:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/02Bismark/GRCh38/

/share/apps/genomics/htslib-1.9/bin/bgzip /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/02Bismark/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fa
