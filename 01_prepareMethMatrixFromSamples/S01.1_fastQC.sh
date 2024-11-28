#$ -N fasterqc
#$ -S /bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=2:00:00
#$ -o /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.out
#$ -e /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.err
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/

#Run fastqc for quality check of raw reads
/share/apps/genomics/FastQC-0.11.9/fastqc -o 01FastQC 00RawFastq/*fastq 


## NB 1h/sample, add cores +++ and run on gz
