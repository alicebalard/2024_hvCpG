#!/bin/bash
#$ -N bssnper2 
#$ -S /bin/bash
#$ -pe smp 1
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=10:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs

BISMARK_OUTDIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/02Bismark"
GENOME_DIR="$BISMARK_OUTDIR/GRCh38"
DEDUPBAM="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/02Bismark/test.bam"
## From before in pipeline ##
#############################

echo "**** Start of sorting test bam file : $(date) ****" 
## The bam file needs to be sorted before
/share/apps/genomics/samtools-1.9/bin/samtools sort -o $DEDUPBAM.sorted.bam $DEDUPBAM

echo "**** End of sorting test bam file : $(date) ****" 

echo "**** Start of bssnper2 on test bam file : $(date) ****" 

BSSNPER2_OUTDIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/03Bssnper2"

# add bssnper2 to my path
export PATH=/home/abalard/bssnper2/:$PATH

bssnper2 $DEDUPBAM.sorted.bam --ref $GENOME_DIR/GCF_000001405.40_GRCh38.p14_genomic.fa.gz --vcf $BBSSNPER2_OUTDIR/test.vcf

echo "**** End of bssnper2 on test bam file : $(date) ****" 

echo "**** Start of sorted 088 bam file : $(date) ****" 
## Test on bigger if it works

DEDUPBAM="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/02Bismark/SRR28532088_pe.deduplicated.bam"

## The bam file needs to be sorted before
/share/apps/genomics/samtools-1.9/bin/samtools sort -o $DEDUPBAM.sorted.bam $DEDUPBAM
echo "**** End of sorted 088 bam file : $(date) ****" 

echo "**** Start of bssbper2 on sorted 088 bam file : $(date) ****" 
bssnper2 $DEDUPBAM.sorted.bam --ref $GENOME_DIR/GCF_000001405.40_GRCh38.p14_genomic.fa.gz --vcf $BBSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf

echo "**** End of bssbper2 on sorted 088 bam file : $(date) ****" 

