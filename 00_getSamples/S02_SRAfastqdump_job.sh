#!/bin/bash
#$ -N fastqdump
#$ -S /bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -pe smp 8
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-48   # <-- Set this to cat temp.samplelist | wc -l

################################################
## DL files from a SRA list: previous step (S01)
DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/00RawFastq/dataset3"
cd $DIR
SAMPLELIST="$DIR/temp.samplelist" ## outputed in S01

SAMPLE=$(sed -n "${SGE_TASK_ID}p" "$SAMPLELIST")

echo "Running on host: $(hostname)"
echo "Start time: $(date)"
echo "Processing sample: $SAMPLE"

# Check if FASTQ file compressed exists
if [ ! -f "${SAMPLE}_1.fastq.gz" ]; then
    echo "Run fasterq-dump for the sample: $SAMPLE"
    /share/apps/genomics/sratoolkit.3.0.2/bin/fasterq-dump -e 8 "$SAMPLE"

    echo "Gzip both forward and reverse samples:"
    gzip ${SAMPLE}_1.fastq ; gzip ${SAMPLE}_2.fastq

    echo "Remove SRA folder"
    rm -rf $SAMPLE
else
    echo "Skipping $SAMPLE: Files already exist."
fi
