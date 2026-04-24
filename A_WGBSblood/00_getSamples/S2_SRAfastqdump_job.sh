#!/bin/bash
#$ -N fastqdump
#$ -S /bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -pe smp 8 # Request N cores per task 
#$ -l h_rt=200:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

DIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/00RawFastq/dataset2"
cd $DIR

echo "Running on host: $(hostname)"
echo "Start time: $(date)"

###############################################
## DL files from a SRA list: previous step (S1)

SAMPLELIST="$DIR/temp.samplelist"
echo "samples list:"; cat $SAMPLELIST

while IFS= read -r SAMPLE; do
    # Check if FASTQ file compressed exists
    if [ ! -f "${SAMPLE}_1.fastq.gz" ]; then
	echo "Run fasterq-dump for the sample:"; cat $SAMPLE
	/share/apps/genomics/sratoolkit.3.0.2/bin/fasterq-dump -e 8 "$SAMPLE"

	echo "Gzip both forward and reverse samples:"
	gzip ${SAMPLE}_1.fastq ; gzip ${SAMPLE}_2.fastq

	echo "Remove SRA folder"
	rm -rf $SAMPLE
    else
   	echo "Skipping $SAMPLE: Files already exist."
    fi
done < "$SAMPLELIST"
