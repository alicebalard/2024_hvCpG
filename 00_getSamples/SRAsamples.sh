#$ -N SRA sample
#$ -S /bin/bash
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -t 1-2
#$ -l h_rt=240:00:00
#$ -o /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.out
#$ -e /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.err
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/00RawFastq

echo "Running on host: $(hostname)"
echo "Start time: $(date)"

## DL files from a SRA list

SAMPLELIST=listTest.txt
## SAMPLELIST=sraAccList1.txt

## Dataset1
/share/apps/genomics/sratoolkit.3.0.2/bin/prefetch --max-size 30GB --option-file $SAMPLELIST

# Get the sample for this task ID
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $SAMPLELIST)

if [ -n "$SAMPLE" ]; then
    echo "Processing $SAMPLE..."
    # Run fasterq-dump for the sample 
    /share/apps/genomics/sratoolkit.3.0.2/bin/fasterq-dump "$SAMPLE" --threads $NSLOTS 
    
    # Check if fasterq-dump was successful
    if [ $? -eq 0 ]; then
        echo "$SAMPLE downloaded and compressed successfully."
    else
        echo "Error downloading or compressing $SAMPLE."
    fi
    ## Gzip both forward and reverse samples
    gzip ${SAMPLE}_1.fastq ; gzip ${SAMPLE}_2.fastq
else
    echo "No sample found for task ID $SGE_TASK_ID"
fi

echo "End time: $(date)"
