#$ -N fasterq_dump_job
#$ -S /bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -t 1-10
#$ -l h_rt=24:00:00
#$ -o /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.out
#$ -e /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.err
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/00RawFastq

echo "Running on host: $(hostname)"
echo "Start time: $(date)"

echo "Running with $NSLOTS threads"

SAMPLELIST=listTest.txt
## SAMPLELIST=sraAccList1.txt

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
