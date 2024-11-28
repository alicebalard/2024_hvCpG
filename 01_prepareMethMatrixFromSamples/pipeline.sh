#$ -N SRA sample
#$ -S /bin/bash
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -t 1-$(ls data/WGBS_human/00RawFastq/*fastq | wc -l)
#$ -l h_rt=240:00:00
#$ -o /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.out
#$ -e /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs/$JOB_NAME_$JOB_ID.err
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/

echo "Running on host: $(hostname)"
echo "Start time: $(date)"

## Run on all the fastq files available. Run by batches to avoid storing insane amounts of big s, then remove fastq

## Get the list of FASTQ files
FILES=(data/WGBS_human/00RawFastq/*fastq)

## Get the current file based on the task ID
CURRENT_FILE=${FILES[$SGE_TASK_ID-1]}

## -i "$CURRENT_FILE" -o "output_${SGE_TASK_ID}.txt"

###########################################
## Trim galore with automatic parameters - at least 20 M reads with length above 100 bp.


##################################################
# Run fastqc for quality check of trimmed reads ##
## /share/apps/genomics/FastQC-0.11.9/fastqc -o 01FastQC $XXX
## NB 1h/sample with 4G, add cores +++ and run on gz


#########################################################################
## Bisulfite conversion efficiency measure → store in matrix (BCREval) ##
#python code/2024_hvCpG/01_prepareMethMatrixFromSamples/BCReval.py -n 8 -i $INPUT -o $INPUT.BCREval.out
#cat $INPUT.BCREval.out | awk '{sum = $(NF-2) + $(NF-1) + $NF; mean = sum / 3; result = 1 - mean; printf "%.4f\n", result}' ## save somewhere this Conversion ratio (CR)
#rm $INPUT.BCREval.out

## Bismark: Align to GRCh38 -> Output BAM

## Bismark: Deduplication

## Bismark: base-level methylation calling (--pairedend and --no_overlap parameters set)

## Remove fastq.gz file

echo "End time: $(date)"
