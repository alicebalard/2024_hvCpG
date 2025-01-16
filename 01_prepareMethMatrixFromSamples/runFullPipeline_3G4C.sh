#!/bin/bash
#$ -N runPipeline
#$ -S /bin/bash
#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -t 1-2
#$ -tc 2
#$ -pe smp 4  # Request N cores per task 
#$ -l h_rt=40:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample 
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you.

## Angie from cs cluster:
# There are some working compute nodes with 40 threads available on the cluster, so your job should be able to start.
# Also, making sure that your memory requirement and runtime limits are within a reasonable range can further reduce waiting times.
# Please note that in parallel environments (like SMP), the total memory requested would be calculated as:
# tmem * number of thread
# In your case, this means you are requesting 12 GB * 40 threads = 480 GB of memory.

echo "Number of cores requested per task: $NSLOTS (must be >4 for trimming)"

echo "Test: 3G 4cores 40h"
## Some tips from CS (<3) team to have a job queuing less long:
# you can try reducing your total requested memory (smp value x tmem setting) to 15G (e.g. for smp=4, set tmem and h_vmem to 3G). This would allow the jobs to run on our lowest specd. 4-core nodes which have 15G ram (thereare over 500 of them in the cluster).
# If your pipeline is effectively using smp and you wish to try another setting like 10 to see if they run faster, it's still worth limiting the memory to 3G which will give you access to > 140 compute hosts.

## We make a temporary folder for all intermediate files that we delete at the end.
mkdir -p /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/TEMP
TEMP_OUTDIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/TEMP"

################################
# Create the files to loop over if it does not exist:
if [ ! -e "$TEMP_OUTDIR/bysample_list_of_files.tmp" ]; then
    ## NB: input files must be named filename_1.fastq.gz and filename_2.fastq.gz
    ls -1 /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/00RawFastq/*_1.fastq.gz > $TEMP_OUTDIR/list_of_files.tmp
    sed 's/_1.fastq.gz$//' $TEMP_OUTDIR/list_of_files.tmp > $TEMP_OUTDIR/bysample_list_of_files.tmp
fi

## Select the correct line of list of files at each iteration
INPUT=$(sed -n "${SGE_TASK_ID}p" $TEMP_OUTDIR/bysample_list_of_files.tmp)

INPUT_1=${INPUT}_1.fastq.gz
INPUT_2=${INPUT}_2.fastq.gz

## The script runs on all the fastq files available, by batches to avoid storing insane amounts of big files,
## then we'll remove fastq

###########################
# Log the start of the job. 
echo "**** Job $JOB_NAME.$JOB_ID started at $(date) ****"
echo "Task ID: $SGE_TASK_ID"
echo "Selected input file: $INPUT"
echo "We work on sample ${INPUT##*/}"

## We'll feed information to this table
DATA_TABLE_OUT="/SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/data_out_01/infoAllsamples_Pipeline.csv"

# Check if the file exists
if [ ! -e "$DATA_TABLE_OUT" ]; then
    # Create the file
    echo "Sample,Nbr_reads,BS_conversion_rate" > $DATA_TABLE_OUT
    echo "File created: $DATA_TABLE_OUT"
else
    echo "File already exists: $DATA_TABLE_OUT"
fi

################################################################################
echo "**** Start of step 1: Trim galore with automatic parameters $(date) ****" 
## ~6h/sample with 1 thread, 5Gb
## ~5h/sample for 4C 10G
## ~4h/sample for 4C 5G

# Output directory
OUTPUT_DIR_01="$TEMP_OUTDIR/01Trimmed_data"

# Create output directory if it does not exist
mkdir -p $OUTPUT_DIR_01

## Need to be BEFORE the script
TRIMMED_1="$OUTPUT_DIR_01/${INPUT##*/}_1_val_1.fq.gz"
TRIMMED_2="$OUTPUT_DIR_01/${INPUT##*/}_2_val_2.fq.gz"
FASTQC_1="$OUTPUT_DIR_01/${INPUT##*/}_1_val_1_fastqc.html"
FASTQC_2="$OUTPUT_DIR_01/${INPUT##*/}_2_val_2_fastqc.html"

if [ ! -f "$TRIMMED_1" ] || [ ! -f "$TRIMMED_2" ] || [ ! -f "$FASTQC_1" ] || [ ! -f "$FASTQC_2" ]; then
    # source python 3.6.4 to run cutadapt 2.4
    source /share/apps/source_files/python/python-3.6.4.source
    # add FastQC to my path
    export PATH=/share/apps/genomics/FastQC-0.11.9/:$PATH
    # add pigz to my path
    export PATH=/share/apps/pigz-2.6/:$PATH
    echo "Run Trim Galore command, with fastQC on the trimmed files"
    /share/apps/genomics/TrimGalore-0.6.7/trim_galore --fastqc --path_to_cutadapt /share/apps/genomics/cutadapt-2.5/bin/cutadapt --paired --trim1 --cores 1 --output_dir $OUTPUT_DIR_01 --no_report_file $INPUT_1 $INPUT_2  ## the help advise for 4 cores, but Ed Martin from CS advise that with pigz, 1 core is better on this system. Or we could use bgzip
else
    echo "Files already trimmed."
fi

echo "**** End of step 1: $(date) ****"

##############################################################################

## Check if the file was already processed for nbr reads and BS rate
if ! grep -q "${INPUT##*/}," "$DATA_TABLE_OUT"; then
    echo "Step 2: Measure the number of reads in the trimmed file → store in matrix"
    nbrReads=$(zcat "$TRIMMED_1" | wc -l | awk '{print $1 / 4}')

    ## ~20min/sample 4C 10G, ~25min/sample 4C 5G
    echo "**** End of step 2: $(date) ****"

    echo "Step 3: Bisulfite conversion efficiency measure (BCREval) → store in matrix" ### Takes ~1.5h with 5G and 10 cores (but no parallelisation here)
    python /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/01_prepareMethMatrixFromSamples/BCReval.py -n 8 -i "$TRIMMED_1" -o "$TEMP_OUTDIR/${INPUT##*/}.BCREval.out"

    ## Save the Conversion ratio (CR)
    CR=$(cat "$TEMP_OUTDIR/${INPUT##*/}.BCREval.out" | awk '{sum = $(NF-2) + $(NF-1) + $NF; mean = sum / 3; result = 1 - mean; printf "%.4f\n", result}')

    # Construct the output string using parameter expansion
    echo "${INPUT##*/},${nbrReads},${CR}" >> "$DATA_TABLE_OUT"

    ## Keep only unique rows, if we re-run the code on the same samples for debugging or so
    sort "$DATA_TABLE_OUT" | uniq > temp
    mv temp "$DATA_TABLE_OUT"

    echo "**** End of step 3: $(date) ****"
else 
    echo "Step 2 (calculate nbr reads) and step 3 (calculate BS conversion rate) already done before" 
fi

####################################################################
echo "**** Start of step 4: Bismark: Align to GRCh38 -> Output BAM"

## Output directory:
BISMARK_OUTDIR="$TEMP_OUTDIR/02Bismark"

# Create output directory if it does not exist
mkdir -p $BISMARK_OUTDIR

## see manual at https://felixkrueger.github.io/Bismark/

BISMARK="/share/apps/genomics/Bismark-0.22.3"
BOWTIE2="/share/apps/genomics/bowtie2-2.4.1/"

# add samtools to my path
export PATH=/share/apps/genomics/samtools-1.9/bin/:$PATH

GENOME_DIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/GRCh38"

## in GRCh38 folder: wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
## rename for Bismark: mv GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna.gz GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fa.gz  

#  Warning message "Chromosomal sequence could not be extracted for":
# It does occur sometimes whenever reads map to the very edges of chromosomes, or scaffolds, and is generally nothing you need to be worried about (losing a handful of reads in a few hundred million is not really a deal breaker...). (Felix Krueger)

echo "**** Start of step 4.1 Genome indexing..."

# Check if the folder exists
if [ ! -e "$GENOME_DIR/Bisulfite_Genome/" ]; then
    $BISMARK/bismark_genome_preparation --path_to_aligner $BOWTIE2 --verbose $GENOME_DIR
    echo "Genome indexed"
else
    echo "Genome already indexed"
fi

echo "**** End of step 4.1: $(date) ****"

echo "**** Start of step 4.2 Genome alignement with $NSLOTS threads..."
### takes 13h/samples with 5 cores 12Gb

# Check if the alignement was done
BAM="$BISMARK_OUTDIR/${INPUT##*/}_pe.bam"
if ! test -f "$BAM" ; then
    $BISMARK/bismark --genome $GENOME_DIR -parallel $NSLOTS --path_to_bowtie2 $BOWTIE2 --output_dir $BISMARK_OUTDIR -1 $TRIMMED_1 -2 $TRIMMED_2
    #NB: Specifying --basename in conjuction with --multicore is currently not supported (but we are aiming to fix this soon). Please lose either --basename or --multicore to proceed
    echo "Alignment done."
else
    echo "Alignement already done."
fi

echo "**** End of step 4.2: $(date) ****"

echo "**** Start of step 4.3 Deduplication of the bam file..."
## 2.5h with 12G 10nodes 72h

# Check if the deduplication was done
DEDUPBAM="$BISMARK_OUTDIR/${INPUT##*/}_pe.deduplicated.bam"
if ! test -f "$DEDUPBAM" ; then
    $BISMARK/deduplicate_bismark --bam $BAM --output_dir $BISMARK_OUTDIR
    echo "Deduplication done."
else
    echo "Deduplication already done."
fi

echo "**** End of step 4.3: $(date) ****"

echo "4.4 Methylation extraction..."
$BISMARK/bismark_methylation_extractor --gzip --ignore_r2 2 -o $BISMARK_OUTDIR -parallel $NSLOTS $DEDUPBAM

## --ignore_r2 <int>        Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing
##                         results only. Since the first couple of bases in Read 2 of BS-Seq experiments
##                         show a severe bias towards non-methylation as a result of end-repairing
##                         sonicated fragments with unmethylated cytosines (see M-bias plot), it is
##                         recommended that the first couple of bp of Read 2 are removed before
##                         starting downstream analysis. Please see the section on M-bias plots in the
##                         Bismark User Guide for more details.
#
## This will produce three methytlation output files:
## CpG_context_SAMPLE_bismark_bt2.txt.gz
## CHG_context_SAMPLE_bismark_bt2.txt.gz
## CHH_context_SAMPLE_bismark_bt2.txt.gz
## as well as a Bismark coverage file: SAMPLE_pe.deduplicated.bismark.cov.gz

## move the needed final file to SAN to keep it for further analyses
cp $BISMARK_OUTDIR/${INPUT##*/}_PE_report.txt /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/01Methcall/
cp $BISMARK_OUTDIR/${INPUT##*/}_pe.deduplicated.bismark.cov.gz /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/01Methcall/

echo "**** End of step 4.4: $(date) ****"
## took 4h30 on 12G 10nodes

### TO DO WHEN DO
### Remove intermediate files (or in scratch?)
#
#echo "In the future, do variant calling with code from Noah Kessler if he shares!"
### Maria: RM single nucleotide polymorphism (SNP)-related probes identified by Zhou et al. (48) that contain SNPs (MAF > 1%) that are within 5 bp of the CpG interrogation site and/or SNPs effecting probe hybridization
#
### NB one VCF can be ~5GB compressed; BCF can be better for compression
#
## Log the end of the job
echo "**** Job $JOB_NAME.$JOB_ID finished at $(date) ****"
