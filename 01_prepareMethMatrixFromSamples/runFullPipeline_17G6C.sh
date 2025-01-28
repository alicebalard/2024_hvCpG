#!/bin/bash
#$ -N runPipeline
#$ -S /bin/bash
#$ -l tmem=17G
#$ -l h_vmem=17G
#$ -t 1-2
#$ -tc 2
#$ -pe smp 6 # Request N cores per task 
#$ -l tscratch=1000G # request scratch space
#$ -l h_rt=40:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

# Start resource monitoring in the background
bash ~/monitor_resources.sh &
echo "Number of cores requested per task: $NSLOTS"

## Runs in ~30h/sample with 6C and 17G

##############
##************
## START UP ##
##**********##
##############

## We make a temporary folder for all intermediate files that we delete at the end.
mkdir -p /scratch0/abalard/TEMP.${JOB_ID}
TEMP_OUTDIR="/scratch0/abalard/TEMP.${JOB_ID}"
cd $TEMP_OUTDIR

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

######################
##********************
## STEP 1. Trimming ##
##********************
######################

echo "**** Start of step 1: Trim galore with automatic parameters $(date) ****" 
## ~6h/sample with 1 thread, 5Gb
## ~5h/sample for 4C 10G
## ~4h/sample for 4C 5G
## 4h30 to 5h30 with 17G 6C
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

############################################################
##********************************************************##
## STEP 2 & 3. Calculate nbr reads and BS conversion rate ##
##********************************************************##
############################################################

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

#####################################################
##*************************************************##
## STEP 4. Bismark alignement and methylation call ##
##*************************************************##
#####################################################

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
## crashes on 2 cores 7Gb
## crashes with 2 cores and 15Gb
## takes 24 to 32h/samples with 4 cores 15Gb
## 18-19h with 17G 6C

# Bismark holds the reference genome in memory and in addition to that runs four parallel instances of Bowtie. The memory usage is dependent on the size of the reference genome. For a large eukaryotic genome (human or mouse) we experienced a typical memory usage of around 12GB. We thus recommend running Bismark on a machine with 5 CPU cores and at least 12 GB of RAM. The memory requirements of Bowtie 2 are somewhat larger (possibly to allow gapped alignments). When running Bismark using Bowtie 2 we therefore recommend a system with at least 5 cores and > 16GB of RAM.

# Check if the alignement was done
BAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.bam"
if ! test -f "$BAM" ; then
    $BISMARK/bismark --genome $GENOME_DIR --parallel $NSLOTS --path_to_bowtie2 $BOWTIE2 --output_dir $BISMARK_OUTDIR -1 $TRIMMED_1 -2 $TRIMMED_2
    echo "Alignment done."
else
    echo "Alignement already done."
fi
echo "**** End of step 4.2: $(date) ****"

echo "**** Start of step 4.3 Deduplication of the bam file..."
## 2.5h with 12G 10nodes 72h
## 1.5-2h with 17G 6C

# Run if deduplication was not done
DEDUPBAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bam"
if ! test -f "$DEDUPBAM" ; then
    $BISMARK/deduplicate_bismark --bam $BAM --output_dir $BISMARK_OUTDIR
    echo "Deduplication done."
else
    echo "Deduplication already done."
fi
echo "**** End of step 4.3: $(date) ****"

echo "4.4 Methylation extraction..."
## 3h in 17G 6C
## took 4h30 on 12G 10nodes

METHCOV="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"

if ! test -f "$METHCOV" ; then
    $BISMARK/bismark_methylation_extractor --gzip --ignore_r2 2 -o $BISMARK_OUTDIR -parallel $NSLOTS --merge_non_CpG --bedGraph $DEDUPBAM
    echo "Methylation called."
else
    echo "Methylation already called."
fi
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
cp $BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_PE_report.txt /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/01Methcall/
cp $BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/01Methcall/

echo "**** End of step 4.4: $(date) ****"

##################################
##******************************##
## STEP 5. bssnper2 SNP calling ##
##******************************##
##################################
### Maria's paper: RM single nucleotide polymorphism (SNP)-related probes identified by Zhou et al. (48) that contain SNPs (MAF > 1%) that are within 5 bp of the CpG interrogation site and/or SNPs effecting probe hybridization

## bssnper2 code is from Noah Kessler
## ~1h with 17G 6C

echo "**** Start of sorting test bam file : $(date) ****" 
## The bam file needs to be sorted before bssnper2
if ! test -f "$DEDUPBAM.sorted.bam" ; then
    samtools sort -o $DEDUPBAM.sorted.bam --threads $NSLOTS $DEDUPBAM
else
    echo "Bam file already sorted."
fi
echo "**** End of sorting test bam file : $(date) ****"

echo "**** Start of bssbper2 SNP call : $(date) ****" 
BSSNPER2_OUTDIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/02SNPcall/" # not a temporary file

## run bssnper2 if not called yet
if ! test -f "$BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf.gz" ; then
    /home/abalard/bssnper2/bssnper2 $DEDUPBAM.sorted.bam --ref $GENOME_DIR/GCF_000001405.40_GRCh38.p14_genomic.fa --vcf $BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf
    ## Compress and index vcf files
    /share/apps/htslib-1.20/bgzip --threads $NSLOTS $BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf
    /share/apps/htslib-1.20/tabix -p vcf --threads $NSLOTS $BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf.gz
else
    echo "SNPs already called"
fi

echo "**** End of bssbper2 : $(date) ****" 

##############
##************
## CLEAN UP ##
##**********##
##############

## What is needed for age and sex determination? figure out before rm the files!

## rm raw data:
# rm -rf /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/00RawFastq/*${INPUT##*/}*

## Rm TEMP directory whenever the job exits, regardless of whether it finished successfully or not
function finish {
    rm -rf $OUTPUT_DIR_01/*${INPUT##*/}*
    rm -rf $BISMARK_OUTDIR/*${INPUT##*/}*
    rm -rf $BSSNPER2_OUTDIR/*${INPUT##*/}*
}
trap finish EXIT ERR INT TERM

## Log the end of the job
echo "**** Job $JOB_NAME.$JOB_ID finished at $(date) ****"
