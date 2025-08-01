#!/bin/bash
#$ -N runPipeline
#$ -S /bin/bash
#$ -l tmem=17G
#$ -l h_vmem=17G
#$ -t 1-40
#$ -pe smp 6 # Request N cores per task 
#$ -l h_rt=240:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

## Runs in ~30h/sample with 6C and 17G

CODEDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG"
DATADIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human"

## Step 1 - Creates a TRIMOM_OUTDIR called TEMP/01Trimmed_data. Output: TRIMMED_1="$TRIMOM_OUTDIR/${INPUT##*/}_1_val_1.fq.gz"; TRIMMED_2="$TRIMOM_OUTDIR/${INPUT##*/}_2_val_2.fq.gz"

## Step 2 & 3: Fill in DATA_TABLE_OUT ($CODEDIR/01_prepareMethMatrixFromSamples/infoAllsamples_Pipeline.csv) with nbr reads and BS conversion rate

## Step 4: Methylation extraction 
## 4.2. Align. output: BAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.bam"
## 4.3. Deduplicated BAM. output: DEDUPBAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bam"
## 4.4. Methylation extraction. output: METHCOV="$DATADIR/01Methcall/$(basename "${DEDUPBAM%.*}").bismark.cov.gz" and $METHCOV.20X.cov.gz

## Step 5. SNP call with bssnper2; Create BSSNPER2_OUTDIR="$DATADIR/02SNPcall/"

##************
## START UP ##
##**********##

## We make a temporary folder for all intermediate files that we delete at the end.
TEMP_OUTDIR="$DATADIR/TEMP"

mkdir -p $TEMP_OUTDIR
cd $TEMP_OUTDIR

# Start resource monitoring in the background
## source ~/monitor_resources.sh & # source allows to use variables define in this script, in particular $TEMP_OUTDIR. Log stored in WGBS_human
## NB: use only for 1-2 samples otherwise it's a massive log file

# Create the files to loop over if it does not exist:
if [ ! -e "$TEMP_OUTDIR/bysample_list_of_files.tmp" ]; then
    ## NB: input files must be named filename_1.fastq.gz and filename_2.fastq.gz
    ls -1 $DATADIR/00RawFastq/*_1.fastq.gz > $TEMP_OUTDIR/list_of_files.tmp
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
DATA_TABLE_OUT="$CODEDIR/01_prepareMethMatrixFromSamples/infoAllsamples_Pipeline.csv"

# Check if the files exists
if [ ! -e "$DATA_TABLE_OUT" ]; then
    # Create the file
    echo "Sample,Nbr_reads,BS_conversion_rate" > $DATA_TABLE_OUT
    echo "File created: $DATA_TABLE_OUT"
else
    echo "File already exists: $DATA_TABLE_OUT"
fi

## Create output directories if they do not exist

TRIMOM_OUTDIR="$TEMP_OUTDIR/01Trimmed_data"
mkdir -p $TRIMOM_OUTDIR

BISMARK_OUTDIR="$TEMP_OUTDIR/02Bismark"
mkdir -p $BISMARK_OUTDIR

#########################################################################
## An if statement to see if the file has already been (partly) processed

## if METHCOV does not exist, then only run steps 1 to 4
METHCOV="$DATADIR/01Methcall/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"

if [ -e "$METHCOV" ]; then
    echo "Methcov has already been processed, jump to step 5: "
else
    echo "Methcov has not been processed, run step 1 to 4: "
    
    ##********************
    ## STEP 1. Trimming ##
    ##********************

    echo "**** Start of step 1: Trim galore with automatic parameters $(date) ****" 
    ## ~6h/sample with 1 thread, 5Gb; ~5h/sample for 4C 10G; ~4h/sample for 4C 5G; 4h30 to 5h30 with 17G 6C

    # source python 3.6.4 to run cutadapt 2.4
    source /share/apps/source_files/python/python-3.6.4.source
    # add FastQC to my path
    export PATH=/share/apps/genomics/FastQC-0.11.9/:$PATH
    # add pigz to my path
    export PATH=/share/apps/pigz-2.6/:$PATH

    ## Check if the file was already processed for nbr reads and BS rate
    if [ ! -e "$TRIMOM_OUTDIR/${INPUT##*/}_2_val_2.fq.gz" ]; then
	echo "Run Trim Galore command without fastqc"
	/share/apps/genomics/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /share/apps/genomics/cutadapt-2.5/bin/cutadapt --paired --trim1 --cores 1 --output_dir $TRIMOM_OUTDIR --no_report_file $INPUT_1 $INPUT_2  ## the help advise for 4 cores, but Ed Martin from CS advise that with pigz, 1 core is better on this system. Or we could use bgzip
    else
	echo "File already trimmed"
    fi

    echo "**** End of step 1: $(date) ****"

    ## output:
    TRIMMED_1="$TRIMOM_OUTDIR/${INPUT##*/}_1_val_1.fq.gz"
    TRIMMED_2="$TRIMOM_OUTDIR/${INPUT##*/}_2_val_2.fq.gz"

    ##********************************************************##
    ## STEP 2 & 3. Calculate nbr reads and BS conversion rate ##
    ##********************************************************##

    ## Check if the file was already processed for nbr reads and BS rate
    if ! grep -q "${INPUT##*/}," "$DATA_TABLE_OUT"; then
	echo "Step 2: Measure the number of reads in the trimmed file → store in matrix"
	nbrReads=$(zcat "$TRIMMED_1" | wc -l | awk '{print $1 / 4}')

	## ~20min/sample 4C 10G, ~25min/sample 4C 5G
	echo "**** End of step 2: $(date) ****"

	echo "Step 3: Bisulfite conversion efficiency measure (BCREval) → store in matrix" ### Takes ~1.5h with 5G and 10 cores (but no parallelisation here)
	python $CODEDIR/01_prepareMethMatrixFromSamples/BCReval.py -n 8 -i "$TRIMMED_1" -o "$TEMP_OUTDIR/${INPUT##*/}.BCREval.out"

	## Save the Conversion ratio (CR)
	CR=$(cat "$TEMP_OUTDIR/${INPUT##*/}.BCREval.out" | awk '{sum = $(NF-2) + $(NF-1) + $NF; mean = sum / 3; result = 1 - mean; printf "%.4f\n", result}')

	# Construct the output string using parameter expansion
	echo "${INPUT##*/},${nbrReads},${CR}" >> "$DATA_TABLE_OUT"

	## rm temporary table
	rm "$TEMP_OUTDIR/${INPUT##*/}.BCREval.out"

	## Keep only unique rows, if we re-run the code on the same samples for debugging or so
	sort "$DATA_TABLE_OUT" | uniq > temp
	mv temp "$DATA_TABLE_OUT"

	echo "**** End of step 3: $(date) ****"
    else 
	echo "Step 2 (calculate nbr reads) and step 3 (calculate BS conversion rate) already done before" 
    fi

    ##*************************************************##
    ## STEP 4. Bismark alignement and methylation call ##
    ##*************************************************##

    echo "**** Start of step 4: Bismark: Align to GRCh38 -> Output BAM"
    ## see manual at https://felixkrueger.github.io/Bismark/

    BISMARK="/share/apps/genomics/Bismark-0.22.3"
    BOWTIE2="/share/apps/genomics/bowtie2-2.4.1/"

    # add samtools to my path
    export PATH=/share/apps/genomics/samtools-1.9/bin/:$PATH

    GENOME_DIR="$DATADIR/GRCh38"

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

    ## NB: this steps generates ~90G/samples of temporary files
    ### takes 13h/samples with 5 cores 12Gb; crashes on 2 cores 7Gb; crashes with 2 cores and 15Gb;  takes 24 to 32h/samples with 4 cores 15Gb; 18-19h with 17G 6C

    # Bismark holds the reference genome in memory and in addition to that runs four parallel instances of Bowtie. The memory usage is dependent on the size of the reference genome. For a large eukaryotic genome (human or mouse) we experienced a typical memory usage of around 12GB. We thus recommend running Bismark on a machine with 5 CPU cores and at least 12 GB of RAM. The memory requirements of Bowtie 2 are somewhat larger (possibly to allow gapped alignments). When running Bismark using Bowtie 2 we therefore recommend a system with at least 5 cores and > 16GB of RAM.

    ## Check if the file was already processed for nbr reads and BS rate
    if [ ! -e "$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.bam" ]; then
	echo "Run Trim Galore command without fastqc"
	$BISMARK/bismark --genome $GENOME_DIR --parallel $NSLOTS --gzip --path_to_bowtie2 $BOWTIE2 --temp_dir $BISMARK_OUTDIR --output_dir $BISMARK_OUTDIR -1 $TRIMMED_1 -2 $TRIMMED_2
    else
	echo "File already trimmed"
    fi

    echo "**** End of step 4.2: $(date) ****"

    ## output:
    BAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.bam"

    ## Rm trimmed files, not useful any longer (~ 30G/sample):
    rm -rf $TRIMOM_OUTDIR/*${INPUT##*/}*

    echo "**** Start of step 4.3 Deduplication of the bam file..."
    ## Tests: 2.5h with 12G 10nodes 72h; 1.5-2h with 17G 6C

    $BISMARK/deduplicate_bismark --bam $BAM --output_dir $BISMARK_OUTDIR

    echo "**** End of step 4.3: $(date) ****"

    ## output:
    DEDUPBAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bam"

    ## Rm bam files, not useful any longer (~40G):
    rm $BAM

    echo "Start of step 4.4: Methylation extraction $(date)"
    ## Tests: 3h in 17G 6C; took 4h30 on 12G 10nodes

    $BISMARK/bismark_methylation_extractor --gzip --bedGraph --ignore_r2 2 -o $BISMARK_OUTDIR --genome_folder $GENOME_DIR --parallel $NSLOTS $DEDUPBAM
    ## NB --bedGraph important option to get cov file

    METHCOV="${BISMARK_OUTDIR}/$(basename "${DEDUPBAM%.*}").bismark.cov.gz"

    ## move the needed final files to SAN to keep it for further analyses
    mv $BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_PE_report.txt $DATADIR/01Methcall/
    mv $METHCOV $DATADIR/01Methcall/

    ## redefine METHCOV
    METHCOV="$DATADIR/01Methcall/$(basename "${DEDUPBAM%.*}").bismark.cov.gz"

    ## Select 20X minimum coverage 
    less $METHCOV | awk '$5 + $6 >= 20' - | gzip > $METHCOV.20X.cov.gz

    echo "**** End of step 4.4: $(date) ****"

fi

##******************************##
## STEP 5. bssnper2 SNP calling ##
##******************************##

DEDUPBAM="$BISMARK_OUTDIR/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bam"

### Maria's paper: RM single nucleotide polymorphism (SNP)-related probes identified by Zhou et al. (48) that contain SNPs (MAF > 1%) that are within 5 bp of the CpG interrogation site and/or SNPs effecting probe hybridization

## bssnper2 code is from Noah Kessler
## ~1h with 17G 6C

echo "**** Start of sorting test bam file : $(date) ****" 

## The bam file needs to be sorted before bssnper2
samtools sort -o $DEDUPBAM.sorted.bam --threads $NSLOTS $DEDUPBAM

echo "**** Start of bssbper2 SNP call : $(date) ****" 
BSSNPER2_OUTDIR="$DATADIR/02SNPcall/" # not a temporary file

if [ ! -e "$BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf.gz" ]; then
    ## run bssnper2
    /home/abalard/bssnper2/bssnper2 $DEDUPBAM.sorted.bam --ref $GENOME_DIR/GCF_000001405.40_GRCh38.p14_genomic.fa --vcf $BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf
    rm $DEDUPBAM
    ## Compress and index vcf files
    /share/apps/htslib-1.20/bgzip --threads $NSLOTS $BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf
    /share/apps/htslib-1.20/tabix -p vcf --threads $NSLOTS $BSSNPER2_OUTDIR/${DEDUPBAM##*/}.vcf.gz
else
    echo "bssnper2 already ran on this sample"
fi

echo "**** End of bssbper2 : $(date) ****" 

##************
## CLEAN UP ##
##**********##

## Rm TEMP directory whenever the job exits, regardless of whether it finished successfully or not
function finish {
    rm -rf $TRIMOM_OUTDIR/*${INPUT##*/}*
    rm -rf $BISMARK_OUTDIR/*${INPUT##*/}* # remove all bam and methylation called files
}

trap finish EXIT ERR INT TERM ## allows for cleanup operations or other final tasks to be performed regardless of how the script terminates

## Log the end of the job
echo "**** Job $JOB_NAME.$JOB_ID finished at $(date) ****"
