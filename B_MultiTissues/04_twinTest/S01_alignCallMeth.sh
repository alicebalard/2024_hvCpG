#!/bin/bash
#$ -N runPipelineTwins
#$ -S /bin/bash
#$ -l tmem=17G
#$ -l h_vmem=17G
#$ -t 1-14          ## <- number of SAMPLES to process (7 MZ pairs = 14 samples). Update!
#$ -tc 50
#$ -pe smp 6
#$ -l h_rt=240:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

## ─────────────────────────────────────────────────────────────────────────────
## WGBS pipeline: FASTQ -> methylation calls (.cov), MZ twin samples.
## Caucasian female twins belonging to the TwinsUK/MuTHER cohort
## WGBS on adipose tissues biopsie (seven MZ pairs, six DZ pairs)
##   - drives on SAMPLES (sample_accession_id), not files
##   - each sample has MULTIPLE runs (Run1/Run2...) and R1/R2, each in its own
##     file_accession_id subfolder -> we merge all runs per sample first
##   - filters to MZ samples from sample_file.csv (sample_alias contains "MZ")
## ─────────────────────────────────────────────────────────────────────────────

CODEDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG"
DATADIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015"
FASTQROOT="$DATADIR/AdiposeWGBS"                       # contains one subfolder per file_accession_id
GENOME_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/GRCh38"

SAMPLE_CSV="$DATADIR/sample_file_adiposeWGBS.csv"
## ZYGOSITY="MZ"       ## <- "MZ" or "DZ": which twins to process (change + rerun for DZ)
ZYGOSITY="DZ" 

BISMARK="/share/apps/genomics/Bismark-0.22.3"
BOWTIE2="/share/apps/genomics/bowtie2-2.4.1/"
export PATH=/share/apps/genomics/samtools-1.9/bin/:$PATH
export PATH=/share/apps/genomics/FastQC-0.11.9/:$PATH
export PATH=/share/apps/pigz-2.6/:$PATH

TEMP_OUTDIR="$DATADIR/TEMP"; mkdir -p "$TEMP_OUTDIR"
TRIMOM_OUTDIR="$TEMP_OUTDIR/01Trimmed_data"; mkdir -p "$TRIMOM_OUTDIR"
BISMARK_OUTDIR="$TEMP_OUTDIR/02Bismark";     mkdir -p "$BISMARK_OUTDIR"
METHCALL_OUTDIR="$DATADIR/01Methcall";       mkdir -p "$METHCALL_OUTDIR"
cd "$TEMP_OUTDIR"

## ── Build the per-sample list ONCE (task 1 territory; safe if it already exists)
## One line per SAMPLE = sample_accession_id, restricted to the chosen zygosity.
SAMPLE_LIST="$TEMP_OUTDIR/samples_${ZYGOSITY}.tmp"
if [ ! -e "$SAMPLE_LIST" ]; then
    ## columns: 1 sample_accession_id, 2 sample_alias, 3 file_name, 4 file_accession_id
    ## keep rows whose alias contains the zygosity, take unique sample_accession_id
    awk -F',' -v z="$ZYGOSITY" 'NR>1 && $2 ~ z {print $1}' "$SAMPLE_CSV" \
        | sort -u > "$SAMPLE_LIST"
    echo "Samples to process ($ZYGOSITY): $(wc -l < "$SAMPLE_LIST")"
fi

## Pick this task's sample
SAMPLE=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST")
[ -z "$SAMPLE" ] && { echo "No sample for task $SGE_TASK_ID"; exit 0; }

echo "**** Job $JOB_NAME.$JOB_ID task $SGE_TASK_ID started $(date) ****"
echo "Sample (sample_accession_id): $SAMPLE"

## ─────────────────────────────────────────────────────────────────────────────
## STEP 0. Gather + MERGE all runs of this sample into one R1 and one R2
## For this sample_accession_id, sample_file.csv lists several file_accession_id
## folders (Run1.R1, Run1.R2, Run2.R1, ...). We concatenate all R1s and all R2s.
## ─────────────────────────────────────────────────────────────────────────────
MERGED_R1="$TEMP_OUTDIR/${SAMPLE}_R1.fastq.gz"
MERGED_R2="$TEMP_OUTDIR/${SAMPLE}_R2.fastq.gz"

if [ ! -e "$MERGED_R1" ] || [ ! -e "$MERGED_R2" ]; then
    # for this sample: get (file_accession_id, file_name) of its R1 and R2 files
    #   R1 rows -> concat into MERGED_R1 ; R2 rows -> concat into MERGED_R2
    # match R1/R2 by the ".R1." / ".R2." token in file_name
    R1_FILES=$(awk -F',' -v s="$SAMPLE" '$1==s && $3 ~ /\.R1\./ {print $4"/"$3}' "$SAMPLE_CSV")
    R2_FILES=$(awk -F',' -v s="$SAMPLE" '$1==s && $3 ~ /\.R2\./ {print $4"/"$3}' "$SAMPLE_CSV")

    echo "  R1 run files:"; echo "$R1_FILES" | sed 's/^/    /'
    echo "  R2 run files:"; echo "$R2_FILES" | sed 's/^/    /'

    # concatenating gzipped FASTQs is valid (gzip streams concatenate)
    : > "$MERGED_R1"; for f in $R1_FILES; do cat "$FASTQROOT/$f" >> "$MERGED_R1"; done
    : > "$MERGED_R2"; for f in $R2_FILES; do cat "$FASTQROOT/$f" >> "$MERGED_R2"; done
fi

INPUT_1="$MERGED_R1"
INPUT_2="$MERGED_R2"

## skip everything if the final .cov already exists (resume support)
METHCOV="$METHCALL_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
if [ -e "$METHCOV" ]; then
    echo "Already processed to methylation calls, nothing to do: $METHCOV"
    exit 0
fi

## ── STEP 1. Trimming (Trim Galore) ───────────────────────────────────────────
echo "**** Step 1 Trim Galore $(date) ****"
source /share/apps/source_files/python/python-3.6.4.source     # for cutadapt
TRIMMED_1="$TRIMOM_OUTDIR/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_2="$TRIMOM_OUTDIR/${SAMPLE}_R2_val_2.fq.gz"
if [ ! -e "$TRIMMED_2" ]; then
    /share/apps/genomics/TrimGalore-0.6.7/trim_galore \
        --path_to_cutadapt /share/apps/genomics/cutadapt-2.5/bin/cutadapt \
        --paired --trim1 --cores 1 --no_report_file \
        --output_dir "$TRIMOM_OUTDIR" "$INPUT_1" "$INPUT_2"
fi
## Trim Galore names outputs after the input basename; normalise to our names
mv -f "$TRIMOM_OUTDIR/${SAMPLE}_R1_val_1.fq.gz" "$TRIMMED_1" 2>/dev/null || true
mv -f "$TRIMOM_OUTDIR/${SAMPLE}_R2_val_2.fq.gz" "$TRIMMED_2" 2>/dev/null || true
echo "**** End step 1 $(date) ****"

## ── STEP 2. Bismark genome indexing (once) ───────────────────────────────────
if [ ! -e "$GENOME_DIR/Bisulfite_Genome/" ]; then
    $BISMARK/bismark_genome_preparation --path_to_aligner $BOWTIE2 --verbose "$GENOME_DIR"
fi

## ── STEP 3. Bismark alignment ────────────────────────────────────────────────
echo "**** Step 3 Bismark align ($NSLOTS threads) $(date) ****"
BAM="$BISMARK_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_pe.bam"
if [ ! -e "$BAM" ]; then
    $BISMARK/bismark --genome "$GENOME_DIR" --parallel $NSLOTS --gzip \
        --path_to_bowtie2 $BOWTIE2 --temp_dir "$BISMARK_OUTDIR" \
        --output_dir "$BISMARK_OUTDIR" -1 "$TRIMMED_1" -2 "$TRIMMED_2"
fi
rm -f "$TRIMMED_1" "$TRIMMED_2" "$MERGED_R1" "$MERGED_R2"   # free space
echo "**** End step 3 $(date) ****"

## ── STEP 4. Deduplicate ──────────────────────────────────────────────────────
echo "**** Step 4 Deduplicate $(date) ****"
$BISMARK/deduplicate_bismark --bam "$BAM" --output_dir "$BISMARK_OUTDIR"
DEDUPBAM="$BISMARK_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
rm -f "$BAM"
echo "**** End step 4 $(date) ****"

## ── STEP 5. Methylation extraction -> .cov ───────────────────────────────────
echo "**** Step 5 Methylation extraction $(date) ****"
$BISMARK/bismark_methylation_extractor --gzip --bedGraph --ignore_r2 2 \
    -o "$BISMARK_OUTDIR" --genome_folder "$GENOME_DIR" --parallel $NSLOTS "$DEDUPBAM"

COV="$BISMARK_OUTDIR/$(basename "${DEDUPBAM%.*}").bismark.cov.gz"
mv -f "$BISMARK_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_PE_report.txt" "$METHCALL_OUTDIR/" 2>/dev/null || true
mv -f "$COV" "$METHCALL_OUTDIR/"
## also keep a 20X-filtered version
zcat "$METHCALL_OUTDIR/$(basename "$COV")" | awk '$5+$6 >= 20' | gzip \
    > "$METHCALL_OUTDIR/$(basename "$COV").20X.cov.gz"
echo "**** End step 5 $(date) ****"

## ── CLEAN UP ─────────────────────────────────────────────────────────────────
function finish {
    rm -rf "$TRIMOM_OUTDIR"/*"${SAMPLE}"* "$BISMARK_OUTDIR"/*"${SAMPLE}"*
    rm -f "$MERGED_R1" "$MERGED_R2"
}
trap finish EXIT ERR INT TERM

echo "**** Job finished $(date) ****"
