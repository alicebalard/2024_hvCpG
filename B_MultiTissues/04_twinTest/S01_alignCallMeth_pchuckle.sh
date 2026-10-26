#!/bin/bash
#$ -N runPipelineTwins
#$ -S /bin/bash
#$ -l tmem=17G
#$ -l h_vmem=17G
#$ -t 1-60          ## generous: adipose=26 samples, blood=?; extra tasks exit harmlessly. Raise if a tissue has >60 samples.
#$ -tc 50
#$ -pe smp 6
#$ -l h_rt=240:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

## FASTQ -> methylation calls (.cov) for TwinsUK WGBS twins, any tissue, all zygosities.
## Usage: qsub -v TISSUE=adipose  S01_alignCallMeth.sh      (default)
##        qsub -v TISSUE=blood    S01_alignCallMeth.sh
set -euo pipefail

TISSUE="${TISSUE:-adipose}"                         # adipose | blood
DATADIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015"
GENOME_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/GRCh38"

case "$TISSUE" in
  adipose) FASTQSUB="AdiposeWGBS"; METHCALL_OUTDIR="$DATADIR/01Methcall/Adipose" ;;
  blood)   FASTQSUB="BloodWGBS";   METHCALL_OUTDIR="$DATADIR/01Methcall/Blood" ;;
  *) echo "TISSUE must be adipose or blood" >&2; exit 1 ;;
esac
FASTQROOT="$DATADIR/$FASTQSUB"
SAMPLE_CSV="$DATADIR/sample_file_${TISSUE}WGBS.csv"   # sample_file_adiposeWGBS.csv / sample_file_bloodWGBS.csv

BISMARK="/share/apps/genomics/Bismark-0.22.3"
BOWTIE2="/share/apps/genomics/bowtie2-2.4.1/"
export PATH=/share/apps/genomics/samtools-1.9/bin/:$PATH
export PATH=/share/apps/genomics/FastQC-0.11.9/:$PATH
export PATH=/share/apps/pigz-2.6/:$PATH

TEMP_OUTDIR="$DATADIR/TEMP_$TISSUE"; mkdir -p "$TEMP_OUTDIR"
TRIMOM_OUTDIR="$TEMP_OUTDIR/01Trimmed_data"; mkdir -p "$TRIMOM_OUTDIR"
BISMARK_OUTDIR="$TEMP_OUTDIR/02Bismark";     mkdir -p "$BISMARK_OUTDIR"
mkdir -p "$METHCALL_OUTDIR"
cd "$TEMP_OUTDIR"

## ALL samples of this tissue (both MZ and DZ) -> one line per sample_accession_id
SAMPLE_LIST="$TEMP_OUTDIR/samples_${TISSUE}.tmp"
if [ ! -e "$SAMPLE_LIST" ]; then
    awk -F',' 'NR>1 {print $1}' "$SAMPLE_CSV" | sort -u > "$SAMPLE_LIST"
    echo "Samples to process ($TISSUE): $(wc -l < "$SAMPLE_LIST")  <-- set -t to at least this"
fi

SAMPLE=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST")
[ -z "$SAMPLE" ] && { echo "No sample for task $SGE_TASK_ID (list exhausted)"; exit 0; }

echo "**** $JOB_NAME.$JOB_ID task $SGE_TASK_ID  tissue=$TISSUE  sample=$SAMPLE  $(date) ****"

## STEP 0. Merge all runs of this sample into one R1 / R2
MERGED_R1="$TEMP_OUTDIR/${SAMPLE}_R1.fastq.gz"
MERGED_R2="$TEMP_OUTDIR/${SAMPLE}_R2.fastq.gz"
if [ ! -e "$MERGED_R1" ] || [ ! -e "$MERGED_R2" ]; then
    R1_FILES=$(awk -F',' -v s="$SAMPLE" '$1==s && $3 ~ /\.R1\./ {print $4"/"$3}' "$SAMPLE_CSV")
    R2_FILES=$(awk -F',' -v s="$SAMPLE" '$1==s && $3 ~ /\.R2\./ {print $4"/"$3}' "$SAMPLE_CSV")
    : > "$MERGED_R1"; for f in $R1_FILES; do cat "$FASTQROOT/$f" >> "$MERGED_R1"; done
    : > "$MERGED_R2"; for f in $R2_FILES; do cat "$FASTQROOT/$f" >> "$MERGED_R2"; done
fi
INPUT_1="$MERGED_R1"; INPUT_2="$MERGED_R2"

## resume: skip if final .cov exists
METHCOV="$METHCALL_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
if [ -e "$METHCOV" ]; then echo "Already done: $METHCOV"; exit 0; fi

## STEP 1. Trim Galore
source /share/apps/source_files/python/python-3.6.4.source
TRIMMED_1="$TRIMOM_OUTDIR/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_2="$TRIMOM_OUTDIR/${SAMPLE}_R2_val_2.fq.gz"
if [ ! -e "$TRIMMED_2" ]; then
    /share/apps/genomics/TrimGalore-0.6.7/trim_galore \
        --path_to_cutadapt /share/apps/genomics/cutadapt-2.5/bin/cutadapt \
        --paired --trim1 --cores 1 --no_report_file \
        --output_dir "$TRIMOM_OUTDIR" "$INPUT_1" "$INPUT_2"
fi
mv -f "$TRIMOM_OUTDIR/${SAMPLE}_R1_val_1.fq.gz" "$TRIMMED_1" 2>/dev/null || true
mv -f "$TRIMOM_OUTDIR/${SAMPLE}_R2_val_2.fq.gz" "$TRIMMED_2" 2>/dev/null || true

## STEP 2. Genome index (once)
if [ ! -e "$GENOME_DIR/Bisulfite_Genome/" ]; then
    $BISMARK/bismark_genome_preparation --path_to_aligner $BOWTIE2 --verbose "$GENOME_DIR"
fi

## STEP 3. Align
BAM="$BISMARK_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_pe.bam"
if [ ! -e "$BAM" ]; then
    $BISMARK/bismark --genome "$GENOME_DIR" --parallel $NSLOTS --gzip \
        --path_to_bowtie2 $BOWTIE2 --temp_dir "$BISMARK_OUTDIR" \
        --output_dir "$BISMARK_OUTDIR" -1 "$TRIMMED_1" -2 "$TRIMMED_2"
fi
rm -f "$TRIMMED_1" "$TRIMMED_2" "$MERGED_R1" "$MERGED_R2"

## STEP 4. Deduplicate
$BISMARK/deduplicate_bismark --bam "$BAM" --output_dir "$BISMARK_OUTDIR"
DEDUPBAM="$BISMARK_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
rm -f "$BAM"

## STEP 5. Methylation extraction -> .cov
$BISMARK/bismark_methylation_extractor --gzip --bedGraph --ignore_r2 2 \
    -o "$BISMARK_OUTDIR" --genome_folder "$GENOME_DIR" --parallel $NSLOTS "$DEDUPBAM"
COV="$BISMARK_OUTDIR/$(basename "${DEDUPBAM%.*}").bismark.cov.gz"
mv -f "$BISMARK_OUTDIR/${SAMPLE}_R1_val_1_bismark_bt2_PE_report.txt" "$METHCALL_OUTDIR/" 2>/dev/null || true
mv -f "$COV" "$METHCALL_OUTDIR/"

function finish { rm -rf "$TRIMOM_OUTDIR"/*"${SAMPLE}"* "$BISMARK_OUTDIR"/*"${SAMPLE}"*; rm -f "$MERGED_R1" "$MERGED_R2"; }
trap finish EXIT ERR INT TERM
echo "**** done $SAMPLE $(date) ****"
