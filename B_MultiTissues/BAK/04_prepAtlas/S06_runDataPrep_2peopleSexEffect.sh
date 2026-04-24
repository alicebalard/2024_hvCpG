#!/bin/bash
#$ -N prepMetadata2sexgroups
#$ -S /bin/bash
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
source /share/apps/source_files/python/python-3.13.0a6.source

DIR_DATA="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/"
DIR_CODE="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/"
META_IN="${DIR_CODE}SupTab1_Loyfer2023_amended.csv"

prep_atlas() {
    # Arguments
    STEP=$1             # e.g., 1 or 2 or 5
    NAME=$2             # name of the analysis (affects file/script naming)
    MIN_DATASETS=$3     # CpG covered in at least MIN_DATASETS
    MIN_SAMPLES=${4:-3} # min samples per group (default 3; pairs use 2)

    # Derived paths
    META_OUT="${META_IN}_${STEP}_${NAME}.csv"
    FILES_OUT="${DIR_DATA}10X_${STEP}_${NAME}"
    CODE_1="${DIR_CODE}prepare_metadata_${STEP}_${NAME}.py"

    # Create output directory
    mkdir -p "$FILES_OUT"

    # Prepare metadata
    python3 "$CODE_1" --meta "$META_IN" --output "$META_OUT" --n_groups 22 ## Maximim 22 groups, as there are only 22 groups with 2 males

    # Run preparation script
    python3 "${DIR_CODE}prepare_beta_matrices.py" \
        --beta_files "${DIR_DATA}betaFiles/GSM*.hg38.beta" \
        --cpg_bed "${DIR_DATA}wgbs_tools/references/hg38/CpG.bed.gz" \
        --output_folder "$FILES_OUT" \
        --meta "$META_OUT" \
        --minCov 10 \
        --min_samples "$MIN_SAMPLES" \
        --min_datasets "$MIN_DATASETS"
}

## ------------------------------------------------------------------
## New analyses: Pairs by sex composition (step 5)
## Each script writes a filtered metadata with exactly 2 samples per
## Analysis group, with sex composition MM, FF, or MF respectively.
## We then run prepare_beta_matrices with --min_samples 2.
## ------------------------------------------------------------------

# 1) Pairs of 2 males
prep_atlas 15 "pairs_MM" 2 2

# 2) Pairs of 2 females
prep_atlas 16 "pairs_FF" 2 2

# 3) Pairs of 1 male + 1 female
prep_atlas 17 "pairs_MF" 2 2
