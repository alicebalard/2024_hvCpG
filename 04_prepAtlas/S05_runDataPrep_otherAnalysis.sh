#!/bin/bash
#$ -N prepMetadataSubAnalyses
#$ -S /bin/bash
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y

## Prepare data for 8 different re-analyses of the Atlas dataset

source /share/apps/source_files/python/python-3.13.0a6.source

DIR_DATA="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/"
DIR_CODE="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/"
META_IN="${DIR_CODE}SupTab1_Loyfer2023_amended.csv"

prep_atlas() {
    # Arguments
    STEP=$1            # e.g., 1 or 2
    NAME=$2            # name of the analysis
    MIN_DATASETS=$3     # CpG covered in at least MIN_DATASET

    # Derived paths
    META_OUT="${META_IN}_${STEP}_$2.csv"
    FILES_OUT="${DIR_DATA}10X_${STEP}_$2"
    CODE_1="${DIR_CODE}prepare_metadata_${STEP}_$2.py"

    # Create output directory
    mkdir -p "$FILES_OUT"

    # Prepare metadata
    python3 "$CODE_1" --meta "$META_IN" --output "$META_OUT"

    # Run preparation script
    python3 "${DIR_CODE}prepare_beta_matrices.py" \
        --beta_files "${DIR_DATA}betaFiles/GSM*.hg38.beta" \
        --cpg_bed "${DIR_DATA}wgbs_tools/references/hg38/CpG.bed.gz" \
        --output_folder "$FILES_OUT" \
        --meta "$META_OUT" \
        --minCov 10 \
        --min_samples 3 \
        --min_datasets "$MIN_DATASETS" 
}

## 1/ Endo, meso, ecto --> byDevLayer
# CODE_1="${DIR_CODE}prepare_metadata_1_byDevLayer.py"
# prep_atlas 1 "byDevLayer" 3

## 2/ Rm samples which have multiple cell types sampled --> rmMultSamples
# prep_atlas 2 "rmMultSamples" 25 ## 25 groups with >=3 samples, calculated before

## 3/ Keep only cells found in Maria datasets --> correspMariaTissues
# prep_atlas 3 "correspMariaTissues" 27

## 4/ male only --> maleOnly
# prep_atlas 4 "maleOnly" 6

## 5/ all but male only --> allButMaleOnly
# prep_atlas 5 "allButMaleOnly" 40

## 6/ female only --> femaleOnly
# prep_atlas 6 "femaleOnly" 15

## 7/ all but female only --> allButfemaleOnly
# prep_atlas 7 "allButfemaleOnly" 31

## 8/ by tissue rather than by cell type --> byTissue
# prep_atlas 8 "byTissue" 22
