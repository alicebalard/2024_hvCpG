#!/bin/bash
#$ -N runalgo_atlas
#$ -S /bin/bash
#$ -pe smp 5 
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=50:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-120 ## to accomodate also tests with more CpGs covered
#$ -tc 30

# ---- Config ----
CHUNK_SIZE=250000 ## what is the size of the chunk per array task
BATCH_SIZE=5000 ## how many CpGs are loaded at once
RSCRIPT="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/02_runAlgo/Atlas/S01.1_runalgov6_atlas_cscluster.R"

echo "**** Job $JOB_NAME.$SGE_TASK_ID started at $(date) ****"

##P0=0.80
##P1=0.65
##MININD=3
##
##for ANALYSIS in "atlas_general" "02_rmMultSamples" "04_maleOnly" "05_femaleOnly6gp" "06_bothsexes6gp" "09_immuneOnly" "10_noImmune" "11_noImmune_sample11gp" "12_endo" "12_2_endo6gp" "13_meso" "13_2_meso6gp" "14_ecto" "18_mesoEndo" "19_endoEcto" "20_mesoEcto"; do
##    echo "[INFO] Running analysis: $ANALYSIS"
##    Rscript $RSCRIPT $ANALYSIS $SGE_TASK_ID $CHUNK_SIZE $BATCH_SIZE $P0 $P1 $MININD
##done

##P0=0.80
##P1=0.9
##MININD=3
##ANALYSIS="atlas_general"
##
##echo "[INFO] Running analysis: $ANALYSIS"
##Rscript $RSCRIPT $ANALYSIS $SGE_TASK_ID $CHUNK_SIZE $BATCH_SIZE $P0 $P1 $MININD

P0=0.55
P1=0.65
MININD=3
ANALYSIS="atlas_general"
echo "[INFO] Running analysis: $ANALYSIS"
Rscript $RSCRIPT $ANALYSIS $SGE_TASK_ID $CHUNK_SIZE $BATCH_SIZE $P0 $P1 $MININD

## to do: pairs (if needed)
##MININD=2 # very few, for specific analysis
##for ANALYSIS in "10X_15_pairs_MM" "10X_16_pairs_FF" "10X_17_pairs_MF"; do
##    echo "[INFO] Running analysis: $ANALYSIS"
##    Rscript $RSCRIPT $ANALYSIS $SGE_TASK_ID $CHUNK_SIZE $BATCH_SIZE $P0 $P1 $MININD
##done

echo "**** Job $JOB_NAME.$SGE_TASK_ID finished at $(date) ****"
