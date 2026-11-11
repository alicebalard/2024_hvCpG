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

# ---- Config ----
CHUNK_SIZE=250000 ## what is the size of the chunk per array task
BATCH_SIZE=5000 ## how many CpGs are loaded at once
RSCRIPT="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/02_runAlgo/Atlas/S01.1_runalgov6_atlas_cscluster.R"

echo "**** Job $JOB_NAME.$SGE_TASK_ID started at $(date) ****"

# ---- pick the run variant here ----

## 2 others done, do the no filter now
DATA_SUFFIX="_noFilter"     # "" SNP-only | "_noSDASM" SNP+SD-ASM | "_noFilter" nothing removed
RES_SUBDIR="noFilter"       # "SNPrm"     | "SNP_SDASMrm"          | "noFilter"

## 3 layers individually, with higher thresholds
P0=0.80
P1=0.65
MININD=3
for ANALYSIS in  "atlas_general" "12_endo" "13_meso" "14_ecto"; do
    echo "[INFO] Running analysis: $ANALYSIS"
    Rscript $RSCRIPT $ANALYSIS $SGE_TASK_ID $CHUNK_SIZE $BATCH_SIZE $P0 $P1 $MININD "$DATA_SUFFIX" "$RES_SUBDIR"
done

#for ANALYSIS in  "atlas_general" "02_rmMultSamples" "12_endo" "13_meso" "14_ecto" "12_2_endo6gp" "13_2_meso6gp"; do
#    echo "[INFO] Running analysis: $ANALYSIS"
#    Rscript $RSCRIPT $ANALYSIS $SGE_TASK_ID $CHUNK_SIZE $BATCH_SIZE $P0 $P1 $MININD "$DATA_SUFFIX" "$RES_SUBDIR"
#done


echo "**** Job $JOB_NAME.$SGE_TASK_ID finished at $(date) ****"
