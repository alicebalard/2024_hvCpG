#!/bin/bash
#$ -N runalgo_atlas
#$ -S /bin/bash
#$ -pe smp 12
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=10:00:00
#$ -wd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/logs
#$ -R y
#$ -t 1-93
#$ -tc 30

## Timing tests:
##batch size time in min time in hour time in days Days for 23036026 CpGs Size of chunks N chunks  Time per chunk in day  n threads  n ram (G)
##1000  	20	     0.3	0.013	     29.9468338	            250000	 92.144104	0.325	              10	5
##1000  	40	     0.7	0.029	     66.8044754	            250000	 92.144104	0.725	              5	        10
##1000  	13	     0.2	0.008	     18.4288208	            250000	 92.144104	0.2	              20	5
##1000  	12	     0.2	0.008	     18.4288208	            250000	 92.144104	0.2	              15	5
##1000  	12	     0.2	0.008	     18.4288208	            250000	 92.144104	0.2	              16	5
##2000  	30	     0.5	0.021	     24.1878273	            250000	 92.144104	0.2625	              15	5
##5000  	7	     0.1	0.004	     18.4288208	            250000	 92.144104	0.2	              5	5
## I chose: 15 threads, 5G each, 93 chunks of 250k CpGs which should each run for 5h.
## The code ~/monitor_resources.sh 5 & was used to check that all resources were used
## 15 doesn't run, let's try 12

CHUNK_SIZE=250000 ## How big are chunks sent to arrays?
BATCH_SIZE=10000 ## How many CpGs are loaded at the same time?

echo "**** Job $JOB_NAME.$SGE_TASK_ID started at $(date) ****"

DATA_DIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/"

Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/S02.1_runalgov6_atlas_cscluster.R "$DATA_DIR" "$SGE_TASK_ID" "$CHUNK_SIZE" "$BATCH_SIZE"

echo "**** Job $JOB_NAME.$SGE_TASK_ID finished at $(date) ****"
