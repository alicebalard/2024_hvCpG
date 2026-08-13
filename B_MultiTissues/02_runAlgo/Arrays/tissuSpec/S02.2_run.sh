#!/usr/bin/env bash
set -euo pipefail

# Config
DATA_DIR="/home/alice/arraysh5"
BATCH_SIZE=10000
RES_DIR="/home/alice/2024_hvCpG/B_MultiTissues/resultsDir_gitIgnored/Arrays/tissue"
NTHREADS=50
R_SCRIPT="/home/alice/2024_hvCpG/B_MultiTissues/02_runAlgo/Arrays/tissuSpec/S02.1_runTissueSpecificAlgo_arrays_COMPONENTS.R"

echo "**** Job started at $(date) ****"
R --vanilla --args "$DATA_DIR" "$BATCH_SIZE" "$RES_DIR" "$NTHREADS" < "$R_SCRIPT"
echo "**** Job finished at $(date) ****"
