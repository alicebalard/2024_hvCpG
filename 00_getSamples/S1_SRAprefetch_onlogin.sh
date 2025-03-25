#!/bin/bash                                                                                                                                         

## To run on "large" and not in a job

echo "Running on host: $(hostname)"
echo "Start time: $(date)"

DIR="/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/00RawFastq/dataset2"

cd $DIR

###########################
## DL files from a SRA list

## dataset 1: 20 samples --> done
## tail -n +2 /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/datasetInfo/SraRunTable_dataset1.csv | cut -f 1 > $DIR/temp.samplelist

## dataset 2: 40 samples
tail -n +2 /SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/datasetInfo/SraRunTable_dataset2.csv | cut -f 1 > $DIR/temp.samplelist

##########################
SAMPLELIST=temp.samplelist

echo "Samples:"
cat $SAMPLELIST

echo "prefetch the files:"
/share/apps/genomics/sratoolkit.3.0.2/bin/prefetch --max-size 30GB -vv -O $DIR --option-file $SAMPLELIST
