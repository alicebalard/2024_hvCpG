###################################
##*******************************##
## STEP 6. Calculate age and sex ##
##*******************************##
###################################

## Step 6. Calculate age and sex
## Based on METHCOV, call a R script to calculate age and sex and add it to:
## AGE_TABLE_OUT="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/04SampleInfos/ageAllsamples_Pipeline.csv"
## SEX_TABLE_OUT="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/04SampleInfos/sexAllsamples_Pipeline.csv"

<<<<<<< HEAD
##METHCOV20=$(less "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/01Methcall/${INPUT##*/}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.20X.cov.gz")

=======
>>>>>>> 817221b90ff0ddae2513e5ea2ba2deb472980c30
##cd /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/03AgeSex/
##
#### Horvath methylation array file:
##HORV="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/03AgeSex/Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.bed"
##
#### Sample name:
##SAMPLE_HORV=$(echo $METHCOV | sed 's|.*/\([^_]*\)_.*|\1|')
##
#### Convert WGBS methylation call aligned to hg38 with Bismark to HorvathMammalMethylChip40:
#### print METHCOV, calculate end as start + 1, intersect with Horvath file
##less $METHCOV | cut -f1,2,3,4 | awk -F'\t' 'BEGIN{OFS="\t"} {$3=$3+1; print}' | /share/apps/genomics/bedtools-2.30.0/bin/bedtools intersect -a - -b $HORV -wb | cut -f1,2,3,4,8 > $SAMPLE_HORV.horv.bed
##
##AGE_TABLE_OUT="$CODEDIR/ageAllsamples_Pipeline.csv"
##SEX_TABLE_OUT="$CODEDIR/sexAllsamples_Pipeline.csv"
##
#### Call a R script to calculate age and fill in a table
##Rscript /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/02_calculateAgeSex/calculateAgeSex.R $SAMPLE_HORV.horv.bed ${INPUT##*/} $METHCOV $AGE_TABLE_OUT $SEX_TABLE_OUT
##
#### Make sure that duplicates are removed
##sort $AGE_TABLE_OUT| uniq > temp ; cat temp > $AGE_TABLE_OUT
##sort $SEX_TABLE_OUT| uniq > temp ; cat temp > $SEX_TABLE_OUT
