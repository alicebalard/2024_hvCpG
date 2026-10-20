## NB. composition is the main bulk-tissue confounder for a twin concordance test.
# The twin samples are bulk tissue, so each measured methylation value is a cell-composition-weighted average of the underlying cell-type methylomes.
# If two twins of a pair differ in cell composition, their bulk methylation differs for reasons unrelated to the epiallele/genetic effect we’re testing, and composition itself is partly heritable and shared-environment-driven, so MZ pairs may be more alike in composition than DZ pairs. That would make MZ within-pair differences smaller than DZ even with no cis-genetic effect on methylation, i.e. it can mimic exactly the signal our MZ<DZ test looks for. Estimating composition lets us either adjust it out or show it isn’t driving the result.

CODEDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/"

# 1. build the hg38 blood reference once

Rscript $CODEDIR/04_twinTest/build_reference_jaffe.R --build hg38 --chain  $CODEDIR/dataIn/hg19ToHg38.over.chain --genome_pkg BSgenome.Hsapiens.UCSC.hg38 --output $CODEDIR/dataIn/jaffe_hg38_reference.tsv

#Jaffe coefficients: 600 CpGs, cell types: CD8T, CD4T, NK, Bcell, Mono, Gran
#Mapped uniquely to hg38: 600 of 600
#Same orientation after conversion: 600 of 600
#Annotated position is the C of a CpG: 599; the G of a CpG: 0; not a CpG in this genome (dropped): 1
#Wrote 599 CpGs to /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues//dataIn/jaffe_hg38_reference.tsv

## TBC from here

# 2. targets = both C positions of each reference CpG (pools strands)
awk -F'\t' 'NR>1{print $2"\t"$3; print $2"\t"($3+1)}' \
    $CODEDIR/dataIn/jaffe_hg38_reference.tsv | sort -u \
  | sed '1i chr\tpos' > $CODEDIR/dataIn/jaffe_ref_targets_hg38.tsv

# 3. extract those CpGs from the blood .cov with the flexible S02
S02_extractTargetMeth.sh $CODEDIR/dataIn/jaffe_ref_targets_hg38.tsv \
                          blood $CODEDIR/dataOut/blood_jaffeRef_meth.tsv

# 4. reshape to the 'long' deconvolution input
awk 'BEGIN{FS=OFS="\t"} NR==1{print "sample","chrom","pos","n_valid","meth";next}
     {split($2,a,"_"); print $1,a[1],a[2],$4,$3}' \
     $CODEDIR/dataOut/blood_jaffeRef_meth.tsv > $CODEDIR/dataOut/blood_ref_long.tsv

# 5. estimate proportions
Rscript $CODEDIR/04_twinTest/ont_cell_proportions.R --input $CODEDIR/dataOut/blood_ref_long.tsv \
  --input_type long --scale percent \
  --reference $CODEDIR/dataIn/jaffe_hg38_reference.tsv \
  --output $CODEDIR/dataOut/blood_cellprops
