
CODEDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/"

# 1. build the hg38 blood reference once
Rscript build_reference_jaffe.R --build hg38 \
  --chain  B_MultiTissues/dataIn/hg19ToHg38.over.chain \
  --genome_pkg BSgenome.Hsapiens.UCSC.hg38 \
  --output B_MultiTissues/dataIn/jaffe_hg38_reference.tsv

# 2. targets = both C positions of each reference CpG (pools strands)
awk -F'\t' 'NR>1{print $2"\t"$3; print $2"\t"($3+1)}' \
    B_MultiTissues/dataIn/jaffe_hg38_reference.tsv | sort -u \
  | sed '1i chr\tpos' > B_MultiTissues/dataIn/jaffe_ref_targets_hg38.tsv

# 3. extract those CpGs from the blood .cov with the flexible S02
S02_extractTargetMeth.sh  B_MultiTissues/dataIn/jaffe_ref_targets_hg38.tsv \
                          blood  B_MultiTissues/dataOut/blood_jaffeRef_meth.tsv

# 4. reshape to the 'long' deconvolution input
awk 'BEGIN{FS=OFS="\t"} NR==1{print "sample","chrom","pos","n_valid","meth";next}
     {split($2,a,"_"); print $1,a[1],a[2],$4,$3}' \
     B_MultiTissues/dataOut/blood_jaffeRef_meth.tsv > B_MultiTissues/dataOut/blood_ref_long.tsv

# 5. estimate proportions
Rscript ont_cell_proportions.R --input B_MultiTissues/dataOut/blood_ref_long.tsv \
  --input_type long --scale percent \
  --reference B_MultiTissues/dataIn/jaffe_hg38_reference.tsv \
  --output B_MultiTissues/dataOut/blood_cellprops
