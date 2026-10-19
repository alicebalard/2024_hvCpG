METHDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015/01Methcall"
DATAINDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataIn"
DATAOUTDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut"
SEQREPORT="$DATAINDIR/sequence_report.tsv"

#### 1. Endoderm
#TARGETS="$DATAOUTDIR/S05_endo_res_annotation.tsv"
#
#for f in "$METHDIR"/*_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz; do
#  sample=$(basename "$f" | sed 's/_R1_val_1.*//')
#
#  zcat "$f" | awk -v s="$sample" -v seqreport="$SEQREPORT" -v targets="$TARGETS" '
#    BEGIN {
#      FS = OFS = "\t"
#      # RefSeq accession (col 8) -> UCSC chromosome (col 11), skip header
#      first = 1
#      while ((getline line < seqreport) > 0) {
#        if (first) { first = 0; continue }
#        split(line, a, "\t")
#        refseq_to_ucsc[a[8]] = a[11]
#      }
#      close(seqreport)
#
#      # targets: chr <tab> pos, WITH header -> skip it
#      first_t = 1
#      while ((getline line < targets) > 0) {
#        if (first_t) { first_t = 0; continue }
#        split(line, a, "\t")
#        target[a[1] "_" a[2]] = 1
#      }
#      close(targets)
#    }
#    {
#      chr = refseq_to_ucsc[$1]
#      if (chr == "") next
#      chr_pos = chr "_" $2
#      if (chr_pos in target) {
#        cov = $5 + $6
#        print s, chr_pos, $4, cov
#      }
#    }
#  ' -
#done > "$DATAOUTDIR/twins_endo_target_meth_all.tsv"
#
#sed -i '1i sample\tchr_pos\tmeth_pct\tcoverage' "$DATAOUTDIR/twins_endo_target_meth_all.tsv"
#
### 2. Mesoderm
TARGETS="$DATAOUTDIR/S05_meso_res_annotation.tsv"

for f in "$METHDIR"/*_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz; do
  sample=$(basename "$f" | sed 's/_R1_val_1.*//')

  zcat "$f" | awk -v s="$sample" -v seqreport="$SEQREPORT" -v targets="$TARGETS" '
    BEGIN {
      FS = OFS = "\t"
      # RefSeq accession (col 8) -> UCSC chromosome (col 11), skip header
      first = 1
      while ((getline line < seqreport) > 0) {
        if (first) { first = 0; continue }
        split(line, a, "\t")
        refseq_to_ucsc[a[8]] = a[11]
      }
      close(seqreport)

      # targets: chr <tab> pos, WITH header -> skip it
      first_t = 1
      while ((getline line < targets) > 0) {
        if (first_t) { first_t = 0; continue }
        split(line, a, "\t")
        target[a[1] "_" a[2]] = 1
      }
      close(targets)
    }
    {
      chr = refseq_to_ucsc[$1]
      if (chr == "") next
      chr_pos = chr "_" $2
      if (chr_pos in target) {
        cov = $5 + $6
        print s, chr_pos, $4, cov
      }
    }
  ' -
done > "$DATAOUTDIR/twins_meso_target_meth_all.tsv"

sed -i '1i sample\tchr_pos\tmeth_pct\tcoverage' "$DATAOUTDIR/twins_meso_target_meth_all.tsv"

