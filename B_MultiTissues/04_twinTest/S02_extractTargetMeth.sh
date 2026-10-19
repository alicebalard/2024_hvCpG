#!/bin/bash
## Extract methylation at an arbitrary CpG list from a tissue's .cov files.
## Usage: S02_extractTargetMeth.sh  <targets.tsv(chr,pos,header)>  <adipose|blood>  <out.tsv>
set -euo pipefail

TARGETS="$1"; TISSUE="$2"; OUT="$3"
DATADIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015"
DATAINDIR="/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataIn"
SEQREPORT="$DATAINDIR/sequence_report.tsv"

case "$TISSUE" in
  adipose) METHDIR="$DATADIR/01Methcall/Adipose" ;;
  blood)   METHDIR="$DATADIR/01Methcall/Blood" ;;
  *) echo "TISSUE must be adipose or blood" >&2; exit 1 ;;
esac

: > "$OUT"
for f in "$METHDIR"/*_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz; do
  [ -e "$f" ] || continue
  sample=$(basename "$f" | sed 's/_R1_val_1.*//')
  zcat "$f" | awk -v s="$sample" -v seqreport="$SEQREPORT" -v targets="$TARGETS" '
    BEGIN{
      FS=OFS="\t"
      first=1; while((getline line < seqreport)>0){ if(first){first=0;continue}
        split(line,a,"\t"); refseq_to_ucsc[a[8]]=a[11] } close(seqreport)
      ft=1; while((getline line < targets)>0){ if(ft){ft=0;continue}
        split(line,a,"\t"); target[a[1]"_"a[2]]=1 } close(targets)
    }
    { chr=refseq_to_ucsc[$1]; if(chr=="") next
      cp=chr"_"$2
      if(cp in target){ cov=$5+$6; print s, cp, $4, cov } }
  ' - >> "$OUT"
done
sed -i '1i sample\tchr_pos\tmeth_pct\tcoverage' "$OUT"
echo "Wrote $OUT  ($(wc -l < "$OUT") lines)"

## to call
# the twin test (S05 candidate CpGs), per tissue
## S02_extractTargetMeth.sh  dataOut/S05_meso_res_annotation.tsv  adipose  dataOut/twins_adipose_meso_meth.tsv
## S02_extractTargetMeth.sh  dataOut/S05_meso_res_annotation.tsv  blood    dataOut/twins_blood_meso_meth.tsv
