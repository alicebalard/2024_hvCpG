#!/bin/bash
# Needs the hg38 RepeatMasker .align and calcDivergenceFromAlign.pl (ships in RepeatMasker/util/)
# 1. get the alignment (verify the URL for your mirror; it is the RepeatMasker .align, not .out)
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.align.gz
gunzip -k hg38.fa.align.gz                                   # -> hg38.fa.align  (large, tens of GB)

# 2. recompute divergence with the CpG modification (Kimura 2-parameter, CpG sites down-weighted).
#    -a writes a new .align whose per-element line carries "Kimura (with divCpGMod) = X"
calcDivergenceFromAlign.pl -s hg38.divsum -a hg38.kimura.align hg38.fa.align

# 3. parse: for every alignment header (coords + repeat), emit it with the following Kimura value
## walks the alignment file; for each element it grabs chr, start, end, strand, repeat name, class/family from the header, then pairs it with the Kimura (with divCpGMod) = X line that follows, writing a 7-column TSV. Takes a few minutes (11Gb input)

awk '
  /#/ && $5 ~ /^chr/ {
    chr=$5; b=$6; e=$7; st="+"; rep="";
    for (i=8; i<=NF; i++) {
      if ($i=="C") st="-";
      if ($i ~ /#/) { rep=$i; break }
    }
    split(rep,a,"#"); name=a[1]; cf=a[2]; have=1; next
  }
  /Kimura/ { if (have) { v=$NF+0;
                         print chr"\t"b"\t"e"\t"st"\t"name"\t"cf"\t"v; have=0 } }
' hg38.kimura.align > hg38_rmsk_kimura.tsv

## 3. Sanity check
wc -l hg38_rmsk_kimura.tsv                  # ~ millions of rows (one per repeat element)
head hg38_rmsk_kimura.tsv                   # chr  start  end  strand  name  class/family  kimura
# divergence should span ~0 (young) to ~40+ (old); quick look at the range:
awk -F'\t' '{print $7}' hg38_rmsk_kimura.tsv | sort -n | sed -n '1p;$p'

## move to local gitignore for R
