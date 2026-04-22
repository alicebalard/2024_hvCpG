## Create a list of SNPs from the 1000 genome project

## 1. Download from https://ftp.ensembl.org/pub/release-115/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz

## 2. Create a list of SNPs in the format I use for my IIHV, only with MAF > 0.005

zcat 1000GENOMES-phase_3.gvf.gz   | awk -F'\t' '
    BEGIN { OFS="\t" }
    $0 !~ /^#/ && $3=="SNV" {
      maf_ok = 0
      n = split($9, a, ";")
      for (i = 1; i <= n; i++) {
        if (a[i] ~ /^(AFR|AMR|EAS|EUR|SAS)=/) {
          split(a[i], b, "=")
          split(b[2], vals, ",")
          for (j in vals) {
            if (vals[j] + 0 >= 0.005) maf_ok = 1
          }
        }
      }
      if (maf_ok) print $1
    }' > snps_maf_gt_0p5pct_chr.txt

zcat 1000GENOMES-phase_3.gvf.gz   | awk -F'\t' '
    BEGIN { OFS="\t" }
    $0 !~ /^#/ && $3=="SNV" {
      maf_ok = 0
      n = split($9, a, ";")
      for (i = 1; i <= n; i++) {
        if (a[i] ~ /^(AFR|AMR|EAS|EUR|SAS)=/) {
          split(a[i], b, "=")
          split(b[2], vals, ",")
          for (j in vals) {
            if (vals[j] + 0 >= 0.005) maf_ok = 1
          }
        }
      }
      if (maf_ok) print $4     
    }' > snps_maf_gt_0p5pct_pos.txt
