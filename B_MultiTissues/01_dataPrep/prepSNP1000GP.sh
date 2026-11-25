## Create a list of SNPs from the 1000 genome project
## Corrected after Hiruni spotted a bug (Nov 2026)

cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore

## 1. Download from https://ftp.ensembl.org/pub/release-115/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz

## 2. Create a list of SNPs in the format I use for my IIHV, only with MAF > 0.01

zcat 1000GENOMES-phase_3.gvf.gz   | awk -F'\t' ' ##  Decompresses and uses tab separator
    BEGIN { OFS="\t" }
    $0 !~ /^#/ && $3=="SNV" {     ## skip headers
      maf_ok = 0  ## Resets the per-line flag
      n = split($9, a, ";")    ## splits the attribute column on ; with 9 columns
      for (i = 1; i <= n; i++) {    ## loops over the columns
        if (a[i] ~ /^(AFR|AMR|EAS|EUR|SAS)=/) {    ## Loops over attributes and keeps only the five 1000G super-population frequency fields
          split(a[i], b, "=")
          split(b[2], vals, ",")   ## splitting on , the value part into vals handles multi-allelic SNVs, one frequency per ALT allele.
	  for (j in vals) {
            f = vals[j] + 0    ## iterate over the array vals, which holds the per-allele frequencies for one population; +0 forces the value to a number
            maf = (f > 0.5) ? 1 - f : f       ## grammar: condition ? value_if_true : value_if_false. Computes the minor-allele frequency: if the ALT frequency f is above 0.5 the ALT is actually the major allele
            if (maf > 0.01) maf_ok = 1         # strict > to match "MAF > 1%"
          }
## Prev: pb if ALT freq > 0.99 (the ALT is the major allele freq)  
##          for (j in vals) {
##            if (vals[j] + 0 >= 0.01) maf_ok = 1   ## For each ALT-allele frequency, +0 coerces the string to a number and flags the site if any ALT frequency in any population is ≥ 0.01
##          }
        }
      }
      if (maf_ok) print "chr" $1 "_" $4  ## Emits chr{chrom}_{pos} (e.g. chr1_10505) for every passing site
    }' > snps_maf_0.01_chr_pos.txt
