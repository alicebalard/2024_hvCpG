## Create a list of SNPs from the 1000 genome project

## 1. Download from https://ftp.ensembl.org/pub/release-115/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz

## 2. Create a list of SNPs in the format I use for my IIHV

zcat /home/alice/Downloads/1000GENOMES-phase_3.gvf.gz | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="SNV" {print $1"_"$4}' > /home/alice/Documents/GIT/2024_hvCpG/gitignore/all_snps_GRCh38_chr_pos.txt

## 3. calculate how many IIHV are SNPS
cd ~/Documents/GIT/2024_hvCpG/gitignore

grep -Ff <(sed '1d' topIntersect90.csv) all_snps_GRCh38_chr_pos.txt | wc -l
#122891
wc -l topIntersect90.csv
#174495 topIntersect90.csv  (NB. 1 extra line)
grep -Ff <(sed '1d' overlapLayers.csv) all_snps_GRCh38_chr_pos.txt | wc -l
#6475342
wc -l overlapLayers.csv
#17474841 overlapLayers.csv (NB. 1 extra line)

