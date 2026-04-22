## Create a list of SNPs from the 1000 genome project

## 1. Download from https://ftp.ensembl.org/pub/release-115/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz

## 2. Create a list of SNPs in the format I use for my IIHV

zcat /home/alice/Documents/GIT/2024_hvCpG/gitignore/1000GENOMES-phase_3.gvf.gz | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="SNV" {print $1}' > /home/alice/Documents/GIT/2024_hvCpG/gitignore/all_snps_GRCh38_chr.txt
zcat /home/alice/Documents/GIT/2024_hvCpG/gitignore/1000GENOMES-phase_3.gvf.gz | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="SNV" {print $4}' > /home/alice/Documents/GIT/2024_hvCpG/gitignore/all_snps_GRCh38_pos.txt
