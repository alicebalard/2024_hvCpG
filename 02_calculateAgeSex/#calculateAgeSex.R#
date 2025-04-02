##############
## A.Balard ##
## Jan 2025 ##

### Retrieve arguments from the bash calling file:
#args <- commandArgs(trailingOnly=TRUE)
#sample_horvbed <- read.table(args[1])
#sample_name <- args[2]
#methcov <- args[3]
#filepath_age <- args[4]
#filepath_sex <- args[5]
#
## Before installing R packages on the CS UCL cluster, run:
## scl enable devtoolset-9 bash
## This uses the good gcc
# library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

###################################
## Run on a list of files (tbc)

## Select only coverage >=20 (to tweak)
WGBS <-read.table(gzfile("/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/01Methcall/SRR28532098_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.20X.cov.gz"))
names(WGBS) <- c("chr", "start", "stop", "percMeth", "mC", "C")

## Convert to beta values
WGBS$beta_values <- WGBS$percMeth / 100

# Map genomic positions to mammalian array probes
# (Requires annotation file from Mammalian Methylation Consortium)
annot <- read.csv("MammalianMethylationArrayAnnotation.csv")


## Conversion WGBS methylation from hg38 to hg19:

## Conversion WGBS into 450k array:





library(tidyverse)
library(dplyr)
library(MammalMethylClock)
## NB: MammalMethylClock is installed in R v4.2.1 (issue on hpc for 4.4.1)
library(wateRmelon)
library(methylKit)
library(GenomicRanges)
## devtools::install_github(repo = "tiagochst/ELMER.data")
library("ELMER.data")
library(minfi)

#########################################################################
## part 1: Age calculation based on universal pan-tissue mammanlian clock
                                        # https://github.com/jazoller96/mammalian-methyl-clocks/tree/main

## 1. prepare the reference file 
refNeeded=F
if (refNeeded){
    download.file("https://raw.githubusercontent.com/shorvath/MammalianMethylationConsortium/refs/tags/v1.0.0/Annotations%2C%20Amin%20Haghani/Mammals/Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv", "Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv")
    Homo_sapiens.hg38.HorvathMammalMethylChip40.v1 <- read.csv("Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv")
    file.remove("Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv")

    ## rm NAs
    Homo_sapiens.hg38.HorvathMammalMethylChip40.v1=na.omit(Homo_sapiens.hg38.HorvathMammalMethylChip40.v1[c("seqnames", "start", "end", "CGid", "strand")])

    ## keep only one position per CGid
    Homo_sapiens.hg38.HorvathMammalMethylChip40.v1 = Homo_sapiens.hg38.HorvathMammalMethylChip40.v1 %>% dplyr::distinct(seqnames, start, end, .keep_all = TRUE)

    ## rm the final ".1" of some CGid
    Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$CGid = gsub("\\.1$", "", Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$CGid)

    ## check: ok
    Homo_sapiens.hg38.HorvathMammalMethylChip40.v1[grep("\\.", Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$CGid),]
    df1 <- data.frame(seqnames = Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$seqnames)
    df1$start = Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$start

    ## If strand "-", the C is the second base 
    df1$start[Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$strand %in% "-"] = Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$end[Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$strand %in% "-"]
    df1$end = df1$start + 1
    df1$CGid = Homo_sapiens.hg38.HorvathMammalMethylChip40.v1$CGid

    ## Give correct naming for chromosomes hg38
    df2 <- read.table("Homo_sapiens.hg38.chromNames.tsv", header = T)
    df2$seqnames <- paste0("chr", df2$Chromosome_name)
    df1$seqnames <- df2$RefSeq_accession[match(df1$seqnames, df2$seqnames)]

    write.table(unique(na.omit(df1)), "/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/03AgeSex/Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.bed", sep = "\t", quote = F, row.names = F, col.names = F)
}

read.table("/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/03AgeSex/Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.bed")


## 2. calculate age 
dfAge <- data.frame(CGid=sample_horvbed$V5,
                    methPerc=sample_horvbed$V4/100)

dfAgeWide <- dfAge %>%
    pivot_wider(names_from = CGid, values_from = methPerc)%>%
    mutate(SID=sample_name, tissue = "tissue") %>% 
    data.frame

ys.output <- predictAge(
    dfAgeWide[!names(dfAgeWide) %in% c("SID", "tissue")],
    dfAgeWide[names(dfAgeWide) %in% c("SID", "tissue")], 
    tissue.names = c("tissue"), species.name = "Homo sapiens")

print("age calculated")

## Save the new data
if (file.exists(filepath_age)) {
    x = read.csv(filepath_age)
    x = rbind(x, ys.output)
    x = unique(x) # keep unique rows
    write.csv(x, filepath_age, quote = F, row.names=F)
} else {
    write.csv(ys.output, filepath_age, quote = F, row.names=F)   
}

print("age saved")

###########################
## part 2: Sex predictor ##
## https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07675-2
## 4345 sex-associated sites

## Load Illumina450K array manifest for hg38
data("hm450.hg38.manifest")

## prepare a hm450.hg38.manifest at the same format than bismark outdir

## issue with chrM (6 probes)
hm450.hg38.manifest = hm450.hg38.manifest[!hm450.hg38.manifest@seqnames %in% "chrM"]

## If strand "-", the C is the second base
hmneg = hm450.hg38.manifest[hm450.hg38.manifest@strand %in% "-"]
hmpos = hm450.hg38.manifest[hm450.hg38.manifest@strand %in% "+"]

pos <-c(hmneg@ranges@start-1, hmpos@ranges@start)

## Give correct naming for chromosomes hg38
df2 <- read.table("/SAN/ghlab/pophistory/Alice/hvCpG_project/data/WGBS_human/03AgeSex/Homo_sapiens.hg38.chromNames.tsv", header = T)
df2$seqnames <- paste0("chr", df2$Chromosome_name)
seqnames <- c(hmneg@seqnames, hmpos@seqnames)
seqnames <- df2$RefSeq_accession[match(as.character(seqnames), df2$seqnames)]

## Create the GRanges object
hm450_gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = pos, end = pos), strand = "*")
names(hm450_gr) <- c(names(hmneg), names(hmpos))

## Read in coverage file
bismark_cov <- read.table(methcov, header=FALSE, 
                          col.names=c("chr", "start", "end", "methylation", "meth_reads", "unmeth_reads"))

## Create a GRange object
bismark_gr <- GRanges(seqnames=bismark_cov$chr, 
                      ranges=IRanges(start=bismark_cov$start, end=bismark_cov$end),
                      methylation=bismark_cov$methylation,
                      coverage=bismark_cov$meth_reads + bismark_cov$unmeth_reads)

## Find overlaps between WGBS data and hm450 array positions:
overlaps <- findOverlaps(hm450_gr, bismark_gr)

## Extract methylation values for 450K array positions
methylation_450k <- bismark_gr$methylation[subjectHits(overlaps)]
coverage_450k <- bismark_gr$coverage[subjectHits(overlaps)]

## Create a MethylSet object
meth_set <- MethylSet(Meth = matrix(methylation_450k * coverage_450k / 100, ncol=1),
                      Unmeth = matrix((100 - methylation_450k) * coverage_450k / 100, ncol=1),
                      colData = DataFrame(Sample = "WGBS_sample"),
                      rowData = hm450_gr[queryHits(overlaps),],
                      annotation = c(array = "IlluminaHumanMethylation450k",
                                     annotation = "ilmn12.hg38"))

estimateSex <- function(betas, do_plot=FALSE){
    betas <- as.matrix(betas)
    single_sample <- FALSE
    if(ncol(betas) == 1) {
        betas <- cbind(betas, betas)
        single_sample <- TRUE
    }
    if (length(grep('_', head(rownames(betas), n = 10L)))==10){
        betas <- epicv2clean(betas)
    }
                                        # predict sex by two PCAs on X and Y chromosomes
    data("sexCoef")
                                        # Z score normalization
    betas <- betas[rownames(betas) %in% sex_coef$IlmnID, ]
    message('Normalize beta values by Z score...')
    autosomes <- sex_coef$IlmnID[!(sex_coef$CHR %in% c('X', 'Y'))]
    auto_betas <- betas[rownames(betas) %in% autosomes, ]
    d_mean <- colMeans(auto_betas, na.rm=TRUE)
    d_sd <- colSds(auto_betas, na.rm=TRUE, useNames=T) ## updated here
    z_beta <- (t(betas) - d_mean) / d_sd
    message('Fishished Zscore normalization.')

                                        # Sex prediction
    pred_XY <- list()
    for(chr in c('X', 'Y')){
        coefs <- sex_coef[sex_coef$pca == chr,]
        miss_probes <- setdiff(coefs$IlmnID, colnames(z_beta))
        if(length(miss_probes) > 0){
            warning('Missing ', length(miss_probes), ' probes!\n', paste(c(miss_probes), collapse=", "))  
            coefs <- coefs[!(coefs$IlmnID %in% miss_probes), ]
        }
        chr_beta <- z_beta[, coefs$IlmnID]
        chr_beta[is.na(chr_beta)] <- 0
        pred_chr <- t(t(chr_beta) - coefs$mean) %*% coefs$coeff
        pred_XY[[chr]] <- pred_chr
    }
    pred_XY <- data.frame(pred_XY)

    pred_XY$'predicted_sex' <- 'Female'
    pred_XY$'predicted_sex'[(pred_XY$X < 0) & (pred_XY$Y > 0)] <- 'Male'
    pred_XY$'predicted_sex'[(pred_XY$X > 0) & (pred_XY$Y > 0)] <- '47,XXY'
    pred_XY$'predicted_sex'[(pred_XY$X < 0) & (pred_XY$Y < 0)] <- '45,XO'
    if(single_sample){
        pred_XY <- pred_XY[1, ]
    }
    
    if(do_plot){
        plot_predicted_sex(pred_XY)
    }else{
        message('You can visualize the predicted results by set "do_plot=TRUE".\n')
    }
    return(pred_XY)
}

est <- estimateSex(betas(meth_set))

print("sex calculated")

## Save the new data
y = data.frame(sample_name=sample_name, predicted_sex=est$predicted_sex)

if (file.exists(filepath_sex)) {
    x = read.csv(filepath_sex)
    x = rbind(x, y)
    x = unique(x)
    write.csv(x, filepath_sex, quote = F, row.names=F)
} else {
    write.csv(y, filepath_sex, quote = F, row.names=F)
}

print("sex saved")
