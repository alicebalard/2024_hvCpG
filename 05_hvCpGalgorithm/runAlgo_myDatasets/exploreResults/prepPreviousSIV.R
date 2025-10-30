###############################################
## Prepare previous putative ME list on hg38 ##
###############################################
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/makeProbes2GenDictionary.R"))

message("Creates:
        \na vector of 1773 SIV from Harris 2012 (HarrisSIV_hg38)
        \none of 1579 ESS from Van Baak 2018 (VanBaakESS_hg38)
        \na GRange object for Kessler 2018 676 SIV regions (KesslerSIV_GRanges_hg38)
        \na GRange object for Gunasekara 2019 9926 corSIV regions (corSIV_GRanges_hg38)
        \na vector of 3644 hvCpG from Derakhshan 2022 (DerakhshanhvCpGs_hg38)
        \na vector for matching mQTL controls (mQTLcontrols_hg38)
        \na vector for 259 Silver 2022 SoCCpGs on 10WGBS (SoCCpGs_hg38)")

#######################################
## Harris2012_1776SIV_10children450k ##
HarrisSIV <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Harris2012_1776SIV_10children450k.xls"), sheet = 3)
HarrisSIV_hg38 <- dico[match(HarrisSIV$Probe, dico$CpG), "chrpos_hg38"]
HarrisSIV_hg38 <- na.omit(HarrisSIV_hg38)
length(HarrisSIV_hg38) # 1773
rm(HarrisSIV)

###########################
## VanBaak2018_ESS_HM450 ##
VanBaakESS <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/VanBaak2018_1580ESS_450k.xlsx"), sheet = 2)
## only ESS hits
VanBaakESS <- VanBaakESS[VanBaakESS$`ESS hit`,]
VanBaakESS_hg38 <- dico[match(VanBaakESS$CG, dico$CpG), "chrpos_hg38"]
VanBaakESS_hg38 <- na.omit(VanBaakESS_hg38)
length(VanBaakESS_hg38) # 1579
rm(VanBaakESS)

###########################################
## Kessler2018_687SIVregions_2WGBS hg19! ##
KesslerSIV <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Kessler2018_supTables.xlsx"), sheet = 2, skip = 1)

KesslerSIV_GRanges <- GRanges(
  seqnames = KesslerSIV$Chromosome,
  ranges = IRanges(start = KesslerSIV$`ME start`, 
                   end = KesslerSIV$`ME end`))

## liftover to hg38, keep uniquely mapping regions
mapped <- liftOver(KesslerSIV_GRanges, chain)
keep <- lengths(mapped) == 1
KesslerSIV_GRanges_hg38 <- unlist(mapped[keep])
length(KesslerSIV_GRanges_hg38)
rm(mapped, keep, KesslerSIV, KesslerSIV_GRanges)

#######################################
## Gunasekara2019_9926CoRSIVs_10WGBS ##
# Load corSIV intervals (already in hg38)
corSIV <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Gunasekara2019_9926CoRSIVs_10WGBS.xls"), sheet = 3)
corSIV <- unique(corSIV$USCS_Coordinates_CoRSIV)
corSIV_split <- tstrsplit(corSIV, "[:-]", fixed = FALSE)

corSIV_GRanges_hg38 <- GRanges(
  seqnames = corSIV_split[[1]],
  ranges = IRanges(start = as.integer(corSIV_split[[2]]), end = as.integer(corSIV_split[[3]])),
  strand = "*")

length(corSIV_GRanges_hg38)
rm(corSIV, corSIV_split)

#######################################
## Derakhshan 2022 (previous hvCpGs) ##

data <- read.table(here("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)
DerakhshanhvCpGs_hg38 <- dico[match(data$hvCpG_name, dico$CpG), "chrpos_hg38"]
mQTLcontrols_hg38 <- dico[match(data$controlCpG_name, dico$CpG), "chrpos_hg38"]

rm(data)

###############################
## Silver2022_SoCCpGs_10WGBS ##
SoCCpGs <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Silver2022_259SoC_hg19.xlsx"), sheet = 6, skip = 2)
SoCCpGs_hg38 <- na.omit(dico[match(SoCCpGs$cpg, dico$CpG), "chrpos_hg38"])
length(SoCCpGs_hg38) #259
