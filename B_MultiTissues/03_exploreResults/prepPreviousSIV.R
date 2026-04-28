###############################################
## Prepare previous putative ME list on hg38 ##
###############################################
source(here("B_MultiTissues/03_exploreResults/makeProbes2GenDictionary.R"))

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
HarrisSIV <- readxl::read_excel(here("B_MultiTissues/dataIn/Harris2012_1776SIV_10children450k.xls"), sheet = 3)
HarrisSIV_hg38 <- dico[match(HarrisSIV$Probe, dico$CpG), "chrpos_hg38"]
HarrisSIV_hg38 <- na.omit(HarrisSIV_hg38)
length(HarrisSIV_hg38) # 1773
rm(HarrisSIV)

###########################
## VanBaak2018_ESS_HM450 ##
VanBaakESS <- readxl::read_excel(here("B_MultiTissues/dataIn/VanBaak2018_1580ESS_450k.xlsx"), sheet = 2)
## only ESS hits
VanBaakESS <- VanBaakESS[VanBaakESS$`ESS hit`,]
VanBaakESS_hg38 <- dico[match(VanBaakESS$CG, dico$CpG), "chrpos_hg38"]
VanBaakESS_hg38 <- na.omit(VanBaakESS_hg38)
length(VanBaakESS_hg38) # 1579
rm(VanBaakESS)

###########################################
## Kessler2018_687SIVregions_2WGBS hg19! ##
KesslerSIV <- readxl::read_excel(here("B_MultiTissues/dataIn/Kessler2018_supTables.xlsx"), sheet = 2, skip = 1)

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
corSIV <- readxl::read_excel(here("B_MultiTissues/dataIn/Gunasekara2019_9926CoRSIVs_10WGBS.xls"), sheet = 3)
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
SoCCpGs <- readxl::read_excel(here("B_MultiTissues/dataIn/Silver2022_259SoC_hg19.xlsx"), sheet = 6, skip = 2)
SoCCpGs_hg38 <- na.omit(dico[match(SoCCpGs$cpg, dico$CpG), "chrpos_hg38"])
length(SoCCpGs_hg38) #259

################################################
## Prepare equivalent GRanges objects for all ##
################################################

HarrisSIV_hg38_GR <- makeGRfromMyCpGPos(HarrisSIV_hg38, "Harris SIV")  
VanBaakESS_hg38_GR <- makeGRfromMyCpGPos(VanBaakESS_hg38, "VanBaak ESS")
KesslerSIV_GRanges_hg38$set <- "Kessler SIV"
corSIV_GRanges_hg38$set <- "Gunasekara corSIV"
DerakhshanhvCpGs_hg38_GR <- makeGRfromMyCpGPos(DerakhshanhvCpGs_hg38, "Derakhshan hvCpG")

putativeME_GR <- c(DerakhshanhvCpGs_hg38_GR, HarrisSIV_hg38_GR, VanBaakESS_hg38_GR, KesslerSIV_GRanges_hg38,
                   corSIV_GRanges_hg38)

putativeME_GR$set <- factor(putativeME_GR$set, 
                            levels = c("Derakhshan hvCpG", "Harris SIV", 
                                       "VanBaak ESS", "Kessler SIV", "Gunasekara corSIV"))
putativeME_GR$genome <- "hg38"

## Save for paleo project (once)
# saveRDS(putativeME_GR, "../../../2025_paleoMethylVar/gitignore/putativeME_GR.RDS")
# 
# putativeME_GR <- readRDS("/path/to/putativeME_GR.RDS")

###############################################################################
## Prepare CpG associated with vmeQTL identified in MZ twins by Jordana Bell ##
###############################################################################

vmeQTL_hg19probes <- readxl::read_xlsx(here("B_MultiTissues/dataIn/vmeQTL_vCpG_359pair_sig_Zhang2025.xlsx"))
vmeQTL_hg38 <- na.omit(dico$chrpos_hg38[match(vmeQTL_hg19probes$vCpG, dico$CpG)]) ; rm(vmeQTL_hg19probes)

vmeQTL_hg38_GR <- makeGRfromMyCpGPos(vmeQTL_hg38, "Zhang vmeQTL")  
vmeQTL_hg38_GR$genome <- "hg38"

previousSIVprepared = "prepared"