###############################################
## Prepare previous putative ME list on hg38 ##
###############################################
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/makeEpicDictionary.R"))

#######################################
## Harris2012_1776SIV_10children450k ##
HarrisSIV <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Harris2012_1776SIV_10children450k.xls"), sheet = 3)
HarrisSIV_hg38 <- dico[match(HarrisSIV$Probe, dico$CpG), "chrpos_hg38"]
HarrisSIV_hg38 <- na.omit(HarrisSIV_hg38)
length(HarrisSIV_hg38) # 1773

###########################
## VanBaak2018_ESS_HM450 ##
VanBaakESS <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/VanBaak2018_1580ESS_450k.xlsx"), sheet = 2)
## only ESS hits
VanBaakESS <- VanBaakESS[VanBaakESS$`ESS hit`,]
VanBaakESS_hg38 <- dico[match(VanBaakESS$CG, dico$CpG), "chrpos_hg38"]
VanBaakESS_hg38 <- na.omit(VanBaakESS_hg38)
length(VanBaakESS_hg38) # 1579

###########################################
## Kessler2018_687SIVregions_2WGBS hg19! ##
KesslerSIV <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Kessler2018_supTables.xlsx"), sheet = 2, skip = 1)

KesslerSIV_GRanges <- GRanges(
  seqnames = KesslerSIV$Chromosome,
  ranges = IRanges(start = KesslerSIV$`ME start`, 
                   end = KesslerSIV$`ME end`),
  strand = "*")




## liftover to hg38, keep uniquely mapping regions
mapped <- liftOver(KesslerSIV_GRanges, hvCpGandControls$chain)
keep <- lengths(mapped) == 1
hg38_unique <- unlist(mapped[keep])

## find the match with Atlas cpg
cpg_46 <- read.table("~/Documents/Project_hvCpG/selected_cpgs_min3_in46_datasets.txt")$V1
# Parse with regex
parsed <- str_match(cpg_46, "(chr[0-9XYM]+)_(\\d+)-(\\d+)")
# Build GRanges
cpg_46_GR <- GRanges(
  seqnames = parsed[,2],
  ranges   = IRanges(start = as.numeric(parsed[,3]),
                     end   = as.numeric(parsed[,4]))
)

overlaps <- findOverlaps(query = KesslerSIV_GRanges, subject = cpg_46_GR)

# Extract the overlapping ranges
CpG_overlapping     <- cpg_46_GR[subjectHits(overlaps)]

KesslerSIV_hg38 <- paste0(CpG_overlapping@seqnames, "_", CpG_overlapping@ranges)
KesslerSIV_hg38 <- na.omit(KesslerSIV_hg38)
length(KesslerSIV_hg38) # 819

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

## find the match with Atlas cpg
overlaps <- findOverlaps(query = corSIV_GRanges_hg38, subject = cpg_46_GR)
CpG_overlapping     <- cpg_46_GR[subjectHits(overlaps)]
corSIV_hg38 <- paste0(CpG_overlapping@seqnames, "_", CpG_overlapping@ranges)

corSIV_hg38 <- na.omit(corSIV_hg38)
length(corSIV_hg38) # 70352

