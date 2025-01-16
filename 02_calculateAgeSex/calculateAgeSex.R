##############
## A.Balard ##
## Jan 2025 ##

#########################################################################
## part 1: Age calculation based on universal pan-tissue mammanlian clock
# https://github.com/jazoller96/mammalian-methyl-clocks/tree/main

## NB: MammalMethylClock is installed in R v4.2.1 (issue on hpc for 4.4.1)

## Convert WGBS methylation call aligned to hg38 with Bismark to HorvathMammalMethylChip40
download.file(
  "https://raw.githubusercontent.com/shorvath/MammalianMethylationConsortium/refs/tags/v1.0.0/Annotations%2C%20Amin%20Haghani/Mammals/Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv",
  "Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv"
)
Homo_sapiens.hg38.HorvathMammalMethylChip40.v1 <- read.csv("Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv")
file.remove("Homo_sapiens.hg38.HorvathMammalMethylChip40.v1.csv")

gr1 <- makeGRangesFromDataFrame(Homo_sapiens.hg38.HorvathMammalMethylChip40.v1, 
                                keep.extra.columns=TRUE, na.rm=TRUE)

df2 <- read.table("Downloads/madeupWGBS.txt")
names(df2) <- c("chr", "Cpos", "strand", "context", "methPerc", "methCount", "totalCount")

## Horvath chips report C and G position, WGBS only C
df2$start <- NA; df2$end <- NA
df2$start[df2$strand %in% "+"] <- df2$Cpos[df2$strand %in% "+"]
df2$end[df2$strand %in% "+"] <- df2$Cpos[df2$strand %in% "+"] +1
df2$start[df2$strand %in% "-"] <- df2$Cpos[df2$strand %in% "-"] -1
df2$end[df2$strand %in% "-"] <- df2$Cpos[df2$strand %in% "-"]

gr2 <- makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE)

intersection <- plyranges::join_overlap_inner(gr1, gr2)

dfAge <- data.frame(CGid=intersection$CGid, methPerc=intersection$methPerc)

## For only one sample. Code for all samples.

## browseVignettes("MammalMethylClock")
dfAgeWide <- dfAge %>%
  pivot_wider(names_from = CGid, values_from = methPerc)%>%
  mutate(SID="sampleID", tissue = "tissue") %>% 
  data.frame


ys.output <- predictAge(dfAgeWide[!names(dfAgeWide) %in% c("SID", "tissue")],
                        dfAgeWide[names(dfAgeWide) %in% c("SID", "tissue")], 
                        tissue.names = c("tissue"), species.name = "Homo sapiens")
head(ys.output)

###########################
## part 2: Sex predictor ##
## https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07675-2
## 4345 sex-associated sites

library("usethis")
install.packages('BiocManager')
BiocManager::install('wateRmelon')
library(wateRmelon)


