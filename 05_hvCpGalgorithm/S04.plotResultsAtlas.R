## July 2025
## A. Balard
## Plot files for scan Atlas

packages <- c("dplyr", "data.table", "readxl", "progress", "ggplot2",
              "tidyr", "scales", "viridis", "ggrastr", "Cairo",
              "genomation", "GenomicRanges", "GenomicFeatures")

lapply(packages, library, character.only = TRUE)

## This code does:
### Histogram of coverage across datasets
### Manhattan plot
### Test enrichment of features for high alpha
### Test for enrichment in other putative MEs for hvCpGs with alpha > threshold

###########################################
## Histogram of coverage across datasets ##
###########################################

t5 <- read.table("../04_prepAtlas/CpG_coverage_freqtable5X.tsv", header = T)
t10 <- read.table("../04_prepAtlas/CpG_coverage_freqtable10X.tsv", header = T)

# Add coverage type label
t5$coverage <- "≥5"
t10$coverage <- "≥10"

# Combine into one data.table
t_combined <- rbind(t5, t10)

# Filter out CpGs with zero dataset coverage if needed
t_combined <- t_combined[t_combined$datasets_covered_in > 0, ]

t_combined <- t_combined %>% group_by(coverage) %>%
  dplyr::arrange(-dplyr::row_number(datasets_covered_in)) %>%
  mutate(nCpGcum = cumsum(num_CpGs))

options(scipen=0)
# Plot
pdf("figures/freqCpGperdataset.pdf", width = 14, height = 4)
ggplot(t_combined, aes(x = as.factor(datasets_covered_in), y = nCpGcum, fill = coverage)) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    breaks = seq(0, 100000000, by = 10000000),  # 10 million steps
    labels = label_number(scale = 1e-6, suffix = "M")
  ) +  scale_fill_manual(values = c("≥5" = "steelblue", "≥10" = "firebrick")) +
  labs(title = "Nbr of CpG covered across X or more datasets",
       x = ">= X datasets",
       y = "Number of CpGs with ≥3 samples covered",
       fill = "Coverage threshold") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(.2,.5),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "black", size = 0.4),
        legend.key = element_rect(fill = "white", color = NA))
dev.off()

####################
## Manhattan plot ##
####################

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- "resultsDir/Atlas10X"

# Get list of relevant RData files
rdata_files <- dir(parent_dir, pattern = "results_Atlas10X_100000CpGs_0_8p0_0_65p1\\.RData$", 
                   recursive = TRUE, full.names = TRUE)

## Check if all batches have ran
length(rdata_files) ## 274

all_cpg_values <- numeric()
pb <- progress_bar$new(total = length(rdata_files), format = "📦 :current/:total [:bar] :percent")

for (file in rdata_files) {
  e <- new.env()
  load(file, envir = e)
  obj <- e[[ls(e)[1]]]
  if (is.matrix(obj)) obj <- obj[, 1]
  all_cpg_values <- c(all_cpg_values, obj)
  pb$tick()
}

# Create data.table from named vector
dt <- data.table(
  name = names(all_cpg_values),
  alpha = as.numeric(all_cpg_values)
)

# Parse "chr_start-end" in name into chr, start_pos, end_pos
dt[, c("chr", "start_end") := tstrsplit(name, "_", fixed = TRUE)]
dt[, c("start_pos", "end_pos") := tstrsplit(start_end, "-", fixed = TRUE)]

# Convert to integer/numeric if not already
dt[, start_pos := as.integer(start_pos)]
dt[, end_pos := as.integer(end_pos)]

## Check chromosomes present:
message("Chromosomes in the dataset:")
table(unique(dt$chr))

# Convert chr from "chrN" to  factor
dt[, chr := sub("chr", "", chr)]
dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y")))]

## Check chromosomes order:
message("Chromosomes in the dataset:")
table(unique(dt$chr))

# Compute cumulative position offsets for Manhattan plot
setorder(dt, chr, start_pos)

offsets <- dt[, .(max_pos = max(start_pos, na.rm = TRUE)), by = chr]
offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]

dt <- merge(dt, offsets[, .(chr, cum_offset)], by = "chr", all.x = TRUE, sort = FALSE)

# Convert to integer/numeric if not already
dt[, cum_offset := as.numeric(cum_offset)]
dt[, pos2 := start_pos + cum_offset]

# Compute chromosome centers for x-axis labeling
df2 <- dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y"), levels=as.character(c(1:22, "X", "Y")))),
             df2, by = "chr", all.x = TRUE, sort = TRUE)
df2 <- na.omit(df2)

plot <- ggplot() +
  geom_point_rast(data = dt, aes(x = pos2, y = alpha),
                  color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr)) +
  xlab("Chromosome") +
  ylab("pr(Z=1) of being a hvCpG") +
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
# CairoPDF("figures/ManhattanAlphaPlot.pdf", width = 15, height = 3)
# print(plot)
# dev.off()

################################################
## Test enrichment of features for high alpha ##
################################################

#############
threshold=0.7
#############

# Filter valid rows
dt_clean <- dt[!is.na(start_pos) & !is.na(end_pos)]

# Create GRanges
gr_cpg <- GRanges(
  seqnames = paste0("chr", dt_clean$chr),
  ranges = IRanges(start = dt_clean$start_pos, end = dt_clean$end_pos),
  alpha = dt_clean$alpha
)

doAnnot = TRUE

if (doAnnot){
  
  # Import bed file
  bed_features <- readTranscriptFeatures("../../../Project_hvCpG/hg38_GENCODE_V47.bed")
  
  # annotate CpGs with high alpha (>threshold):
  anno_result_highalpha <- annotateWithGeneParts(
    target = gr_cpg[!is.na(mcols(gr_cpg)$alpha) & mcols(gr_cpg)$alpha > threshold],
    feature = bed_features)
  
  # annotate CpGs with low alpha (<=threshold):
  anno_result_lowalpha <- annotateWithGeneParts(
    target = gr_cpg[!is.na(mcols(gr_cpg)$alpha) & mcols(gr_cpg)$alpha <= threshold],
    feature = bed_features)
  
  anno_result_lowalpha@precedence
  anno_result_highalpha@precedence
  
  # Percentages from annotations
  low_anno <- anno_result_lowalpha@precedence
  high_anno <- anno_result_highalpha@precedence
  ## Convert percentages to counts
  N_low <- nrow(dt_clean[alpha <= threshold])
  N_high <- nrow(dt_clean[alpha > threshold])
  
  # Reconstruct counts from percentages
  count_low <- round(low_anno / 100 * N_low)
  count_high <- round(high_anno / 100 * N_high)
  
  ## Build contingency table
  contingency <- rbind(
    LowAlpha = count_low,
    HighAlpha = count_high
  )
  print(contingency)
  
  ##  Perform chi-squared test
  chisq.test(contingency)
  
  ################
  ## Donut plot ##
  df_donut <- as.data.frame(contingency) %>%
    tibble::rownames_to_column("AlphaGroup") %>%
    tidyr::pivot_longer(-AlphaGroup, names_to = "Region", values_to = "Count") %>%
    group_by(AlphaGroup) %>%
    mutate(Percent = Count / sum(Count) * 100)
  
  # Add angle positions within each ring
  df_donut_label <- df_donut %>%
    group_by(AlphaGroup) %>%
    arrange(Region) %>%
    mutate(
      ymax = cumsum(Percent),
      ymin = lag(ymax, default = 0),
      ymid = (ymin + ymax) / 2
    )
  
  pdf("figures/donutFeatures.pdf", width = 4, height = 4)
  # Plot with percent labels inside rings
  ggplot(df_donut_label, aes(x = AlphaGroup, y = Percent, fill = Region)) +
    geom_col(width = 1, color = "white") +
    geom_text(
      aes(y = ymin, label = paste0(round(Percent), "%")),
      color = "white", size = 5, fontface = "bold"
    ) +
    scale_x_discrete(limits = c(" ", "LowAlpha", "HighAlpha")) +
    scale_fill_viridis_d() +
    coord_polar("y") +
    theme_void(base_size = 14) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    ggtitle(paste0(
      "Inner: proba to be hvCpG <= ", threshold,
      "\nOuter: proba to be hvCpG > ", threshold
    ))
  dev.off()
  
  rm(anno_result_highalpha, anno_result_lowalpha)
}

table(dt_clean$alpha >= 0.7)

# FALSE     TRUE 
# 23505671  3252148

dt_clean[dt_clean$alpha >= 0.7,"name"]


## TBC




#################################################################################
## Test for enrichment in other putative MEs for hvCpGs with alpha > threshold ##
#################################################################################

#############################################
## Needed to perform liftover hg19 to hg38 ##
library(rtracklayer)
library(GenomicRanges)

# Download the chain file
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_gz <- "dataPrev/hg19ToHg38.over.chain.gz"
chain_file <- "dataPrev/hg19ToHg38.over.chain"
if (!file.exists(chain_file)) {
  download.file(chain_url, chain_gz)
  R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
}
chain <- import.chain(chain_file)

## Manifest illumina450k to check arrays
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

################################
## To get alpha in each group ##
getMEdt <- function(gr_cpg, GRanges_hg38, MEgroup){
  MEgroup = unlist(strsplit(deparse(substitute(GRanges_hg38)), "_"))[1]
  hits = findOverlaps(gr_cpg, GRanges_hg38)
  cpg_in_MEs = gr_cpg[queryHits(hits)]
  print(length(cpg_in_MEs))
  cpg_in_MEs_dt <- data.table(
    chr = as.character(seqnames(cpg_in_MEs)),
    start_pos = start(cpg_in_MEs),
    end_pos = end(cpg_in_MEs),
    alpha = mcols(cpg_in_MEs)$alpha,
    ME = MEgroup
  )
  return(cpg_in_MEs)
}

## E.g. getMEdt(gr_cpg, VanBaakESS_GRanges_hg38)

#######################################
## Harris2012_1776SIV_10children450k ##
HarrisSIV <- readxl::read_excel("dataPrev/Harris2012_1776SIV_10children450k.xls", sheet = 3)
HarrisSIV_GRanges <- GRanges(
  seqnames = anno450k[match(HarrisSIV$Probe, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(HarrisSIV$Probe, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(HarrisSIV$Probe, anno450k$Name),"pos"],
                                  anno450k[match(HarrisSIV$Probe, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(HarrisSIV$Probe, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(HarrisSIV$Probe, anno450k$Name),"pos"] + 1,
                                anno450k[match(HarrisSIV$Probe, anno450k$Name),"pos"])),
  strand = anno450k[match(HarrisSIV$Probe, anno450k$Name),"strand"])

HarrisSIV_GRanges_hg38 <- unlist(liftOver(HarrisSIV_GRanges, chain))
## Select only the ones tested for alpha in Atlas
HarrisSIV_GRanges_hg38 <- getMEdt(gr_cpg, HarrisSIV_GRanges_hg38)

###########################
## VanBaak2018_ESS_HM450 ##
VanBaakESS <- readxl::read_excel("dataPrev/VanBaak2018_1580ESS_450k.xlsx", sheet = 2)
## only ESS hits
VanBaakESS <- VanBaakESS[VanBaakESS$`ESS hit`,]

VanBaakESS_GRanges <- GRanges(
  seqnames = anno450k[match(VanBaakESS$CG, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(VanBaakESS$CG, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(VanBaakESS$CG, anno450k$Name),"pos"],
                                  anno450k[match(VanBaakESS$CG, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(VanBaakESS$CG, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(VanBaakESS$CG, anno450k$Name),"pos"] + 1,
                                anno450k[match(VanBaakESS$CG, anno450k$Name),"pos"])),
  strand = anno450k[match(VanBaakESS$CG, anno450k$Name),"strand"])

VanBaakESS_GRanges_hg38 <- unlist(liftOver(VanBaakESS_GRanges, chain))
## Select only the ones tested for alpha in Atlas
VanBaakESS_GRanges_hg38 <- getMEdt(gr_cpg, VanBaakESS_GRanges_hg38)

###########################################
## Kessler2018_687SIVregions_2WGBS hg19! ##
KesslerSIV <- readxl::read_excel("dataPrev/Kessler2018_supTables.xlsx", sheet = 2, skip = 1)
KesslerSIV_GRanges <- GRanges(
  seqnames = KesslerSIV$Chromosome,
  ranges = IRanges(start = KesslerSIV$`ME start`, 
                   end = KesslerSIV$`ME end`),
  strand = "*")

KesslerSIV_GRanges_hg38 <- unlist(liftOver(KesslerSIV_GRanges, chain))
## Select only the ones tested for alpha in Atlas
KesslerSIV_GRanges_hg38 <- getMEdt(gr_cpg, KesslerSIV_GRanges_hg38)

#######################################
## Gunasekara2019_9926CoRSIVs_10WGBS ##
# Load corSIV intervals (already in hg38)
corSIV <- readxl::read_excel("dataPrev/Gunasekara2019_9926CoRSIVs_10WGBS.xls", sheet = 3)
corSIV <- unique(corSIV$USCS_Coordinates_CoRSIV)
corSIV_split <- tstrsplit(corSIV, "[:-]", fixed = FALSE)
corSIV_GRanges_hg38 <- GRanges(
  seqnames = corSIV_split[[1]],
  ranges = IRanges(start = as.integer(corSIV_split[[2]]), end = as.integer(corSIV_split[[3]])),
  strand = "*")
## Select only the ones tested for alpha in Atlas
corSIV_GRanges_hg38 <- getMEdt(gr_cpg, corSIV_GRanges_hg38)

#######################################
## Silver2022_SoCCpGs_10WGBS ##
arrayRef <- readxl::read_excel("dataPrev/Silver2022_259SoC_hg19.xlsx", sheet = 3, skip = 2)
SoCCpGs <- readxl::read_excel("dataPrev/Silver2022_259SoC_hg19.xlsx", sheet = 6, skip = 2)

SoCCpGs_GRanges <- GRanges(
  seqnames = paste0("chr", arrayRef$chr[match(SoCCpGs$cpg, arrayRef$cpg)]),
  ranges = IRanges(start = arrayRef$loc[match(SoCCpGs$cpg, arrayRef$cpg)],
                   end = arrayRef$loc[match(SoCCpGs$cpg, arrayRef$cpg)] + 1),
  strand = "*")

SoCCpGs_GRanges_hg38 <- unlist(liftOver(SoCCpGs_GRanges, chain))
## Select only the ones tested for alpha in Atlas
SoCCpGs_GRanges_hg38 <- getMEdt(gr_cpg, SoCCpGs_GRanges_hg38)

###########################
## Derakhshan2022_hvCpGs ##
DerakhshanhvCpGs <- readxl::read_excel("dataPrev/Derakhshan2022_4143hvCpGs_450k.xlsx", sheet = 6, skip = 3)

DerakhshanhvCpGs_GRanges <- GRanges(
  seqnames = anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"],
                                  anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"] + 1,
                                anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"])),
  strand = anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"strand"])

DerakhshanhvCpGs_GRanges_hg38 <- unlist(liftOver(DerakhshanhvCpGs_GRanges, chain))
## Select only the ones tested for alpha in Atlas
DerakhshanhvCpGs_GRanges_hg38 <- getMEdt(gr_cpg, DerakhshanhvCpGs_GRanges_hg38)

##########################################
## Matching genetic controls to hvCpGs  ##
mQTLcontrols <- read.table("../03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

mQTLcontrols_GRanges <- GRanges(
  seqnames = anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"],
                                  anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"] + 1,
                                anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"])),
  strand = anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"strand"])

mQTLcontrols_GRanges_hg38 <- unlist(liftOver(mQTLcontrols_GRanges, chain))
## Select only the ones tested for alpha in Atlas
mQTLcontrols_GRanges_hg38 <- getMEdt(gr_cpg, mQTLcontrols_GRanges_hg38)

####################
## Check overlaps ##

# Function to convert any GRanges to individual base positions as "chr:pos"
gr_to_pos <- function(gr) {
  unique(paste(gr@seqnames, gr@ranges))
}

sets <- list(
  HarrisSIV = gr_to_pos(HarrisSIV_GRanges_hg38),
  VanBaakESS = gr_to_pos(VanBaakESS_GRanges_hg38),
  KesslerSIV = gr_to_pos(KesslerSIV_GRanges_hg38),
  CoRSIV = gr_to_pos(corSIV_GRanges_hg38),
  SoCCpGs = gr_to_pos(SoCCpGs_GRanges_hg38),
  hvCpG = gr_to_pos(DerakhshanhvCpGs_GRanges_hg38),
  mQTLcontrols = gr_to_pos(mQTLcontrols_GRanges_hg38)
)

## To correct adding the CpGs sequenced in Atlas
library(UpSetR)
library(gridGraphics)
library(grid)

# Create the plot in a base graphics device
pdf(NULL)  # draw to null device to avoid displaying
upset(fromList(sets), nsets = 7, order.by = "freq")
grid_plot <- grid.grab()  # Capture as a grid object
dev.off()

# Now save the captured grid object to a real PDF
pdf("figures/upsetPreviousME.pdf", width = 12, height = 5)
grid.draw(grid_plot)
dev.off()

#############
## Combine ##

makedtMEset <- function(gr = gr_cpg, name = "all CpGs in atlas"){
  data.table(
    chr = as.character(seqnames(gr)),
    start_pos = start(gr),
    end_pos = end(gr),
    alpha = mcols(gr)$alpha,
    ME = name)
}

## Our dataset of CpG tested is in gr_cpg
MEsetdt <- na.omit(rbind(
  HarrisSIV = makedtMEset(HarrisSIV_GRanges_hg38, "HarrisSIV"),
  VanBaakESS = makedtMEset(VanBaakESS_GRanges_hg38, "VanBaakESS"),
  KesslerSIV = makedtMEset(KesslerSIV_GRanges_hg38, "KesslerSIV"),
  CoRSIV = makedtMEset(corSIV_GRanges_hg38, "CoRSIV"),
  SoCCpGs = makedtMEset(SoCCpGs_GRanges_hg38, "SoCCpGs"),
  DerakhshanhvCpGs = makedtMEset(DerakhshanhvCpGs_GRanges_hg38, "DerakhshanhvCpGs"),
  mQTLcontrols = makedtMEset(mQTLcontrols_GRanges_hg38, "mQTLcontrols")))

p <- ggplot(MEsetdt, aes(x = ME, y = alpha)) +
  geom_jitter(data = MEsetdt,
              aes(fill=ME), pch=21, size = 3, alpha = .1)+ 
  geom_violin(aes(col=ME))+
  geom_boxplot(aes(col=ME), width = .1) + 
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Probability of being a hvCpG")

pdf("figures/boxplot_otherMEs.pdf", width = 10, height = 4)
p
dev.off()

###################################################################
## Compare values for Derakshan hvCpGs from arrays vs from Atlas ##

makedtMEset(DerakhshanhvCpGs_GRanges_hg38, "DerakhshanhvCpGs")

Rred3array <- read.table("results_MariasarraysREDUCED_3samples_15datasets_6906CpGs_0_8p0_0_65p1.tsv",
                                                  header = T, sep = " ")
Rred3array_GRanges <- GRanges(
  seqnames = anno450k[match(Rred3array$CpG, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(Rred3array$CpG, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(Rred3array$CpG, anno450k$Name),"pos"],
                                  anno450k[match(Rred3array$CpG, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(Rred3array$CpG, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(Rred3array$CpG, anno450k$Name),"pos"] + 1,
                                anno450k[match(Rred3array$CpG, anno450k$Name),"pos"])),
  strand = anno450k[match(Rred3array$CpG, anno450k$Name),"strand"],
  alpha = Rred3array$alpha,
  ishvCpG = Rred3array$ishvCpG)

Rred3array_GRanges_hg38 <- unlist(liftOver(Rred3array_GRanges, chain))

# Find overlapping CpGs
hits <- findOverlaps(Rred3array_GRanges_hg38, DerakhshanhvCpGs_GRanges_hg38)

# Extract alpha values
alpha_dt <- data.table(
  alpha_x = Rred3array_GRanges_hg38$alpha[queryHits(hits)],
  alpha_y = DerakhshanhvCpGs_GRanges_hg38$alpha[subjectHits(hits)]
)

# Remove rows with NA values if needed
alpha_dt <- na.omit(alpha_dt)

cor_val <- cor(alpha_dt$alpha_x, alpha_dt$alpha_y)
ggplot(alpha_dt, aes(x = alpha_x, y = alpha_y)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 0.05, y = 0.95, hjust = 0, label = paste("r =", round(cor_val, 2))) +
  labs(
    x = "Alpha (Rred3array)",
    y = "Alpha (Derakhshan)"
  ) +
  theme_minimal(base_size = 14)

#######################################################################
## Investigate variation of purified cell types hvCpGs in Atlas data ##














## Statistical testing: does alpha differ between cpg types?
# kruskal.test(alpha ~ ME, data = MEsetdt)
# Kruskal-Wallis rank sum test
# 
# data:  alpha by ME
# Kruskal-Wallis chi-squared = 37.045, df = 5, p-value = 0.0000005868

## Pairwise testing:
pairwise.wilcox.test(MEsetdt$alpha, MEsetdt$ME,
                     p.adjust.method = "BH")  # Benjamini-Hochberg correction



