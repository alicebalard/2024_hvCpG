## July 2025
## A. Balard
## Plot files for scan Atlas

packages <- c("dplyr", "data.table", "readxl", "progress", "ggplot2",
              "tidyr", "scales", "viridis", "ggrastr", "Cairo",
              "genomation", "GenomicRanges", "GenomicFeatures", "boot")

lapply(packages, library, character.only = TRUE)

## This code does:
### I. Histogram of coverage across datasets
### II. Load data & Manhattan plot
### III. Test enrichment of features for high alpha
### IV. Test for enrichment in other putative MEs for hvCpGs with alpha > threshold

##############################################
## I. Histogram of coverage across datasets ##
##############################################

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

t_combined[t_combined$coverage %in% "≥10" &  t_combined$datasets_covered_in %in% 46,"num_CpGs"] /
  sum(t_combined[t_combined$coverage %in% "≥10","num_CpGs"])

## 84% of all CpGs (23/27.5M) are covered in 46 cell types
rm(t_combined, t5, t10)

####################################
## II. Load data & Manhattan plot ##
####################################

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- "resultsDir/Atlas10X"

# Get list of relevant RData files
rdata_files <- dir(parent_dir, pattern = "results_Atlas10X_100000CpGs_0_8p0_0_65p1\\.RData$", 
                   recursive = TRUE, full.names = TRUE)

## Check if all batches have ran
length(rdata_files) ## 287

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

################################## *** NEW *** #####################################
## Select only positions that are covered in all cell types, in at least 3 people ##
cpgs46M <- read.table("~/Documents/10X/selected_cpgs_min3_in46_datasets.txt")

dt <- dt[dt$name %in% cpgs46M$V1,] ## 7th August: 22.6M
rm(cpgs46M) 
#######################################################################

# Parse "chr_start-end" in name into chr, start_pos, end_pos. NB: takes a couple of minutes
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
# df2 <- dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
# df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y"), levels=as.character(c(1:22, "X", "Y")))),
#              df2, by = "chr", all.x = TRUE, sort = TRUE)
# df2 <- na.omit(df2)
# 
# plot <- ggplot() +
#   geom_point_rast(data = dt, aes(x = pos2, y = alpha),
#                   color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
#   theme_classic() + theme(legend.position = "none") +
#   scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr)) +
#   xlab("Chromosome") +
#   ylab("pr(Z=1) of being a hvCpG") +
#   theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
# CairoPDF("figures/ManhattanAlphaPlot.pdf", width = 15, height = 3)
# print(plot)
# dev.off()

#####################################################
## III. Test enrichment of features for high alpha ##
#####################################################

#############
threshold=0.7
#############

# Filter valid rows
dt_clean <- dt[!is.na(start_pos) & !is.na(end_pos)]
rm(dt)

# Create GRanges
gr_cpg <- GRanges(
  seqnames = paste0("chr", dt_clean$chr),
  ranges = IRanges(start = dt_clean$start_pos, end = dt_clean$end_pos),
  alpha = dt_clean$alpha
)

doAnnot = F

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

#######################################################################
## Investigate variation of purified cell types hvCpGs in Atlas data ##

table(dt_clean$alpha >= 0.7)
# FALSE     TRUE 
# 20708740  1872287 

cpg_names_all <- rhdf5::h5read("/home/alice/Documents/10X/all_scaled_matrix.h5", "cpg_names")
metadata <- read.table("/home/alice/Documents/10X/sample_metadata.tsv", sep ="\t", header = TRUE)
medsd_lambdas <- read.table("/home/alice/Documents/10X/all_medsd_lambda.tsv", sep = "\t", header = TRUE)

source_M_batchCpG <- function(cpg_indices) {
  # Read a block of CpGs (rows = samples, cols = CpGs)
  M <- rhdf5::h5read(
    "/home/alice/Documents/10X/all_scaled_matrix.h5",
    "scaled_matrix",
    index = list(NULL, cpg_indices)
  )
  rownames(M) <- metadata$sample
  colnames(M) <- cpg_names_all[cpg_indices]
  return(M)
}

myrun_batch <- function(batch_size = 10000, 
                        cpg_ids) {
  cpgPos_all <- match(cpg_ids, cpg_names_all)
  # Prepare result table
  dfres <- data.table(dataset = unique(metadata$dataset), hv = 0, nonhv = 0)
  setkey(dfres, dataset)
  meta_split <- split(metadata, metadata$dataset)
  dataset_samples <- lapply(meta_split, function(df) df$sample)
  medsd_dt <- data.table(medsd_lambdas)
  setkey(medsd_dt, dataset)
  # Loop through CpGs in batches
  batches <- split(cpgPos_all, ceiling(seq_along(cpgPos_all) / batch_size))
  batchtimer <- 1
  for (batch in batches) {
    print(paste0("Processing batch ", batchtimer, "..."))
    M_block <- source_M_batchCpG(batch)  # Samples x CpGs matrix
    for (dataset_name in names(dataset_samples)) {
      samples <- dataset_samples[[dataset_name]]
      idx <- match(samples, rownames(M_block))
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) next
      M_ds <- M_block[idx, , drop = FALSE]
      # Calculate SD across samples for each CpG in this dataset
      sds <- apply(M_ds, 2, sd, na.rm = TRUE)
      threshold <- medsd_dt[dataset == dataset_name, median_sd] *
        medsd_dt[dataset == dataset_name, lambda]
      dfres[dataset == dataset_name, hv := hv + sum(sds >= threshold, na.rm = TRUE)]
      dfres[dataset == dataset_name, nonhv := nonhv + sum(sds < threshold, na.rm = TRUE)]
    }
    batchtimer <- batchtimer + 1
  }
  return(dfres)
}

dfres <- myrun_batch(cpg_ids =  sample(unlist(dt_clean[dt_clean$alpha >= 0.7, "name"]), 100))

# system.time(dfres <- myrun_batch(
#   batch_size = 10000,
#   cpg_ids = unlist(dt_clean[dt_clean$alpha >= 0.7, "name"]))) # should take ~ 5 to 10h

dfres <- dfres %>% mutate(propHv = hv / (hv + nonhv))

pdf("figures/propPurCellVar.pdf", width = 10, height = 7)
ggplot(dfres, aes(x=propHv, y = dataset))+
  geom_point() +
  theme_minimal(base_size = 14) + ylab("") +
  xlab("Proportion of purified cell hvCpG in the top 5% variable CpGs in this dataset\n(ran on x random hvCpGs)")
dev.off()

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
getMEdt <- function(gr_cpg, GRanges_hg38, MEgroup_name = NULL) {
  hits <- findOverlaps(gr_cpg, GRanges_hg38)
  # Subset only overlapping CpGs
  cpg_in_regions <- gr_cpg[queryHits(hits)]
  # Add ME region info from GRanges_hg38 if wanted
  if (!is.null(MEgroup_name)) {
    mcols(cpg_in_regions)$ME <- MEgroup_name
  }
  return(cpg_in_regions)  # Still a GRanges object, one row per CpG/alpha
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
length(HarrisSIV_GRanges_hg38) # 1402

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
length(VanBaakESS_GRanges_hg38) #1269

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
length(KesslerSIV_GRanges_hg38) ## 3994 positions

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
length(corSIV_GRanges_hg38) # 70477 pos

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
length(SoCCpGs_GRanges_hg38) # 203 pos

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
length(DerakhshanhvCpGs_GRanges_hg38) # 3430

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
length(mQTLcontrols_GRanges_hg38) # 2929

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

## Check that covered in 46 cell types
table(paste0(MEsetdt$chr, "_", MEsetdt$start_pos, "-", MEsetdt$end_pos) %in%
        dt_clean$name)

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

###################################################
## SD in hvCpG and controls per dataset in atlas ##

library(boot)
library(data.table)
library(matrixStats)

myrun_batch_boot <- function(batch_size = 10000, cpg_ids, B = 1000) {
  cpgPos_all <- match(cpg_ids, cpg_names_all)
  
  # Prepare result tables
  datasets <- unique(metadata$dataset)
  dfres <- data.table(dataset = datasets)
  boot_res_list <- list()
  
  # Preprocess
  setkey(dfres, dataset)
  meta_split <- split(metadata, metadata$dataset)
  dataset_samples <- lapply(meta_split, function(df) df$sample)
  medsd_dt <- data.table(medsd_lambdas)
  setkey(medsd_dt, dataset)
  
  # Batching
  batches <- split(cpgPos_all, ceiling(seq_along(cpgPos_all) / batch_size))
  batchtimer <- 1
  
  for (batch in batches) {
    message(paste0("Processing batch ", batchtimer, " of ", length(batches)))
    M_block <- source_M_batchCpG(batch)  # Samples x CpGs matrix
    
    for (dataset_name in names(dataset_samples)) {
      samples <- dataset_samples[[dataset_name]]
      idx <- match(samples, rownames(M_block))
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) next
      
      M_ds <- M_block[idx, , drop = FALSE]
      
      # SDs for each CpG in this dataset
      sds <- apply(M_ds, 2, sd, na.rm = TRUE)
      
      # Bootstrapped median SD per dataset (pooled per dataset, across CpGs)
      boot_median <- function(data, indices) {
        median(data[indices], na.rm = TRUE)
      }
      
      boot_obj <- boot(sds, statistic = boot_median, R = B)
      ci <- tryCatch(boot.ci(boot_obj, type = "perc"), error = function(e) NULL)
      if (!is.null(ci) && !is.null(ci$percent)) {
        lower_ci <- ci$percent[4]
        upper_ci <- ci$percent[5]
      } else {
        lower_ci <- NA
        upper_ci <- NA
      }
      
      boot_res_list[[dataset_name]] <- data.frame(
        dataset = dataset_name,
        median_sd = median(sds, na.rm = TRUE),
        lower = lower_ci,
        upper = upper_ci
      )
    }
    batchtimer <- batchtimer + 1
  }
  
  boot_df <- do.call(rbind, boot_res_list)
  final_df <- merge(dfres, boot_df, by = "dataset", all = TRUE)
  return(final_df)
}

df_hvstats_controls <- myrun_batch_boot(
  batch_size = 10000,
  cpg_ids = paste0(mQTLcontrols_GRanges_hg38@seqnames, "_", mQTLcontrols_GRanges_hg38@ranges),
  B = 1000)

df_hvstats_hvCpGs <- myrun_batch_boot(
  batch_size = 1000,
  cpg_ids = paste0(DerakhshanhvCpGs_GRanges_hg38@seqnames, "_", DerakhshanhvCpGs_GRanges_hg38@ranges),
  B = 1000)

## Newly detected hvCpG as a comparison
df_hvstats_newhvCpGs07 <- myrun_batch_boot(
  batch_size = 1000,
  cpg_ids = sample(dt_clean$name[dt_clean$alpha > 0.7 & !is.na(dt_clean$alpha)], 3000),
  B = 1000)

df <- rbind(df_hvstats_controls %>% mutate(type = "mQTLcontrol"),
            df_hvstats_hvCpGs %>% mutate(type = "hvCpG (Derakhshan)"),
            df_hvstats_newhvCpGs07 %>% mutate(type = "3000 purCells-hvCpGs proba 70%+"))

df$type <- factor(df$type, levels = levels(factor(df$type))[c(1,3,2,4,5)])

ggplot(df, aes(x = dataset, y = median_sd, fill = type, col = type)) +
  geom_pointrange(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 1),
    shape = 21, size = 0.4, stroke = 0.5
  ) +
  geom_point(
    pch = 21, size = 2,
    position = position_dodge(width = 1)
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "black", size = 0.4),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(title = NULL), colour = guide_legend(NULL)) +
  xlab("") +
  ylab("Median SD ± 95% CI")

## Only hvCpGs of Maria and controls
df <- rbind(df_hvstats_controls %>% mutate(type = "mQTLcontrol"),
            df_hvstats_hvCpGs %>% mutate(type = "hvCpG (Derakhshan)"))

ggplot(df,
       aes(x = dataset, y = median_sd, fill = type, col = type)) +
  geom_pointrange(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 1),
    shape = 21, size = 0.4, stroke = 0.5
  ) +
  geom_point(
    pch = 21, size = 2,
    position = position_dodge(width = 1)
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "black", size = 0.4),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(title = NULL), colour = guide_legend(NULL)) +
  xlab("") + coord_cartesian(ylim = c(0,2))+
  ylab("Median SD ± 95% CI")

### Check one example of high alpha CpG, and see the profile
unlog <- function(x) {
  odds <- 2^x
  beta <- odds / (1 + odds)
  beta
}

source("../05_hvCpGalgorithm/hvCpG_algorithm_detection_v4scan.R")

library(patchwork)

beta_to_log2M <- function(beta, epsilon = 1e-6) {
  # Clip beta values for numerical stability
  beta <- pmin(pmax(beta, epsilon), 1 - epsilon)
  log2(beta / (1 - beta))
}

exploreAlgo <- function(x,title){
  alpha_trace <<- list()
  
  ## Prepare data in the environment:
  prep = prepData("testLocalPC")
  metadata = prep$metadata
  medsd_lambdas = prep$medsd_lambdas
  cpg_names_all = prep$cpg_names_all
  source_M_1CpG = prep$source_M_1CpG
  
  runAndSave(
    analysis = "testLocalPC",
    cpgPos_vec = x,
    resultDir = "~/Documents/Project_hvCpG/RESULT/",
    NCORES = 1,
    p0 = 0.80,
    p1 = 0.65,
    overwrite = TRUE
  ) 
  load("~/Documents/Project_hvCpG/RESULT/results_testLocalPC_1CpGs_0_8p0_0_65p1.RData")
  message("alpha:")
  res = results_testLocalPC_1CpGs_0_8p0_0_65p1
  print(res)
  
  cpgRaw = source_M_batchCpG(cpg_indices = x)
  cpgRaw <- as.data.frame(cpgRaw)
  cpgRaw$sample <- rownames(cpgRaw)
  
  # Reshape to long format: sample, CpG, value
  long_df <- reshape2::melt(cpgRaw, id.vars = "sample", variable.name = "CpG",
                            value.name = "value")
  # Merge with metadata to get dataset info
  long_df <- left_join(long_df, metadata, by = "sample")
  
  # Now plot: dataset on x, value on y
  p0 <- ggplot(long_df, aes(x = dataset, y = unlog(value))) +
    geom_point(position = position_jitter(width = 0.01), alpha = 0.7, size = 2) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Raw methylation") + xlab("") +
    theme(axis.text.x = element_text(size =6)) +
    facet_wrap(~ CpG, scales = "free_y", ncol = 3)
  
  p1 <- ggplot(long_df, aes(x = dataset, y = value)) +
    geom_point(position = position_jitter(width = 0.01), alpha = 0.7, size = 2) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Logit methylation") + xlab("") +
    theme(axis.text.x = element_text(size =6)) +
    facet_wrap(~ CpG, scales = "free_y", ncol = 3)
  
  sd_long_df <- long_df %>% group_by(dataset, CpG) %>% 
    summarise(sdMethyl = sd(value, na.rm = T))
  
  sd_long_df <- merge(sd_long_df, medsd_lambdas)
  sd_long_df$p95 <- sd_long_df$lambda * sd_long_df$median_sd
  
  sd_long_df <- sd_long_df[sd_long_df$dataset %in% long_df$dataset,]
  
  p2 <- ggplot(sd_long_df, aes(x = dataset)) +
    geom_segment(aes(
      x = dataset,
      xend = dataset,
      y = median_sd,
      yend = sdMethyl,
      color = sdMethyl > p95
    ),
    arrow = arrow(length = unit(0.15, "cm")),
    position = position_jitter(width = 0.2)) +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green")) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    xlab("") +
    facet_wrap(~ CpG, scales = "free_y", ncol = 3) +
    guides(color = "none") +
    ggtitle(label = "Up = sd j > median sd j,k", subtitle = "green = top5% var")
  
  alpha_df <- do.call(rbind, lapply(alpha_trace, as.data.frame))
  p3 <- ggplot(alpha_df, aes(x = alpha, y = loglik)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    labs(title = "Alpha Optimization Trace",
         x = expression(alpha),
         y = "Log-likelihood")
  
  clip_and_logit <- function(beta, epsilon = 0.1) {
    beta <- pmin(pmax(beta, epsilon), 1 - epsilon)
    log2(beta / (1 - beta))
  }
  
  p4 <- ggplot(long_df, aes(x = dataset, y = clip_and_logit(unlog(value)))) +
    geom_point(position = position_jitter(width = 0.01), alpha = 0.7, size = 2) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Clipped THEN logit methylation") + xlab("") +
    theme(axis.text.x = element_text(size =6)) +
    facet_wrap(~ CpG, scales = "free_y", ncol = 3)
  
  layout <- (p0 / p1) | (p2 / p3/p4)
  print(layout)
  print(table(sd_long_df$sdMethyl > sd_long_df$thr))
}

# a hvCpG:
exploreAlgo(x = match(paste0(DerakhshanhvCpGs_GRanges_hg38@seqnames, "_", DerakhshanhvCpGs_GRanges_hg38@ranges)[4],
                      cpg_names_all), title = "a hvCpG from Maria")

# a high alpha
exploreAlgo(x = match(dt_clean$name[dt_clean$alpha > 0.9 & !is.na(dt_clean$alpha)] %>% head(2) %>% tail(1),
                      cpg_names_all), title = "a CpG with high alpha")

# a low alpha
exploreAlgo(x = match(dt_clean$name[dt_clean$alpha < 0.1 & !is.na(dt_clean$alpha)] %>% head(3) %>% tail(1),
                      cpg_names_all), title = "a CpG with low alpha")



exploreAlgo(x = 57175, title = "a hvCpG from Maria")



# a high alpha that should be low 
exploreAlgo(x = match(dt_clean$name[dt_clean$alpha > 0.9 & !is.na(dt_clean$alpha)] %>% head(2) %>% tail(1),
                      cpg_names_all), title = "a CpG with high alpha")

## To log to match lambda 

## Find variable CpGs in data
findVar <- function(x){
  
  layout <- (p1 / p2 )
  print(layout)
}

findVar(85300)

match(dt_clean$name, cpg_names_all) %>% head
cpg_names_all[244] #  [1] 13344 20890 20896 20940 48812* 48813 49062 57175** 57177 85300

match(paste0(DerakhshanhvCpGs_GRanges_hg38@seqnames, "_", 
             DerakhshanhvCpGs_GRanges_hg38@ranges)[1:10],
      cpg_names_all)

# ➤ A CpG with low variation across datasets but high alpha
# In few datasets, the data may still slightly favor the hvCpG model (e.g., by chance).
# 
# Even a small deviation from expected variability under the null could make the alternative statistically better, just because there’s so little data to contradict it.
# 
# With only 2–3 datasets, the log-likelihood difference needed to swing toward hvCpG can be small.
# 
# ➤ A CpG with lots of variability but low alpha
# Even if variability is high, if it is consistent across datasets, the model could conclude that it’s intrinsically noisy, not heterogeneous across datasets.
# 
# That would support a CpG model, not a hvCpG one.
# 
# Alpha will stay low because the more datasets you have, the stronger evidence you need to say "this CpG is inconsistent across datasets."


beta_to_log2M(0)
