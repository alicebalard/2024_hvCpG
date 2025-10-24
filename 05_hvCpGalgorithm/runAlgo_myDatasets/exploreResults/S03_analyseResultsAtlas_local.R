#################################################
## Plot results of algorithm ran on atlas data ##
#################################################
library(here)

source(here("05_hvCpGalgorithm", "quiet_library.R"))

## Load array results
resArray <- readRDS(here("05_hvCpGalgorithm/dataOut/resArray.RDS"))

## Add Maria's results
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R"))
hvCpGandControls <- prephvCpGandControls(codeDir = here())

## This code does:
### I. Histogram of coverage across datasets
### II. Load data & Manhattan plot
### III. Test enrichment of features for high alpha
### IV. Test for enrichment in other putative MEs for hvCpGs with alpha > threshold

## Data in WGBS atlas:
## from the CS cluster: sample_groups <- h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/all_matrix_noscale.h5","sample_groups")
sample_groups <- readRDS(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/sample_groups.RDS"))

ggplot(data.frame(table(sample_groups)), aes(x = Freq)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of number of samples per dataset",
    x = "Number of samples",
    y = "Count of datasets"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 1))

table(sample_groups)

##############################################
## I. Histogram of coverage across datasets ##
##############################################

t5 <- read.table(here("04_prepAtlas/CpG_coverage_freqtable5X.tsv"), header = T)
t10 <- read.table(here("04_prepAtlas/CpG_coverage_freqtable10X.tsv"), header = T)

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
pdf(here("05_hvCpGalgorithm", "figures", "freqCpGperdataset.pdf"), width = 14, height = 4)
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
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA))
dev.off()

t_combined[t_combined$coverage %in% "≥10" &  t_combined$datasets_covered_in %in% 46,"num_CpGs"] /
  sum(t_combined[t_combined$coverage %in% "≥10","num_CpGs"])

## 84% of all CpGs (23/27.5M) are covered in 46 cell types
rm(t_combined, t5, t10)

####################################
## II. Load data & Manhattan plot ##
####################################

prepAtlasdt <- function(){
  # Define parent folder containing all "Atlas_batchXXX" folders
  parent_dir <- here("05_hvCpGalgorithm/resultsDir/Atlas/Atlas10X/")
  
  # Get list of relevant RData files
  rdata_files <- dir(parent_dir, pattern = "results_Atlas10X_[0-9]+CpGs_0_8p0_0_65p1\\.RData$", 
                     recursive = TRUE, full.names = TRUE)
  
  ## Check if all batches have ran
  length(rdata_files) == 231
  
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
  
  rm(e, pb, all_cpg_values, obj, file, parent_dir, rdata_files)
  
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
  dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y", "M")))]
  
  ## Check chromosomes order:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  ## Mark group membership in dt
  dt[, group := NA_character_]
  dt[name %in% hvCpGandControls$DerakhshanhvCpGs_names, group := "hvCpG_Derakhshan"]
  dt[name %in% hvCpGandControls$mQTLcontrols_names, group := "mQTLcontrols"]
  
  # Compute cumulative position offsets for Manhattan plot
  setorder(dt, chr, start_pos)
  
  offsets <- dt[, .(max_pos = max(start_pos, na.rm = TRUE)), by = chr]
  offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]
  
  dt <- merge(dt, offsets[, .(chr, cum_offset)], by = "chr", all.x = TRUE, sort = FALSE)
  
  # Convert to integer/numeric if not already
  dt[, cum_offset := as.numeric(cum_offset)]
  dt[, pos2 := start_pos + cum_offset]
  
  return(dt)
}

system.time(Atlas_dt <- prepAtlasdt())

if (exists("doIprepAtlas") && isTRUE(doIprepAtlas)) {
  stop("stop here to only prepare atlas_dt")
}

# Compute chromosome centers for x-axis labeling
df2 <- Atlas_dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y", "M"), levels=as.character(c(1:22, "X", "Y", "M")))),
             df2, by = "chr", all.x = TRUE, sort = TRUE)
df2 <- na.omit(df2)

# Compute chromosome boundaries
df_bounds <- Atlas_dt[, .(min_pos = min(pos2, na.rm = TRUE), 
                          max_pos = max(pos2, na.rm = TRUE)), by = chr]

# Midpoints between chromosomes = where to draw dotted lines
df_bounds[, next_start := data.table::shift(min_pos, n = 1, type = "lead")]
vlines <- df_bounds[!is.na(next_start), .(xintercept = (max_pos + next_start)/2)]

plot <- ggplot() +
  # background cloud
  geom_point_rast(data = Atlas_dt[is.na(group)], 
                  aes(x = pos2, y = alpha),
                  color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
  # hvCpG highlights
  geom_point(data = Atlas_dt[group == "hvCpG_Derakhshan"],
             aes(x = pos2, y = alpha),
             color = "#DC3220", size = 1, alpha = 0.7) +
  # mQTL controls highlights
  geom_point(data = Atlas_dt[group == "mQTLcontrols"],
             aes(x = pos2, y = alpha),
             color = "#005AB5", size = 1, alpha = 0.7) +
  # Add dotted separators
  geom_vline(data = vlines, aes(xintercept = xintercept),
             linetype = 3, color = "grey60") +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
CairoPDF(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_previoushvCpGplotted_atlas.pdf"), width = 15, height = 3)
print(plot)
dev.off()

## Without the layer of previous hvCpG and controls plotted:
plot <- ggplot() +
  geom_point_rast(data = Atlas_dt, 
                  aes(x = pos2, y = alpha),
                  color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
   # Add dotted separators
  geom_vline(data = vlines, aes(xintercept = xintercept),
             linetype = 3, color = "grey60") +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
CairoPDF(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_atlas.pdf"), width = 15, height = 3)
print(plot)
dev.off()

##################################
## Save names high alpha points ##
##################################
table(Atlas_dt$alpha >= 0.7)
# FALSE     TRUE 
# 22335423   700603 

## Map on arrays
dico <- readRDS(here("05_hvCpGalgorithm/dataOut/anno_combined_probes_hg38.rds"))

matches <- match(
  x = unlist(Atlas_dt[Atlas_dt$alpha >= 0.7, "name"]),
  table = dico$chrpos_hg38
)

highAlphaPos <- dico[na.omit(matches), ]

table(highAlphaPos$array)
# 450k    450k and EPIC      EPIC 
# 340          4483          3431 

# 4483+340 = 4823 on the 450k array

saveRDS(highAlphaPos, here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/fetalSIV/highAlphaPos_atlas0.7.RDS"))

##########################################
## What are the gaps in Manhattan plot? ##
# Compute the gap between consecutive CpGs on the same chromosome
Atlas_dt[, gap := start_pos - data.table::shift(end_pos), by = chr]

# Identify large gaps (>= 500k bp)
gaps_dt <- Atlas_dt[gap >= 500000, .(
  chr,
  gap_start = data.table::shift(end_pos),
  gap_end = start_pos,
  gap_size = gap
)]

# Drop first NA (since shift introduces one per chromosome)
gaps_dt[!is.na(gap_size)]
## Load the results at https://genome.ucsc.edu/cgi-bin/hgTracks

####################################
## Mitochondrial DNAm variability ##
####################################
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09541-9?utm_source=chatgpt.com
## Near absence of 5mC, so expected

ggplot() +
  geom_point(data = Atlas_dt[Atlas_dt$chr == "M",], 
             aes(x = pos2, y = alpha),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2[df2$chr == "M","center"],
                     labels = as.character(df2[df2$chr == "M","chr"]), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,.1)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

###################################
## Y chromosome DNAm variability ##
###################################
ggplot() +
  geom_point(data = Atlas_dt[Atlas_dt$chr == "Y",], 
             aes(x = start_pos, y = alpha),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  geom_hline(yintercept = .7, linetype = 3)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

Atlas_dt[Atlas_dt$chr == "Y" & Atlas_dt$alpha > 0.7,]
# High alpha in 3 regions: chrY:5,043,848-6,534,238, chrY:10,107,290-11,747,410, chrY:56,822,399-56,841,336

#####################################################
## Find regions of high alpha using sliding window ##
#####################################################
# 
# high_alpha_thresh <- 0.7
# Atlas_dt[, high_alpha := alpha > high_alpha_thresh]
# 
# # Function to compute high-alpha fraction in each window
# find_high_alpha_windows <- function(dt, window_size, step_size) {
#   res <- list()
#   for (ch in unique(dt$chr)) {
#     chr_dt <- dt[chr == ch]
#     starts <- seq(min(chr_dt$start_pos), max(chr_dt$end_pos), by = step_size)
#     for (s in starts) {
#       e <- s + window_size - 1
#       w <- chr_dt[start_pos >= s & end_pos <= e]
#       if (nrow(w) > 0) {
#         frac <- sum(w$high_alpha) / nrow(w)
#         res[[length(res)+1]] <- data.table(chr = ch, window_start = s, window_end = e, frac_high_alpha = frac)
#       }
#     }
#   }
#   rbindlist(res)
# }
# 
# # system.time(high_alpha_windows <- find_high_alpha_windows(Atlas_dt, window_size = 100000, step_size = 50000))
## 100000, 50000: 10 minutes to run
# saveRDS(object = high_alpha_windows, 
#         file = here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/high_alpha_windows.RDS"))
# 
# # Keep windows where >50% of CpGs have high alpha
# hotspots <- high_alpha_windows[frac_high_alpha > 0.5]
# hotspots
# chr window_start window_end frac_high_alpha
# <char>        <num>      <num>           <num>
#   1:      1       517452     617451       1.0000000
# 2:      3     92911806   93011805       0.6666667
# 3:      6     32496889   32596888       0.5670498
# 4:      6     32546889   32646888       0.5853659
# 5:      6     32596889   32696888       0.5498721
# 6:      9     61210596   61310595       0.6666667
# 7:     15     21850602   21950601       0.6000000
# 8:     16     37260230   37360229       1.0000000
# 9:     16     37310230   37410229       1.0000000


# subAtlasdt <- Atlas_dt[Atlas_dt$chr ==3 & 
#                          Atlas_dt$start_pos >= 92911806 & 
#                          Atlas_dt$end_pos <= 93011805,]









###################################################################
## Calculate proba hvCpG minus matching control: is it always +? ##

data <- read.table(here("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)

x = hvCpGandControls$dictionary$hg38[match(data$hvCpG_name, hvCpGandControls$dictionary$illu450k)]
y = hvCpGandControls$dictionary$hg38[match(data$controlCpG_name, hvCpGandControls$dictionary$illu450k)]

# Build mapping from hvCpG -> control
pairs <- data.frame(
  hvCpG = x,
  control = y,
  stringsAsFactors = FALSE
)

# Merge hvCpG alphas
res <- Atlas_dt[!is.na(group)]

hv_alpha <- res[, c("name", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- res[, c("name", "alpha")]
colnames(ctrl_alpha) <- c("control", "alpha_control")

# Join everything
merged <- pairs %>%
  left_join(hv_alpha, by = "hvCpG") %>%
  left_join(ctrl_alpha, by = "control") %>%
  mutate(diffAlpha=alpha_hvCpG-alpha_control)

merged <- merged %>%
  mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
  filter(!is.na(diffAlpha))

pdf(here("05_hvCpGalgorithm/figures/DifferenceOfProbabilityForhvCpG-matching_controlInAtlas.pdf"),
    width = 4, height = 5)
ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control) in atlas")+
  ylab("Difference of probability")
dev.off()

ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control) in atlas")+
  ylab("Difference of probability") +
  facet_grid(.~chr)

## How different from 50/50 expected? Binomial test

mean_diff <- mean(merged$diffAlpha, na.rm = TRUE)
median_diff <- median(merged$diffAlpha, na.rm = TRUE)
summary(merged$diffAlpha)
# If mean/median ≈ 0.02 (like your later ME plot), it means hvCpGs are slightly enriched for higher hypervariability probability compared to matched controls.

## How different from 50/50 expectation?
n_pos <- sum(merged$diffAlpha > 0, na.rm = TRUE)
n_neg <- sum(merged$diffAlpha < 0, na.rm = TRUE)
n_total <- n_pos + n_neg

binom.test(n_pos, n_total, p = 0.5, alternative = "greater")

## Effect size (proportion above zero)
n_pos / n_total

#####################################################
## III. Test enrichment of features for high alpha ##
#####################################################

##############
threshold=0.7#
##############

# Filter valid rows
dt_clean <- Atlas_dt[!is.na(start_pos) & !is.na(end_pos)]
rm(Atlas_dt)

# Create GRanges
gr_cpg <- GRanges(
  seqnames = paste0("chr", dt_clean$chr),
  ranges = IRanges(start = dt_clean$start_pos, end = dt_clean$end_pos),
  alpha = dt_clean$alpha
)

# Import bed file
bed_features <- genomation::readTranscriptFeatures(
  "~/Documents/Project_hvCpG/hg38_GENCODE_V47.bed")

# annotate CpGs with high alpha (>threshold):
anno_result_highalpha <- genomation::annotateWithGeneParts(
  target = gr_cpg[!is.na(mcols(gr_cpg)$alpha) & mcols(gr_cpg)$alpha > threshold],
  feature = bed_features)

# annotate CpGs with low alpha (<=threshold):
anno_result_lowalpha <- genomation::annotateWithGeneParts(
  target = gr_cpg[!is.na(mcols(gr_cpg)$alpha) & mcols(gr_cpg)$alpha <= threshold],
  feature = bed_features)

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
#            promoter    exon   intron intergenic
# LowAlpha   3112722 1179458 12071239    4987877
# HighAlpha    53647   21547   343090     166446

##  Perform chi-squared test
chisq.test(contingency)
# Pearson's Chi-squared test
# data:  contingency
# X-squared = 21801, df = 3, p-value < 0.00000000000000022

# promoters and exons are depleted, and intergenic and introns are enriched, in hypervariable CpGs

################
## Barplot ##
df_plot <- as.data.frame(contingency) %>%
  tibble::rownames_to_column("AlphaGroup") %>%
  tidyr::pivot_longer(-AlphaGroup, names_to = "Region", values_to = "Count") %>%
  group_by(AlphaGroup) %>%
  mutate(Percent = Count / sum(Count) * 100)

pdf(here("05_hvCpGalgorithm/figures/barplotFeaturesLowHighAlpha.pdf"), width = 5, height = 4)
ggplot(df_plot, aes(x=Region, y=Percent, fill = AlphaGroup))+
  geom_bar(position="dodge", stat="identity") +
  theme_minimal(base_size = 14)+
  scale_fill_manual(
    values = c("red", "skyblue"),
    name = "p(hv)",      # optional
    labels = c(">70%", "<=70%")   # new names for legend keys
  ) +
  theme(axis.title.x = element_blank())
dev.off()

rm(anno_result_highalpha, anno_result_lowalpha)

#################################################################################
## Test for enrichment in other putative MEs for hvCpGs with alpha > threshold ##
#################################################################################

#######################################
## Harris2012_1776SIV_10children450k ##
HarrisSIV <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Harris2012_1776SIV_10children450k.xls"), sheet = 3)
HarrisSIV_hg38 <- hvCpGandControls$dictionary[
  match(HarrisSIV$Probe, hvCpGandControls$dictionary$illu450k), "hg38"]
HarrisSIV_hg38 <- na.omit(HarrisSIV_hg38)
length(HarrisSIV_hg38) # 1773

###########################
## VanBaak2018_ESS_HM450 ##
VanBaakESS <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/VanBaak2018_1580ESS_450k.xlsx"), sheet = 2)
## only ESS hits
VanBaakESS <- VanBaakESS[VanBaakESS$`ESS hit`,]
VanBaakESS_hg38 <- hvCpGandControls$dictionary[
  match(VanBaakESS$CG, hvCpGandControls$dictionary$illu450k), "hg38"]
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

#######################################
## Silver2022_SoCCpGs_10WGBS ##
arrayRef <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Silver2022_259SoC_hg19.xlsx"), sheet = 3, skip = 2)
SoCCpGs <- readxl::read_excel(here("05_hvCpGalgorithm/dataPrev/Silver2022_259SoC_hg19.xlsx"), sheet = 6, skip = 2)

SoCCpGs_GRanges <- GRanges(
  seqnames = paste0("chr", arrayRef$chr[match(SoCCpGs$cpg, arrayRef$cpg)]),
  ranges = IRanges(start = arrayRef$loc[match(SoCCpGs$cpg, arrayRef$cpg)],
                   end = arrayRef$loc[match(SoCCpGs$cpg, arrayRef$cpg)] + 1),
  strand = "*")

## liftover to hg38, keep uniquely mapping regions
mapped <- liftOver(SoCCpGs_GRanges, hvCpGandControls$chain)
keep <- lengths(mapped) == 1
hg38_unique <- unlist(mapped[keep])

## find the match with Atlas cpg
overlaps <- findOverlaps(query = hg38_unique, subject = cpg_46_GR)
CpG_overlapping     <- cpg_46_GR[subjectHits(overlaps)]
SoCCpGs_hg38 <- paste0(CpG_overlapping@seqnames, "_", CpG_overlapping@ranges)

SoCCpGs_hg38 <- na.omit(SoCCpGs_hg38)
length(SoCCpGs_hg38) #177

#################################
## Check overlaps with upset plot 
sets <- list(
  HarrisSIV = HarrisSIV_hg38,
  VanBaakESS = VanBaakESS_hg38,
  KesslerSIV = KesslerSIV_hg38,
  CoRSIV = corSIV_hg38,
  SoCCpGs = SoCCpGs_hg38,
  hvCpG = hvCpGandControls$DerakhshanhvCpGs_names,
  mQTLcontrols = hvCpGandControls$mQTLcontrols_names
)

# Create the plot in a base graphics device
pdf(NULL)  # draw to null device to avoid displaying
upset(fromList(sets), nsets = 7, order.by = "freq")
grid_plot <- grid.grab()  # Capture as a grid object
dev.off()

# Now save the captured grid object to a real PDF
pdf(here("05_hvCpGalgorithm/figures/upsetPreviousME.pdf"), width = 12, height = 5)
grid.draw(grid_plot)
dev.off()

#######################################################
## Check alpha for the different MEs: are they high? ##

MEsetdt <- rbindlist(lapply(names(sets), function(nm) {
  dt_clean[.(sets[[nm]]), on = .(name), .(name, alpha)][, ME := nm]
}))

MEsetdt <- na.omit(MEsetdt) ## 33478 so far (2 sept)

## Control as baseline
MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

p1 <- ggplot(MEsetdt, aes(x = ME, y = alpha)) +
  geom_jitter(data = MEsetdt,
              aes(fill=ME), pch=21, size = 3, alpha = .1)+
  geom_violin(aes(col=ME))+
  geom_boxplot(aes(col=ME), width = .1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Probability of being a hvCpG")

p1

## Statistical comparisons of alpha between MEs

# 1️⃣ Fit the model with mQTLcontrols as baseline
fit <- lm(alpha ~ ME, data = MEsetdt)

# 2️⃣ Get estimated marginal means and contrasts vs baseline
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

# 3️⃣ Prepare for plotting
contrasts <- contrasts %>%
  mutate(ME = contrast,  # rename for clarity
         lower = estimate - 1.96*SE,
         upper = estimate + 1.96*SE)

# 4️⃣ Plot
p2 <- ggplot(contrasts, aes(x = ME, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    y = "Difference in alpha vs mQTLcontrols",
    x = "ME group",
    title = "Comparison of ME groups to mQTLcontrols",
    subtitle = "lm with multiple comparison correction (Sidak)"
  ) +
  theme_minimal()

pdf(here("05_hvCpGalgorithm/figures/alphaComparisonBetweenMEtypes.pdf"), width = 13, height = 4)
cowplot::plot_grid(p1,p2, rel_widths = c(1, .8))
dev.off()

#################################
### Regions of hypervariation ###
#################################
plot <- ggplot() +
  geom_point_rast(data = dt_clean, aes(x = pos2, y = alpha), color = "black",
                  size = 0.01, alpha = 0.01, raster.dpi = 72) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
CairoPDF(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_atlas_2.pdf"), width = 15, height = 3)
print(plot)
dev.off()

## TBC










###############
## FIND PAX8
# Chromosome 2, NC_000002.12 (113215997..113278921,

# 1. Convert array and atlas to GRanges object
gr_atlas <- GRanges(
  seqnames = paste0("chr", Atlas_dt$chr),
  ranges = IRanges(start = Atlas_dt$start_pos, end = Atlas_dt$end_pos),
  mcols = Atlas_dt
)

gr_array <- GRanges(
  seqnames = resArray$chr,
  ranges = IRanges(start = resArray$pos, end = resArray$pos +1),
  mcols = resArray
)

# 2. PAX8 as a GRanges
gr_PAX8 <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(start = 113215997, end = 113278921)
)

# 3. Find overlaps between CpGs and the gene region
hits <- findOverlaps(gr_atlas, gr_PAX8)
hits2 <- findOverlaps(gr_array, gr_PAX8)

# 4. Extract matching rows from original df
df_hits <- gr_atlas[queryHits(hits), ]
df_hits2 <- gr_array[queryHits(hits2), ]

ggplot(data.frame(df_hits), aes(x=start)) + 
  geom_point(aes(y=mcols.alpha)) +
  # geom_point(data.frame(df_hits2), aes(y=mcols.alpha_array_all, fill=mcols.group), size = 3, pch =21) +
  theme_minimal(base_size = 14) +
  # Add the rectangle
  annotate(
    "rect",
    xmin = 113278394,
    xmax = 113279523,
    ymin = -Inf,      # spans entire y-range
    ymax = Inf,
    alpha = 0.2,
    fill = "blue"
  ) +
  # Add a label for the promoter
  annotate(
    "text",
    x = (113278394 + 113279523) / 2,
    y = max(df_hits$mcols.alpha, na.rm = TRUE)-.2,
    label = "promoter",
    vjust = -1,
    size = 5,
    fontface = "bold",
    color = "blue"
  )
  