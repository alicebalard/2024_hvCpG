#################################################
## Plot results of algorithm ran on atlas data ##
#################################################
library(here)

source(here("05_hvCpGalgorithm", "quiet_library.R"))

## Load array results
resArray <- readRDS(here("05_hvCpGalgorithm/dataOut/resArray.RDS"))

## This code does:
### I. Histogram of coverage across datasets
### II. Load data & Manhattan plot
### III. Test enrichment of features for high alpha
### IV. Test for enrichment in other putative MEs for hvCpGs with alpha > threshold

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

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- here("05_hvCpGalgorithm/resultsDir/Atlas10X/")

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
dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y")))]

## Check chromosomes order:
message("Chromosomes in the dataset:")
table(unique(dt$chr))

## Add Maria's results
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R"))
hvCpGandControls <- prephvCpGandControls(codeDir = here())

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

# Compute chromosome centers for x-axis labeling
df2 <- dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y"), levels=as.character(c(1:22, "X", "Y")))),
             df2, by = "chr", all.x = TRUE, sort = TRUE)
df2 <- na.omit(df2)

plot <- ggplot() +
  # background cloud
  geom_point_rast(data = dt[is.na(group)], 
                  aes(x = pos2, y = alpha),
                  color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
  # hvCpG highlights
  geom_point(data = dt[group == "hvCpG_Derakhshan"],
             aes(x = pos2, y = alpha),
             color = "#DC3220", size = 1, alpha = 0.7) +
  # mQTL controls highlights
  geom_point(data = dt[group == "mQTLcontrols"],
             aes(x = pos2, y = alpha),
             color = "#005AB5", size = 1, alpha = 0.7) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
CairoPDF(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_atlas.pdf"), width = 15, height = 3)
print(plot)
dev.off()

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
res <- dt[!is.na(group)]

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

#####################################################
## III. Test enrichment of features for high alpha ##
#####################################################

##############
threshold=0.75
##############

# Filter valid rows
dt_clean <- dt[!is.na(start_pos) & !is.na(end_pos)]
rm(dt)

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

##  Perform chi-squared test
chisq.test(contingency)

# HighAlpha CpGs are depleted in promoters/exons and enriched in intergenic regions (with a mild enrichment in introns).

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
  scale_fill_viridis_d() +
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
cpg_46 <- read.table("~/Documents/Project_hvCpG/10X/selected_cpgs_min3_in46_datasets.txt")$V1
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
length(KesslerSIV_hg38) # 940

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
length(corSIV_hg38) # 71320

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
length(SoCCpGs_hg38) #204

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

tail(dt_clean)
tail(resArray)
