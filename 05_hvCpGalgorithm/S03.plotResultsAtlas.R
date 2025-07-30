## July 2025
## A. Balard
## Plot files for scan Atlas

library(dplyr)
library(data.table)
library(progress)
library(readxl)
library(ggplot2)
library(tidyr)
library(scales)
library(viridis)
library(ggrastr)     # For rasterized geom_point
library(Cairo)       # For high-quality PDF
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)

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

options(scipen=999)
# Plot
pdf("figures/freqCpGperdataset.pdf", width = 15, height = 4)
ggplot(t_combined, aes(x = as.factor(datasets_covered_in), y = num_CpGs, fill = coverage)) +
  geom_col(position = "dodge") +
  scale_y_log10() +
  scale_fill_manual(values = c("≥5" = "steelblue", "≥10" = "firebrick")) +
  labs(
    title = "CpG coverage across datasets",
    x = "Number of datasets with ≥3 samples covered",
    y = "Number of CpGs (log scale)",
    fill = "Coverage threshold"
  ) +
  theme_minimal(base_size = 14)
dev.off()

####################
## Manhattan plot ##
####################

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- "resultsDir"

# Get list of relevant RData files
rdata_files <- dir(parent_dir, pattern = "results_Atlas_100000CpGs_0_8p0_0_65p1\\.RData$", 
                   recursive = TRUE, full.names = TRUE)

## Check if all batches have ran
length(rdata_files)

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
dt[, chr := factor(chr, levels = as.character(c(1:22, "X")))]

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
df2 <- merge(data.frame(chr = factor(c(1:22, "X"), levels=as.character(c(1:22, "X")))),
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
CairoPDF("figures/ManhattanAlphaPlot.pdf", width = 15, height = 3)
print(plot)
dev.off()

################################################
## Test enrichment of features for high alpha ##
################################################

#############
threshold=0.9
#############

# Filter valid rows
dt_clean <- dt[!is.na(start_pos) & !is.na(end_pos)]

# Create GRanges
gr_cpg <- GRanges(
  seqnames = paste0("chr", dt_clean$chr),
  ranges = IRanges(start = dt_clean$start_pos, end = dt_clean$end_pos),
  alpha = dt_clean$alpha
)

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

#############################################################
## Test for enrichment in corSIV of hvCpG with alpha > 0.9 ##
#############################################################

# Load corSIV intervals
corSIV <- readxl::read_excel("13059_2019_1708_MOESM1_ESM.xls", sheet = 3)
corSIV <- unique(corSIV$USCS_Coordinates_CoRSIV)

# Parse corSIV to data.table with numeric chromosome, start, end
corSIV_split <- tstrsplit(corSIV, "[:-]", fixed = FALSE)
corSIV_dt <- data.table(
  chr = sub("chr", "", corSIV_split[[1]]),
  start_pos = as.integer(corSIV_split[[2]]),
  end_pos = as.integer(corSIV_split[[3]])
)
 
# Key the tables for foverlaps
setkey(corSIV_dt, chr, start_pos, end_pos)

dt_subset <- dt[, .(chr = as.character(chr), start_pos, end_pos)]
setkey(dt_subset, chr, start_pos, end_pos)
dt_subset[, row_id := .I]

# Find overlaps
overlaps <- foverlaps(dt_subset, corSIV_dt, nomatch = 0L)
overlap_idx <- overlaps[, row_id]

# Mark points in corSIV intervals
dt[, in_corSIV := FALSE]
dt[overlap_idx, in_corSIV := TRUE]

# Create binary column for alpha > threshold
dt[, alpha_high := alpha > threshold]

# Count combinations of (in_corSIV, alpha_high)
tbl <- dt[, .N, by = .(in_corSIV, alpha_high)]

# Format as contingency table
# Row 1: alpha > threshold; Row 2: alpha <= threshold
# Column 1: in corSIV; Column 2: not in corSIV

A <- tbl[in_corSIV == TRUE & alpha_high == TRUE, N]
B <- tbl[in_corSIV == FALSE & alpha_high == TRUE, N]
C <- tbl[in_corSIV == TRUE & alpha_high == FALSE, N]
D <- tbl[in_corSIV == FALSE & alpha_high == FALSE, N]

# Replace missing counts with zero
A <- ifelse(length(A) == 0, 0, A)
B <- ifelse(length(B) == 0, 0, B)
C <- ifelse(length(C) == 0, 0, C)
D <- ifelse(length(D) == 0, 0, D)

# Build matrix
contingency_matrix <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE,
                             dimnames = list(`Alpha > threshold` = c("Yes", "No"),
                                             `In corSIV` = c("Yes", "No")))

# Perform Fisher's Exact Test
fisher_result <- fisher.test(contingency_matrix)

# View result
print(contingency_matrix)

print(fisher_result)
# 
# Fisher's Exact Test for Count Data
# 
# data:  contingency_matrix
# p-value = 7.389e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.064985 1.155738
# sample estimates:
# odds ratio 
#   1.109572 
# 🔍 Fisher’s Exact Test Output
# p-value: 7.389e-07 → Highly significant
# → Strong evidence that the distribution of alpha > threshold differs between corSIV and non-corSIV regions.
# 
# Odds ratio: 1.11
# → Sites with alpha > threshold are ~11% more likely to occur in corSIV regions than expected by chance.
# 
# 95% CI: [1.065, 1.156]
# → This interval does not include 1, further supporting statistical significance.

# CCl: We observe a statistically significant enrichment of high-alpha CpG sites (alpha > threshold) in corSIV regions (Fisher’s exact test: OR = 1.11, p = 7.4 × 10⁻⁷), suggesting a non-random association between regulatory variability and corSIV domains."

p <- ggplot(dt, aes(x = in_corSIV, y = alpha)) +
 geom_violin(aes(fill = in_corSIV), trim = FALSE) + 
  geom_boxplot(width = 0.2) + theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#FC4E07", "#00AFBB"))+
  xlab("CpG in corSIV (Gunasekara et al. 2019)") +
  ylab("Probability of being a hvCpG") + theme(legend.position = "none")

pdf("figures/boxplot_corSIV.pdf", width = 5, height = 5)
p
dev.off()
