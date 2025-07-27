## July 2025
## A. Balard
## Plot files for scan Atlas

library(data.table)
library(progress)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(readxl)
library(Cairo)

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- "resultsDir"

# Get list of relevant RData files
rdata_files <- dir(parent_dir, pattern = "results_Atlas_100000CpGs_0_8p0_0_65p1\\.RData$", 
                   recursive = TRUE, full.names = TRUE)

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
dt[, `:=`(start_pos = as.integer(start_pos), end_pos = as.integer(end_pos))]

# Convert chr from "chrN" to integer factor
dt[, chr := as.integer(sub("chr", "", chr))]
dt[, chr := factor(chr, levels = 1:22)]

# Compute cumulative position offsets for Manhattan plot
setorder(dt, chr, start_pos)
offsets <- dt[, .(max_pos = max(start_pos, na.rm = TRUE)), by = chr]
offsets[, cum_offset := shift(cumsum(max_pos), fill = 0)]
dt <- merge(dt, offsets[, .(chr, cum_offset)], by = "chr", all.x = TRUE, sort = FALSE)
dt[, pos2 := start_pos + cum_offset]

# Compute chromosome centers for x-axis labeling
df2 <- dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
df2 <- merge(data.frame(chr = factor(1:22, levels=1:22)), df2, by = "chr", all.x = TRUE, sort = TRUE)

# Load corSIV intervals
corSIV <- readxl::read_excel("13059_2019_1708_MOESM1_ESM.xls", sheet = 3)
corSIV <- unique(corSIV$USCS_Coordinates_CoRSIV)

# Parse corSIV to data.table with numeric chromosome, start, end
corSIV_split <- tstrsplit(corSIV, "[:-]", fixed = FALSE)
corSIV_dt <- data.table(
  chr = as.integer(sub("chr", "", corSIV_split[[1]])),
  start_pos = as.integer(corSIV_split[[2]]),
  end_pos = as.integer(corSIV_split[[3]])
)

# Key the tables for foverlaps
setkey(corSIV_dt, chr, start_pos, end_pos)

dt_subset <- dt[, .(chr = as.integer(as.character(chr)), start_pos, end_pos)]
setkey(dt_subset, chr, start_pos, end_pos)
dt_subset[, row_id := .I]

# Find overlaps
overlaps <- foverlaps(dt_subset, corSIV_dt, nomatch = 0L)
overlap_idx <- overlaps[, row_id]

# Mark points in corSIV intervals
dt[, in_corSIV := FALSE]
dt[overlap_idx, in_corSIV := TRUE]

# Plot with base points black, highlighted points red and bigger
plot <- ggplot() +
  geom_point_rast(data = dt[in_corSIV == FALSE], aes(x = pos2, y = alpha),
                  color = "black", size = 0.3, alpha = 0.005, raster.dpi = 72) +
  geom_point_rast(data = dt[in_corSIV == TRUE], aes(x = pos2, y = alpha),
                  color = "red", size = 1, alpha = 0.1, raster.dpi = 72) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center,
                     labels = c(1:15, "", 17, "", 19, "", 21, "")) +
  xlab("Chromosome")

# Save plot
CairoPNG("ManhattanAlphaPlot.png", width = 1600, height = 600)
print(plot)
dev.off()

###########################################################
## Test for enrichment of hvCpG with alpha > p in corSIV ##
###########################################################

threshold=0.7

# Create binary column for alpha > threshold
dt[, alpha_high := alpha > threshold]

# Count combinations of (in_corSIV, alpha_high)
tbl <- dt[, .N, by = .(in_corSIV, alpha_high)]

# Format as contingency table
# Row 1: alpha > 0.6; Row 2: alpha <= 0.6
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
#                In corSIV
# Alpha > 0.6   Yes      No
#         Yes  2689  854700
#         No  17381 6130075

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
# → Strong evidence that the distribution of alpha > 0.6 differs between corSIV and non-corSIV regions.
# 
# Odds ratio: 1.11
# → Sites with alpha > 0.6 are ~11% more likely to occur in corSIV regions than expected by chance.
# 
# 95% CI: [1.065, 1.156]
# → This interval does not include 1, further supporting statistical significance.

# CCl: We observe a statistically significant enrichment of high-alpha CpG sites (alpha > 0.6) in corSIV regions (Fisher’s exact test: OR = 1.11, p = 7.4 × 10⁻⁷), suggesting a non-random association between regulatory variability and corSIV domains."