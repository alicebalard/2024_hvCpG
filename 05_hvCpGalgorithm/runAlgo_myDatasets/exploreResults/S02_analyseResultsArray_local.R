setwd("~/Documents/GIT/2024_hvCpG/")

source("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R")
source("05_hvCpGalgorithm/quiet_library.R")
hvCpGandControls <- prephvCpGandControls(codeDir = "~/Documents/GIT/2024_hvCpG/")

load("05_hvCpGalgorithm/resultsDir/Arrays/results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1.RData")

results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1 <- as.data.frame(results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1)
results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1$chrpos <- hvCpGandControls$dictionary$hg38[
  match(rownames(results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1), hvCpGandControls$dictionary$illu450k)]

## Indicate the hvCpG of Maria and controls
results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1$group <- NA
results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1$group[
  results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1$chrpos %in% 
    hvCpGandControls$DerakhshanhvCpGs_names] <- "hvCpG_Derakhshan"
results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1$group[
  results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1$chrpos %in% 
    hvCpGandControls$mQTLcontrols_names] <- "mQTLcontrols"

res = results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1

# Parse chromosome and position
res <- res %>%
  mutate(
    chr = str_extract(chrpos, "^chr[0-9XY]+"),
    pos = as.numeric(str_extract(chrpos, "(?<=_)[0-9]+"))
  )

# Order chromosomes
chr_order <- paste0("chr", c(1:22, "X", "Y"))
res$chr <- factor(res$chr, levels = chr_order)

# Compute cumulative position for genome-wide x-axis
chr_sizes <- res %>%
  group_by(chr) %>%
  summarise(max_pos = max(pos)) %>%
  mutate(cum_start = lag(cumsum(max_pos), default = 0))

res <- res %>%
  left_join(chr_sizes, by = "chr") %>%
  mutate(cum_pos = pos + cum_start)

# Compute midpoints for chromosome labels
chr_mid <- res %>%
  group_by(chr) %>%
  summarise(mid = (min(cum_pos) + max(cum_pos)) / 2)

# Assign alternating black/grey per chromosome
chr_colors <- data.frame(
  chr = chr_order,
  point_col = rep(c("black", "grey60"), length.out = length(chr_order))
)

# Merge with res
res <- res %>%
  left_join(chr_colors, by = "chr")

# Plot
pdf("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_array.pdf", width = 15, height = 3)
## colorblind friendly
ggplot() +
  geom_point(data = res, aes(x = cum_pos, y = alpha, col = point_col),
             alpha = 0.1, size = 1) +
  # Highlight hvCpG
  geom_point(data = res[res$group %in% "hvCpG_Derakhshan", ],
             aes(x = cum_pos, y = alpha),
             col = "#DC3220", alpha = 0.8) +
  # Highlight mQTL controls
  geom_point(data = res[res$group %in% "mQTLcontrols", ],
             aes(x = cum_pos, y = alpha),
             col = "#005AB5", alpha = 0.8) +
  scale_color_identity() +
  scale_x_continuous(breaks = chr_mid$mid,
                     labels = gsub("chr", "", chr_mid$chr),
                     expand = c(0, 0)) + ## rm padding
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")
dev.off()

###################################################################
## Calculate proba hvCpG minus matching control: is it always +? ##

data <- read.table(file.path(codeDir = "~/Documents/GIT/2024_hvCpG/", "03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)

x = hvCpGandControls$dictionary$hg38[match(data$hvCpG_name, hvCpGandControls$dictionary$illu450k)]
y = hvCpGandControls$dictionary$hg38[match(data$controlCpG_name, hvCpGandControls$dictionary$illu450k)]

# Build mapping from hvCpG -> control
pairs <- data.frame(
  hvCpG = x,
  control = y,
  stringsAsFactors = FALSE
)

# Merge hvCpG alphas
hv_alpha <- res[, c("chrpos", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- res[, c("chrpos", "alpha")]
colnames(ctrl_alpha) <- c("control", "alpha_control")

# Join everything
merged <- pairs %>%
  left_join(hv_alpha, by = "hvCpG") %>%
  left_join(ctrl_alpha, by = "control") %>%
  mutate(diffAlpha=alpha_hvCpG-alpha_control)

pdf("05_hvCpGalgorithm/figures/DifferenceOfProbabilityForhvCpG-matching_controlInArray.pdf", width = 4, height = 5)
ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control) in array")+
  ylab("Difference of probability")
dev.off()
