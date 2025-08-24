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

# Alternate background colors
bg <- chr_sizes %>%
  mutate(bg_color = rep(c("white", "grey95"), length.out = n()))

# Plot
ggplot(res, aes(x = cum_pos, y = alpha, color = chr)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_point(data = res[res$group %in% "hvCpG_Derakhshan",],
             aes(x = cum_pos, y = alpha), pch=21, col = "red")+
  geom_point(data = res[res$group %in% "mQTLcontrols",],
             aes(x = cum_pos, y = alpha), pch=21, col = "black")+
  scale_fill_identity() +
  scale_color_manual(values = rep(c("steelblue", "darkorange"), length.out = length(chr_order))) +
  scale_x_continuous(breaks = chr_mid$mid, labels = gsub("chr", "", chr_mid$chr)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = "Probability alpha of being a hvCpG")
