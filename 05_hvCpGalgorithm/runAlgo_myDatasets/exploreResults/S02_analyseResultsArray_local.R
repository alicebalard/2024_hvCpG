## This script does:
## Load full results on array ##
## Calculate proba hvCpG minus matching control: is it always +? ##
## Load full results on array with only 3 individuals/ds ##

## Prepare
library(here)

source(here("05_hvCpGalgorithm/quiet_library.R"))
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R"))

hvCpGandControls <- prephvCpGandControls(codeDir = "~/Documents/GIT/2024_hvCpG/")

################################
## Load full results on array ##
################################

load(here("05_hvCpGalgorithm/resultsDir/Arrays/results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1.RData"))

resArrayAll <- as.data.frame(results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1)
rm(results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1)

prepareChrDataset <- function(res){
  res$chrpos <- hvCpGandControls$dictionary$hg38[
    match(rownames(res), hvCpGandControls$dictionary$illu450k)]
  
  ## Indicate the hvCpG of Maria and controls
  res$group <- NA
  res$group[
    res$chrpos %in% 
      hvCpGandControls$DerakhshanhvCpGs_names] <- "hvCpG_Derakhshan"
  res$group[
    res$chrpos %in% 
      hvCpGandControls$mQTLcontrols_names] <- "mQTLcontrols"
  
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
  
  # Assign alternating black/grey per chromosome
  chr_colors <- data.frame(
    chr = chr_order,
    point_col = rep(c("black", "grey60"), length.out = length(chr_order))
  )
  
  # Merge with res
  res <- res %>%
    left_join(chr_colors, by = "chr")
  
  # Remove values with no position
  res <- res[!is.na(res$chrpos),]
  
  return(res)
}

resArrayAll <- prepareChrDataset(resArrayAll)

# Plot
# Compute midpoints for chromosome labels
chr_mid <- resArrayAll %>%
  group_by(chr) %>%
  summarise(mid = (min(cum_pos) + max(cum_pos)) / 2)
pdf(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_array.pdf"), width = 15, height = 3)
## colorblind friendly
ggplot() +
  geom_point(data = resArrayAll, 
             aes(x = cum_pos, y = alpha, col = "black"),
             alpha = 0.05, size = 1) +
  # Highlight hvCpG
  geom_point(data = resArrayAll[resArrayAll$group %in% "hvCpG_Derakhshan", ],
             aes(x = cum_pos, y = alpha),
             col = "#DC3220", alpha = 0.8) +
  # Highlight mQTL controls
  geom_point(data = resArrayAll[resArrayAll$group %in% "mQTLcontrols", ],
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
###################################################################

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
hv_alpha <- resArrayAll[, c("chrpos", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- resArrayAll[, c("chrpos", "alpha")]
colnames(ctrl_alpha) <- c("control", "alpha_control")

# Join everything
merged <- pairs %>%
  left_join(hv_alpha, by = "hvCpG") %>%
  left_join(ctrl_alpha, by = "control") %>%
  mutate(diffAlpha=alpha_hvCpG-alpha_control)

pdf(here("05_hvCpGalgorithm/figures/DifferenceOfProbabilityForhvCpG-matching_controlInArray.pdf"), width = 4, height = 5)
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

###########################################################
## Load full results on array with only 3 individuals/ds ##
###########################################################

load(here("05_hvCpGalgorithm/resultsDir/Arrays/results_Arrays_3indperds_394240CpGs_0_8p0_0_65p1.RData"))

resArray3ind <- as.data.frame(results_Arrays_3indperds_394240CpGs_0_8p0_0_65p1)
rm(results_Arrays_3indperds_394240CpGs_0_8p0_0_65p1)

resArray3ind <- prepareChrDataset(resArray3ind)

names(resArray3ind)[names(resArray3ind) %in% "alpha"] <- "alpha_array_3ind"
resCompArray <- dplyr::left_join(resArray3ind, resArrayAll)
names(resCompArray)[names(resCompArray) %in% "alpha"] <- "alpha_array_all"

p1 <- ggplot(resCompArray, 
       aes(x=alpha_array_all, y=alpha_array_3ind, fill = group, col = group)) +
  geom_point(data = resCompArray[is.na(resCompArray$group),],
             pch = 21, alpha = 0.05) +
  geom_point(data = resCompArray[!is.na(resCompArray$group),],
             pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"),
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  theme(legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Probability of being hypervariable",
       x = "P(hv) using full array datasets",
       y = "P(hv) using reduced (3 ind/ds) array datasets")

pdf(here("05_hvCpGalgorithm/figures/arrayfullvsred3ind_myalgo"), width = 9, height = 7)
p1
dev.off()

## What cutoff to get the same number of sites than Maria?
poscutoff = 0.94
negcutoff = 0.5

nrow(resCompArray[resCompArray$alpha_array_all > poscutoff,])# == 3535

## True positive: detected with full AND reduced array / all detected with full array
nrow(resCompArray[resCompArray$alpha_array_all > poscutoff &
               resCompArray$alpha_array_3ind > poscutoff,]) /
  nrow(resCompArray[resCompArray$alpha_array_all > poscutoff,]) * 100

## False positive = % CpGs detected only using reduced datasets
nrow(resCompArray[resCompArray$alpha_array_all < negcutoff &
                    resCompArray$alpha_array_3ind > poscutoff,]) /
  nrow(resCompArray[resCompArray$alpha_array_3ind > poscutoff,]) * 100

## False negative = % CpGs detected only using full datasets
nrow(resCompArray[resCompArray$alpha_array_all > poscutoff &
                    resCompArray$alpha_array_3ind < negcutoff,]) /
  nrow(resCompArray[resCompArray$alpha_array_all > poscutoff,]) * 100

## TBC

# ##############################################
# ## --- Test 5: batch correction effect? --- ##
# ##############################################
# 
# ## Compare results array with either raw data uncorrected or corrected
# load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Arrays_noCorrectionInRaw/results_Arrays_noCorrectionInRaw__406334CpGs_0_8p0_0_65p1.RData")
# 
# resArrayNoCor <- results_Arrays_noCorrectionInRaw__406334CpGs_0_8p0_0_65p1
# rm(results_Arrays_noCorrectionInRaw__406334CpGs_0_8p0_0_65p1)
# 
# resArrayNoCor <- resArrayNoCor %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "cpgprobe") %>%
#   dplyr::rename(alpha_array_nocor = alpha)
# 
# resArrayNoCor$chrpos = hvCpGandControls$dictionary$hg38[
#   match(resArrayNoCor$cpgprobe,
#         hvCpGandControls$dictionary$illu450k)]
# 
# resCommon_Array_Atlas_rawArray <- full_join(res_Alpha_Atlas, resArrayNoCor)
# 
# p1 <- ggplot(resCommon_Array_Atlas_rawArray, 
#              aes(x=alpha_array_all, y=alpha_array_nocor, fill = group, col = group)) +
#   geom_point(data = resCommon_Array_Atlas_rawArray[is.na(resCommon_Array_Atlas_rawArray$group),],
#              pch = 21, alpha = 0.05) +
#   geom_point(data = resCommon_Array_Atlas_rawArray[!is.na(resCommon_Array_Atlas_rawArray$group),],
#              pch = 21, alpha = 0.4) +
#   geom_smooth(method = "lm", fill = "black") +
#   scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
#                     labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
#   scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
#   theme_minimal(base_size = 14) +
#   guides(fill = guide_legend(position = "inside"))+
#   theme(legend.position.inside = c(0.3,0.8),
#         legend.box = "horizontal", legend.title = element_blank(),
#         legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
#         legend.key = element_rect(fill = "white", color = NA)) +
#   labs(title = "Probability of being hypervariable",
#        x = "P(hv) considering array data after full correction",
#        y = "P(hv) considering array data without correction")
# 
# p2 <- ggplot(resCommon_Array_Atlas_rawArray, 
#              aes(x=alpha_atlas, y=alpha_array_nocor, fill = group, col = group)) +
#   geom_point(data = resCommon_Array_Atlas_rawArray[is.na(resCommon_Array_Atlas_rawArray$group),],
#              pch = 21, alpha = 0.05) +
#   geom_point(data = resCommon_Array_Atlas_rawArray[!is.na(resCommon_Array_Atlas_rawArray$group),],
#              pch = 21, alpha = 0.4) +
#   geom_smooth(method = "lm", fill = "black") +
#   scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
#                     labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
#   scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
#   theme_minimal(base_size = 14) +
#   guides(fill = guide_legend(position = "inside"))+
#   theme(legend.position.inside = c(0.3,0.8),
#         legend.box = "horizontal", legend.title = element_blank(),
#         legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
#         legend.key = element_rect(fill = "white", color = NA)) +
#   labs(title = "Probability of being hypervariable",
#        x = "P(hv) considering atlas WGBS data",
#        y = "P(hv) considering array data without correction")
# 
# # --- Turn off legends inside plots ---
# p1_clean <- p1 + theme(legend.position = "none")
# p2_clean <- p2 + theme(legend.position = "none")
# 
# # --- Extract one legend (e.g. from p1) ---
# legend <- cowplot::get_legend(
#   p1 + theme(legend.position = "bottom",
#              legend.box = "horizontal",
#              legend.title = element_blank(),
#              legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
#              legend.key = element_rect(fill = "white", color = NA))
# )
# 
# # --- Arrange plots with legend as 4th panel ---
# pdf(here("05_hvCpGalgorithm/figures/test5_compArrayVsnocor.pdf"), width = 10, height = 7)
# cowplot::plot_grid(p1_clean, p2_clean, legend,
#                    ncol = 2, rel_heights = c(.9, .3))  # grid layout: 2 cols × 2 rows
# dev.off()

## rm junk
rm(x,y, pairs, merged, chr_mid, hv_alpha, data, ctrl_alpha, resArray3ind, resArrayAll)

saveRDS(resCompArray, here("05_hvCpGalgorithm/dataOut/resArray.RDS"))
