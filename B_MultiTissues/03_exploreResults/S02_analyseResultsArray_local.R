## This script does:
## Load full results on array ##
## Calculate proba hvCpG minus matching control: is it always +? ##
## Load full results on array with only 3 individuals/ds ##

## Prepare
library(here)
## Load libraries
if (!exists("libLoaded")) {
  source(here("B_MultiTissues", "quiet_library.R"))}

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Add previous MEs including Maria's results
## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}

################################
## Load full results on array ##
################################

load(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1.RData"))

resArrayAll <- as.data.frame(results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1)
rm(results_arrayAll_algov5_394240CpGs_0_8p0_0_65p1)

prepareChrDataset <- function(res){
  res$chrpos <- dico$chrpos_hg38[
    match(rownames(res), dico$CpG)]
  
  ## Indicate the hvCpG of Maria and controls
  res$group <- NA
  res$group[
    res$chrpos %in% DerakhshanhvCpGs_hg38] <- "hvCpG_Derakhshan"
  res$group[
    res$chrpos %in% mQTLcontrols_hg38] <- "mQTLcontrols"
  
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
## colorblind friendly
p1_manhattanArray <- ggplot() +
  geom_point(data = subset(resArrayAll, is.na(group)),
    aes(x = cum_pos, y = alpha), color = "gray", alpha = .5, size = .8) +
  geom_point(data = subset(resArrayAll, !is.na(group)),
    aes(x = cum_pos, y = alpha, col = group), alpha = .8, size = 1) +
  scale_color_manual(values = c("hvCpG_Derakhshan" = "#DC3220",
      "mQTLcontrols" = "#005AB5")) +
  scale_x_continuous(breaks = chr_mid$mid,
    labels = gsub("chr", "", chr_mid$chr), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.box.margin = margin(10, 0, 0, 0),
        legend.margin = margin(-5, 0, 0, 0),
        plot.margin = margin(t = 1, r = 5, b = 5, l = 5)) 

###################################################################
## Calculate proba hvCpG minus matching control: is it always +? ##
###################################################################

x = DerakhshanhvCpGs_hg38
y = mQTLcontrols_hg38

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

p2_DiffProbhvCpG_matchingcontrol_Array <- ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control)", subtitle =  "in array")+
  ylab("Difference of probability") +
  coord_cartesian(ylim = c(-1,1))

################################################################
## Load full results on array with only 2 or 3 individuals/ds ##
################################################################

makePlotNrob <- function(resCompArray, N){
  mycor <- cor(resCompArray$alpha_array_all, resCompArray$alpha_array_reduce, use = "complete.obs")
  
  ggplot(resCompArray,
         aes(x=alpha_array_all, y=alpha_array_reduce)) +
    geom_point(data = resCompArray[is.na(resCompArray$group),], aes(col = group),
               alpha = 0.05) +
    geom_point(data = resCompArray[!is.na(resCompArray$group),], aes(col = group),
               alpha = 0.4) +
    geom_smooth(method = "lm", fill = "grey", col = "grey") +
    scale_color_manual(values = c("#DC3220", "#005AB5", "grey"),
                       labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank()) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable",
         x = "P(hv) using full array datasets",
         y = paste0("P(hv) using reduced (", N, " ind/ds) array datasets")) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))
}

load(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/results_Arrays_3indperds_394240CpGs_0_8p0_0_65p1.RData"))

resArray3ind <- as.data.frame(results_Arrays_3indperds_394240CpGs_0_8p0_0_65p1)
rm(results_Arrays_3indperds_394240CpGs_0_8p0_0_65p1)
resArray3ind <- prepareChrDataset(resArray3ind)
names(resArray3ind)[names(resArray3ind) %in% "alpha"] <- "alpha_array_reduce"
resCompArray_allvs3 <- dplyr::left_join(resArray3ind, resArrayAll)
names(resCompArray_allvs3)[names(resCompArray_allvs3) %in% "alpha"] <- "alpha_array_all"
p3ind <- makePlotNrob(resCompArray_allvs3, 3)

resArray2ind <- as.data.frame(readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/results_Arrays_2indperds_394240CpGs_0_8p0_0_65p1.rds")))
resArray2ind <- prepareChrDataset(resArray2ind)
names(resArray2ind)[names(resArray2ind) %in% "alpha"] <- "alpha_array_reduce"
resCompArray_allvs2 <- dplyr::left_join(resArray2ind, resArrayAll)
names(resCompArray_allvs2)[names(resCompArray_allvs2) %in% "alpha"] <- "alpha_array_all"

p2ind <- makePlotNrob(resCompArray_allvs2, 2)

## From script S01:
venn_obj <- readRDS(here("B_MultiTissues/dataOut/figures/arrayCutoffLowPower2or3ind.RDS"))

plot_venn3 <- makeVennArrayReduced(df_circles = venn_obj$df_circles, v = venn_obj$v, counts = venn_obj$counts,
                                   fmt_fn = function(n, tot) paste0(n, "\n(", round(100*n/tot, 1), "%)"))
plot_venn3

###############################
## Make figure of array test ##
###############################

figure2 <- cowplot::plot_grid(
  cowplot::plot_grid(p1_manhattanArray, p2_DiffProbhvCpG_matchingcontrol_Array, ncol = 2, rel_widths = c(2,1), labels = c("A", "B")), 
  cowplot::plot_grid(plot_venn3,
                     p2ind + theme(legend.position = "none"),
                     p3ind + theme(legend.position = "none"),
                     cowplot::get_legend(p3ind), 
                     ncol = 4, rel_widths = c(1,1,1,.3), labels = c("C", "D", "E")), nrow = 2)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/Figure2.png"),
  plot = figure2, width = 18, height = 12,
  dpi = 300,        # 300 DPI = standard publication quality
  bg = "white"
)

#############################################################
## What cutoff to get the same number of sites than Maria? ##
#############################################################
table(resCompArray$group)
# hvCpG_Derakhshan     mQTLcontrols 
# 3535             3359

poscutoff = 0.94
negcutoff = 0.5

nrow(resCompArray[resCompArray$alpha_array_all > poscutoff,])# == 3674

## True positive: detected with full AND reduced array / all detected with full array
nrow(resCompArray[resCompArray$alpha_array_all > poscutoff &
               resCompArray$alpha_array_reduce > poscutoff,]) /
  nrow(resCompArray[resCompArray$alpha_array_all > poscutoff,]) * 100 #  49.21067

## False positive = % CpGs detected only using reduced datasets
nrow(resCompArray[resCompArray$alpha_array_all < negcutoff &
                    resCompArray$alpha_array_reduce > poscutoff,]) /
  nrow(resCompArray[resCompArray$alpha_array_reduce > poscutoff,]) * 100 # 0.4022526

## False negative = % CpGs detected only using full datasets
nrow(resCompArray[resCompArray$alpha_array_all > poscutoff &
                    resCompArray$alpha_array_reduce < negcutoff,]) /
  nrow(resCompArray[resCompArray$alpha_array_all > poscutoff,]) * 100 # 3.320631

## rm junk
rm(x,y, pairs, merged, chr_mid, hv_alpha, data, ctrl_alpha, resArray3ind, resArrayAll)