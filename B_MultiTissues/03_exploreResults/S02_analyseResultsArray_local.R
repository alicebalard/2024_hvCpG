## This script does:
## Load full results on array ##
## Calculate proba hvCpG minus matching control: is it always +? ##
## Load full results on array with only 3 individuals/ds ##

#####################################################################
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
#####################################################################

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

## Save for next scripts
saveRDS(resArrayAll, here("B_MultiTissues/dataOut/resArray.RDS"))

# Plot
# Compute midpoints for chromosome labels
chr_mid <- resArrayAll %>%
  group_by(chr) %>%
  summarise(mid = (min(cum_pos) + max(cum_pos)) / 2)
## colorblind friendly
p1_manhattanArray <- ggplot() +
  geom_point(data = subset(resArrayAll, is.na(group)),
             aes(x = cum_pos, y = alpha), color = "gray", alpha = .4, size = .8) +
  geom_point(data = subset(resArrayAll, !is.na(group)),
             aes(x = cum_pos, y = alpha, col = group), alpha = .6, size = 1) +
  scale_color_manual(values = c("hvCpG_Derakhshan" = "#DC3220",
                                "mQTLcontrols" = "#005AB5"),
                     labels = c("hvCpG_Derakhshan" = "Derakhshan hvCpG", "mQTLcontrols" = "mQTL controls")) +
  scale_x_continuous(breaks = chr_mid$mid,
                     labels = gsub("chr", "", chr_mid$chr), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(x = "Chromosome", y = "Pr(hv)")+
  theme(legend.position = "inside",
        legend.position.inside      = c(.8, 1.2),   # above the plot area
        legend.justification.inside = c(1, 1),
        plot.margin = margin(t = 40, r = 5, b = 5, l = 5),  # space for legend above
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linewidth = 0.5, linetype = "blank", colour = "black"))

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
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.3)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 11)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  ylab("Difference of Pr(hv)") +
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
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation:\n r = %.2f\n", mycor)) +
    labs(title = element_blank(),
         x = "P(hv) using full datasets",
         y = paste0("Pr(hv) using reduced (", N, " ind/ds) datasets")) +
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
x <- readRDS(here("B_MultiTissues/dataOut/figures/arrayCutoffLowPower2or3ind.RDS"))

library(sf)

plot_venn3 <- ggVennDiagram(x, label = "both", label_alpha = 0, 
                            label_color = "white", category.names = 
                              c("Full \ndatasets","2 ind per \ndataset", "3 ind per dataset")) +
  scale_fill_gradient(low = "grey", high = "black") +
  theme(legend.position = "none")+
  coord_sf(xlim = c(-5, 10), ylim = c(-10, 5))

plot_venn3 # "Cutoff algorithm"

## What cutoff to get the same number of sites than Maria? with Bayesian approach ##
table(resCompArray_allvs3$group)
# hvCpG_Derakhshan     mQTLcontrols 
# 3535             3359

top3535 <- resCompArray_allvs3 |>
  dplyr::slice_max(order_by = alpha_array_all, n = 3535, with_ties = FALSE)

top3535_2 <- resCompArray_allvs2 |>
  dplyr::slice_max(order_by = alpha_array_all, n = 3535, with_ties = FALSE)

# sanity check
table(top3535$chrpos %in% top3535_2$chrpos)

min(top3535$alpha_array_all)

pos3ind <- top3535$chrpos[top3535$alpha_array_reduce >= min(top3535$alpha_array_all)]
pos2ind <- top3535_2$chrpos[top3535_2$alpha_array_reduce >= min(top3535_2$alpha_array_all)]

y <- list(
  "Full array"     = top3535_2$chrpos,
  "Array 2 ind/ds" = pos2ind,
  "Array 3 ind/ds" = pos3ind
)

plot_venn3_Bayes <- ggVennDiagram(y, label = "both", label_alpha = 0, 
                                  label_color = "white", category.names = 
                                    c("Full \ndatasets","2 ind per \ndataset", "3 ind per dataset")) +
  scale_fill_gradient(low = "grey", high = "black") +
  theme(legend.position = "none")+
  coord_sf(xlim = c(-5, 10), ylim = c(-10, 5))

plot_venn3_Bayes

###############################
## Make figure of array test ##
###############################

lab <- list(size = 14, x = 0.01, y = 0.99, hjust = 0, vjust = 1)
mg  <- theme(plot.margin = margin(15, 5, 5, 5))

row1 <- cowplot::plot_grid(
  p1_manhattanArray + theme(plot.margin = margin(40, 5, 5, 5)), #+
  # ggtitle("Detection of IIHV with both methods (red = cutoff, Pr(hv) in y = Bayesian)"),
  p2_DiffProbhvCpG_matchingcontrol_Array + theme(plot.margin = margin(50, 5, 5, 5)),
  ncol = 2, rel_widths = c(4, 1),
  labels = c("A. Detection of IIHV with both methods (red = cutoff, Pr(hv) in y = Bayesian)", 
             "B. Difference between Pr(hv) of \nDerakhshan hvCpGs and \nPr(hv) of matching controls"),
  label_size = lab$size, label_x = lab$x, label_y = lab$y,
  hjust = lab$hjust, vjust = lab$vjust)

row2_1 <- cowplot::plot_grid(
  plot_venn3      + theme_void(base_size = 10) + theme(legend.position = "none", plot.margin = margin(40, 5, 5, 5)),
  plot_venn3_Bayes + theme_void(base_size = 10) + theme(legend.position = "none", plot.margin = margin(40, 5, 5, 5)),
  labels = c("C. Detection of highly variable CpGs \nwith reduced datasets (cutoff)", 
             "D. Detection of highly variable CpGs \nwith reduced datasets (Bayesian)"), nrow = 1,
  label_size = lab$size, label_x = lab$x, label_y = lab$y,
  hjust = lab$hjust, vjust = lab$vjust)

row2_2 <- cowplot::plot_grid(
  p2ind + theme_minimal(base_size = 11) + theme(legend.position = "none") + mg,
  p3ind + theme_minimal(base_size = 11) + theme(legend.position = "none") + mg,
  labels = c("E. Bayesian: full vs 2 ind/ds", "F. Bayesian: full vs 3 ind/ds"), nrow = 1,
  label_size = lab$size, label_x = lab$x, label_y = lab$y,
  hjust = lab$hjust, vjust = lab$vjust)

row2 <- cowplot::plot_grid(row2_1, row2_2)

figure2 <- cowplot::plot_grid(row1, row2, nrow = 2)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/Figure2.png"),
  plot = figure2, width = 16, height = 8,
  dpi = 300, bg = "white")

rm(x,y, pairs, merged, chr_mid, hv_alpha, data, ctrl_alpha, resArray3ind, resArrayAll)