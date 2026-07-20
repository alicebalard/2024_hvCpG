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

## with sum per dataset:
pathSum <- "B_MultiTissues/resultsDir_gitIgnored/Arrays/prev_sumInd/"

## with mean per dataset:
pathMean <- "B_MultiTissues/resultsDir_gitIgnored/Arrays/"

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

resArrayAll_sum <- as.data.frame(
  readRDS(here(paste0(pathSum, "results_Arrays_all_406036CpGs_0_8p0_0_65p1.rds"))))
resArrayAll_sum <- prepareChrDataset(resArrayAll_sum)

resArrayAll_sum_strict <- as.data.frame(
  readRDS(here(paste0(pathSum, "results_Arrays_all_406036CpGs_0_8p0_0_9p1.rds"))))
resArrayAll_sum_strict <- prepareChrDataset(resArrayAll_sum_strict)

resArrayAll_mean <- as.data.frame(
  readRDS(here(paste0(pathMean, "results_Arrays_all_406036CpGs_0_8p0_0_65p1.rds"))))
resArrayAll_mean <- prepareChrDataset(resArrayAll_mean)

resArrayAll_mean_strict <- as.data.frame(
  readRDS(here(paste0(pathMean, "results_Arrays_all_406036CpGs_0_8p0_0_9p1.rds"))))
resArrayAll_mean_strict <- prepareChrDataset(resArrayAll_mean_strict)

######################################
## NB: save the best for later scripts
## Save for next scripts
saveRDS(resArrayAll_mean_strict, here("B_MultiTissues/dataOut/resArraymean0.9p10.8p0.RDS"))

################################################################################
## Compare with previous buggy version (fix 13th July 2026)
resArrayAll_prev <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/prev_sumInd/resArray_beforep0p1bugcorrection.RDS"))
resArrayAll_prev <- resArrayAll_prev[c("alpha", "chrpos")]
names(resArrayAll_prev) <- c("alpha_beforeBugp0p1", "chrpos")

compare <- merge(resArrayAll_prev, resArrayAll_sum, all =T)
ggplot(compare, aes(x = alpha_beforeBugp0p1, y = alpha)) +
  geom_abline(slope = 1) +
  geom_point(alpha=.01)+
  theme_bw()
################################################################################

################################################################################
## Compare with different p0 p1

# separate object just for the p0/p1 scatter comparison - keep resArrayAllStrict intact for makeScript2Fig()
resArrayAll_sum_strict_cmp <- resArrayAll_sum_strict[c("alpha", "chrpos")]
names(resArrayAll_sum_strict_cmp) <- c("alpha_stricterp0p1", "chrpos")

# library(data.table)
setDT(resArrayAll_sum_strict_cmp); setDT(resArrayAll_sum)
setkey(resArrayAll_sum_strict_cmp, chrpos)
setkey(resArrayAll_sum, chrpos)

compare <- resArrayAll_sum_strict_cmp[resArrayAll_sum, nomatch = 0]
ggplot(compare, aes(x = alpha, y = alpha_stricterp0p1, colour = group)) +
  geom_point(data = compare[is.na(compare$group),], alpha = .1) + 
  geom_point(data = compare[!is.na(compare$group),], alpha = .1) + 
  geom_abline(slope = 1) + theme_bw()
################################################################################

################################################################################
## Compare between SUM over individuals or MEAN over individuals  
a = resArrayAll_sum_strict
b = resArrayAll_mean_strict
names(b)[names(b) == "alpha"] <- "alpha_mean"

comp <- dplyr::left_join(a, b[,c("chrpos", "alpha_mean")], by = "chrpos")

ggplot(comp, aes(x = alpha, y = alpha_mean)) +
  geom_point(size = 3, colour = "black", alpha = .01) +
  theme_minimal(base_size = 12)

cor(comp$alpha, comp$alpha_mean)

# rm(a,b,comp)
################################################################################


makeScript2Fig <- function(resArray, path, p0p1 = "0_8p0_0_65p1"){
  # Plot
  # Compute midpoints for chromosome labels
  chr_mid <- resArray %>%
    group_by(chr) %>%
    summarise(mid = (min(cum_pos) + max(cum_pos)) / 2)
  ## colorblind friendly
  p1_manhattanArray <- ggplot() +
    geom_point(data = subset(resArray, is.na(group)),
               aes(x = cum_pos, y = alpha), color = "gray", alpha = .4, size = .8) +
    geom_point(data = subset(resArray, !is.na(group)),
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
          legend.position.inside      = c(0.9, 1.15),   # above the plot area
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
  hv_alpha <- resArray[, c("chrpos", "alpha")]
  colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")
  
  # Merge control alphas
  ctrl_alpha <- resArray[, c("chrpos", "alpha")]
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
                 alpha = 0.01) +
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
  
  resArray3ind <- as.data.frame(readRDS(here(paste0(
    path, "results_Arrays_3ind_406036CpGs_", p0p1, ".rds"))))
  resArray3ind <- prepareChrDataset(resArray3ind)
  names(resArray3ind)[names(resArray3ind) %in% "alpha"] <- "alpha_array_reduce"
  
  resCompArray_allvs3 <- dplyr::left_join(
    resArray3ind,
    resArray[, c("chrpos", "alpha")],
    by = "chrpos"
  )
  names(resCompArray_allvs3)[names(resCompArray_allvs3) %in% "alpha"] <- "alpha_array_all"
  
  p3ind <- makePlotNrob(resCompArray_allvs3, 3)
  
  resArray2ind <- as.data.frame(readRDS(here(paste0(
    path, "results_Arrays_2ind_406036CpGs_", p0p1, ".rds"))))
  resArray2ind <- prepareChrDataset(resArray2ind)
  names(resArray2ind)[names(resArray2ind) %in% "alpha"] <- "alpha_array_reduce"
  
  resCompArray_allvs2 <- dplyr::left_join(
    resArray2ind,
    resArray[, c("chrpos", "alpha")],
    by = "chrpos"
  )
  names(resCompArray_allvs2)[names(resCompArray_allvs2) %in% "alpha"] <- "alpha_array_all"
  
  p2ind <- makePlotNrob(resCompArray_allvs2, 2)
  
  ## From script S01:
  x <- readRDS(here("B_MultiTissues/dataOut/arrayCutoffLowPower2or3ind.RDS"))
  
  plot_venn3 <- ggVennDiagram(x, label = "both", label_alpha = 0,
                              label_color = "white",
                              category.names = c("Full \ndatasets",
                                                 "2 ind per \ndataset",
                                                 "3 ind per dataset")) +
    scale_fill_gradient(low = "grey", high = "black") +
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(-5, 10), ylim = c(-10, 5), clip = "off")
  
  plot_venn3 # "Cutoff algorithm"
  
  message("What cutoff to get the same number of sites than Maria? with Bayesian approach")
  print(table(resCompArray_allvs3$group))
  
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
    coord_cartesian(xlim = c(-5, 10), ylim = c(-10, 5), clip = "off")
  
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
    labels = c("A. Detection of highly variable CpGs with both methods (red = cutoff, Pr(hv) in y = Bayesian)", 
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
  
  return(figure2)
}

fig2_p080p165 <- makeScript2Fig(resArrayAll_sum, pathSum, p0p1 = "0_8p0_0_65p1")

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/testOnArray_p080p165.png"),
  plot = fig2_p080p165, width = 18, height = 10,
  dpi = 300, bg = "white")

fig2_p080p190 <- makeScript2Fig(resArrayAll_sum_strict, pathSum, p0p1 = "0_8p0_0_9p1")

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/testOnArray_p080p190.png"),
  plot = fig2_p080p190, width = 18, height = 10,
  dpi = 300, bg = "white")

################################################################################
## Evolution of top-3535 CpG recovery as individuals/dataset increases        ##
## (p0=80%, p1=90%): 2, 3, 4, 5, 6, 10 ind/ds vs full datasets                ##
################################################################################

Ns <- c(2, 3, 4, 5, 6, 10, 15, 20)

## Baseline: top 3535 CpGs using the full datasets (same cutoff as before)
top3535_full <- resArrayAll_sum_strict |> dplyr::slice_max(order_by = alpha, n = 3535, with_ties = FALSE)
cutoff_full  <- min(top3535_full$alpha)

recovery <- lapply(Ns, function(N) {
  f <- here(paste0(pathSum, "results_Arrays_", N,
                   "ind_406036CpGs_0_8p0_0_9p1.rds"))
  if (!file.exists(f)) {
    message("SKIP N=", N, " - file not found (not computed yet?): ", f)
    return(NULL)
  }
  
  resN <- as.data.frame(readRDS(f)) |> prepareChrDataset()
  names(resN)[names(resN) == "alpha"] <- "alpha_reduced"
  
  comp <- dplyr::left_join(resN, resArrayAll_sum_strict[, c("chrpos", "alpha")], by = "chrpos")
  
  ## Of the top 3535 CpGs under full data, how many are still >= the same
  ## cutoff when using only N individuals/dataset?
  recovered <- sum(comp$chrpos %in% top3535_full$chrpos & comp$alpha_reduced >= cutoff_full)
  
  data.frame(N = N, n_recovered = recovered, pct_recovered = recovered / nrow(top3535_full))
}) |> dplyr::bind_rows()

## Add the full-dataset reference point (100% recovery, by definition)
recovery <- dplyr::bind_rows(recovery, data.frame(N = NA, n_recovered = nrow(top3535_full), pct_recovered = 1))
recovery$N_label <- ifelse(is.na(recovery$N), "Full", as.character(recovery$N))
recovery$N_label <- factor(recovery$N_label, levels = c(as.character(sort(Ns)), "Full"))

print(recovery)

plotRecovery <- ggplot(recovery, aes(x = N_label, y = pct_recovered, group = 1)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "black") +
  geom_text(aes(label = scales::percent(pct_recovered, accuracy = 1)), vjust = -1, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  theme_minimal(base_size = 13) +
  labs(x = "Individuals per dataset", y = "% of top 3,535 CpGs (full data) still recovered",
       title = "Recovery of top hvCpG hits as individuals/dataset increases",
       subtitle = "p0=80%, p1=90% - cutoff fixed at the full-dataset top-3,535 threshold")

# is the correlation trend smoother?
corByN <- lapply(Ns, function(N) {
  f <- here(paste0(pathSum, "results_Arrays_", N,
                   "ind_406036CpGs_0_8p0_0_9p1.rds"))
  if (!file.exists(f)) return(NULL)
  resN <- as.data.frame(readRDS(f)) |> prepareChrDataset()
  names(resN)[names(resN) == "alpha"] <- "alpha_reduced"
  comp <- dplyr::left_join(resN, resArrayAll_sum_strict[, c("chrpos", "alpha")], by = "chrpos")
  data.frame(N = N, r = cor(comp$alpha, comp$alpha_reduced, use = "complete.obs"))
}) |> dplyr::bind_rows()

print(corByN)

plotCorByN <- ggplot(corByN, aes(x = N, y = r)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "black") +
  geom_text(aes(label = scales::percent(r, accuracy = 1)), vjust = -1, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(.5, 1.05)) +
  theme_minimal(base_size = 12) +
  labs(x = "Individuals per dataset", y = "Correlation (Pearson r) of Pr(hv) with full data",
       title = "Agreement with full-data Pr(hv) improves with sample size",
       subtitle = "p0=80%, p1=90%")

plotCorByN

ggsave(here("B_MultiTissues/dataOut/figures/script02/plotcorByN_vs_Nind_p080p190.png"),
       plotCorByN, width = 8, height = 6, dpi = 300, bg = "white")







sd0 <- 0.015
sd1 <- 0.05    # ratio = 3.33 - the actual gap the model is choosing between

Ns <- c(10, 20, 50, 100, 150, 200, 300, 400)
ratios <- c(1.0, 1.5, 2.5, 3.33)   # 3.33 = exactly matches sd1 (genuinely "hv"-level spread)
nReps <- 200

set.seed(1)
marginSim <- expand.grid(N = Ns, trueSDratio = ratios, rep = seq_len(nReps)) %>%
  rowwise() %>%
  mutate(margin = {
    x <- rnorm(N, 0, sd0 * trueSDratio)
    sum(dnorm(x, mean(x), sd1, log = TRUE)) - sum(dnorm(x, mean(x), sd0, log = TRUE))
  }) %>%
  ungroup() %>%
  group_by(N, trueSDratio) %>%
  summarise(mean_margin = mean(margin), se = sd(margin) / sqrt(n()), .groups = "drop")

ggplot(marginSim, aes(x = N, y = mean_margin, colour = factor(trueSDratio))) +
  geom_ribbon(aes(ymin = mean_margin - 2*se, ymax = mean_margin + 2*se, fill = factor(trueSDratio)),
              alpha = .15, colour = NA) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  theme_minimal(base_size = 13) +
  labs(x = "N individuals", y = "Mean log-lik margin (hv - typical) over 200 reps",
       colour = "True SD\n(× typical)", fill = "True SD\n(× typical)",
       title = "N amplifies whatever the data actually supports, in either direction",
       subtitle = paste0("sd0=", sd0, ", sd1=", sd1, " (", round(sd1/sd0, 2),
                         "x gap) - only elevations near that gap trend positive with N"))

## idea:
# tissueRes <- readRDS(here("path/to/results_bytissue_Array_tissue0.5_....rds"))  # CpG x dataset matrix
# 
# supportRule <- rowMeans(tissueRes > log(0.5), na.rm = TRUE) >= (2/3)  # adjust threshold as appropriate




################################################################################
## hvCpGs missed by the previous cutoff method but now confidently flagged    ##
## by the Bayesian method (p0=80%, p1=90%)                                    ##
################################################################################

alphaThresh <- 0.99   # "close or equal to 1"

## "Previously detected by cutoff" = Derakhshan hvCpG_Derakhshan group
## "Not detected by cutoff" = background OR matched mQTL controls (both were
## NOT flagged by the old method) - now newly confident under the Bayesian model
resArrayAll_sum_strict <- resArrayAll_sum_strict %>%
  mutate(discoveryGroup = dplyr::case_when(
    group == "hvCpG_Derakhshan" ~ "known (cutoff-detected)",
    (is.na(group) | group == "mQTLcontrols") & alpha >= alphaThresh ~ "new discovery (alpha>=0.99)",
    TRUE ~ "other"
  ))

table(resArrayAll_sum_strict$discoveryGroup)
# known (cutoff-detected) new discovery (alpha>=0.99)     other 
# 3644                       10288                        392084 

newDiscovery <- resArrayAll_sum_strict %>% dplyr::filter(discoveryGroup == "new discovery (alpha>=0.99)")
knownCutoff  <- resArrayAll_sum_strict %>% dplyr::filter(discoveryGroup == "known (cutoff-detected)")

## prepareChrDataset() drops the original cg IDs - recover them via dico for the h5 lookup
newDiscovery$cpg <- dico$CpG[match(newDiscovery$chrpos, dico$chrpos_hg38)]
knownCutoff$cpg  <- dico$CpG[match(knownCutoff$chrpos, dico$chrpos_hg38)]

message(nrow(newDiscovery), " CpGs newly confident (alpha>=", alphaThresh,
        ") that the previous cutoff method never flagged.")
# 10288 CpGs newly confident (alpha>=0.99) that the previous cutoff method never flagged.

################################################################################
## Are they backed by evidence as broad/strong as the known hvCpGs, or        ##
##    weaker/narrower? Per-dataset sd0-vs-sd1 test generalised                ##
################################################################################

h5file     <- "/home/alice/arraysh5/all_matrix_noscale.h5"
metaFile   <- "/home/alice/arraysh5/sample_metadata.tsv"
lambdaFile <- "/home/alice/arraysh5/all_medsd_lambda.tsv"

metadata  <- read.table(metaFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ds_params <- read.table(lambdaFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(sd0 = pmax(median_sd, 1e-4), sd1 = pmax(lambda * median_sd, 1e-4))

cpg_names_all <- rhdf5::h5read(h5file, "cpg_names")
samples       <- rhdf5::h5read(h5file, "samples")

computePerDatasetSummary <- function(cpgs) {
  rowIdx <- match(cpgs, cpg_names_all)
  keep <- !is.na(rowIdx)
  if (any(!keep)) message(sum(!keep), " CpGs not found in h5, dropped.")
  cpgs <- cpgs[keep]; rowIdx <- rowIdx[keep]
  
  rawMat <- rhdf5::h5read(h5file, "matrix", index = list(rowIdx, NULL), native = TRUE)
  rownames(rawMat) <- cpgs
  colnames(rawMat) <- samples
  
  lapply(cpgs, function(cpg) {
    d <- data.frame(sample = samples, beta = as.numeric(rawMat[cpg, ])) %>%
      inner_join(metadata, by = "sample") %>%
      filter(!is.na(beta)) %>%
      inner_join(ds_params[, c("dataset", "sd0", "sd1")], by = "dataset") %>%
      group_by(dataset, sd0, sd1) %>%
      summarise(n = n(), sd_obs = sd(beta),
                ll0 = sum(dnorm(beta, mean(beta), sd0[1], log = TRUE)),
                ll1 = sum(dnorm(beta, mean(beta), sd1[1], log = TRUE)),
                .groups = "drop") %>%
      mutate(z_hv = ll1 > ll0, margin = ll1 - ll0)
    
    data.frame(cpg = cpg,
               n_datasets = nrow(d),
               prop_datasets_hv = mean(d$z_hv),
               prop_individuals_hv = sum(d$n[d$z_hv]) / sum(d$n),
               total_ll_margin = sum(d$margin))
  }) %>% bind_rows()
}

## Known-cutoff set can be large (~3600) - sample for a fair, fast comparison
set.seed(1)
knownCutoff_sample <- knownCutoff %>% slice_sample(n = min(300, nrow(knownCutoff)))

newDiscovery_sample <- newDiscovery %>% slice_sample(n = min(300, nrow(newDiscovery)))

message("Computing per-dataset support for ", nrow(newDiscovery_sample), " new-discovery CpGs...")
newDiscovery_support <- computePerDatasetSummary(newDiscovery_sample$cpg)

message("Computing per-dataset support for ", nrow(knownCutoff_sample), " known cutoff-detected CpGs (sampled)...")
knownCutoff_support <- computePerDatasetSummary(knownCutoff_sample$cpg)

supportComparison <- bind_rows(
  newDiscovery_support %>% mutate(discoveryGroup = "new discovery"),
  knownCutoff_support  %>% mutate(discoveryGroup = "known (cutoff-detected)")
)

plotSupportCompare <- supportComparison %>%
  pivot_longer(c(n_datasets, prop_datasets_hv, prop_individuals_hv, total_ll_margin),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         n_datasets = "N datasets contributing",
                         prop_datasets_hv = "Prop. datasets 'looks hv'",
                         prop_individuals_hv = "Prop. individuals 'looks hv' (weighted)",
                         total_ll_margin = "Total log-lik margin (hv vs typical)"))

ggplot(plotSupportCompare, aes(x = discoveryGroup, y = value, fill = discoveryGroup)) +
  geom_violin()+
  geom_boxplot(outlier.size = .5, alpha = .7, width = .3) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = NULL,
       title = "Do new-discovery CpGs have evidence as strong/broad as known hvCpGs?",
       subtitle = "p0=80%, p1=90%")

# ggsave(here("B_MultiTissues/dataOut/figures/script02/newDiscovery_vs_known_support.png"),
#        width = 10, height = 7, dpi = 300, bg = "white")

## Use the known-cutoff group's own 5th percentile as the "normal support" floor
weakThresh <- quantile(knownCutoff_support$prop_datasets_hv, .05, na.rm = TRUE)
message("Known hvCpGs' 5th percentile of prop_datasets_hv = ", round(weakThresh, 2))

newDiscovery_support <- newDiscovery_support %>%
  mutate(supportTier = ifelse(prop_datasets_hv < weakThresh, "narrow support", "broad support"))

table(newDiscovery_support$supportTier)
prop.table(table(newDiscovery_support$supportTier))

# pull the 10 lowest prop_datasets_hv new discoveries and plot their per-dataset 
# spread directly, so we can see whether each is a genuine "driven by a few strongly-spread
# large datasets" case (defensible) or something noisier:
narrowSupport <- newDiscovery_support %>%
  dplyr::filter(supportTier == "narrow support") %>%
  arrange(prop_datasets_hv) %>%
  slice_head(n = 10) %>%
  left_join(newDiscovery %>% dplyr::select(cpg, chrpos, alpha), by = "cpg")

print(narrowSupport %>% dplyr::select(cpg, chrpos, alpha, n_datasets, prop_datasets_hv, prop_individuals_hv, total_ll_margin))

rowIdx <- match(narrowSupport$cpg, cpg_names_all)
rawMat <- rhdf5::h5read(h5file, "matrix", index = list(rowIdx, NULL), native = TRUE)
rownames(rawMat) <- narrowSupport$cpg
colnames(rawMat) <- samples

narrowSupport <- narrowSupport %>%
  mutate(cpg_label = sprintf("%s\nalpha=%.2f | prop_datasets_hv=%.2f", cpg, alpha, prop_datasets_hv))

rawLongNarrow <- as.data.frame(rawMat) %>%
  tibble::rownames_to_column("cpg") %>%
  pivot_longer(-cpg, names_to = "sample", values_to = "beta") %>%
  dplyr::filter(!is.na(beta)) %>%
  left_join(narrowSupport %>% dplyr::select(cpg, cpg_label, prop_datasets_hv), by = "cpg") %>%
  mutate(dataset = metadata$dataset[match(sample, metadata$sample)])

plotNarrowSupport <- ggplot(rawLongNarrow, aes(x = dataset, y = beta)) +
  geom_violin() +
  geom_boxplot(outlier.size = .4, width = .2) +
  geom_jitter(width = .15, size = .3, alpha = .4) +
  facet_wrap(~ reorder(cpg_label, prop_datasets_hv), scales = "free_x", ncol = 2) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Dataset (source tissue / cell type)", y = "Beta value",
       title = "Raw data for the 10 weakest-supported new-discovery CpGs",
       subtitle = "alpha≈1, but only a minority of datasets actually look hypervariable")

ggsave(here("B_MultiTissues/dataOut/figures/script02/narrowSupportNewDiscoveryRaw.png"),
       plotNarrowSupport, width = 13, height = 14, dpi = 300, bg = "white")







################################################################################
## Which mQTL controls flipped from "not flagged" (buggy code, p0=80%, p1=65%)
## to "flagged" (fixed code, p0=80%, p1=90%) -- and why?
##
## The comparison you asked for conflates two separate changes:
##   1. the bug fix (marginalisation bug, fixed 13 July)
##   2. the parameter change (p1: 65% -> 90%)
## Since we also have the FIXED code run at the OLD p1=65% on disk
## (results_Arrays_all_406036CpGs_0_8p0_0_65p1.rds), we can add it as a third
## reference point and split "why" into "how much is the bug fix" vs
## "how much is the stricter p1", instead of only seeing the combined jump.
################################################################################

## Adjust to wherever you're running this (LSHTM server vs ing-p5 vs local copy)
h5file <- "/home/alice/arraysh5/all_matrix_noscale.h5"

###############################################################
## 1. Three alpha estimates for the same CpGs                ##
###############################################################

## a) buggy code, p0=80%, p1=65%  (what you had before the fix)
prev <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/resArray_beforep0p1bugcorrection.RDS"))
prev <- prev[c("alpha_buggy_65", "chrpos")]
names(prev) <- c("alpha_buggy_65", "chrpos")

## b) fixed code, SAME params p0=80%, p1=65% -> isolates the bug-fix effect alone
fixed65 <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/results_Arrays_all_406036CpGs_0_8p0_0_65p1.rds"))
fixed65 <- data.frame(cpg = rownames(fixed65), alpha_fixed_65 = as.numeric(fixed65[, "alpha"]))
fixed65$chrpos <- dico$chrpos_hg38[match(fixed65$cpg, dico$CpG)]

## c) fixed code, stricter params p0=80%, p1=90% -> adds the parameter-change effect
fixed90 <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Arrays/results_Arrays_all_406036CpGs_0_8p0_0_9p1.rds"))
fixed90 <- data.frame(cpg = rownames(fixed90), alpha_fixed_90 = as.numeric(fixed90[, "alpha"]))
fixed90$chrpos <- dico$chrpos_hg38[match(fixed90$cpg, dico$CpG)]

merged <- prev %>%
  inner_join(fixed65[, c("chrpos", "alpha_fixed_65")], by = "chrpos") %>%
  inner_join(fixed90, by = "chrpos") %>%                 # brings in cpg (cg ID) + alpha_fixed_90
  dplyr::filter(chrpos %in% mQTLcontrols_hg38)

################################################################
## 2. Controls that flipped: low before, high now             ##
################################################################

flipped <- merged %>%
  filter(alpha_buggy_65 < 0.70, alpha_fixed_90 > 0.90) %>%    # adjust thresholds as needed
  mutate(jump = alpha_fixed_90 - alpha_buggy_65) %>%
  arrange(desc(jump)) %>%
  slice_head(n = 10)

flipped %>%
  mutate(due_to_bugfix = alpha_fixed_65 - alpha_buggy_65,      # buggy p1=65 -> fixed p1=65
         due_to_p1     = alpha_fixed_90 - alpha_fixed_65) %>%  # fixed p1=65 -> fixed p1=90
  dplyr::select(cpg, chrpos, alpha_buggy_65, alpha_fixed_65, alpha_fixed_90, due_to_bugfix, due_to_p1) 

###############################################################
## 3. Pull raw data for these flipped CpGs (same as S14)     ##
###############################################################

cpg_names_all <- rhdf5::h5read(h5file, "cpg_names")
samples       <- rhdf5::h5read(h5file, "samples")
sample_groups <- rhdf5::h5read(h5file, "sample_groups")

rowIdx <- match(flipped$cpg, cpg_names_all)
stopifnot(!anyNA(rowIdx))

rawMat <- rhdf5::h5read(h5file, "matrix", index = list(rowIdx, NULL), native = TRUE)
rownames(rawMat) <- flipped$cpg
colnames(rawMat) <- samples

flipped <- flipped %>%
  mutate(cpg_label = sprintf("%s\nbuggy p1=.65: %.2f | fixed p1=.65: %.2f | fixed p1=.90: %.2f",
                             cpg, alpha_buggy_65, alpha_fixed_65, alpha_fixed_90))

rawLong <- as.data.frame(rawMat) %>%
  tibble::rownames_to_column("cpg") %>%
  pivot_longer(-cpg, names_to = "sample", values_to = "beta") %>%
  filter(!is.na(beta)) %>%
  left_join(flipped[, c("cpg", "cpg_label", "jump")], by = "cpg") %>%
  mutate(dataset = sample_groups[match(sample, samples)])

#######################
## 4. Plot           ##
#######################

plotFlipped <- ggplot(rawLong, aes(x = dataset, y = beta)) +
  geom_boxplot(outlier.size = .4, width = .6) +
  geom_jitter(width = .15, size = .3, alpha = .4) +
  facet_wrap(~ reorder(cpg_label, -jump), scales = "free_x", ncol = 2) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Dataset (source tissue / cell type)", y = "Beta value",
       title = "Controls flipped from low (buggy, p0=80%,p1=65%) to high (fixed, p0=80%,p1=90%) alpha",
       subtitle = "Facet labels show alpha under all three settings, to separate the bug-fix effect from the p1 effect")

ggsave(here("B_MultiTissues/dataOut/figures/script02/flippedControlsRaw.png"),
       plotFlipped, width = 13, height = 14, dpi = 300, bg = "white")

################################################################################
## For the flipped CpGs, test PER DATASET whether the observed
## spread looks more like the algorithm's "typical" (sd0) or "hypervariable"
## (sd1) Gaussian -- using the SAME dataset-level parameters hyperVarMeth
## itself uses (all_medsd_lambda.tsv), not an eyeballed boxplot.
##
## NB: p1 is the assumed SENSITIVITY (if truly hv, ~p1 of datasets should
## look spread) -- it is not a rule that alpha can only be high if >=p1 of
## datasets look spread. alpha comes from a marginal likelihood across
## datasets, so a few large, strongly-spread datasets can outweigh many
## small, tight-looking ones. This script reports BOTH the simple
## proportion-of-datasets metric AND an individual-weighted version, so you
## can see whether a high alpha is backed by a broad, consistent pattern or
## driven by a minority of high-N datasets.
################################################################################

## Adjust to wherever you're running this
h5file     <- "/home/alice/arraysh5/all_matrix_noscale.h5"
metaFile   <- "/home/alice/arraysh5/sample_metadata.tsv"
lambdaFile <- "/home/alice/arraysh5/all_medsd_lambda.tsv"

metadata <- read.table(metaFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ds_params <- read.table(lambdaFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(sd0 = pmax(median_sd, 1e-4),
         sd1 = pmax(lambda * median_sd, 1e-4))

###########################################
## 1. Raw values for the CpGs to check   ##
###########################################

cpg_names_all <- rhdf5::h5read(h5file, "cpg_names")
samples       <- rhdf5::h5read(h5file, "samples")

rowIdx <- match(flipped$cpg, cpg_names_all)
stopifnot(!anyNA(rowIdx))

rawMat <- rhdf5::h5read(h5file, "matrix", index = list(rowIdx, NULL), native = TRUE)
rownames(rawMat) <- flipped$cpg
colnames(rawMat) <- samples

##########################################################################
## 2. Per (CpG, dataset): observed spread vs that dataset's own sd0/sd1 ##
##########################################################################

perDatasetCall <- lapply(flipped$cpg, function(cpg) {
  data.frame(sample = samples, beta = as.numeric(rawMat[cpg, ])) %>%
    inner_join(metadata, by = "sample") %>%
    filter(!is.na(beta)) %>%
    inner_join(ds_params[, c("dataset", "sd0", "sd1")], by = "dataset") %>%
    group_by(dataset, sd0, sd1) %>%
    summarise(n = n(), mu = mean(beta), sd_obs = sd(beta),
              ll0 = sum(dnorm(beta, mean(beta), sd0[1], log = TRUE)),
              ll1 = sum(dnorm(beta, mean(beta), sd1[1], log = TRUE)),
              .groups = "drop") %>%
    mutate(cpg = cpg,
           z_hv = ll1 > ll0,                       # which Gaussian better explains THIS dataset
           sd_ratio_vs_typical = sd_obs / sd0)      # "more variable than average" for that dataset
}) %>% bind_rows()

## sd_ratio_vs_typical > 1   -> more variable than a typical CpG in that dataset
## sd_ratio_vs_typical > lambda (i.e. sd_obs > sd1) -> exceeds that dataset's own hv threshold
print(perDatasetCall %>% arrange(cpg, sd_ratio_vs_typical), n = 200)

###########################################################
## 3. Per-CpG: how consistently does it look hv, vs p1?  ##
###########################################################

summaryByCpG <- perDatasetCall %>%
  group_by(cpg) %>%
  dplyr::summarise(
    n_datasets           = n(),
    prop_datasets_hv      = mean(z_hv),                # simple proportion of datasets
    prop_individuals_hv   = sum(n[z_hv]) / sum(n),      # weighted by how many people back each call
    .groups = "drop"
  ) %>%
  left_join(flipped %>% 
              dplyr::select(cpg, alpha_buggy_65, alpha_fixed_65, alpha_fixed_90), by = "cpg")

print(summaryByCpG, n = 20)

################################################################
## 4. Plot: does the alpha=1 call hold up dataset-by-dataset? ##
################################################################

plotDF <- summaryByCpG %>%
  pivot_longer(c(prop_datasets_hv, prop_individuals_hv), names_to = "metric", values_to = "prop")

plotCheck <- ggplot(plotDF, aes(x = reorder(cpg, -prop), y = prop, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.90, linetype = "dashed", colour = "firebrick") +
  geom_hline(yintercept = 0.20, linetype = "dashed", colour = "steelblue") +
  annotate("text", x = Inf, y = 0.90, label = "p1 = 90% (expected if truly hv)",
           hjust = 1.05, vjust = -0.5, size = 3, colour = "firebrick") +
  annotate("text", x = Inf, y = 0.20, label = "1-p0 = 20% (expected false-positive rate)",
           hjust = 1.05, vjust = -0.5, size = 3, colour = "steelblue") +
  coord_flip() +
  scale_fill_manual(values = c(prop_datasets_hv = "grey40", prop_individuals_hv = "grey70"),
                    labels = c("prop. of datasets", "prop. of individuals (weighted)")) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank()) +
  labs(x = NULL, y = "Proportion classified 'looks hypervariable'",
       title = "Do the alpha≈1 calls hold up dataset-by-dataset?",
       subtitle = "Compared against the sensitivity/specificity assumed by the model (p1=90%, p0=80%)")

ggsave(here("B_MultiTissues/dataOut/figures/script02/perDatasetVariabilityCheck.png"),
       plotCheck, width = 9, height = 6, dpi = 300, bg = "white")

