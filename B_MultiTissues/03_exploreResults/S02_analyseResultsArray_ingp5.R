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

## with the new logic:
pathNew <- "B_MultiTissues/resultsDir_gitIgnored/Arrays/"

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

################################################################################
## old
resArrayAll_sum <- as.data.frame(
  readRDS(here(paste0(pathSum, "results_Arrays_all_406036CpGs_0_8p0_0_65p1.rds"))))
resArrayAll_sum <- prepareChrDataset(resArrayAll_sum)

resArrayAll_sum_strict <- as.data.frame(
  readRDS(here(paste0(pathSum, "results_Arrays_all_406036CpGs_0_8p0_0_9p1.rds"))))
resArrayAll_sum_strict <- prepareChrDataset(resArrayAll_sum_strict)

################################################################################
## New
resArray_0_8p0_0_8p1 <- as.data.frame(
  readRDS(here(paste0(pathNew, "results_Arrays_all_406036CpGs_0_8p0_0_8p1_0_01a0.rds"))))
resArray_0_8p0_0_8p1 <- prepareChrDataset(resArray_0_8p0_0_8p1)

resArray_0_8p0_0_65p1 <- as.data.frame(
  readRDS(here(paste0(pathNew, "results_Arrays_all_406036CpGs_0_8p0_0_65p1_0_01a0.rds"))))
resArray_0_8p0_0_65p1 <- prepareChrDataset(resArray_0_8p0_0_65p1)

resArray_0_9p0_0_9p1 <- as.data.frame(
  readRDS(here(paste0(pathNew, "results_Arrays_all_406036CpGs_0_9p0_0_9p1_0_01a0.rds"))))
resArray_0_9p0_0_9p1 <- prepareChrDataset(resArray_0_9p0_0_9p1)

resArray_0_65p0_0_8p1 <- as.data.frame(
  readRDS(here(paste0(pathNew, "results_Arrays_all_406036CpGs_0_65p0_0_8p1_0_01a0.rds"))))
resArray_0_65p0_0_8p1 <- prepareChrDataset(resArray_0_65p0_0_8p1)

resArray_0_65p0_0_65p1 <- as.data.frame(
  readRDS(here(paste0(pathNew, "results_Arrays_all_406036CpGs_0_65p0_0_65p1_0_01a0.rds"))))
resArray_0_65p0_0_65p1 <- prepareChrDataset(resArray_0_65p0_0_65p1)

################################################################################

getPlotComp <- function(A, B, column1 = "alpha", column2) {
  setDT(A); setDT(B)
  merged <- merge(A, B, all = FALSE)
  
  ggplot(merged, aes(x = .data[[column1]], y = .data[[column2]])) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth() +
    theme_minimal()
}

getPlotScores <- function(resArray){
  
  set.seed(1234)
  pa <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArrayAll_sum[sample(1:nrow(resArrayAll_sum),100000),], 
                    column2 = "post_hv")
  
  pb <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArrayAll_sum[sample(1:nrow(resArrayAll_sum),100000),], 
                    column2 = "logBF")
  
  pc <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArrayAll_sum[sample(1:nrow(resArrayAll_sum),100000),], 
                    column2 = "logBF_per_ds")
  
  pd <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArrayAll_sum[sample(1:nrow(resArrayAll_sum),100000),], 
                    column2 = "lambda_hat")
  
  pe <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArray[sample(1:nrow(resArray),100000),], 
                    column1 = "logBF", column2 = "logBF_per_ds")
  
  pf <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArray[sample(1:nrow(resArray),100000),], 
                    column1 = "logBF", column2 = "n_hv_ds")
  
  pg <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArray[sample(1:nrow(resArray),100000),], 
                    column1 = "logBF_per_ds", column2 = "n_hv_ds")
  
  ph <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArrayAll_sum[sample(1:nrow(resArrayAll_sum),100000),], 
                    column2 = "n_hv_ds")
  
  plot <- cowplot::plot_grid(pa, pb, pc, pd, pe, pf, pg, ph)
  return(plot)
}

getPlotScores(resArray = resArray_0_65p0_0_65p1)
getPlotScores(resArray = resArray_0_9p0_0_9p1)
getPlotScores(resArray = resArray_0_8p0_0_65p1)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/plotDifferentScores_0_8p0_0_65p1.png"),
  plot = getPlotScores(resArray = resArray_0_8p0_0_65p1), width = 18, height = 10,
  dpi = 300, bg = "white")

################################################################################
## Compare different p0 p1                                                    ##  
################################################################################

makePlotCompareP0P1 <- function(path1, path2, score = "logBF"){
  extract_params <- function(path) {
    pat <- "_(\\d+_\\d+)p0_(\\d+_\\d+)p1_(\\d+_\\d+)a0"
    m   <- regexec(pat, path, perl = TRUE)
    groups <- regmatches(path, m)[[1]]
    nums <- gsub("_", ".", groups[2:4])
    vals <- as.numeric(nums)
    setNames(vals, c("p0", "p1", "a0"))
  }
  
  listA <- extract_params(path1)
  listB <- extract_params(path2)
  
  x <- paste0("logBF_", listA[1], "p0_", listA[2], "p1")
  A <- as.data.frame(readRDS(here(paste0(pathNew, path1))))
  A <- prepareChrDataset(A)
  names(A)[names(A) %in% score] <- x
  
  y <- paste0(score, listB[1], "p0_", listB[2], "p1")
  B <- as.data.frame(readRDS(here(paste0(pathNew, path2))))
  B <- prepareChrDataset(B)
  B <- B[c(score, "chrpos")]
  names(B)[names(B) %in% score] <- y
  
  setDT(A); setDT(B)
  setkey(A, chrpos); setkey(B, chrpos)
  
  compare <- A[B, nomatch = 0]
  ggplot(compare, aes(x = .data[[x]], y = .data[[y]], colour = group)) +
    geom_point(data = compare[compare$group %in% "mQTLcontrols",], alpha = .1) +
    geom_point(data = compare[compare$group %in% "hvCpG_Derakhshan",], alpha = .1) +
    geom_abline(slope = 1) + theme_bw()
}

p1 <- makePlotCompareP0P1(path1 = "results_Arrays_all_406036CpGs_0_8p0_0_65p1_0_01a0.rds",
                          path2 = "results_Arrays_all_406036CpGs_0_8p0_0_8p1_0_01a0.rds",
                          score = "logBF")

p2 <- makePlotCompareP0P1(path1 = "results_Arrays_all_406036CpGs_0_65p0_0_8p1_0_01a0.rds",
                          path2 = "results_Arrays_all_406036CpGs_0_8p0_0_8p1_0_01a0.rds",
                          score = "logBF")

p3 <- makePlotCompareP0P1(path1 = "results_Arrays_all_406036CpGs_0_65p0_0_8p1_0_01a0.rds",
                          path2 = "results_Arrays_all_406036CpGs_0_65p0_0_65p1_0_01a0.rds",
                          score = "logBF")

p4 <- makePlotCompareP0P1(path1 = "results_Arrays_all_406036CpGs_0_8p0_0_8p1_0_01a0.rds",
                          path2 = "results_Arrays_all_406036CpGs_0_9p0_0_9p1_0_01a0.rds",
                          score = "logBF")

legend <- get_legend( p4 + theme(legend.box.margin = margin(0, 0, 0, 12)))

plotp0p1 <- cowplot::plot_grid(p1 + theme(legend.position = "none"),
                               p2+ theme(legend.position = "none"),
                               p3+ theme(legend.position = "none"),
                               p4+ theme(legend.position = "none"),
                               legend)
## Ranking is perfectly preserved whatever the p0 and p1
ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/plotp0p1.png"),
  plot = plotp0p1, width = 18, height = 10,
  dpi = 300, bg = "white")

######################################
## NB: save the best for later scripts
## Save for next scripts
saveRDS(resArray_0_8p0_0_65p1, here("B_MultiTissues/dataOut/resArray_0_8p0_0_65p1.RDS"))

################################################################################
makeScript2Fig <- function(resArray, path, p0p1, score = "logBF", shift = TRUE){
  
  if (shift){
    message("We shift the score so it starts at zero")
    resArray[[score]] <- resArray[[score]] + abs(min(resArray[[score]]))
  }
  
  # Plot
  # Compute midpoints for chromosome labels
  chr_mid <- resArray %>%
    group_by(chr) %>%
    summarise(mid = (min(cum_pos) + max(cum_pos)) / 2)
  ## colorblind friendly
  p1_manhattanArray <- ggplot() +
    geom_point(data = subset(resArray, is.na(group)),
               aes(x = cum_pos, y = .data[[score]]), color = "gray", alpha = .4, size = .8) +
    geom_point(data = subset(resArray, !is.na(group)),
               aes(x = cum_pos, y = .data[[score]], col = group), alpha = .6, size = 1) +
    scale_color_manual(values = c("hvCpG_Derakhshan" = "#DC3220",
                                  "mQTLcontrols" = "#005AB5"),
                       labels = c("hvCpG_Derakhshan" = "Derakhshan hvCpG",
                                  "mQTLcontrols" = "mQTL controls")) +
    scale_x_continuous(breaks = chr_mid$mid,
                       labels = gsub("chr", "", chr_mid$chr), expand = c(0, 0)) +
    theme_minimal(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(x = "Chromosome", y = score)+
    theme(legend.position = "inside",
          legend.position.inside      = c(0.9, 1.15),   # above the plot area
          legend.justification.inside = c(1, 1),
          plot.margin = margin(t = 40, r = 5, b = 5, l = 5),  # space for legend above
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linewidth = 0.5, linetype = "blank", colour = "black"))
  
  if (shift){
    p1_manhattanArray <- 
      p1_manhattanArray + 
      ylab(paste0(score, " shifted to start at zero")) 
  } else {
    p1_manhattanArray <- 
      p1_manhattanArray + 
      ylab(paste0(score))
  }
  
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
  
  # Merge hvCpG logBFs
  hv_score <- resArray[, c("chrpos", score)]
  colnames(hv_score) <- c("hvCpG", "score_hvCpG")
  
  # Merge control logBFs
  ctrl_score <- resArray[, c("chrpos", score)]
  colnames(ctrl_score) <- c("control", "score_control")
  
  # Join everything
  merged <- pairs %>%
    left_join(hv_score, by = "hvCpG") %>%
    left_join(ctrl_score, by = "control") %>%
    dplyr::mutate(diffscore = score_hvCpG-score_control)
  
  p2_DiffProbhvCpG_matchingcontrol_Array <-
    ggplot(merged, aes(x="", y=diffscore))+
    geom_jitter(data=merged[merged$diffscore>=0,], col="black", alpha=.3)+
    geom_jitter(data=merged[merged$diffscore<0,], fill="yellow",col="black",
                pch=21, alpha=.5)+
    geom_violin(width=.5, fill = "grey", alpha=.8) +
    geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
    theme_minimal(base_size = 11)+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  if (shift){
    p2_DiffProbhvCpG_matchingcontrol_Array <- 
      p2_DiffProbhvCpG_matchingcontrol_Array + 
      ylab(paste0("Difference of ", score, " shifted to start at zero")) 
  } else {
    p2_DiffProbhvCpG_matchingcontrol_Array <- 
      p2_DiffProbhvCpG_matchingcontrol_Array + 
      ylab(paste0("Difference of ", score))
  }
  
  ##
  
  ###############################
  ## Make figure of array test ##
  ###############################
  
  lab <- list(size = 14, x = 0.01, y = 0.99, hjust = 0, vjust = 1)
  mg  <- theme(plot.margin = margin(15, 5, 5, 5))
  
  row1 <- cowplot::plot_grid(
    p1_manhattanArray + theme(plot.margin = margin(40, 5, 5, 5)), 
    p2_DiffProbhvCpG_matchingcontrol_Array + theme(plot.margin = margin(50, 5, 5, 5)),
    ncol = 2, rel_widths = c(4, 1),
    labels = c("A. Detection of highly variable CpGs with both methods",
               "B. Difference of score between\nDerakhshan hvCpGs and\nmatching controls"),
    label_size = lab$size, label_x = lab$x, label_y = lab$y,
    hjust = lab$hjust, vjust = lab$vjust)
  row1
}

fig2 <- makeScript2Fig(resArray_0_8p0_0_65p1, path = pathNew)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/Fig2_newAug26_resArray_0_8p0_0_65p1.png"),
  plot = fig2, width = 18, height = 5,
  dpi = 300, bg = "white")


## TBC when sod is not down



################################################################################
makeScript2Fig <- function(resArray, path, p0p1){
  
  ################################################################
  ## Load full results on array with only 2 or 3 individuals/ds ##
  ################################################################
  
  # makePlotNrob <- function(resCompArray, N){
  #   mycor <- cor(resCompArray$alpha_array_all, 
  # resCompArray$alpha_array_reduce, use = "complete.obs")
  #   
  #   ggplot(resCompArray,
  #          aes(x=alpha_array_all, y=alpha_array_reduce)) +
  #     geom_point(data = resCompArray[is.na(resCompArray$group),], aes(col = group),
  #                alpha = 0.01) +
  #     geom_point(data = resCompArray[!is.na(resCompArray$group),], aes(col = group),
  #                alpha = 0.4) +
  #     geom_smooth(method = "lm", fill = "grey", col = "grey") +
  #     scale_color_manual(values = c("#DC3220", "#005AB5", "grey"),
  #                        labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  #     theme_minimal(base_size = 14) +
  #     theme(legend.title = element_blank()) +
  #     annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation:\n r = %.2f\n", mycor)) +
  #     labs(title = element_blank(),
  #          x = "P(hv) using full datasets",
  #          y = paste0("Pr(hv) using reduced (", N, " ind/ds) datasets")) +
  #     coord_cartesian(xlim = c(0,1), ylim = c(0,1))
  # }
  # 
  # resArray3ind <- as.data.frame(readRDS(here(paste0(
  #   path, "results_Arrays_3ind_406036CpGs_", p0p1, ".rds"))))
  # resArray3ind <- prepareChrDataset(resArray3ind)
  # names(resArray3ind)[names(resArray3ind) %in% "alpha"] <- "alpha_array_reduce"
  # 
  # resCompArray_allvs3 <- dplyr::left_join(
  #   resArray3ind,
  #   resArray[, c("chrpos", "alpha")],
  #   by = "chrpos"
  # )
  # names(resCompArray_allvs3)[names(resCompArray_allvs3) %in% "alpha"] <- "alpha_array_all"
  # 
  # p3ind <- makePlotNrob(resCompArray_allvs3, 3)
  # 
  # resArray2ind <- as.data.frame(readRDS(here(paste0(
  #   path, "results_Arrays_2ind_406036CpGs_", p0p1, ".rds"))))
  # resArray2ind <- prepareChrDataset(resArray2ind)
  # names(resArray2ind)[names(resArray2ind) %in% "alpha"] <- "alpha_array_reduce"
  # 
  # resCompArray_allvs2 <- dplyr::left_join(
  #   resArray2ind,
  #   resArray[, c("chrpos", "alpha")],
  #   by = "chrpos"
  # )
  # names(resCompArray_allvs2)[names(resCompArray_allvs2) %in% "alpha"] <- "alpha_array_all"
  # 
  # p2ind <- makePlotNrob(resCompArray_allvs2, 2)
  # 
  # ## From script S01:
  # x <- readRDS(here("B_MultiTissues/dataOut/arrayCutoffLowPower2or3ind.RDS"))
  # 
  # plot_venn3 <- ggVennDiagram(x, label = "both", label_alpha = 0,
  #                             label_color = "white",
  #                             category.names = c("Full \ndatasets",
  #                                                "2 ind per \ndataset",
  #                                                "3 ind per dataset")) +
  #   scale_fill_gradient(low = "grey", high = "black") +
  #   theme(legend.position = "none") +
  #   coord_cartesian(xlim = c(-5, 10), ylim = c(-10, 5), clip = "off")
  # 
  # plot_venn3 # "Cutoff algorithm"
  # 
  # message("What cutoff to get the same number of sites than Maria? with Bayesian approach")
  # print(table(resCompArray_allvs3$group))
  # 
  # top3535 <- resCompArray_allvs3 |>
  #   dplyr::slice_max(order_by = alpha_array_all, n = 3535, with_ties = FALSE)
  # 
  # top3535_2 <- resCompArray_allvs2 |>
  #   dplyr::slice_max(order_by = alpha_array_all, n = 3535, with_ties = FALSE)
  # 
  # # sanity check
  # table(top3535$chrpos %in% top3535_2$chrpos)
  # 
  # min(top3535$alpha_array_all)
  # 
  # pos3ind <- top3535$chrpos[top3535$alpha_array_reduce >= min(top3535$alpha_array_all)]
  # pos2ind <- top3535_2$chrpos[top3535_2$alpha_array_reduce >= min(top3535_2$alpha_array_all)]
  # 
  # y <- list(
  #   "Full array"     = top3535_2$chrpos,
  #   "Array 2 ind/ds" = pos2ind,
  #   "Array 3 ind/ds" = pos3ind
  # )
  # 
  # plot_venn3_Bayes <- ggVennDiagram(y, label = "both", label_alpha = 0, 
  #                                   label_color = "white", category.names = 
  #                                     c("Full \ndatasets","2 ind per \ndataset", "3 ind per dataset")) +
  #   scale_fill_gradient(low = "grey", high = "black") +
  #   theme(legend.position = "none")+
  #   coord_cartesian(xlim = c(-5, 10), ylim = c(-10, 5), clip = "off")
  # 
  # plot_venn3_Bayes
  # 
  # ###############################
  # ## Make figure of array test ##
  # ###############################
  # 
  # lab <- list(size = 14, x = 0.01, y = 0.99, hjust = 0, vjust = 1)
  # mg  <- theme(plot.margin = margin(15, 5, 5, 5))
  # 
  
  # 
  # row2_1 <- cowplot::plot_grid(
  #   plot_venn3      + theme_void(base_size = 10) + theme(legend.position = "none", plot.margin = margin(40, 5, 5, 5)),
  #   plot_venn3_Bayes + theme_void(base_size = 10) + theme(legend.position = "none", plot.margin = margin(40, 5, 5, 5)),
  #   labels = c("C. Detection of highly variable CpGs \nwith reduced datasets (cutoff)", 
  #              "D. Detection of highly variable CpGs \nwith reduced datasets (Bayesian)"), nrow = 1,
  #   label_size = lab$size, label_x = lab$x, label_y = lab$y,
  #   hjust = lab$hjust, vjust = lab$vjust)
  # 
  # row2_2 <- cowplot::plot_grid(
  #   p2ind + theme_minimal(base_size = 11) + theme(legend.position = "none") + mg,
  #   p3ind + theme_minimal(base_size = 11) + theme(legend.position = "none") + mg,
  #   labels = c("E. Bayesian: full vs 2 ind/ds", "F. Bayesian: full vs 3 ind/ds"), nrow = 1,
  #   label_size = lab$size, label_x = lab$x, label_y = lab$y,
  #   hjust = lab$hjust, vjust = lab$vjust)
  # 
  # row2 <- cowplot::plot_grid(row2_1, row2_2)
  # 
  # figure2 <- cowplot::plot_grid(row1, row2, nrow = 2)
  # 
  # return(figure2)
  return(p1_manhattanArray)
}

################################################################################
## Evolution of top-3535 CpG recovery as individuals/dataset increases        ##
## (p0=80%, p1=90%): 2, 3, 4, 5, 6, 10 ind/ds vs full datasets                ##
################################################################################

Ns <- c(2, 3, 4, 5, 6, 10, 15, 20)

## Baseline: top 3535 CpGs using the full datasets (same cutoff as before)
top3535_full <- resArrayAll_mean_strict |> dplyr::slice_max(order_by = alpha, n = 3535, with_ties = FALSE)
cutoff_full  <- min(top3535_full$alpha)

recovery <- lapply(Ns, function(N) {
  f <- here(paste0(pathMean, "results_Arrays_", N,
                   "ind_406036CpGs_0_8p0_0_9p1.rds"))
  if (!file.exists(f)) {
    message("SKIP N=", N, " - file not found (not computed yet?): ", f)
    return(NULL)
  }
  
  resN <- as.data.frame(readRDS(f)) |> prepareChrDataset()
  names(resN)[names(resN) == "alpha"] <- "alpha_reduced"
  
  comp <- dplyr::left_join(resN, resArrayAll_mean_strict[, c("chrpos", "alpha")], by = "chrpos")
  
  ## Of the top top CpGs under full data, how many are still >= the same
  ## cutoff when using only N individuals/dataset?
  recovered <- sum(comp$chrpos %in% top3535_full$chrpos & comp$alpha_reduced >= cutoff_full)
  
  data.frame(N = N, n_recovered = recovered, pct_recovered = recovered / nrow(top3535_full))
}) |> dplyr::bind_rows()

## Add the full-dataset reference point (100% recovery, by definition)
recovery <- dplyr::bind_rows(recovery, data.frame(N = NA, n_recovered = nrow(top3535_full), pct_recovered = 1))
recovery$N_label <- ifelse(is.na(recovery$N), "Full", as.character(recovery$N))
recovery$N_label <- factor(recovery$N_label, levels = c(as.character(sort(Ns)), "Full"))

print(recovery)

recovery <- recovery[!is.na(recovery$N),]

plotRecovery <- ggplot(recovery, aes(x = N, y = pct_recovered, group = 1)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "black") +
  geom_text(aes(label = scales::percent(pct_recovered, accuracy = .1)), vjust = -1, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  scale_x_continuous(breaks = c(2,3,4,5,10,15,20), labels = c(2,3,4,5,10,15,20)) +
  theme_minimal(base_size = 13) +
  labs(x = "Individuals per dataset", y = "% of top 3,535 CpGs (full data) still recovered",
       title = "Recovery of top hvCpG hits as individuals/dataset increases",
       subtitle = "p0=80%, p1=90% - cutoff fixed at the full-dataset top-3,535 threshold")

plotRecovery

# is the correlation trend smoother?
corByN <- lapply(Ns, function(N) {
  f <- here(paste0(pathMean, "results_Arrays_", N,
                   "ind_406036CpGs_0_8p0_0_9p1.rds"))
  if (!file.exists(f)) return(NULL)
  resN <- as.data.frame(readRDS(f)) |> prepareChrDataset()
  names(resN)[names(resN) == "alpha"] <- "alpha_reduced"
  comp <- dplyr::left_join(resN, resArrayAll_mean_strict[, c("chrpos", "alpha")], by = "chrpos")
  data.frame(N = N, r = cor(comp$alpha, comp$alpha_reduced, use = "complete.obs"))
}) |> dplyr::bind_rows()

print(corByN)

plotCorByN <- ggplot(corByN, aes(x = N, y = r)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "black") +
  geom_text(aes(label = scales::percent(r, accuracy = .1)), vjust = -1, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(.5, 1.05)) +
  scale_x_continuous(breaks = c(2,3,4,5,10,15,20), labels = c(2,3,4,5,10,15,20)) +
  theme_minimal(base_size = 12) +
  labs(x = "Individuals per dataset", y = "Correlation (Pearson r) of Pr(hv) with full data",
       title = "Agreement with full-data Pr(hv) improves with sample size",
       subtitle = "p0=80%, p1=90%")

plotCorByN

ggsave(here("gitignore/improveAccuracyWithMoreNperds.png"),
       plotRecovery / plotCorByN, width = 9, height = 10, dpi = 300, bg = "white")