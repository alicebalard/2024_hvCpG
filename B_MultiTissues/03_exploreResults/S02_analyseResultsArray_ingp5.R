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
      chr = stringr::str_extract(chrpos, "^chr[0-9XY]+"),
      pos = as.numeric(stringr::str_extract(chrpos, "(?<=_)[0-9]+"))
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
    dplyr::mutate(cum_pos = pos + cum_start)
  
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

######################################
## NB: save one for later scripts
## Save for next scripts
saveRDS(resArray_0_8p0_0_65p1, here("B_MultiTissues/dataOut/resArray_0_8p0_0_65p1.RDS"))

resArray3ind_0_8p0_0_65p1 <- as.data.frame(
  readRDS(here(paste0(pathNew, "results_Arrays_3ind_406036CpGs_0_8p0_0_65p1_0_01a0.rds"))))
resArray3ind_0_8p0_0_65p1 <- prepareChrDataset(resArray3ind_0_8p0_0_65p1)
saveRDS(resArray3ind_0_8p0_0_65p1, here("B_MultiTissues/dataOut/resArray3ind_0_8p0_0_65p1.RDS"))

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
  
  resArray$prophv <- resArray$n_hv_ds / resArray$n_ds
  
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
                    column1 = "logBF", column2 = "prophv")
  
  pg <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArray[sample(1:nrow(resArray),100000),], 
                    column1 = "logBF_per_ds", column2 = "prophv")
  
  ph <- getPlotComp(resArray[sample(1:nrow(resArray),100000),], 
                    resArrayAll_sum[sample(1:nrow(resArrayAll_sum),100000),], 
                    column2 = "prophv")
  
  plot <- cowplot::plot_grid(pa, pb, pc, pd, pe, pf, pg, ph)
  return(plot)
}

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/plotDifferentScores_0_8p0_0_65p1.png"),
  plot = getPlotScores(resArray = resArray_0_8p0_0_65p1), width = 18, height = 10,
  dpi = 300, bg = "white")

################################################################################
## Compare p0/p1 settings — one metric decides: sensitivity at fixed specificity
################################################################################

## 1. one by the other
slim <- function(x) {
  x <- as.data.table(x)[, .(chrpos, logBF_per_ds)]
  setkey(x, chrpos)
  x
}

a <- slim(resArray_0_65p0_0_65p1); setnames(a, "logBF_per_ds", "lp_65")
b <- slim(resArray_0_9p0_0_9p1);   setnames(b, "logBF_per_ds", "lp_90")

m <- a[b, nomatch = 0]   # inner join on the shared key — CpGs present in both
rm(a, b); gc()

cor_val <- m[, cor(lp_65, lp_90, use = "complete.obs")]   # on ALL rows, not the subsample

plot65by90 <- ggplot(m, aes(lp_65, lp_90)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Pearson's r = %.3f", cor_val)) +
  labs(x = "Hypervariability score (p0=0.65, p1=0.65)",
       y = "Hypervariability score (p0=0.9,  p1=0.9)") +
  theme_minimal(base_size = 13)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/plot65by90.png"),
  plot = plot65by90,
  width = 6, height = 6,
  dpi = 300, bg = "white")

## 2. more detailed choice (ROC curve)

# ── load one settings file, label the two control groups, return scored data ──
load_scored <- function(file, base, score = "logBF_per_ds") {
  dat <- prepareChrDataset(as.data.frame(readRDS(file.path(base, file))))
  dat$group <- data.table::fifelse(dat$chrpos %in% DerakhshanhvCpGs_hg38, "pos",
                                   data.table::fifelse(dat$chrpos %in% mQTLcontrols_hg38,     "neg", NA))
  dat <- dat[!is.na(dat$group) & is.finite(dat[[score]]), ]
  data.table::data.table(group = dat$group, score = dat[[score]])
}

# ── everything about one setting: AUC, pAUC(FPR<0.1), sensitivity@spec ─────────
eval_setting <- function(d, spec = 0.95) {
  r <- pROC::roc(d$group == "pos", d$score, direction = "<", quiet = TRUE)
  data.table::data.table(
    auc      = as.numeric(pROC::auc(r)),
    pauc_lowFPR = as.numeric(pROC::auc(r, partial.auc = c(0.9, 1),
                                       partial.auc.focus = "specificity", partial.auc.correct = TRUE)),
    sens_at_spec = as.numeric(pROC::coords(r, x = spec, input = "specificity",
                                           ret = "sensitivity")),
    n_pos = sum(d$group == "pos"), n_neg = sum(d$group == "neg"),
    roc = list(r))                              # keep for optional DeLong test
}

# ── settings to compare ───────────────────────────────────────────────────────
base <- here("B_MultiTissues/resultsDir_gitIgnored/Arrays")
settings <- c(
  "p0=0.65,p1=0.8" = "results_Arrays_all_406036CpGs_0_65p0_0_8p1_0_01a0.rds",
  "p0=0.8,p1=0.65" = "results_Arrays_all_406036CpGs_0_8p0_0_65p1_0_01a0.rds",
  "p0=0.8,p1=0.8"  = "results_Arrays_all_406036CpGs_0_8p0_0_8p1_0_01a0.rds",
  "p0=0.9,p1=0.9"  = "results_Arrays_all_406036CpGs_0_9p0_0_9p1_0_01a0.rds")

SCORE <- "logBF_per_ds"; SPEC <- 0.95     # operating point: 5% mQTL false positives

# ── the whole comparison in one table ─────────────────────────────────────────
data_list <- lapply(settings, load_scored, base = base, score = SCORE)
res <- data.table::rbindlist(lapply(data_list, eval_setting, spec = SPEC),
                             idcol = "setting")
res_show <- res[order(-sens_at_spec), .(setting, auc, pauc_lowFPR, sens_at_spec, n_pos, n_neg)]
print(res_show)

cat(sprintf("\n>> Choose: %s  (sensitivity %.1f%% at %d%% specificity)\n",
            res_show$setting[1], 100*res_show$sens_at_spec[1], round(100*SPEC)))
# >> Choose: p0=0.8,p1=0.65  (sensitivity 51.3% at 95% specificity)

# the density panels (visual sanity check only — not needed to decide)
plot_density <- function(d, title) ggplot(d, aes(score, colour = group, fill = group)) +
  geom_density(alpha = .25, linewidth = .9) + geom_vline(xintercept = 0, linetype = 2) +
  scale_colour_manual(values = c(pos="#D55E00", neg="#0072B2")) +
  scale_fill_manual(values   = c(pos="#D55E00", neg="#0072B2")) +
  labs(x = SCORE, title = title) + theme_bw(12) + theme(legend.position = "none")

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/plotp0p1.png"),
  plot = cowplot::plot_grid(plotlist = Map(plot_density, data_list, names(settings))),
  width = 18, height = 10,
  dpi = 300, bg = "white")

################################################################################
makeScript2Fig <- function(resArray, path, p0p1, score = "logBF_per_ds",
                           label = "logBF per ds", shift = TRUE){
  
  setDT(resArray)
  if (shift){
    message("We shift the score so it starts at zero (for the Manhattan and violin only)")
    resArray[, score_shifted := get(score) + abs(min(get(score), na.rm = TRUE))]
  } else {
    resArray[, score_shifted := get(score)]
  }
  
  # Plot
  band_cols <- c("0" = "grey20", "1" = "grey55")
  
  # --- chromosome order (by genomic position), midpoints, and parity: one source of truth ---
  setDT(resArray)
  
  chr_tab <- resArray[!is.na(cum_pos),
                      .(minp = min(cum_pos), maxp = max(cum_pos)), by = chr]
  setorder(chr_tab, minp)                         # genomic order, gap-free
  chr_tab[, `:=`(mid    = (minp + maxp) / 2,
                 parity = factor(seq_len(.N) %% 2))]
  
  # attach parity to the points
  resArray <- merge(resArray, chr_tab[, .(chr, parity)], by = "chr", all.x = TRUE, sort = FALSE)
  
  ## colorblind friendly
  p1_manhattanArray <- ggplot() +
    # background cloud: alternating grey per chromosome
    geom_point(data = resArray[is.na(group)],
               aes(x = cum_pos, y = .data[["score_shifted"]], colour = parity),
               alpha = .4, size = .8) +
    scale_colour_manual(values = band_cols, na.translate = FALSE, guide = "none") +
    # highlighted sets keep their own colour scale
    ggnewscale::new_scale_colour() +
    geom_point(data = resArray[!is.na(group)],
               aes(x = cum_pos, y = .data[["score_shifted"]], colour = group), alpha = .6, size = 1) +
    scale_colour_manual(values = c("hvCpG_Derakhshan" = "#DC3220",
                                   "mQTLcontrols"     = "#005AB5"),
                        labels = c("hvCpG_Derakhshan" = "Derakhshan hvCpG",
                                   "mQTLcontrols"     = "mQTL controls")) +
    scale_x_continuous(breaks = chr_tab$mid,
                       labels = gsub("chr", "", chr_tab$chr), expand = c(0, 0)) +
    theme_minimal(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(x = "Chromosome", y = score) +
    theme(legend.position            = "inside",
          legend.position.inside      = c(0.9, 1.15),
          legend.justification.inside = c(1, 1),
          plot.margin  = margin(t = 40, r = 5, b = 5, l = 5),
          legend.title = element_blank(),
          legend.text  = element_text(size = 14),
          legend.background = element_rect(linewidth = 0.5, linetype = "blank", colour = "black"))
  
  if (shift){
    p1_manhattanArray <-
      p1_manhattanArray +
      ylab(paste0(label, " shifted to start at zero"))
  } else {
    p1_manhattanArray <-
      p1_manhattanArray +
      ylab(paste0(label))
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
  hv_score <- resArray[, .(hvCpG = chrpos, score_hvCpG = score_shifted)]
  
  # Merge control logBFs
  ctrl_score <- resArray[, .(control = chrpos, score_control = score_shifted)]
  
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
      ylab(paste0("Difference of ", label, " shifted to start at zero"))
  } else {
    p2_DiffProbhvCpG_matchingcontrol_Array <-
      p2_DiffProbhvCpG_matchingcontrol_Array +
      ylab(paste0("Difference of ", label))
  }
  
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
  
  ################################################################
  ## Load full results on array with only 2 or 3 individuals/ds ##
  ################################################################
  
  makePlotNrob <- function(resCompArray, N){
    mycor <- cor(resCompArray$logBF_per_ds_array_all,
                 resCompArray$logBF_per_ds_array_reduce, use = "complete.obs")
    
    xr <- range(resCompArray$logBF_per_ds_array_all,    na.rm = TRUE)
    yr <- range(resCompArray$logBF_per_ds_array_reduce, na.rm = TRUE)
    
    ggplot(resCompArray,
           aes(x = logBF_per_ds_array_all, y = logBF_per_ds_array_reduce)) +
      geom_point(data = resCompArray[is.na(resCompArray$group),], aes(col = group), alpha = 0.01) +
      geom_point(data = resCompArray[!is.na(resCompArray$group),], aes(col = group), alpha = 0.4) +
      geom_smooth(method = "lm", fill = "grey", col = "grey") +
      scale_color_manual(values = c("#DC3220", "#005AB5", "grey"),
                         labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
      # anchor annotation to data corner, not to 0.5/0.9
      annotate("text", x = xr[1], y = yr[2], hjust = 0, vjust = 1,
               label = sprintf("Pearson r = %.2f", mycor)) +
      coord_cartesian(xlim = xr, ylim = yr) +          # explicit, data-driven — nothing clipped
      theme_minimal(base_size = 14) +
      theme(legend.title = element_blank()) +
      labs(x = "Hypervariability score, full datasets",
           y = sprintf("Hypervariability score, reduced (%d ind/ds)", N))
  }
  
  resArray3ind <- as.data.frame(readRDS(here(paste0(
    path, "results_Arrays_3ind_406036CpGs_", p0p1, "_0_01a0.rds"))))
  resArray3ind <- prepareChrDataset(resArray3ind)
  names(resArray3ind)[names(resArray3ind) %in% "logBF_per_ds"] <- "logBF_per_ds_array_reduce"
  
  resCompArray_allvs3 <- dplyr::left_join(
    resArray3ind,
    resArray[, c("chrpos", "logBF_per_ds")],
    by = "chrpos"
  )
  names(resCompArray_allvs3)[names(resCompArray_allvs3) %in% "logBF_per_ds"] <- "logBF_per_ds_array_all"
  
  p3ind <- makePlotNrob(resCompArray_allvs3, 3)
  
  resArray2ind <- as.data.frame(readRDS(here(paste0(
    path, "results_Arrays_2ind_406036CpGs_", p0p1, "_0_01a0.rds"))))
  resArray2ind <- prepareChrDataset(resArray2ind)
  names(resArray2ind)[names(resArray2ind) %in% "logBF_per_ds"] <- "logBF_per_ds_array_reduce"
  
  resCompArray_allvs2 <- dplyr::left_join(
    resArray2ind,
    resArray[, c("chrpos", "logBF_per_ds")],
    by = "chrpos"
  )
  names(resCompArray_allvs2)[names(resCompArray_allvs2) %in% "logBF_per_ds"] <- "logBF_per_ds_array_all"
  
  p2ind <- makePlotNrob(resCompArray_allvs2, 2)
  
  ###############################################
  ## Compare top-K in full vs reduced analyses ##
  ###############################################
  
  ### Cutoff method from script S01:
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
  
  ## K = full-data cutoff hits (already computed): length(x$`Full array`)
  K <- length(x$`Full array`)
  
  topK <- function(dt, score_col, k = K) {
    dt  <- as.data.frame(dt)
    ok  <- !is.na(dt[[score_col]])
    ord <- order(dt[[score_col]][ok], decreasing = TRUE)
    dt$chrpos[ok][ord][seq_len(min(k, sum(ok)))]
  }
  
  # ---- Bayesian: top-K per analysis ----
  full_top <- topK(resCompArray_allvs3, "logBF_per_ds_array_all")
  red3_top <- topK(resCompArray_allvs3, "logBF_per_ds_array_reduce")
  red2_top <- topK(resCompArray_allvs2, "logBF_per_ds_array_reduce")
  
  # ---- Cutoff: from the saved Venn list x ----
  cut_full <- x$`Full array`
  cut_2    <- x$`Array 2 ind/ds`
  cut_3    <- x$`Array 3 ind/ds`
  
  recovery <- data.frame(
    method   = c("Bayesian","Bayesian","Cutoff","Cutoff"),
    n_ind    = c(2, 3, 2, 3),
    n_full   = c(length(full_top), length(full_top), length(cut_full), length(cut_full)),
    n_recovered = c(
      length(intersect(full_top, red2_top)),
      length(intersect(full_top, red3_top)),
      length(intersect(cut_full, cut_2)),
      length(intersect(cut_full, cut_3))
    )
  )
  recovery$pct_recovered <- round(100 * recovery$n_recovered / recovery$n_full, 1)
  print(recovery)
  row2_2 <- cowplot::plot_grid(
    p2ind + theme_minimal(base_size = 11) + theme(legend.position = "none") + mg,
    p3ind + theme_minimal(base_size = 11) + theme(legend.position = "none") + mg,
    labels = c("C. Bayesian: full vs 2 ind/ds", "D. Bayesian: full vs 3 ind/ds"), nrow = 1,
    label_size = lab$size, label_x = lab$x, label_y = lab$y,
    hjust = lab$hjust, vjust = lab$vjust)
  
  row2 <- cowplot::plot_grid(row2_2)
  
  figure2 <- cowplot::plot_grid(row1, row2, nrow = 2)
  
  return(figure2)
}

fig2 <- makeScript2Fig(resArray_0_8p0_0_65p1, path = pathNew, p0p1 = "0_8p0_0_65p1")
# method n_ind n_full n_recovered pct_recovered
# 1 Bayesian     2   4377        1273          29.1
# 2 Bayesian     3   4377        1644          37.6
# 3   Cutoff     2   4377         143           3.3
# 4   Cutoff     3   4377         693          15.8

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script02/Fig2_newAug26_resArray_0_8p0_0_65p1.png"),
  plot = fig2, width = 18, height = 10,
  dpi = 300, bg = "white")

top5pc <- quantile(resArray_0_8p0_0_65p1$logBF_per_ds, prob=1-5/100)
top <- resArray_0_8p0_0_65p1[
  resArray_0_8p0_0_65p1$logBF_per_ds >= top5pc,]

topNeg <- top[top$group %in% "mQTLcontrols",]
min(topNeg$n_hv_ds / topNeg$n_ds)
nrow(topNeg) # 383

topPos <- top[top$group %in% "hvCpG_Derakhshan",]
min(topPos$n_hv_ds / topPos$n_ds)
nrow(topPos)  # 3027

## --> 10 times more hvCpGs than mQTLcontrols in the top 5% logBF_per_ds
## both detected in more than 92% of the datasets

################################################################################
## Evolution of top-4377 CpG recovery as individuals/dataset increases        ##
## (p0=80%, p1=90%): 2, 3, 4, 5, 6, 10 ind/ds vs full datasets                ##
################################################################################

Ns <- c(2, 3, 4, 5, 10, 15, 20)
K  <- 4377

## Baseline: top-K CpGs under the full datasets
top_full <- resArray_0_8p0_0_65p1 |>
  dplyr::slice_max(order_by = logBF_per_ds, n = K, with_ties = FALSE)
full_hits <- top_full$chrpos

recovery <- lapply(Ns, function(N) {
  f <- here(paste0(pathNew, "results_Arrays_", N,
                   "ind_406036CpGs_0_8p0_0_65p1_0_01a0.rds"))
  if (!file.exists(f)) { message("SKIP N=", N, " - file not found: ", f); return(NULL) }
  
  resN <- as.data.frame(readRDS(f)) |> prepareChrDataset()
  
  ## top-K by the REDUCED run's OWN score (rank, not full-data threshold)
  topN_reduced <- resN |>
    dplyr::filter(!is.na(logBF_per_ds)) |>
    dplyr::slice_max(order_by = logBF_per_ds, n = K, with_ties = FALSE)
  
  recovered <- sum(full_hits %in% topN_reduced$chrpos)   # overlap of the two top-K sets
  data.frame(N = N, n_recovered = recovered, pct_recovered = recovered / K)
}) |> dplyr::bind_rows()

## Add the full-dataset reference point (100% recovery, by definition)
recovery <- dplyr::bind_rows(recovery, data.frame(N = NA, n_recovered = nrow(top4377_full),
                                                  pct_recovered = 1))
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
  labs(x = "Individuals per dataset", y = "% of top 4377 CpGs (full data) still recovered",
       title = "Recovery of top hvCpG hits as individuals/dataset increases",
       subtitle = "p0=80%, p1=65% - cutoff fixed at the full-dataset top-4377 threshold")

plotRecovery

# is the correlation trend smoother?
corByN <- lapply(Ns, function(N) {
  f <- here(paste0(pathNew, "results_Arrays_", N,
                   "ind_406036CpGs_0_8p0_0_65p1_0_01a0.rds"))
  if (!file.exists(f)) return(NULL)
  resN <- as.data.frame(readRDS(f)) |> prepareChrDataset()
  names(resN)[names(resN) == "logBF_per_ds"] <- "logBF_per_ds_reduced"
  comp <- dplyr::left_join(resN, resArray_0_8p0_0_65p1[, c("chrpos", "logBF_per_ds")], by = "chrpos")
  data.frame(N = N, r = cor(comp$logBF_per_ds, comp$logBF_per_ds_reduced, use = "complete.obs"))
}) |> dplyr::bind_rows()

print(corByN)

plotCorByN <- ggplot(corByN, aes(x = N, y = r)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "black") +
  geom_text(aes(label = scales::percent(r, accuracy = .1)), vjust = -1, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(.5, 1.05)) +
  scale_x_continuous(breaks = c(2,3,4,5,10,15,20), labels = c(2,3,4,5,10,15,20)) +
  theme_minimal(base_size = 12) +
  labs(x = "Individuals per dataset", y = "Correlation (Pearson r) of logBF per ds with full data",
       title = "Agreement with full-data hypervariability score improves with sample size",
       subtitle = "p0=80%, p1=65%")

plotCorByN

ggsave(here("B_MultiTissues/dataOut/figures/script02/improveAccuracyWithMoreNperds.png"),
       plotRecovery / plotCorByN, width = 9, height = 12, dpi = 300, bg = "white")
