################################################################################ 
# Identify hypervariable CpGs that are hypervariable only in immune cells
# Hyp: different results using arrays on tissue vs WGBS of purified cells
# because of the different datasets composition, esp. the ratio of immune cells
################################################################################
library(here)

## Load libraries
if (!exists("libLoaded")) {
  source(here("05_hvCpGalgorithm", "quiet_library.R"))}

## Load functions
if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))}

## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))}

## Hypothesis testing:
## Prediction 1:
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_11_vs_9_immuneEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_11_vs_9_immuneEffect",
    xlab = "Pr(hv) calculated on WGBS atlas datasets, immune cells removed (11 gp)",
    ylab = "Pr(hv) calculated on WGBS atlas datasets, only immune cells (11 gp)")
}

if (!exists("Z_inner_immvsnoimm")){
  system.time(Z_inner_immvsnoimm <- makeZ_inner(
    X = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha")) ## takes 2 minutes
  print("Z_inner_immvsnoimm created")
}

stable <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]
cellUniversal <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
immune <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
notimmune <-  Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]

message(paste0("We find ", 
               length(stable), " (", round(length(stable)/length(Z_inner_immvsnoimm$name)*100), "%) stable CpGs, ",
               length(cellUniversal), " (", round(length(cellUniversal)/length(Z_inner_immvsnoimm$name)*100), "%) cell universal hvCpGs, ", 
               length(immune), " (", round(length(immune)/length(Z_inner_immvsnoimm$name)*100), "%) immune cell hvCpGs and ",
               length(notimmune), " (", round(length(notimmune)/length(Z_inner_immvsnoimm$name)*100), "%) hvCpGs undetected in immune cells, out of a total of ",
               length(Z_inner_immvsnoimm$name), " CpGs investigated"))

# We find 20268230 (84%) stable CpGs, 1190281 (5%) cell universal hvCpGs, 1507272 (6%) 
# immune cell hvCpGs and 1123906 (5%) hvCpGs undetected in immune cells, out of a total of 24089689 CpGs investigated

## Prediction 1: Using only the immune cells in WGBS Atlas should give a result closer to what is observed using the array

testPred1=FALSE # change if needs to rerun

if (testPred1){
  Z_inner_all <- makeZ_inner(
    X = resCompArray,
    Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
    whichAlphaX = "alpha_array_all",
    whichAlphaY = "alpha")
  
  Z_inner_noimm <- makeZ_inner(
    X = resCompArray,
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups")),
    whichAlphaX = "alpha_array_all",
    whichAlphaY = "alpha")
  
  Z_inner_onlyimm <- makeZ_inner(
    X = resCompArray,
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha_array_all",
    whichAlphaY = "alpha")
  
  c_noimm <- cor.test(Z_inner_noimm$alpha_X, Z_inner_noimm$alpha_Y)
  c_all <- cor.test(Z_inner_all$alpha_X, Z_inner_all$alpha_Y)
  c_onlyimm <- cor.test(Z_inner_onlyimm$alpha_X, Z_inner_onlyimm$alpha_Y)
  
  slope_noimm <- lm(data = Z_inner_noimm, alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
  slope_all <- lm(data = Z_inner_all,  alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
  slope_onlyimm <- lm(data = Z_inner_onlyimm,  alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
  
  p <- data.frame(x=seq(0,1,0.01)) %>%
    dplyr::mutate(noImmun = slope_noimm * x,
                  all = slope_all * x,
                  onlyimm = slope_onlyimm * x) %>%
    pivot_longer(cols = c(noImmun, all, onlyimm)) %>%
    ggplot(aes(x = x, y = value, group = name, col = name)) +
    geom_abline(slope = 1, linetype = 3) +
    geom_line(linewidth = 2) +
    scale_color_viridis_d() +
    theme_minimal(base_size = 14) +
    labs(x = "Pr(hv) array datasets", y = "Pr(hv) WGBS datasets", col = "WGBS datasets") +
    ylim(c(0,1))
  
  ggplot2::ggsave(
    filename = here::here("05_hvCpGalgorithm/figures/correlations/correlation_prediction2.pdf"),
    plot = p, width = 9, height = 7
  )
  
  # Combine datasets
  Z_inner_all$group <- "all"
  Z_inner_noimm$group <- "noImm"
  Z_inner_onlyimm$group <- "onlyImm"
  
  combined <- rbind(Z_inner_all, Z_inner_noimm, Z_inner_onlyimm)
  
  # Fit model with interaction
  model <- lm(alpha_Y ~ alpha_X * group, data = combined)
  summary(model)
  
  library(emmeans)
  emm <- emtrends(model, ~ group, var = "alpha_X")
  pairs(emm)
  
  # contrast        estimate      SE     df t.ratio p.value
  # all - noImm     -0.00505 0.00166 988312  -3.042  0.0066
  # all - onlyImm   -0.10615 0.00169 988312 -62.974  <.0001
  # noImm - onlyImm -0.10110 0.00163 988312 -62.005  <.0001
  
  ## Difference with controls in both cases?
  data <- read.table(here("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)
  
  # X: no immune; Y: only immune
  X <- readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups"))
  Y <- readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly"))
  
  x = dico$chrpos_hg38[match(data$hvCpG_name, dico$CpG)]
  y = dico$chrpos_hg38[match(data$controlCpG_name, dico$CpG)]
  
  # Build mapping from hvCpG -> control
  pairs <- data.frame(
    hvCpG = x,
    control = y,
    stringsAsFactors = FALSE
  )
  
  # Merge hvCpG alphas
  res <- X[!is.na(group)]
  
  hv_alpha <- res[, c("name", "alpha")]
  colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")
  
  # Merge control alphas
  ctrl_alpha <- res[, c("name", "alpha")]
  colnames(ctrl_alpha) <- c("control", "alpha_control")
  
  # Join everything
  merged <- pairs %>%
    left_join(hv_alpha, by = "hvCpG") %>%
    left_join(ctrl_alpha, by = "control") %>%
    mutate(diffAlpha=alpha_hvCpG-alpha_control)
  
  merged1 <- merged %>%
    mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
    filter(!is.na(diffAlpha))
  
  p1 <- ggplot(merged1, aes(x="diff", y=diffAlpha))+
    geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
    geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
    geom_violin(width=.5, fill = "grey", alpha=.8) +
    geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
    theme_minimal(base_size = 14)+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
    ggtitle("P(hvCpG) minus P(matching control)", subtitle = "no immune cells")+
    ylab("Difference of probability") +
    coord_cartesian(ylim = c(-1,1))
  
  # Merge hvCpG alphas
  res <- Y[!is.na(group)]
  
  hv_alpha <- res[, c("name", "alpha")]
  colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")
  
  # Merge control alphas
  ctrl_alpha <- res[, c("name", "alpha")]
  colnames(ctrl_alpha) <- c("control", "alpha_control")
  
  # Join everything
  merged <- pairs %>%
    left_join(hv_alpha, by = "hvCpG") %>%
    left_join(ctrl_alpha, by = "control") %>%
    mutate(diffAlpha=alpha_hvCpG-alpha_control)
  
  merged2 <- merged %>%
    mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
    filter(!is.na(diffAlpha))
  
  p2 <- ggplot(merged2, aes(x="diff", y=diffAlpha))+
    geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
    geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
    geom_violin(width=.5, fill = "grey", alpha=.8) +
    geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
    theme_minimal(base_size = 14)+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
    ggtitle("P(hvCpG) minus P(matching control)", subtitle =  "immune cells only")+
    ylab("Difference of probability") +
    coord_cartesian(ylim = c(-1,1))
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  wilcox.test(merged1$diffAlpha, merged2$diffAlpha)
  
  mean(merged1$diffAlpha); median(merged1$diffAlpha)
  mean(merged2$diffAlpha); median(merged$diffAlpha, na.rm = T)
}

##############################################################
## Are these hvCpG detected also in individual germ layers? ##
##############################################################

cellUniversal_GR <- makeGRfromMyCpGPos(cellUniversal, "cellUniversal")
immune_GR <- makeGRfromMyCpGPos(immune, "immune")
notimmune_GR <- makeGRfromMyCpGPos(notimmune, "notimmune")
stable_GR <- makeGRfromMyCpGPos(stable, "stable")

makeVenn=FALSE
if (makeVenn){  
  gr_list <- list(cellUniversal = cellUniversal_GR, 
                  immune = immune_GR, 
                  notimmune = notimmune_GR, 
                  stable = stable_GR, 
                  DerakhshanhvCpGs = DerakhshanhvCpGs_hg38_GR) ## test previous hvCpGs
  
  lapply(names(gr_list), function(name){
    x <- gr_list[[name]]
    gr_names <- paste0(as.character(seqnames(x)), "_", start(x))
    p1 <- plotMyVenn(0.5, endo = endo[endo$name %in% gr_names,],
                     meso = meso[meso$name %in% gr_names,],
                     ecto = ecto[ecto$name %in% gr_names,],
                     all = allLayers[allLayers$name %in% gr_names,])
    p2 <- plotMyVenn(0.75, endo = endo[endo$name %in% gr_names,],
                     meso = meso[meso$name %in% gr_names,],
                     ecto = ecto[ecto$name %in% gr_names,],
                     all = allLayers[allLayers$name %in% gr_names,])
    p3 <- plotMyVenn(0.9, endo = endo[endo$name %in% gr_names,],
                     meso = meso[meso$name %in% gr_names,],
                     ecto = ecto[ecto$name %in% gr_names,],
                     all = allLayers[allLayers$name %in% gr_names,])
    pdf(file = paste0("../../05_hvCpGalgorithm/figures/vennGermLayers/Venn_", 
                      name, "_GR.pdf"), width = 13, height = 6)
    print(cowplot::plot_grid(p1 + theme(legend.position = "none"),
                             p2 + theme(legend.position = "none"),
                             p3 + theme(legend.position = "none"), 
                             nrow = 1, labels = name)) 
    dev.off()
  })
}

################################################################################
# Prediction 2. The putative ME are those CpG sites which variability is high no
# matter which collection of cells one uses (“cell-universal” hvCpGs)

###################################################################################
## Find target CpGs in each quadrants: putative MEs but also vmeQTL and MHC CpGs ##

# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
chr6pos = Z_inner_immvsnoimm$name[grep("chr6", Z_inner_immvsnoimm$name)]
chr6pos = as.integer(sub(".*_", "", chr6pos))
MHCpos <- paste0("chr6_", chr6pos[chr6pos >= 28510120 & chr6pos <= 33480577])
length(MHCpos) # 56491
MHCpos_GR <- makeGRfromMyCpGPos(MHCpos, "MHC")

testEnrichment_ResandPlot <- function(targetCpGs){
  hits <- GenomicRanges::findOverlaps(cellUniversal_GR, targetCpGs)
  me_cellUniversal <- targetCpGs[unique(subjectHits(hits))]
  
  hits <- GenomicRanges::findOverlaps(immune_GR, targetCpGs)
  me_immune <- targetCpGs[unique(subjectHits(hits))]
  
  hits <- GenomicRanges::findOverlaps(notimmune_GR, targetCpGs)
  me_notimmune <- targetCpGs[unique(subjectHits(hits))]
  
  hits <- GenomicRanges::findOverlaps(stable_GR, targetCpGs)
  me_stable <- targetCpGs[unique(subjectHits(hits))]
  
  ###########################################
  # Compute enrichment for one ME set
  quads <- list(
    stable        = stable_GR,
    immune        = immune_GR,
    notimmune     = notimmune_GR,
    cellUniversal = cellUniversal_GR
  )
  
  # ---- Run it (ME sets in targetCpGs$set will be tested separately)
  res_quadrants <- test_enrichment_quadrants(quads, targetCpGs, me_col = "set")
  
  # Add log2OR & significant or not
  res_plot2 <- res_quadrants %>%
    mutate(log2OR = log2(odds_ratio),
           signif = p_adj_BH < 0.05) %>%
    group_by(CpG_set) %>%
    ungroup()
  
  res_plot2$signif <- factor(res_plot2$signif, levels = c(TRUE, FALSE))
  
  res_plot2$quadrant <- factor(res_plot2$quadrant, 
                               levels = c("stable", "notimmune", "immune", "cellUniversal"))
  
  if("Derakhshan hvCpG" %in% res_plot2$CpG_set){
    res_plot2$CpG_set <- factor(res_plot2$CpG_set, 
                                levels = c("Derakhshan hvCpG", "Harris SIV", 
                                           "VanBaak ESS", "Kessler SIV", "Gunasekara corSIV"))
  }
  
  plot <- ggplot(res_plot2, aes(x = quadrant, y = log2OR, fill = signif)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("black", "grey")) +
    labs(
      x = "Quadrant",
      y = expression(log[2]~"(odds ratio)"),
      title = "CpG sets enrichment by quadrant vs other three",
      subtitle = "numbers = CpGs in quadrant (% of CpGs in quadrant that these CpGs represent)") +
    facet_wrap(~ CpG_set, scales = "free_x", nrow = 1) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )+
    geom_text(
      aes(label = paste0(
        inME, "\n(", 
        ifelse(round(pct_inME, 2)>0.01, round(pct_inME, 2), round(pct_inME, 3)), "%)"
      ),vjust = ifelse(log2OR >= 0, -0.3, 1.3)  # uses log2OR inside aes
      ),size = 2.8
    )
  
  return(list(res=res_plot2, plot = plot))
}

testPutativeMEs <- testEnrichment_ResandPlot(targetCpGs = putativeME_GR)
pdf(here("05_hvCpGalgorithm/figures/quadrantEnrichment/quadrantsEnrich_putativeMEs.pdf"), width = 12, height = 9)
testPutativeMEs$plot
dev.off()

testPutativeMHC <- testEnrichment_ResandPlot(targetCpGs = MHCpos_GR)
pdf(here("05_hvCpGalgorithm/figures/quadrantEnrichment/quadrantsEnrich_MHC.pdf"), width = 6, height = 9)
testPutativeMHC$plot
dev.off()

testPutativevmeQTL <- testEnrichment_ResandPlot(targetCpGs = vmeQTL_hg38_GR)
pdf(here("05_hvCpGalgorithm/figures/quadrantEnrichment/quadrantsEnrich_vmeQTL.pdf"), width = 6, height = 9)
testPutativevmeQTL$plot
dev.off()

## Prediction 3. CpGs that are only variable using immune cells only (“immune variable hvCpGs”)
# or using non-immune cells capture variability acquired later in life

X = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups"))
Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly"))

plotGeneQuadrant <- function(gene_name, gene_chr, gene_start, gene_end){
  
  X_gene <- X[chr == gene_chr & pos >= gene_start & pos <= gene_end]
  Y_gene <- Y[chr == gene_chr & pos >= gene_start & pos <= gene_end]
  
  X_gene <- X_gene[, !c("cum_offset", "pos2"), with=FALSE]
  Y_gene <- Y_gene[, !c("cum_offset", "pos2"), with=FALSE]
  
  names(X_gene)[names(X_gene) %in% "alpha"] = "alpha_no_immune_cells"
  names(Y_gene)[names(Y_gene) %in% "alpha"] = "alpha_immune_cells_only"
  
  geneDT <- full_join(X_gene, Y_gene)
  
  plot <- ggplot(geneDT, aes(x = alpha_no_immune_cells, y = alpha_immune_cells_only))+
    geom_vline(xintercept = .5, linetype = 2)+
    geom_hline(yintercept = .5, linetype = 2)+
    geom_point() +
    ggtitle(gene_name)+
    theme_minimal() +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))
  
  longDT <- geneDT %>% pivot_longer(cols = starts_with("alpha"), names_to = "analysis", values_to = "alpha")
  plot2 <- ggplot(longDT, aes(x=pos, y=alpha, fill = analysis)) +
    geom_point(pch=21) +
    theme_linedraw()
  
  return(list(plot=plot, wideDT=geneDT, longDT=longDT, plot2 = plot2))
}

## https://pmc.ncbi.nlm.nih.gov/articles/PMC11603180/
p_NOD2 <- plotGeneQuadrant("NOD2", "16", 50693606L, 50733075L)
p_NOD2$plot; p_NOD2$plot2

plotGeneQuadrant("TNFRSF1A", "12", 6328771L, 6342076L)



# ####################
# ## Annotate genes ##
# ####################
# print("Test enrichment by categories:")
# 
# ## 1. Keep CpGs in regions where at least 5 CpGs are in 50bp distance to each other
# ## 2. annotate with associated genes (in gene body or +/- 10kb from TSS)
# ## 3. run GO term enrichment with clusterProfiler::enrichGO
# 
# ## Create universe
# universe <- annotateCpGs_txdb(
#   clusterCpGs(Z_inner_immvsnoimm$name, max_gap = 50, min_size = 5),
#   tss_window = 10000
# )
# 
# print(paste0("Gene universe contains ", length(universe), " genes"))
# ## "Gene universe contains 30190 genes"
# 
# ## Annotate the different CpGs categories 
# resAnnot_cellUniversal <- CpG_GO_pipeline(cellUniversal, universe = universe)
# resAnnot_immune <- CpG_GO_pipeline(immune, universe = universe)
# resAnnot_notimmune <- CpG_GO_pipeline(notimmune, universe = universe)
# 
# resAnnot <- list(resAnnot_cellUniversal=resAnnot_cellUniversal,
#                  resAnnot_immune=resAnnot_immune,
#                  resAnnot_notimmune=resAnnot_notimmune)
# 
# # 1) Flatten the nested list: groups -> ontologies -> result rows
# # resAnnot names: resAnnot_cellUniversal, resAnnot_immune, resAnnot_notimmune
# grp_label_map <- c(
#   resAnnot_cellUniversal = "cellUniversal",
#   resAnnot_immune        = "immune",
#   resAnnot_notimmune     = "not immune"
# )
# 
# df_all <- purrr::imap(resAnnot, function(ont_list, grp_name) {
#   purrr::imap_dfr(ont_list, function(er, ont_name) {
#     if (is.null(er) || nrow(er@result) == 0) return(tibble())
#     as_tibble(er@result) %>%
#       mutate(group_raw = grp_name,  ontology = ont_name)
#   })
# }) %>% bind_rows()
# 
# # 2) Harmonize labels / factors
# df_all <- df_all %>%
#   mutate(
#     group    = grp_label_map[group_raw],
#     group    = factor(group, levels = c("cellUniversal", "immune", "not immune")),
#     ontology = factor(ontology, levels = c("BP", "MF", "CC"))
#   )
# 
# # 3) Keep only terms that are significant AND have counts AND more than 5 genes
# df_sig <- df_all %>%
#   filter(!is.na(p.adjust) & p.adjust < 0.01) %>%
#   filter(!is.na(Description) & !is.na(FoldEnrichment) & !is.na(Count) & Count > 10 & FoldEnrichment > 2)
# 
# # 4) Optional: reorder Description within each ontology by (high FE, low FDR)
# df_sig <- df_sig %>%
#   group_by(ontology) %>%
#   mutate(Description = fct_reorder2(Description, FoldEnrichment, -p.adjust, .fun = max, .desc = TRUE)) %>%
#   ungroup()
# 
# # 5) Plot: X = group, Y = Description, facet by ontology, size = FE, color = p.adjust
# p <- ggplot(df_sig, aes(x = group, y = Description)) +
#   geom_point(aes(size = FoldEnrichment, color = p.adjust), alpha = 0.9) +
#   scale_size_continuous(name = "FoldEnrichment", range = c(1.5, 8)) +
#   # smaller FDR should look darker; direction = -1 handles that
#   scale_color_viridis_c(name = "FDR (p.adjust)", option = "plasma", direction = -1) +
#   facet_grid(group ~ ontology, scales = "free", space = "free_y") +
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = NULL,
#        title = "GO enrichment across groups and ontologies (FDR < 0.05)") +
#   theme(
#     legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#     legend.background     = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#     legend.key            = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#     legend.position       = "top",
#     axis.text.y           = element_text(size = 8),
#     axis.text.x           = element_text(size = 8, angle = 45, hjust = 1),
#     strip.text            = element_text(face = "bold")
#   )
# 
# print(p)
# 
# p <- ggplot(df_sig[df_sig$ontology %in% "BP",], aes(x = group, y = Description)) +
#   geom_point(aes(size = FoldEnrichment, color = p.adjust), alpha = 0.9) +
#   scale_size_continuous(name = "FoldEnrichment", range = c(1.5, 8)) +
#   # smaller FDR should look darker; direction = -1 handles that
#   scale_color_viridis_c(name = "FDR (p.adjust)", option = "plasma", direction = -1) +
#   facet_grid(group ~ ontology, scales = "free", space = "free_y") +
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = NULL,
#        title = "GO enrichment across groups and ontologies (FDR < 0.05)") +
#   theme(
#     legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#     legend.background     = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#     legend.key            = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#     legend.position       = "top",
#     axis.text.y           = element_text(size = 8),
#     axis.text.x           = element_text(size = 8, angle = 45, hjust = 1),
#     strip.text            = element_text(face = "bold")
#   )
# 
# pdf(here("05_hvCpGalgorithm/figures/GOenrichment.pdf"), width = 10, height = 7)
# p
# dev.off()
# 
# ## Nothing very clear...