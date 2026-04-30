####################################################################
## Plot secondary analyses results of algorithm ran on atlas data ##
####################################################################
library(here)

## Load libraries
if (!exists("libLoaded")) {
source(here("B_MultiTissues", "quiet_library.R"))}

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Load array results for comparison
if (!exists("resCompArray")) {
  source(here("B_MultiTissues/03_exploreResults/S02_analyseResultsArray_local.R"))}

## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}

#####################################
## one run per developmental layer ##
#####################################
reload = FALSE
if (reload){
  endo = readRDS(here::here("gitignore/fullres_10X_12_endo"))
  meso = readRDS(here::here("gitignore/fullres_10X_13_meso"))
  ecto = readRDS(here::here("gitignore/fullres_10X_14_ecto"))
  allLayers = readRDS(here::here("gitignore/fullres_Atlas10X"))
  
  WGBS_Array_datasets <- read.csv(here("05_hvCpGalgorithm/figures/WGBS_Array_datasets.csv"))
  
  table(WGBS_Array_datasets[WGBS_Array_datasets$assay %in% "atlas", "Germ.layer"])
  # ectoderm endoderm mesoderm 
  # 6       21       19 
  
  ## Ran with 6 groups only to match between layers
  endo_6gp = readRDS(here::here("gitignore/fullres_10X_12.2_endo6gp"))
  meso_6gp = readRDS(here::here("gitignore/fullres_10X_13.2_meso6gp"))
}

######################################################################
## The results are comparable when the number of groups is reduced: ##
######################################################################

## Endoderm
if (!file.exists(here::here("05_hvCpGalgorithm/figures/correlations/correlation_endoFullvsReduced6gp.pdf"))){
  ## Use data table to handle large data
  setDT(endo_6gp)
  setDT(endo)
  
  x <- endo_6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
  y <- endo[, .(name, alpha_endo = alpha)]; setkey(y, name)
  
  m <- x[y, nomatch = 0]   # keeps matched names only
  mycor <- cor(m$alpha_6gp, m$alpha_endo, use = "complete.obs")
  set.seed(1234)
  p <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_endo))+
    geom_point(pch = 21, alpha = 0.1) +
    theme_minimal(base_size = 14) +
    ylim(c(0,1)) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable in WGBS atlas endoderm cell types",
         subtitle = "(100k random CpG plotted)",
         x = "Pr(hv) calculated on a subset of cell types (N=6)",
         y = "Pr(hv) calculated on all cell types (N=21)")
  
  ggplot2::ggsave(
    filename = here::here(paste0("05_hvCpGalgorithm/figures/correlations/correlation_endoFullvsReduced6gp.pdf")),
    plot = p, width = 8, height = 8
  )
}

## Mesoderm
if (!file.exists(here::here("05_hvCpGalgorithm/figures/correlations/correlation_mesoFullvsReduced6gp.pdf"))){
  setDT(meso_6gp)
  setDT(meso)
  
  x <- meso_6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
  y <- meso[, .(name, alpha_meso = alpha)]; setkey(y, name)
  
  m <- x[y, nomatch = 0]   # keeps matched names only
  mycor <- cor(m$alpha_6gp, m$alpha_meso, use = "complete.obs")
  
  set.seed(1234)
  p <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_meso))+
    geom_point(pch = 21, alpha = 0.1) +
    theme_minimal(base_size = 14) +
    ylim(c(0,1)) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable in WGBS atlas mesoderm cell types",
         subtitle = "(100k random CpG plotted)",
         x = "Pr(hv) calculated on a subset of cell types (N=6)", 
         y = "Pr(hv) calculated on all cell types (N=19)")
  
  ggplot2::ggsave(
    filename = here::here(paste0("05_hvCpGalgorithm/figures/correlations/correlation_mesoFullvsReduced6gp.pdf")),
    plot = p, width = 8, height = 8
  )
}

#################################################################################
## Create a table with all CpG sites, pr(hv) for each germ layer, and SNP info ##
#################################################################################
## SNP info: Bash code in B_MultiTissues/03_exploreResults/prepSNP1000GP.sh
# exclude CpGs that overlap with common SNPs in 1000 Genomes, with a MAF‑threshold of 0.5%

redo3layertab = TRUE
if (redo3layertab){
  endoGR <- GRanges(seqnames = endo$chr,
                    ranges = IRanges(start = endo$pos, end = endo$pos),
                    alpha_endo = endo$alpha)
  
  ectoGR <- GRanges(seqnames = ecto$chr,
                    ranges = IRanges(start = ecto$pos, end = ecto$pos),
                    alpha_ecto = ecto$alpha)
  
  mesoGR <- GRanges(seqnames = meso$chr,
                    ranges = IRanges(start = meso$pos, end = meso$pos),
                    alpha_meso = meso$alpha)
  
  meso6gpGR <- GRanges(seqnames = meso_6gp$chr,
                       ranges = IRanges(start = meso_6gp$pos, end = meso_6gp$pos),
                       alpha_meso6gp = meso_6gp$alpha)
  
  endo6gpGR <- GRanges(seqnames = endo_6gp$chr,
                       ranges = IRanges(start = endo_6gp$pos, end = endo_6gp$pos),
                       alpha_endo6gp = endo_6gp$alpha)
  
  allLayersGR <- GRanges(seqnames = allLayers$chr,
                         ranges = IRanges(start = allLayers$pos, end = allLayers$pos),
                         alpha_allLayers = allLayers$alpha)
  
  ## 1. Create union of all unique CpG positions
  table3layers <- union(union(union(union(allLayersGR, union(ectoGR, mesoGR)), endoGR), meso6gpGR), endo6gpGR)
  
  ## 2. Use findOverlaps to map alpha values back
  # we want endoGR[i] -> table3layers[endoHits[[i]]] for each i, etc.
  
  endoHits <- findOverlaps(endoGR, table3layers, select = "first")
  ectoHits <- findOverlaps(ectoGR, table3layers, select = "first")
  mesoHits <- findOverlaps(mesoGR, table3layers, select = "first")
  allLayersHits <- findOverlaps(allLayersGR, table3layers, select = "first")
  meso6gpHits <- findOverlaps(meso6gpGR, table3layers, select = "first")
  endo6gpHits <- findOverlaps(endo6gpGR, table3layers, select = "first")
  
  # initialize columns with NA
  mcols(table3layers)$alpha_endo <- NA_real_
  mcols(table3layers)$alpha_ecto <- NA_real_
  mcols(table3layers)$alpha_meso <- NA_real_
  mcols(table3layers)$alpha_allLayers <- NA_real_
  mcols(table3layers)$alpha_meso6gp <- NA_real_
  mcols(table3layers)$alpha_endo6gp <- NA_real_
  
  # copy only the hits
  mcols(table3layers)$alpha_endo[endoHits]   <- mcols(endoGR)$alpha_endo
  mcols(table3layers)$alpha_ecto[ectoHits]   <- mcols(ectoGR)$alpha_ecto
  mcols(table3layers)$alpha_meso[mesoHits]   <- mcols(mesoGR)$alpha_meso
  mcols(table3layers)$alpha_allLayers[allLayersHits]   <- mcols(allLayersGR)$alpha_allLayers
  mcols(table3layers)$alpha_meso6gp[meso6gpHits]   <- mcols(meso6gpGR)$alpha_meso6gp
  mcols(table3layers)$alpha_endo6gp[endo6gpHits]   <- mcols(endo6gpGR)$alpha_endo6gp
  
  ## Add SNP info from 1000 genomes project (prepared in prepSNP1000GP.sh)
  snp_set_chr <- readLines(here("gitignore/snps_maf_gt_0p5pct_chr.txt"))
  snp_set_pos <- as.numeric(readLines(here("gitignore/snps_maf_gt_0p5pct_pos.txt")))
  
  snps_GRanges <- GRanges(
    seqnames = snp_set_chr,
    ranges   = IRanges(start = snp_set_pos, end = snp_set_pos)
  )
  
  hits <- GenomicRanges::findOverlaps(table3layers, snps_GRanges, select = "first")
  
  isSNP <- !is.na(hits)   # TRUE wherever there is a hit
  mcols(table3layers)$isSNP <- isSNP
  
  rm(endo, ecto, meso, endo_6gp, meso_6gp, allLayers, 
     endoGR, ectoGR, mesoGR, allLayersGR, endo6gpGR,  meso6gpGR,
     endoHits, ectoHits, mesoHits, allLayersHits, endo6gpHits, meso6gpHits,
     snp_set_chr, snp_set_pos,snps_GRanges, isSNP)
  
  ## Add chr_pos column to identify positions
  table3layers$chr_pos <- paste0("chr", table3layers@seqnames, "_", table3layers@ranges@start)
  
  ### SAVED ###
  save(table3layers, file = here("gitignore/fullTable3layers.Rda"))
}

if (!exists("table3layers")){
  load(file = here("gitignore/fullTable3layers.Rda"))  
}

########################
## How many are SNPs? ##
########################
a = sum(table3layers[!is.na(table3layers$alpha_endo) & 
                       !is.na(table3layers$alpha_meso) & 
                       !is.na(table3layers$alpha_ecto), ]$isSNP) 

b = length(table3layers[!is.na(table3layers$alpha_endo) & 
                          !is.na(table3layers$alpha_meso) & 
                          !is.na(table3layers$alpha_ecto), ])
message(paste0(round(a/b,2)*100, "% CpGs sequenced in the 3 germ layers are SNPs (", a, " out of ", b, ")"))

a = sum(table3layers[!is.na(table3layers$alpha_endo) & table3layers$alpha_endo >= .9 & 
                       !is.na(table3layers$alpha_meso) & table3layers$alpha_meso >= .9 & 
                       !is.na(table3layers$alpha_ecto) & table3layers$alpha_ecto >= .9, ]$isSNP) 

b = length(table3layers[!is.na(table3layers$alpha_endo) & table3layers$alpha_endo >= .9 & 
                          !is.na(table3layers$alpha_meso) & table3layers$alpha_meso >= .9 & 
                          !is.na(table3layers$alpha_ecto) & table3layers$alpha_ecto >= .9, ])
message(paste0(round(a/b,2)*100, "% CpGs sequenced in the 3 germ layers with pr(hv) >=90% are SNPs (", a, " out of ", b, ")"))

# 9% CpGs sequenced in the 3 germ layers are SNPs (1559985 out of 18195276)
# 52% CpGs sequenced in the 3 germ layers with pr(hv) >=90% are SNPs (94627 out of 182710)
## --> Variability could simply come from SNPs!! We need to exclude them.

#######################################
## Save intersection for alpha > 90% ##
#######################################

top90SNPrm <- table3layers[!is.na(table3layers$alpha_endo) & table3layers$alpha_endo >= .9 & 
                             !is.na(table3layers$alpha_meso) & table3layers$alpha_meso >= .9 & 
                             !is.na(table3layers$alpha_ecto) & table3layers$alpha_ecto >= .9 &
                             !table3layers$isSNP, ]$chr_pos

# makeOverlapAndAll90pc <- function(table3layers, which){
exactly_one <- table3layers[((!is.na(table3layers$alpha_endo) & table3layers$alpha_endo >= .9) | 
                               (!is.na(table3layers$alpha_meso) & table3layers$alpha_meso >= .9) | 
                               (!is.na(table3layers$alpha_ecto) & table3layers$alpha_ecto >= .9)) &
                              !table3layers$isSNP, ]$chr_pos
  
  
  message(paste0("Sites with alpha>=0.9 in exactly 1 group: ", length(exactly_one)))
  message(paste0("Sites with alpha>=0.9 in exactly 2 groups: ", length(exactly_two)))
  message(paste0("Total sites in 1 or 2 groups: ", length(sites_1or2)))
  
  
  
  
  if (which == "all"){
    tab <- table3layers[!is.na(table3layers$alpha_endo) & table3layers$alpha_endo >= .9 & 
                          !is.na(table3layers$alpha_meso) & table3layers$alpha_meso >= .9 & 
                          !is.na(table3layers$alpha_ecto) & table3layers$alpha_ecto >= .9 &
                          !table3layers$isSNP, ]
  } else if (which == "6gp"){
    
    
    
    
    # 1) Overlap among all names (unfiltered)
    sets_unfilt <- lapply(list(endo, ecto, meso, allLayers), function(df) df$name)
    overlap <- Reduce(intersect, sets_unfilt)
    
    message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
    # There are 17474840 overlapping CpGs between these groups (unfiltered).
    
    # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
    sets <- lapply(list(endo, ecto, meso, allLayers),
                   function(df) df$name[df$alpha >= 0.9 & df$name %in% overlap])
    names(sets) <- c("endo", "ecto", "meso", "allLayers")
    
    # 3) Select the intersection
    topIntersect90 <- Reduce(intersect, sets)
    
    # Exactly in 1 group
    only_endo <- setdiff(sets$endo, union(sets$ecto, sets$meso))
    only_ecto <- setdiff(sets$ecto, union(sets$endo, sets$meso))
    only_meso <- setdiff(sets$meso, union(sets$endo, sets$ecto))
    exactly_one <- unique(c(only_endo, only_ecto, only_meso))
    
    # Exactly in 2 groups
    endo_ecto <- intersect(sets$endo, sets$ecto)
    endo_ecto <- setdiff(endo_ecto, sets$meso)
    endo_meso <- intersect(sets$endo, sets$meso)
    endo_meso <- setdiff(endo_meso, sets$ecto)
    ecto_meso <- intersect(sets$ecto, sets$meso)
    ecto_meso <- setdiff(ecto_meso, sets$endo)
    exactly_two <- unique(c(endo_ecto, endo_meso, ecto_meso))
    
    # Combined: sites in exactly 1 or 2 groups
    sites_1or2 <- unique(c(exactly_one, exactly_two))
    
    message(paste0("Sites with alpha>=0.9 in exactly 1 group: ", length(exactly_one)))
    message(paste0("Sites with alpha>=0.9 in exactly 2 groups: ", length(exactly_two)))
    message(paste0("Total sites in 1 or 2 groups: ", length(sites_1or2)))
    
    if (which == "all"){
      ## Save
      saveRDS(overlap, file = here("gitignore/overlapLayers.RDS"))
      saveRDS(sites_1or2, file = here("gitignore/overlap.1or2.RDS"))
      saveRDS(topIntersect90, file = here("gitignore/topIntersect90.RDS"))
    } else if (which == "6gp"){
      ## Save
      saveRDS(overlap, file = here("gitignore/overlapLayers_6gp.RDS"))
      saveRDS(sites_1or2, file = here("gitignore/overlap.1or2_6gp.RDS"))
      saveRDS(topIntersect90, file = here("gitignore/topIntersect90_6gp.RDS"))
    }
    
  }
  
  makeOverlapAndAll90pc(endo, ecto, meso, allLayers, which="all")
  
  # makeOverlapAndAll90pc <- function(endo, ecto, meso, allLayers, which){
  #   
  #   # 1) Overlap among all names (unfiltered)
  #   sets_unfilt <- lapply(list(endo, ecto, meso, allLayers), function(df) df$name)
  #   overlap <- Reduce(intersect, sets_unfilt)
  #   
  #   message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
  #   # There are 17474840 overlapping CpGs between these groups (unfiltered).
  #   
  #   # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
  #   sets <- lapply(list(endo, ecto, meso, allLayers),
  #                  function(df) df$name[df$alpha >= 0.9 & df$name %in% overlap])
  #   names(sets) <- c("endo", "ecto", "meso", "allLayers")
  #   
  #   # 3) Select the intersection
  #   topIntersect90 <- Reduce(intersect, sets)
  #   
  #   # Exactly in 1 group
  #   only_endo <- setdiff(sets$endo, union(sets$ecto, sets$meso))
  #   only_ecto <- setdiff(sets$ecto, union(sets$endo, sets$meso))
  #   only_meso <- setdiff(sets$meso, union(sets$endo, sets$ecto))
  #   exactly_one <- unique(c(only_endo, only_ecto, only_meso))
  #   
  #   # Exactly in 2 groups
  #   endo_ecto <- intersect(sets$endo, sets$ecto)
  #   endo_ecto <- setdiff(endo_ecto, sets$meso)
  #   endo_meso <- intersect(sets$endo, sets$meso)
  #   endo_meso <- setdiff(endo_meso, sets$ecto)
  #   ecto_meso <- intersect(sets$ecto, sets$meso)
  #   ecto_meso <- setdiff(ecto_meso, sets$endo)
  #   exactly_two <- unique(c(endo_ecto, endo_meso, ecto_meso))
  #   
  #   # Combined: sites in exactly 1 or 2 groups
  #   sites_1or2 <- unique(c(exactly_one, exactly_two))
  #   
  #   message(paste0("Sites with alpha>=0.9 in exactly 1 group: ", length(exactly_one)))
  #   message(paste0("Sites with alpha>=0.9 in exactly 2 groups: ", length(exactly_two)))
  #   message(paste0("Total sites in 1 or 2 groups: ", length(sites_1or2)))
  #   
  #   if (which == "all"){
  #     ## Save
  #     saveRDS(overlap, file = here("gitignore/overlapLayers.RDS"))
  #     saveRDS(sites_1or2, file = here("gitignore/overlap.1or2.RDS"))
  #     saveRDS(topIntersect90, file = here("gitignore/topIntersect90.RDS"))
  #   } else if (which == "6gp"){
  #     ## Save
  #     saveRDS(overlap, file = here("gitignore/overlapLayers_6gp.RDS"))
  #     saveRDS(sites_1or2, file = here("gitignore/overlap.1or2_6gp.RDS"))
  #     saveRDS(topIntersect90, file = here("gitignore/topIntersect90_6gp.RDS"))
  #   }
  #   
  # }
  # 
  # makeOverlapAndAll90pc(endo, ecto, meso, allLayers, which="all")
  
  ## Check: is topIntersect90 the same if we use only 6 groups per germ layer?
  makeOverlapAndAll90pc(endo = endo_6gp, ecto, meso = meso_6gp, allLayers, which="6gp")
  
  topIntersect90_6gp <- readRDS(here("gitignore/topIntersect90_6gp.RDS"))
  topIntersect90 <- readRDS(here("gitignore/topIntersect90.RDS"))
  
  length(topIntersect90) # 174494
  length(topIntersect90_6gp) # 147516
  length(intersect(topIntersect90, topIntersect90_6gp)) # 124086
  
  length(intersect(topIntersect90, topIntersect90_6gp))/length(topIntersect90)*100
  ## 71% of these sites are found also when only 6 groups are considered per germ layer
  
  ##############################################
  ## What is the overlap for different alpha? ##
  ##############################################
  
  p1 <- plotMyVenn(0.5, endo = endo, meso = meso, ecto = ecto, all = allLayers)
  p2 <- plotMyVenn(0.75, endo = endo, meso = meso, ecto = ecto, all = allLayers)
  p3 <- plotMyVenn(0.90, endo = endo, meso = meso, ecto = ecto, all = allLayers)
  
  pdf(here("05_hvCpGalgorithm/figures/vennGermLayers/Venn_top90.pdf"), width = 16, height = 8)
  cowplot::plot_grid(p1 + theme(legend.position = "none"),
                     p2 + theme(legend.position = "none"),
                     p3 + theme(legend.position = "none"), 
                     nrow = 1, 
                     labels = "Ectoderm (N=6 cell types), mesoderm (N=19), endoderm (N=21), and all combined")
  dev.off()
  
  ##################################################################
  ## Check: is it similar if we use only 6 groups per germ layer? ##
  ##################################################################
  p1 <- plotMyVenn(0.5, endo = endo_6gp, meso = meso_6gp, ecto = ecto, all = allLayers)
  p2 <- plotMyVenn(0.75, endo = endo_6gp, meso = meso_6gp, ecto = ecto, all = allLayers)
  p3 <- plotMyVenn(0.90, endo = endo_6gp, meso = meso_6gp, ecto = ecto, all = allLayers)
  
  pdf(here("05_hvCpGalgorithm/figures/vennGermLayers/Venn_top90_6gp.pdf"), width = 16, height = 8)
  cowplot::plot_grid(p1 + theme(legend.position = "none"),
                     p2 + theme(legend.position = "none"),
                     p3 + theme(legend.position = "none"), 
                     nrow = 1, 
                     labels = "Ectoderm (N=6 cell types), mesoderm (N=6), endoderm (N=6), and all datasets combined (N=46)")
  dev.off()
  
  ####################################################################################
  ## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
  ####################################################################################
  if (!exists("topIntersect90")){
    topIntersect90 <- readRDS(overlap, file = here("gitignore/topIntersect90.RDS"))}
  
  
  
  total <- allLayers$name[allLayers$name %in% endo$name]
  total <- total[total %in% meso$name]
  total <- total[total %in% ecto$name]
  
  top_cpgs <- intersect(
    intersect(allLayers$name[allLayers$alpha > 0.9], endo$name[endo$alpha > 0.9]),
    intersect(meso$name[meso$alpha > 0.9],ecto$name[ecto$alpha > 0.9]))
  print(paste("Found", length(top_cpgs), "overlapping high-alpha CpGs")) # 174494
  
  listGR <- list(top90 = makeGRfromMyCpGPos(vec = top_cpgs, setname = "topCpGs"),
                 allButTop90 = makeGRfromMyCpGPos(total[!total %in% top_cpgs], "allButTop90"))
  
  # ---- Run it (ME sets in putativeME_GR$set will be tested separately)
  res_quadrants <- test_enrichment_quadrants(listGR, putativeME_GR, me_col = "set")
  
  # Order quadrants within each facet by log2OR
  res_plot2 <- res_quadrants %>%
    mutate(
      log2OR = log2(odds_ratio),
      signif  = p_adj_BH < 0.05
    ) %>%
    dplyr::group_by(CpG_set) %>%
    mutate(quadrant_ord = reorder(quadrant, log2OR)) %>%
    ungroup()
  
  plot <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("grey", "black")) +
    labs(
      x = NULL,
      y = expression(log[2]~"(odds ratio)"),
      title = "ME enrichment by group (vs other group)",
      subtitle = "2x2 Fisher's exact test "
    ) +
    facet_wrap(~ CpG_set, scales = "free_x", nrow = 1) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "top",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )
  
  pdf(here("05_hvCpGalgorithm/figures/topCpGsEnrichME.pdf"), width = 6, height = 6)
  plot
  dev.off()
  
  ############################################################
  ## How many of each putative ME is actually in the top90? ##
  ############################################################
  
  # Overlaps of ME ranges with CpG sites (top90 and allButTop90)
  hits_top <- findOverlaps(putativeME_GR, listGR$top90, ignore.strand = TRUE)
  hits_all <- findOverlaps(putativeME_GR, listGR$allButTop90, ignore.strand = TRUE)
  
  # Logical flags per ME range (query index)
  overlap_top <- logical(length(putativeME_GR))
  overlap_top[unique(queryHits(hits_top))] <- TRUE
  
  overlap_all <- logical(length(putativeME_GR))
  overlap_all[unique(queryHits(hits_all))] <- TRUE
  
  # Build a small data.frame with one row per ME range
  df_me <- as.data.frame(mcols(putativeME_GR)) |>
    mutate(
      in_top90       = overlap_top,
      in_allButTop90 = overlap_all
    )
  
  # Summaries per set
  summary_df <- df_me |>
    group_by(set) |>
    summarise(
      n_total = n(),
      n_in_top90 = sum(in_top90),
      pc_in_top90 = n_in_top90/n_total*100,
      n_in_allButTop90 = sum(in_allButTop90),
      pc_in_allButTop90 = n_in_allButTop90/n_total*100,
      n_in_both = sum(in_top90 & in_allButTop90),   # non 0 if region rather than CpG
      pc_in_both = n_in_both/n_total*100,
      n_in_neither = n_total - n_in_top90 - n_in_allButTop90 + n_in_both,
      pc_in_neither = n_in_neither/n_total*100)
  
  ## Format pretty
  summary_df %>%
    mutate(across(starts_with("pc_"), ~ scales::percent(.x / 100, accuracy = 0.1))) %>%
    gt() %>%
    fmt_number(columns = starts_with("n_"), decimals = 0) %>%
    cols_label(
      set = "Putative ME set",
      n_total = "Total ME sites or regions",
      n_in_top90 = "N overlapping top 90 CpGs",
      pc_in_top90 = "%",
      n_in_allButTop90 = "N overlapping remaining CpGs",
      pc_in_allButTop90 = "%",
      n_in_both = "N overlapping both groups",
      pc_in_both = "%",
      n_in_neither = "N not covered in all 3 layers analyses",
      pc_in_neither = "%"
    ) %>% 
    tab_style(
      style = cell_fill(color = "lightblue"),
      locations = cells_column_labels()
    ) %>%
    tab_options(
      table.font.size = 13,
      data_row.padding = px(3)
    )
  
  ## Screenshot saved in figures/topCpGsEnrichME_table.png
  
  ## Are Harris MEs high in some layers? Done on Meso + Endo
  
  
  
  
  
  #### Same but with only 6 groups per germ layer
  
  # #####################################
  # ## one run per developmental layer ##
  # #####################################
  # 
  # endo = readRDS(here::here("gitignore/fullres_10X_12.2_endo6gp"))
  # meso = readRDS(here::here("gitignore/fullres_10X_13_meso"))
  # ecto = readRDS(here::here("gitignore/fullres_10X_14_ecto"))
  # allLayers = readRDS(here::here("gitignore/fullres_Atlas10X"))
  # 
  # WGBS_Array_datasets <- read.csv(here("05_hvCpGalgorithm/figures/WGBS_Array_datasets.csv"))
  # 
  # table(WGBS_Array_datasets[WGBS_Array_datasets$assay %in% "atlas", "Germ.layer"])
  # # ectoderm endoderm mesoderm 
  # # 6       21       19 
  # 
  # ##############################
  # ## new candidate locus Matt ##
  # ##############################
  # 
  # ## Matt's data are in hg19
  # dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))
  # 
  # head(dataMatt)
  # 
  # # --- Download and import hg19 → hg38 chain file ---
  # chain_dir <- here("B_MultiTissues/dataIn")
  # chain_gz <- file.path(chain_dir, "hg19ToHg38.over.chain.gz")
  # chain_file <- file.path(chain_dir, "hg19ToHg38.over.chain")
  # 
  # if (!file.exists(chain_file)) {
  #   message("⬇️  Downloading UCSC liftOver chain file...")
  #   dir.create(chain_dir, showWarnings = FALSE, recursive = TRUE)
  #   download.file(
  #     url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  #     destfile = chain_gz,
  #     quiet = TRUE
  #   )
  #   R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
  # }
  # chain <- import.chain(chain_file)
  # 
  # # --- Liftover (hg19 → hg38) ---
  # dataMatt_gr <- GRanges(
  #   seqnames = dataMatt$chromosome,
  #   ranges = IRanges(start = dataMatt$start, end = dataMatt$end),
  #   alpha_endo = NA,  alpha_meso = NA,  alpha_ecto = NA,  alpha_all = NA,
  #   hg19_chr = dataMatt$chromosome, hg19_start = dataMatt$start, hg19_end = dataMatt$end,
  #   `%change` = dataMatt$`%change`, padj = dataMatt$padj,
  #   TE_chromosome = dataMatt$TE_chromosome, TE_start = dataMatt$TE_start, 
  #   TE_end = dataMatt$TE_end, TE_family = dataMatt$TE_family,TE_type = dataMatt$TE_type
  # )
  # 
  # mapped <- liftOver(dataMatt_gr, chain)
  # 
  # # Keep one-to-one mappings only
  # keep <- lengths(mapped) == 1
  # dataMatt_hg38_gr <- unlist(mapped[keep])
  # 
  # endo_gr <- GRanges(seqnames = paste0("chr", endo$chr),
  #                    ranges = IRanges(start = endo$pos, end = endo$pos), 
  #                    alpha = endo$alpha)
  # ecto_gr <- GRanges(seqnames = paste0("chr", ecto$chr),
  #                    ranges = IRanges(start = ecto$pos, end = ecto$pos), 
  #                    alpha = ecto$alpha)
  # meso_gr <- GRanges(seqnames = paste0("chr", meso$chr),
  #                    ranges = IRanges(start = meso$pos, end = meso$pos), 
  #                    alpha = meso$alpha)
  # all_gr <- GRanges(seqnames = paste0("chr", allLayers$chr),
  #                   ranges = IRanges(start = allLayers$pos, end = allLayers$pos), 
  #                   alpha = allLayers$alpha)
  # 
  # # Find overlapping ranges
  # overlaps <- findOverlaps(dataMatt_hg38_gr, endo_gr)
  # dataMatt_hg38_gr[queryHits(overlaps),]$alpha_endo <- endo_gr[subjectHits(overlaps),]$alpha
  # 
  # # Find overlapping ranges
  # overlaps <- findOverlaps(dataMatt_hg38_gr, meso_gr)
  # dataMatt_hg38_gr[queryHits(overlaps),]$alpha_meso <- meso_gr[subjectHits(overlaps),]$alpha
  # 
  # # Find overlapping ranges
  # overlaps <- findOverlaps(dataMatt_hg38_gr, ecto_gr)
  # dataMatt_hg38_gr[queryHits(overlaps),]$alpha_ecto <- ecto_gr[subjectHits(overlaps),]$alpha
  # 
  # # Find overlapping ranges
  # overlaps <- findOverlaps(dataMatt_hg38_gr, all_gr)
  # dataMatt_hg38_gr[queryHits(overlaps),]$alpha_all <- all_gr[subjectHits(overlaps),]$alpha
  # 
  # hist(dataMatt_hg38_gr$alpha_endo, breaks = 100)
  # hist(dataMatt_hg38_gr$alpha_ecto, breaks = 100)
  # hist(dataMatt_hg38_gr$alpha_meso, breaks = 100)
  # hist(dataMatt_hg38_gr$alpha_all, breaks = 100)
  # 
  # dataMatt_hg38_gr$is50pc <- (dataMatt_hg38_gr$alpha_endo > .5 & dataMatt_hg38_gr$alpha_meso > .5 & 
  #                               dataMatt_hg38_gr$alpha_ecto > .5)
  # dataMatt_hg38_gr$is40pc <- (dataMatt_hg38_gr$alpha_endo > .4 & dataMatt_hg38_gr$alpha_meso > .4 & 
  #                               dataMatt_hg38_gr$alpha_ecto > .4)
  # 
  # write.csv(dataMatt_hg38_gr, here("gitignore/DEGCAGS_intersect_repeats_Alice_withalpha.csv"), 
  #           row.names = F, quote = F)
  # 
  # #######################################
  # ## Save intersection for alpha > 90% ##
  # #######################################
  # 
  # # 1) Overlap among all names (unfiltered)
  # sets_unfilt <- lapply(list(endo, ecto, meso, allLayers), function(df) df$name)
  # overlap <- Reduce(intersect, sets_unfilt)
  # 
  # message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
  # # There are 17474840 overlapping CpGs between these groups (unfiltered).
  # 
  # # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
  # sets <- lapply(list(endo, ecto, meso, allLayers),
  #                function(df) df$name[df$alpha >= 0.9 & df$name %in% overlap])
  # names(sets) <- c("endo", "ecto", "meso", "allLayers")
  # 
  # # 3) Select the intersection
  # topIntersect90 <- Reduce(intersect, sets)
  # 
  # # Exactly in 1 group
  # only_endo <- setdiff(sets$endo, union(sets$ecto, sets$meso))
  # only_ecto <- setdiff(sets$ecto, union(sets$endo, sets$meso))
  # only_meso <- setdiff(sets$meso, union(sets$endo, sets$ecto))
  # exactly_one <- unique(c(only_endo, only_ecto, only_meso))
  # 
  # # Exactly in 2 groups
  # endo_ecto <- intersect(sets$endo, sets$ecto)
  # endo_ecto <- setdiff(endo_ecto, sets$meso)
  # endo_meso <- intersect(sets$endo, sets$meso)
  # endo_meso <- setdiff(endo_meso, sets$ecto)
  # ecto_meso <- intersect(sets$ecto, sets$meso)
  # ecto_meso <- setdiff(ecto_meso, sets$endo)
  # exactly_two <- unique(c(endo_ecto, endo_meso, ecto_meso))
  # 
  # # Combined: sites in exactly 1 or 2 groups
  # sites_1or2 <- unique(c(exactly_one, exactly_two))
  # 
  # message(paste0("Sites with alpha>=0.9 in exactly 1 group: ", length(exactly_one)))
  # message(paste0("Sites with alpha>=0.9 in exactly 2 groups: ", length(exactly_two)))
  # message(paste0("Total sites in 1 or 2 groups: ", length(sites_1or2)))
  # 
  # ## Save
  # saveRDS(overlap, file = here("gitignore/overlapLayers.RDS"))
  # saveRDS(sites_1or2, file = here("gitignore/overlap.1or2.RDS"))
  # saveRDS(topIntersect90, file = here("gitignore/topIntersect90.RDS"))
  # 
  # ##############################################
  # ## What is the overlap for different alpha? ##
  # ##############################################
  # 
  # # Compute overlap across any number of groups and plot a Venn diagram
  # plotMyVenn <- function(cutoff, ...) {
  #   groups <- list(...)                         # list of data.frames (each has name, alpha)
  #   
  #   # 1) Overlap among all names (unfiltered)
  #   sets_unfilt <- lapply(groups, function(df) df$name)
  #   overlap <- Reduce(intersect, sets_unfilt)
  #   
  #   message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
  #   
  #   # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
  #   sets <- lapply(groups, function(df) df$name[df$alpha >= cutoff & df$name %in% overlap])
  #   
  #   # 3) Optional: add names to sets if you passed named arguments
  #   if (!is.null(dots <- match.call(expand.dots = FALSE)$...) && length(names(dots))) {
  #     nm <- names(dots)
  #     if (any(nzchar(nm))) names(sets) <- ifelse(nzchar(nm), nm, paste0("group", seq_along(sets)))
  #   }
  #   
  #   # 4) Plot
  #   p <- ggVennDiagram(sets, label_alpha = 0) +
  #     scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
  #     ggtitle(paste0("Pr(hv) ≥ ", cutoff), 
  #             subtitle = paste0("Sequenced CpGs N = ", length(overlap)))
  #   
  #   return(p)
  # }
  # 
  # p1 <- plotMyVenn(0.5, endo = endo, meso = meso, ecto = ecto, all = allLayers)
  # p2 <- plotMyVenn(0.75, endo = endo, meso = meso, ecto = ecto, all = allLayers)
  # p3 <- plotMyVenn(0.90, endo = endo, meso = meso, ecto = ecto, all = allLayers)
  # 
  # pdf(here("05_hvCpGalgorithm/figures/vennGermLayers/Venn_top90.pdf"), width = 16, height = 8)
  # cowplot::plot_grid(p1 + theme(legend.position = "none"),
  #                    p2 + theme(legend.position = "none"),
  #                    p3 + theme(legend.position = "none"), 
  #                    nrow = 1, 
  #                    labels = "Ectoderm (N=6 cell types), mesoderm (N=19), endoderm (N=21), and all combined")
  # dev.off()
  # 
  # ####################################################################################
  # ## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
  # ####################################################################################
  # total <- allLayers$name[allLayers$name %in% endo$name]
  # total <- total[total %in% meso$name]
  # total <- total[total %in% ecto$name]
  # 
  # top_cpgs <- intersect(
  #   intersect(allLayers$name[allLayers$alpha > 0.9], endo$name[endo$alpha > 0.9]),
  #   intersect(meso$name[meso$alpha > 0.9],ecto$name[ecto$alpha > 0.9]))
  # print(paste("Found", length(top_cpgs), "overlapping high-alpha CpGs")) # 174494
  # 
  # listGR <- list(top90 = makeGRfromMyCpGPos(vec = top_cpgs, setname = "topCpGs"),
  #                allButTop90 = makeGRfromMyCpGPos(total[!total %in% top_cpgs], "allButTop90"))
  # 
  # # ---- Run it (ME sets in putativeME_GR$set will be tested separately)
  # res_quadrants <- test_enrichment_quadrants(listGR, putativeME_GR, me_col = "set")
  # 
  # # Order quadrants within each facet by log2OR
  # res_plot2 <- res_quadrants %>%
  #   mutate(
  #     log2OR = log2(odds_ratio),
  #     signif  = p_adj_BH < 0.05
  #   ) %>%
  #   group_by(ME_set) %>%
  #   mutate(quadrant_ord = reorder(quadrant, log2OR)) %>%
  #   ungroup()
  # 
  # plot <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
  #   geom_col(width = 0.8) +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  #   scale_fill_manual(values = c("grey", "black")) +
  #   labs(
  #     x = NULL,
  #     y = expression(log[2]~"(odds ratio)"),
  #     title = "ME enrichment by group (vs other group)",
  #     subtitle = "Bars ordered by effect within ME set; dashed line = OR = 1"
  #   ) +
  #   facet_wrap(~ ME_set, scales = "free_x", nrow = 1) +
  #   theme_classic(base_size = 10) +
  #   theme(
  #     axis.text.x = element_text(angle = 30, hjust = 1),
  #     legend.position = "right",
  #     strip.background = element_rect(fill = "white"),
  #     strip.text = element_text(face = "bold")
  #   )
  # 
  # pdf(here("05_hvCpGalgorithm/figures/topCpGsEnrichME.pdf"), width = 9, height = 6)
  # plot
  # dev.off()
  # 
  # ############################################################
  # ## How many of each putative ME is actually in the top90? ##
  # ############################################################
  # 
  # # Overlaps of ME ranges with CpG sites (top90 and allButTop90)
  # hits_top <- findOverlaps(putativeME_GR, listGR$top90, ignore.strand = TRUE)
  # hits_all <- findOverlaps(putativeME_GR, listGR$allButTop90, ignore.strand = TRUE)
  # 
  # # Logical flags per ME range (query index)
  # overlap_top <- logical(length(putativeME_GR))
  # overlap_top[unique(queryHits(hits_top))] <- TRUE
  # 
  # overlap_all <- logical(length(putativeME_GR))
  # overlap_all[unique(queryHits(hits_all))] <- TRUE
  # 
  # # Build a small data.frame with one row per ME range
  # df_me <- as.data.frame(mcols(putativeME_GR)) |>
  #   mutate(
  #     in_top90       = overlap_top,
  #     in_allButTop90 = overlap_all
  #   )
  # 
  # # Summaries per set
  # summary_df <- df_me |>
  #   group_by(set) |>
  #   summarise(
  #     n_total = n(),
  #     n_in_top90 = sum(in_top90),
  #     pc_in_top90 = n_in_top90/n_total*100,
  #     n_in_allButTop90 = sum(in_allButTop90),
  #     pc_in_allButTop90 = n_in_allButTop90/n_total*100,
  #     n_in_both = sum(in_top90 & in_allButTop90),   # non 0 if region rather than CpG
  #     pc_in_both = n_in_both/n_total*100,
  #     n_in_neither = n_total - n_in_top90 - n_in_allButTop90 + n_in_both,
  #     pc_in_neither = n_in_neither/n_total*100)
  # 
  # ## Format pretty
  # summary_df %>%
  #   mutate(across(starts_with("pc_"), ~ scales::percent(.x / 100, accuracy = 0.1))) %>%
  #   gt() %>%
  #   fmt_number(columns = starts_with("n_"), decimals = 0) %>%
  #   cols_label(
  #     set = "Metastable Epiallele Set",
  #     n_total = "Total MEs",
  #     n_in_top90 = "In Top90",
  #     pc_in_top90 = "%",
  #     n_in_allButTop90 = "In Rest",
  #     pc_in_allButTop90 = "%",
  #     n_in_both = "In Both",
  #     pc_in_both = "%",
  #     n_in_neither = "In Neither",
  #     pc_in_neither = "%"
  #   ) %>%
  #   tab_style(
  #     style = cell_fill(color = "lightblue"),
  #     locations = cells_column_labels()
  #   ) %>%
  #   tab_options(
  #     table.font.size = 11,
  #     data_row.padding = px(3)
  #   )
  