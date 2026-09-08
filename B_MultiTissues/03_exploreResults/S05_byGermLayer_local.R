#####################################################################
# S05 — Per-germ-layer hypervariability categories
#
# Four mutually-exclusive categories, each a set of chr_pos:
#   top1pc3layers : top 1% of logBF_per_ds in endo AND in meso AND in ecto (stricter than in S04)
#   top1pcEndoOnly : top 1% in endo  AND bottom 50% in meso AND bottom 50% in ecto
#   top1pcMesoOnly : top 1% in meso  AND bottom 50% in endo AND bottom 50% in ecto
#   top1pcEctoOnly : top 1% in ecto  AND bottom 50% in endo AND bottom 50% in meso
#
# Then run the SAME enrichment tests S04 ran on top99q, once per category:
#   (1) TE enrichment (overall + per repClass)
#   (2) SD-ASM enrichment (overall + per classification)
#   (3) GO enrichment (CpG-density-controlled)
#   (4) protocadherin-locus enrichment
#####################################################################

#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Add previous MEs including Maria's results
## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}
#####################################################################

variant  <- "SNP_SDASMrm"
prep_dir <- here("gitignore/resultsAtlasPrepared", variant)

## Object created in S04 (GRanges of CpGs covered in all 3 layers, with per-layer scores)
load(here("gitignore/table3layers_coveredIn3_26_08_26.Rda"))
dt <- as.data.table(table3layers_coveredIn3)           # one row per covered-in-3 CpG
setnames(dt, "chr_pos", "chr_pos", skip_absent = TRUE) # ensure chr_pos exists

## ── 1. Define the four categories (all from per-layer top-1% / bottom-50% flags) ──
q_top <- function(x) quantile(x, 0.99, na.rm = TRUE)   # top 1% within a layer
q_bot <- function(x) quantile(x, 0.50, na.rm = TRUE)   # bottom 50% cutoff (median)

## define top and bottom threshold in terms of logBF_per_ds
thr <- list(
  endo = c(top = unname(q_top(dt$logBF_per_ds_endo)), bot = unname(q_bot(dt$logBF_per_ds_endo))),
  meso = c(top = unname(q_top(dt$logBF_per_ds_meso)), bot = unname(q_bot(dt$logBF_per_ds_meso))),
  ecto = c(top = unname(q_top(dt$logBF_per_ds_ecto)), bot = unname(q_bot(dt$logBF_per_ds_ecto))),
  endo6gp = c(top = unname(q_top(dt$logBF_per_ds_endo6gp)), bot = unname(q_bot(dt$logBF_per_ds_endo6gp))),
  meso6gp = c(top = unname(q_top(dt$logBF_per_ds_meso6gp)), bot = unname(q_bot(dt$logBF_per_ds_meso6gp))),
  ecto6gp = c(top = unname(q_top(dt$logBF_per_ds_ecto)), bot = unname(q_bot(dt$logBF_per_ds_ecto)))
)

dt[, `:=`(
  hv_endo = logBF_per_ds_endo >= thr$endo["top"],   
  hv_meso = logBF_per_ds_meso >= thr$meso["top"],
  hv_ecto = logBF_per_ds_ecto >= thr$ecto["top"],
  lo_endo = logBF_per_ds_endo <  thr$endo["bot"],
  lo_meso = logBF_per_ds_meso <  thr$meso["bot"],
  lo_ecto = logBF_per_ds_ecto <  thr$ecto["bot"],
  hv_endo6gp = logBF_per_ds_endo6gp >= thr$endo6gp["top"],   
  hv_meso6gp = logBF_per_ds_meso6gp >= thr$meso6gp["top"],
  hv_ecto6gp = logBF_per_ds_ecto >= thr$ecto6gp["top"],
  lo_endo6gp = logBF_per_ds_endo6gp <  thr$endo6gp["bot"],
  lo_meso6gp = logBF_per_ds_meso6gp <  thr$meso6gp["bot"],
  lo_ecto6gp = logBF_per_ds_ecto <  thr$ecto6gp["bot"]
)]

categories <- list(
  top1pc3layers = dt[hv_endo & hv_meso & hv_ecto, chr_pos],   # top 1% in ALL three layers
  top1pcEndoOnly = dt[hv_endo & lo_meso & lo_ecto, chr_pos],   # top 1% in 1 layer, bottom 50% in the others
  top1pcMesoOnly = dt[hv_meso & lo_endo & lo_ecto, chr_pos],
  top1pcEctoOnly = dt[hv_ecto & lo_endo & lo_meso, chr_pos],
  constitutive = dt[lo_ecto & lo_endo & lo_meso, chr_pos],
  top1pc3layers6gp = dt[hv_endo6gp & hv_meso6gp & hv_ecto6gp, chr_pos],
  endo6gp_only = dt[hv_endo6gp & lo_meso6gp & lo_ecto6gp, chr_pos],
  meso6gp_only = dt[hv_meso6gp & lo_endo6gp & lo_ecto6gp, chr_pos],
  ecto6gp_only = dt[hv_ecto6gp & lo_endo6gp & lo_meso6gp, chr_pos]
)

message("Category sizes:")
print(sapply(categories, length))
# top1pc3layers    top1pcEndoOnly    top1pcMesoOnly    top1pcEctoOnly constitutive 
# 60424           18           63          297      6723417     
# top1pc3layers6gp endo6gp_only meso6gp_only  ecto6gp_only 
#      26        47166        51527           48122 

## =============================================================================
## Test overlap between the 1% in each layer
## =============================================================================

if (!file.exists(here("B_MultiTissues/dataOut/figures/script05/venn_scatter_1_50pc.pdf"))){
  cat1pc <- list(
    endo1pc = dt[(hv_endo), chr_pos],
    meso1pc = dt[(hv_meso), chr_pos],
    ecto1pc = dt[(hv_ecto), chr_pos])
  
  print(sapply(cat1pc, length))
  # endo1pc meso1pc ecto1pc 
  # 202467  202467  202467 
  
  p_venn1pc <- ggVennDiagram(cat1pc,
                             category.names = c("Endo 1%", "Meso 1%", "Ecto 1%"),
                             label = "both", label_alpha = 0) +
    scale_fill_gradient(low = "grey95", high = "#2166AC") +
    scale_colour_manual(values = rep("grey30", 3)) +
    labs(title = "Overlap of top-1% hvCpGs across germ layers") +
    theme(legend.position = "none")
  
  q50 <- function(x) quantile(x, 0.50, na.rm = TRUE)
  
  cat50pc <- list(
    endo50pc = dt[logBF_per_ds_endo >= q50(logBF_per_ds_endo), chr_pos],
    meso50pc = dt[logBF_per_ds_meso >= q50(logBF_per_ds_meso), chr_pos],
    ecto50pc = dt[logBF_per_ds_ecto >= q50(logBF_per_ds_ecto), chr_pos])
  
  p_venn50pc <- ggVennDiagram(cat50pc, 
                              category.names = c("Endo 50%", "Meso 50%", "Ecto 50%"),
                              label = "both", label_alpha = 0) +
    scale_fill_gradient(low = "grey95", high = "#B2182B") +   # different hue from the 1% Venn
    scale_colour_manual(values = rep("grey30", 3)) +
    labs(title = "Overlap of top-50% hvCpGs across germ layers") +
    theme(legend.position = "none")
  
  p_venn50pc
  
  # endo vs meso, with the thresholds that define the categories
  p_scatter_1 <- ggplot(dt, aes(logBF_per_ds_endo, logBF_per_ds_meso)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(trans = "log10", name = "CpGs") +
    # meso thresholds (y): top 1% and bottom 50%
    geom_hline(yintercept = thr$meso["top"], colour = "#B2182B", linetype = 2) +
    geom_hline(yintercept = thr$meso["bot"], colour = "grey40",  linetype = 3) +
    # endo thresholds (x)
    geom_vline(xintercept = thr$endo["top"], colour = "#B2182B", linetype = 2) +
    geom_vline(xintercept = thr$endo["bot"], colour = "grey40",  linetype = 3) +
    annotate("rect", xmin = thr$endo["top"], xmax = Inf,
             ymin = -Inf, ymax = thr$meso["bot"],
             fill = NA, colour = "black", linewidth = 0.6) +   # the "endo-specific" corner
    labs(x = "Hypervariability score endoderm", y = "Hypervariability score mesoderm",
         title = "Hypervariability scores are strongly correlated") +
    theme_minimal(base_size = 12)
  
  p_scatter_2 <- ggplot(dt, aes(logBF_per_ds_endo, logBF_per_ds_ecto)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(trans = "log10", name = "CpGs") +
    # ecto thresholds (y): top 1% and bottom 50%
    geom_hline(yintercept = thr$ecto["top"], colour = "#B2182B", linetype = 2) +
    geom_hline(yintercept = thr$ecto["bot"], colour = "grey40",  linetype = 3) +
    # endo thresholds (x)
    geom_vline(xintercept = thr$endo["top"], colour = "#B2182B", linetype = 2) +
    geom_vline(xintercept = thr$endo["bot"], colour = "grey40",  linetype = 3) +
    annotate("rect", xmin = thr$endo["top"], xmax = Inf,
             ymin = -Inf, ymax = thr$ecto["bot"],
             fill = NA, colour = "black", linewidth = 0.6) +   # the "endo-specific" corner
    labs(x = "Hypervariability score endoderm", y = "Hypervariability score ectoderm",
         title = "Hypervariability scores are strongly correlated") +
    theme_minimal(base_size = 12)
  
  p_scatter_3 <- ggplot(dt, aes(logBF_per_ds_ecto, logBF_per_ds_meso)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(trans = "log10", name = "CpGs") +
    # meso thresholds (y): top 1% and bottom 50%
    geom_hline(yintercept = thr$meso["top"], colour = "#B2182B", linetype = 2) +
    geom_hline(yintercept = thr$meso["bot"], colour = "grey40",  linetype = 3) +
    # ecto thresholds (x)
    geom_vline(xintercept = thr$ecto["top"], colour = "#B2182B", linetype = 2) +
    geom_vline(xintercept = thr$ecto["bot"], colour = "grey40",  linetype = 3) +
    annotate("rect", xmin = thr$ecto["top"], xmax = Inf,
             ymin = -Inf, ymax = thr$meso["bot"],
             fill = NA, colour = "black", linewidth = 0.6) +   # the "ecto-specific" corner
    labs(x = "Hypervariability score ectoderm", y = "Hypervariability score mesoderm",
         title = "Hypervariability scores are strongly correlated") +
    theme_minimal(base_size = 12)
  
  ggsave(here("B_MultiTissues/dataOut/figures/script05/venn_scatter_1_50pc.pdf"),
         (p_venn1pc | p_venn50pc) / 
           (p_scatter_1 | p_scatter_2 | p_scatter_3), width = 18, height = 10)
}

## =============================================================================
## Test power of detection if reduced N of group (Ecto has only 6 groups)
## --> low (1% only, change easily)
## =============================================================================

# full = *_only (3-layer design)
# 6gp  = *6gp_only (6-group design)
# shared = intersection between full and 6gp sets

calc_stats <- function(full_vec, gp6_vec, label) {
  full_set  <- unique(full_vec)
  gp6_set   <- unique(gp6_vec)
  shared_set<- intersect(full_set, gp6_set)
  
  n_full   <- length(full_set)
  n_6gp    <- length(gp6_set)
  n_shared <- length(shared_set)
  
  pct_shared_of_full <- if (n_full > 0) n_shared / n_full * 100 else NA
  pct_full_of_6gp    <- if (n_6gp > 0) n_full / n_6gp * 100 else NA
  
  message(sprintf(
    "%s: full=%d, 6gp=%d, shared=%d (%.0f%% of full, %.1f%% of 6gp)",
    label, n_full, n_6gp, n_shared, pct_shared_of_full, pct_full_of_6gp
  ))
}

calc_stats(categories$top1pcEndoOnly,    categories$endo6gp_only,    "Endo_specific")
calc_stats(categories$top1pcMesoOnly,    categories$meso6gp_only,    "Meso_specific")
calc_stats(categories$top1pcEctoOnly,    categories$ecto6gp_only,    "Ecto_specific")
calc_stats(categories$top1pc3layers,    categories$top1pc3layers6gp,    "top1pc3layers")
# Endo_specific: full=18, 6gp=47166, shared=0 (0% of full, 0.0% of 6gp)
# Meso_specific: full=63, 6gp=51527, shared=0 (0% of full, 0.1% of 6gp)
# Ecto_specific: full=297, 6gp=48122, shared=73 (25% of full, 0.6% of 6gp)
# top1pc3layers: full=60424, 6gp=26, shared=4 (0% of full, 232400.0% of 6gp)

## Remove useless categories
categories <- categories[
  names(categories) %in% 
    c("top1pc3layers", "top1pcEndoOnly", "top1pcMesoOnly", "top1pcEctoOnly", "constitutive")]

## =============================================================================
## TE ENRICHMENT — per category, overall + per repClass
## =============================================================================

if (!file.exists(here("B_MultiTissues/dataOut/figures/script05/TE_enrichment_byCategory.png"))){
  
  ## Shared background = all covered-in-3 CpGs NOT in the focal category
  ## (built per test below, so each category is compared against everything else)
  all_cpg <- dt$chr_pos
  
  ## ── GRanges helpers (reuse S04's makeGRfromMyCpGPos) ─────────────────────────
  cat_gr <- lapply(names(categories), function(nm)
    makeGRfromMyCpGPos(categories[[nm]], nm))
  names(cat_gr) <- names(categories)
  bg_all_gr <- makeGRfromMyCpGPos(all_cpg, "all")
  
  if (!file.exists(here(paste0("gitignore/S05_TE_", variant, ".rds")))){
    library(AnnotationHub)
    ah <- AnnotationHub()
    rmskhg38 <- ah[["AH111333"]]
    te_classes <- c("LINE","SINE","LTR","DNA","RC","Retroposon")
    te_regions <- rmskhg38[mcols(rmskhg38)$repClass %in% te_classes]
    te_by_class <- split(te_regions, mcols(te_regions)$repClass)
    
    # helper: pull scalar fields from fisher_test_te's list output into consistent columns
    te_fields <- function(x, label, category) {
      data.table(
        category   = category,
        te_class   = label,
        fg_in      = x$contingency[1,1],
        fg_out     = x$contingency[1,2],
        bg_in      = x$contingency[2,1],
        bg_out     = x$contingency[2,2],
        odds_ratio = unname(x$odds_ratio),
        pvalue     = x$pvalue,
        conf_low   = x$conf_low,
        conf_high  = x$conf_high
      )
    }
    
    te_res <- rbindlist(lapply(names(categories), function(cat) {
      fg <- cat_gr[[cat]]
      bg <- makeGRfromMyCpGPos(setdiff(all_cpg, categories[[cat]]), "bg")
      overall  <- te_fields(fisher_test_te(reduce(te_regions), fg, bg, "TE_all"), "TE_all", cat)
      perclass <- rbindlist(lapply(names(te_by_class), function(cl)
        te_fields(fisher_test_te(reduce(te_by_class[[cl]]), fg, bg, cl), cl, cat)))
      rbind(overall, perclass)
    }))
    te_res[, p.adj := p.adjust(pvalue, "BH")]
    saveRDS(te_res, here(paste0("gitignore/S05_TE_", variant, ".rds")))
  } else te_res <- readRDS(here(paste0("gitignore/S05_TE_", variant, ".rds")))
  
  print(te_res[order(category, -odds_ratio)])
  
  ## ---- Plot: OR per TE class, faceted by category (S04-style forest) ----
  
  # order TE classes by the top1pc3layers OR so facets share a sensible y-order
  te_plot <- copy(te_res)
  te_plot[, `:=`(sig = ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."))]
  
  ord <- te_plot[category == "top1pc3layers"][order(odds_ratio), te_class]
  te_plot[, te_class := factor(te_class, levels = unique(ord))]
  te_plot[, category := factor(category,
                               levels = c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly"))]
  
  ## Add N
  te_plot[, category := factor(category,
                               levels = c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly"),
                               labels = paste0(c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly"),
                                               "\n(N=", sapply(categories, length)[c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly")], ")"))]
  
  TEplot_S05 <- ggplot(te_plot[te_plot$sig %in% "FDR < 0.05",],
                       aes(odds_ratio, te_class, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 3) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), height = 0.25, orientation = "y") +
    geom_point(size = 2.5) +
    facet_wrap(category ~ ., ncol = 1) +
    scale_x_log10() +
    scale_colour_manual(values = c("FDR < 0.05" = "#DC3220", "n.s." = "grey60")) +
    labs(x = "Odds ratio (category vs rest, log scale)", y = NULL, colour = NULL,
         title = "TE-class enrichment by germ-layer hypervariability category") +
    theme_minimal(base_size = 12)
  
  TEplot_S05
  ggsave(here("B_MultiTissues/dataOut/figures/script05/TE_enrichment_byCategory.png"),
         TEplot_S05, width = 16, height = 5, dpi = 300, bg = "white")
  
  ## Add per sub categories of TE
  
  ## ---- Per TE FAMILY (finer than class) ----
  if (!file.exists(here(paste0("gitignore/S05_TEfam_", variant, ".rds")))){
    fam_tab  <- table(mcols(te_regions)$repFamily)
    fam_keep <- names(fam_tab)[fam_tab >= 1000 & !grepl("\\?$", names(fam_tab))]
    fam_by   <- split(te_regions, mcols(te_regions)$repFamily)
    
    te_fam <- rbindlist(lapply(names(categories), function(cat) {
      fg <- cat_gr[[cat]]
      bg <- makeGRfromMyCpGPos(setdiff(all_cpg, categories[[cat]]), "bg")
      rbindlist(lapply(fam_keep, function(f)
        te_fields(fisher_test_te(reduce(fam_by[[f]]), fg, bg, f), f, cat)))
    }))
    te_fam[, p.adj := p.adjust(pvalue, "BH")]
    saveRDS(te_fam, here(paste0("gitignore/S05_TEfam_", variant, ".rds")))
  } else te_fam <- readRDS(here(paste0("gitignore/S05_TEfam_", variant, ".rds")))
  
  print(te_fam[order(category, -odds_ratio)])
  
  ## ---- Plot: OR per TE class, faceted by category (S04-style forest) ----
  
  # order TE classes by the top1pc3layers OR so facets share a sensible y-order
  te_plot <- copy(te_fam)
  te_plot[, `:=`(sig = ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."))]
  
  ord <- te_plot[category == "top1pc3layers"][order(odds_ratio), te_class]
  te_plot[, te_class := factor(te_class, levels = unique(ord))]
  te_plot[, category := factor(category,
                               levels = c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly"))]
  
  ## Add N
  te_plot[, category := factor(category,
                               levels = c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly"),
                               labels = paste0(c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly"),
                                               "\n(N=", sapply(categories, length)[c("top1pc3layers","top1pcEndoOnly","top1pcMesoOnly","top1pcEctoOnly")], ")"))]
  
  TEplot_fam_S05 <- ggplot(te_plot[te_plot$sig %in% "FDR < 0.05",],
                           aes(odds_ratio, te_class, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 3) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), height = 0.25, orientation = "y") +
    geom_point(size = 2.5) +
    facet_wrap(category ~ ., ncol = 1) +
    scale_x_log10() +
    scale_colour_manual(values = c("FDR < 0.05" = "#DC3220", "n.s." = "grey60")) +
    labs(x = "Odds ratio (category vs rest, log scale)", y = NULL, colour = NULL,
         title = "TE-class enrichment by germ-layer hypervariability category") +
    theme_minimal(base_size = 12)
  
  TEplot_fam_S05
  
  ggsave(here("B_MultiTissues/dataOut/figures/script05/TE_enrichment_byFamily.png"),
         TEplot_fam_S05, width = 16, height = 7, dpi = 300, bg = "white")
}

## =============================================================================
## Compare our results with previous MEs
## =============================================================================

top99q_S04     <- readRDS(here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))  # S04 set
top1pc3layers_S05  <- categories$top1pc3layers                                              # S05 set

n_over <- length(intersect(top99q_S04, top1pc3layers_S05))
data.table(
  n_top99q_S04   = length(top99q_S04),
  n_top1pc3layers_S05 = length(top1pc3layers_S05),
  n_overlap      = n_over,
  pc_of_top99q   = round(100 * n_over / length(top99q_S04), 1),
  pc_of_top1pc3layers = round(100 * n_over / length(top1pc3layers_S05), 1)
)

# n_top99q_S04 n_top1pc3layers_S05 n_overlap pc_of_top99q pc_of_top1pc3layers
#    202467           60424     60424         29.8             100

# use top1pc3layers_S05 instead of the top 1% of script 04

## Focal set = top1pc3layers_S05 (was top99q)
if (!exists("listGR")){
  listGR <- list(top1pc3layers    = makeGRfromMyCpGPos(vec = top1pc3layers_S05, setname = "top1pc3layers"),
                 allButtop1pc3layers = makeGRfromMyCpGPos(
                   setdiff(table3layers_coveredIn3$chr_pos, top1pc3layers_S05), "allButtop1pc3layers"))
}

# Fix chromosome names in geomMeanGR (1 -> chr1)
if (sum(grepl("chr", seqlevels(table3layers_coveredIn3))) == 0){
  seqlevels(table3layers_coveredIn3) <- paste0("chr", seqlevels(table3layers_coveredIn3))
}


## Use the GR object with analyses in the 3 layers
## Focal "our hvCpG set" for this figure = top1pc3layers from S05 (top 1% in ALL three layers)
top1pc3layers_S05 <- categories$top1pc3layers

sets <- list(
  mQTLcontrols     = makeGRfromMyCpGPos(vec = mQTLcontrols_hg38, setname = "mQTLcontrols"),
  HarrisSIV        = HarrisSIV_hg38_GR,
  VanBaakSIV       = VanBaakSIV_hg38_GR,
  VanBaakESS       = VanBaakESS_hg38_GR,
  KesslerSIV       = KesslerSIV_GRanges_hg38,
  GunasekaraCorSIV = corSIV_GRanges_hg38,
  DerakhshanhvCpGs = DerakhshanhvCpGs_hg38_GR
)

## Associate a colour to a group
group_cols <- c(
  "background"         = "#999999",
  "mQTLcontrols"       = "#000000",
  "HarrisSIV"          = RColorBrewer::brewer.pal(8, "Set2")[1],
  "KesslerSIV"         = RColorBrewer::brewer.pal(8, "Set2")[2],
  "DerakhshanhvCpGs"   = RColorBrewer::brewer.pal(8, "Set2")[3],
  "GunasekaraCorSIV"   = RColorBrewer::brewer.pal(8, "Set2")[4],
  "VanBaakESS"         = RColorBrewer::brewer.pal(8, "Set2")[5],
  "top1pc3layers"          = "orange",
  "top99q"             = RColorBrewer::brewer.pal(8, "Set2")[6],
  "VanBaakSIV"         = RColorBrewer::brewer.pal(8, "Set2")[7]
)

# ── shared spacing refinements ───────────────────────────────────────────────
# Added AFTER each theme_minimal()/theme_classic() so it is not overwritten.
# Top margin reserves room for the cowplot panel labels; axis-title margins
# push titles off the tick text.
spacing <- theme(
  plot.margin  = margin(t = 24, r = 10, b = 12, l = 10),
  axis.title.y = element_text(margin = margin(r = 10)),
  axis.title.x = element_text(margin = margin(t = 10))
)

# ── ME overlap ────────────────────────────────────────────
MEsetdt <- make_MEsetdt(sets, GR = table3layers_coveredIn3)

MEsetdt <- na.omit(MEsetdt)
nrow(MEsetdt) ## 52.244

# Set controls as baseline
MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

## Statistical comparisons of alpha between MEs
fit <- lm(logBF_per_ds ~ ME, data = MEsetdt)
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

contrasts <- contrasts %>%
  mutate(ME = contrast,
         ME_name = sub(" - mQTLcontrols$", "", contrast),   # match colour to the ME group being compared
         lower = estimate - 1.96 * SE,
         upper = estimate + 1.96 * SE)

# after fit/emmeans/contrasts are computed on the original MEsetdt,
# reorder ME by mean score FOR THE PLOT
MEsetdt_plot <- copy(MEsetdt)
MEsetdt_plot[, ME := forcats::fct_reorder(ME, logBF_per_ds, .fun = median, na.rm = TRUE)]

# recompute N labels on the reordered data (levels now in mean order)
n_labels <- MEsetdt_plot[, .(n = .N), by = ME]
y_top    <- max(MEsetdt_plot$logBF_per_ds, na.rm = TRUE)

pMElogBF_per_ds <- ggplot(MEsetdt_plot, aes(x = ME, y = logBF_per_ds)) +
  geom_jitter(aes(colour = ME), size = 3, alpha = .2) +
  geom_violin(aes(colour = ME)) +
  geom_boxplot(aes(colour = ME), width = .1) +
  geom_text(data = n_labels,
            aes(x = ME, y = y_top, label = format(n, big.mark = ",")),
            vjust = -0.4, size = 5.5, fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = group_cols, name = "CpG set") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Hypervariability score")

contrasts_plot <- contrasts %>%
  mutate(ME_name = forcats::fct_reorder(ME_name, estimate, .desc = TRUE))

pcontrast <- ggplot(contrasts_plot, aes(x = ME_name, y = estimate, colour = ME_name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_colour_manual(values = group_cols, name = "CpG set") +
  coord_flip() +
  labs(y = "Difference in hypervariability score vs mQTLcontrols", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

## pdecay - legend inside the plot, bottom-left corner
pdecay <- plot_decay_curve(MEsetdt) +
  scale_colour_manual(values = group_cols, name = "CpG set")

# ── Save key objects for S07 ──────────────────────────────────────────────────
saveRDS(MEsetdt, here(paste0("gitignore/MEsetdt_", variant, ".rds")))

####################################################################################
## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
####################################################################################

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

plot_top1pc3layersCpGsEnrichME <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("black", "grey")) +
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
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

print(plot_top1pc3layersCpGsEnrichME)

################################################################################
## Load SIV plots calculated in fetalSIV folder script (in ing-p5)            ##
##                                                                            ##
## top1pc3layers is an EXACT SUBSET of top99q, and the fetal object stores    ##
## per-CpG values (interlayer_r, interindividual_var) — computed per CpG,     ##
## independent of group. So we do NOT re-run the fetal script: we relabel the ##
## top99q rows that belong to top1pc3layers as a new "top1pc3layers" group.   ##
## Panels D (interlayer_corr) and E (CpG_summary) are direct relabels;        ##
## panel F (binned_summary_boot) is re-bootstrapped locally for top1pc3layers.##
################################################################################

plots <- readRDS(here("gitignore/intercorrelationSIVfetal_sepSIV.rds"))

## top1pc3layers CpGs from S05 (chr_pos) -> EPIC CpG ids via the same dico
top1pc3layers_S05 <- categories$top1pc3layers
## `dico` maps CpG <-> chrpos_hg38 (same object the fetal script used).
## If not in memory, load it here (adjust path to wherever the fetal dico lives):
# dico <- readRDS(here("gitignore/EPIC_dico_hg38.rds"))
top1pc3layers_CpG <- dico$CpG[dico$chrpos_hg38 %in% top1pc3layers_S05]
message(length(top1pc3layers_CpG), " top1pc3layers CpGs mapped to EPIC ids")

## helper: spin off an "top1pc3layers" group by subsetting the top99q rows
add_top1pc3layers <- function(df) {
  sub <- df[df$group == "top99q" & df$CpG %in% top1pc3layers_CpG, ]
  sub$group <- "top1pc3layers"
  dplyr::bind_rows(df, sub)     # keep top99q AND add top1pc3layers alongside
}

plots$interlayer_corr <- add_top1pc3layers(plots$interlayer_corr)
plots$CpG_summary     <- add_top1pc3layers(plots$CpG_summary)

## sanity: top1pc3layers should be more systemically correlated than top99q (subset effect)
message("mean interlayer_r  top99q: ",
        round(mean(plots$interlayer_corr$interlayer_r[plots$interlayer_corr$group=="top99q"], na.rm=TRUE),3),
        " | top1pc3layers: ",
        round(mean(plots$interlayer_corr$interlayer_r[plots$interlayer_corr$group=="top1pc3layers"], na.rm=TRUE),3))

## panel F: re-bootstrap ONLY the new top1pc3layers group, same bins/bootstrap as fetal script
boot_median_ci <- function(x, nboot = 1000, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(c(median = NA, low = NA, high = NA))
  bootfun <- function(data, idx) median(data[idx], na.rm = TRUE)
  b  <- boot(x, statistic = bootfun, R = nboot)
  ci <- boot.ci(b, type = "perc", conf = conf)
  if (!is.null(ci) && "percent" %in% names(ci))
    c(median = median(x, na.rm = TRUE), low = ci$percent[4], high = ci$percent[5])
  else c(median = median(x, na.rm = TRUE), low = NA, high = NA)
}

binned_top1pc3layers <- plots$CpG_summary %>%
  dplyr::filter(group == "top1pc3layers") %>%
  mutate(bin = cut(interindividual_var,
                   breaks = seq(0, max(interindividual_var, na.rm = TRUE) + 0.1, by = 0.1),
                   include.lowest = TRUE)) %>%
  group_by(group, bin) %>%
  summarise(boot_res = list(boot_median_ci(interlayer_r)), .groups = "drop") %>%
  mutate(median_r = sapply(boot_res, `[[`, "median"),
         low      = sapply(boot_res, `[[`, "low"),
         high     = sapply(boot_res, `[[`, "high")) %>%
  dplyr::select(-boot_res)

plots$binned_summary_boot <- dplyr::bind_rows(plots$binned_summary_boot, binned_top1pc3layers)

## ── panels ───────────────────────────────────────────────────────────────────
# per-group N for panel D
nD <- as.data.table(plots$interlayer_corr)[, .(n = .N), by = group]
yD_top <- max(plots$interlayer_corr$interlayer_r, na.rm = TRUE)

## Order levels
plots$interlayer_corr <- plots$interlayer_corr %>%
  mutate(group = fct_reorder(group, interlayer_r, .fun = median, .desc = FALSE))

pinterlayer_corr <- ggplot(plots$interlayer_corr,
                           aes(x = group, y = interlayer_r, group = group, fill = group)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.3, fill = "white") +
  geom_text(data = nD,
            aes(x = group, y = yD_top, label = format(n, big.mark = ",")),
            vjust = -0.4, size = 5.5, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = group_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_minimal(base_size = 14) +
  labs(y = "Mean inter-germ layer correlation\n(Pearson's r)")

## Order levels
plots$CpG_summary <- plots$CpG_summary %>%
  mutate(group = fct_reorder(group, interindividual_var, .fun = median, .desc = FALSE))

pinterindividual_var <- ggplot(plots$CpG_summary, aes(x = interindividual_var, fill = group)) +
  geom_density(alpha = .8) +
  scale_fill_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(x = "Interindividual variation")

plots$binned_summary_boot <- plots$binned_summary_boot %>%
  mutate(group = fct_reorder(group, median_r, .fun = median, .na_rm = TRUE, .desc = FALSE))

pbinned <- ggplot(plots$binned_summary_boot,
                  aes(x = bin, y = median_r, color = group, fill = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.2,
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(x = "Interindividual variation",
       y = "Inter-germ layer correlation \n(median ± bootstrap CI)")

upperRow <- plot_grid(
  plot_grid(
    pMElogBF_per_ds + spacing + theme(axis.title.x = element_blank()),
    pcontrast + theme_minimal(base_size = 18) + spacing +
      theme(plot.title = element_text(size=16), legend.position = "none"),
    nrow = 2,
    labels = c("A. Distribution of the hypervariability score for each CpG set",
               "B. Comparison of previous CPG sets groups to mQTLcontrols"),
    label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1),
  pdecay + theme_minimal(base_size = 20) + spacing +
    theme(legend.position = "inside",
          legend.position.inside = c(.7, .6),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3)),
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("", "C. Decay curve of hypervariability score per percentile"),
  label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1
)

SIV_plot <- plot_grid(
  pinterlayer_corr + theme_minimal(base_size = 20) + spacing +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          axis.title.x = element_blank(), legend.position = "none"),
  plot_grid(pinterindividual_var + theme_minimal(base_size = 14) + spacing +
              labs(fill = "CpG set"),
            pbinned + theme_minimal(base_size = 16) + spacing +
              theme(axis.text.x = element_text(angle = 20, hjust = 1))+
              labs(colour = "CpG set", fill = "CpG set"),
            ncol = 1, align = "v",
            labels = c("E. Densities of interindividual variation per CpG within the fetal data, by set",
                       "F. Inter-germ-layer correlation per interindividual variation, binned"),
            label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1),
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("D. Mean inter-germ-layer correlation for each CpG set", ""),
  label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1
)

final_plot <- plot_grid(
  upperRow,
  SIV_plot,
  ncol = 1,
  rel_heights = c(1, 1)
) + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/script05/CompareWithpreviousMEs.png"),
  plot = final_plot, width = 26, height = 20,  dpi = 300, bg = "white")

##############################################################
## How many of each putative ME is actually in all Layers ? ##
##############################################################

if(!file.exists(here("B_MultiTissues/dataOut/figures/script05/topCpGsEnrichME_table2.png"))){
  
  # Universe of covered CpGs, each already labelled top99q vs not
  # (both are single-CpG GRanges built from your covered-in-3 sites)
  top_gr  <- listGR$top1pc3layers      
  rest_gr <- listGR$allButtop1pc3layers
  
  # For each ME set, find which COVERED CpGs fall inside that set's regions,
  # then classify those CpGs as top99q or not.
  me_sets <- split(putativeME_GR, putativeME_GR$set)
  
  summary_df <- rbindlist(lapply(names(me_sets), function(s) {
    regions <- reduce(me_sets[[s]], ignore.strand = TRUE)   # collapse overlapping regions
    
    # covered CpGs inside this set's regions
    top_in  <- sum(overlapsAny(top_gr,  regions, ignore.strand = TRUE))
    rest_in <- sum(overlapsAny(rest_gr, regions, ignore.strand = TRUE))
    n_cov   <- top_in + rest_in                              # total covered CpGs in the set
    
    data.table(
      set           = s,
      n_cpg_covered = n_cov,                                 # covered CpGs the set overlaps
      n_top99q      = top_in,
      pc_top99q     = if (n_cov) 100 * top_in  / n_cov else NA,
      n_rest        = rest_in,
      pc_rest       = if (n_cov) 100 * rest_in / n_cov else NA
    )
  }))
  
  summary_df$set <- as.character(summary_df$set)
  summary_df$set[summary_df$set %in% "Gunasekara corSIV"] <- "Gunasekara CoRSIV"
  
  ## add fold-enrichment vs the genome-wide top99q rate (baseline ≈ 1%)
  baseline <- length(top_gr) / (length(top_gr) + length(rest_gr))   # ~0.01
  summary_df[, fold := (pc_top99q / 100) / baseline]
  
  ## Format pretty
  summary_df %>%
    mutate(
      across(starts_with("pc_"), ~ scales::percent(.x / 100, accuracy = 0.1)),
      fold = scales::number(fold, accuracy = 0.1, suffix = "×")
    ) %>%
    gt() %>%
    fmt_number(columns = starts_with("n_"), decimals = 0) %>%
    cols_label(
      set           = "Previously published CpG set",
      n_cpg_covered = "Covered CpGs in set",
      n_top99q      = "N in all layers hvCpGs",
      pc_top99q     = "%",
      n_rest        = "N in rest (non-all layers hvCpGs)",
      pc_rest       = "%",
      fold          = "Fold enrichment"
    ) %>%
    tab_style(
      style = cell_fill(color = "lightblue"),
      locations = cells_column_labels()
    ) %>%
    tab_options(
      table.font.size = 13,
      data_row.padding = px(3)
    )
  ## Screenshot saved in figures/script05/topCpGsEnrichME_table2.png
}

## =============================================================================
## GO ENRICHMENT — CpG-density-controlled, per category
## =============================================================================
if (!file.exists(here(paste0("gitignore/S05_GO_", variant, ".rds")))){
  totalSites <- all_cpg
  minimum_CpG_per_cluster <- 1
  universe <- annotateCpGs_txdb(
    clusterCpGs(totalSites, max_gap = 50, min_size = minimum_CpG_per_cluster),
    tss_window = 10000)
  
  go_res <- lapply(names(categories), function(cat) {
    cpgs <- categories[[cat]]
    if (length(cpgs) < 20) { message("SKIP GO ", cat, " - too few CpGs (", length(cpgs), ")"); return(NULL) }
    CpG_GO_pipeline_lengthControlled(
      cpgs, universe = universe, min_size = minimum_CpG_per_cluster,
      control_method = "cpg_count", all_sites = totalSites)
  })
  names(go_res) <- names(categories)
  saveRDS(go_res, here(paste0("gitignore/S05_GO_", variant, ".rds")))
} else go_res <- readRDS(here(paste0("gitignore/S05_GO_", variant, ".rds")))

# significant BP terms per category
lapply(go_res[names(go_res) %in% c("top1pcEctoOnly", "top1pcMesoOnly", "top1pcEndoOnly")],
       function(r) if (!is.null(r))
         r$BP@result[r$BP@result$p.adjust < 0.05, ])

# significant MF terms per category
lapply(go_res[names(go_res) %in% c("top1pcEctoOnly", "top1pcMesoOnly", "top1pcEndoOnly")],
       function(r) if (!is.null(r))
         r$MF@result[r$MF@result$p.adjust < 0.05, ])

# significant CC terms per category
lapply(go_res[names(go_res) %in% c("top1pcEctoOnly", "top1pcMesoOnly", "top1pcEndoOnly")],
       function(r) if (!is.null(r))
         r$CC@result[r$CC@result$p.adjust < 0.05, ])

## Nothing significant for layers except for ecto (the least conclusive dataset)

################################################################################
## Extract raw methylation for our 4 main categories                         ##
## (top1pc3layers, top1pcEndoOnly, top1pcMesoOnly, constitutive)              ##
## CpGpos stored in object "categories"                                       ##
################################################################################

## We have 4 samples with endoderm and mesoderm cells
message("Category sizes:")
print(sapply(categories, length))

# ── Write CpG lists for python extraction ─────────────────────────────────────

set.seed(1234)
## top1pc3layers CpGs — subsample 1000 as positive control
top1pc3layers_sample <- sample(categories$top1pc3layers, min(1000, length(categories$top1pc3layers)))

## constitutive_sample — subsample 1000 as negative control
constitutive_sample <- sample(categories$constitutive, min(1000, length(categories$constitutive)))

## Neutral (variance-unfiltered) background
neutral_sample <- sample(dt$chr_pos, min(2000, nrow(dt)))

## Exclude other categories
neutral_sample <- neutral_sample[!neutral_sample %in% 
                 c(categories$top1pcEctoOnly, categories$top1pcMesoOnly,
                   categories$top1pcEndoOnly, categories$top1pc3layers,
                   categories$constitutive)]

neutral_sample <- sample(neutral_sample, 1000)

all_cpgs_to_extract <- unique(c(
  categories$top1pcEndoOnly,
  categories$top1pcMesoOnly,
  categories$top1pcEctoOnly,
  top1pc3layers_sample,
  constitutive_sample,
  neutral_sample
))

writeLines(all_cpgs_to_extract,
           here("B_MultiTissues/dataOut/CpG2extractS05.txt"))
message(sprintf("Written: CpG2extractS05.txt (%d CpGs total)",
                length(all_cpgs_to_extract)))
# Written: CpG2extractS05.txt (3378 CpGs total)

## In pchuckle (after git pull):
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/CpG2extractS05.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/CpG2extractS05.tsv \
# --minCov    10

## !! Transfer to local gitignore/CpG2extractS05.tsv

meth_CpG2extractS05 <- fread(here("gitignore/CpG2extractS05.tsv"))

# There are DUPLICATES for some samples, let's remove them to avoid confusion at a later stage
# (we keep one replicate per individual)
meth_CpG2extractS05_dedup <- unique(
  meth_CpG2extractS05,
  by = c("cpg_site", "patient_id", "source_tissue_celltype")
)

##############################################
## Add categories to meth_CpG2extractS05_dedup
cpg_category <- rbindlist(lapply(names(categories), function(cat)
  data.table(cpg_site = categories[[cat]], category = cat)))

## Add neutral ones
cpg_category <- rbind(cpg_category,
                      data.table(cpg_site = neutral_sample, category = "background"))

meth_CpG2extractS05_dedup <- merge(meth_CpG2extractS05_dedup, cpg_category, by = "cpg_site", all.x = TRUE)

############################
## Add germ_layer to allmeth
loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])

meth_CpG2extractS05_dedup <- merge(meth_CpG2extractS05_dedup, tissue_to_layer,
                                   by = "source_tissue_celltype", all.x = TRUE)

###########################################
## How many patients have multiple tissues?

multi_patients <- meth_CpG2extractS05_dedup[, .(n = uniqueN(source_tissue_celltype)),
                                            by = patient_id][n > 1, patient_id]
message(sprintf("%d patients with >1 tissue", length(multi_patients)))
# 32 patients with >1 tissue

## keep only these patients
meth_multi <- meth_CpG2extractS05_dedup[patient_id %in% multi_patients]

# ── Descriptive table of multi-tissue patients ────────────────────────────────
patient_table <- meth_multi[, .(
  n_tissues   = uniqueN(source_tissue_celltype),
  germ_layers = paste(sort(unique(germ_layer)), collapse = "+"),
  n_blood     = uniqueN(source_tissue_celltype[grepl("Blood", source_tissue_celltype, ignore.case=TRUE)]),
  n_nonblood  = uniqueN(source_tissue_celltype[!grepl("Blood", source_tissue_celltype, ignore.case=TRUE)])
), by = patient_id]

summary_table <- patient_table[, .(
  n_patients     = .N,
  median_tissues = as.numeric(median(n_tissues)),
  range_tissues  = sprintf("%d-%d", min(n_tissues), max(n_tissues)),
  n_blood_only   = sum(n_nonblood == 0),
  pct_blood_only = round(100 * sum(n_nonblood == 0) / .N)
), by = germ_layers][order(germ_layers)]
print(summary_table)
#      germ_layers n_patients median_tissues range_tissues n_blood_only pct_blood_only
# 1:        Ecto          6            2.0           2-2            0              0
# 2:        Endo         13            2.0           2-3            0              0
# 3:   Endo+Meso          4            2.5           2-3            0              0
# 4:        Meso          9            5.0           2-7            6             67

### patients with >1 tissue WITHIN the same germ layer
same_layer_multi <- meth_multi[, .(n_tissues_in_layer = uniqueN(source_tissue_celltype)),
                               by = .(patient_id, germ_layer)][n_tissues_in_layer > 1]
message(sprintf("Patients with >1 tissue within the SAME germ layer: %d",
                uniqueN(same_layer_multi$patient_id)))
# Patients with >1 tissue within the SAME germ layer: 30

print(same_layer_multi[, .(n_patients = uniqueN(patient_id)), by = germ_layer])
# germ_layer n_patients
# 1:       Meso          9
# 2:       Ecto          6
# 3:       Endo         15

################################################################################
## Calculate inter-cell types correlation by CpG category                     ##
################################################################################

setDT(meth_multi)

################################################################################
## 1. 4 patients with BOTH Endo and Meso tissues -------------------------------
#    (count distinct layers among {Endo, Meso} per patient; keep those with 2)
################################################################################

patients_EM <- meth_multi[germ_layer %in% c("Endo", "Meso"),
                          .(n_layers = uniqueN(germ_layer)),
                          by = patient_id][n_layers == 2, patient_id]

message(sprintf("Patients with both Endo and Meso: %d", length(patients_EM)))
# Patients with both Endo and Meso: 4

## 2. one value per (CpG, patient, layer) = RANDOMLY chosen tissue -------------
set.seed(1234)
em <- meth_multi[
  patient_id %in% patients_EM & germ_layer %in% c("Endo", "Meso"),
  .(methylation = methylation[sample(.N, 1)]),        # pick one of the patient's tissues at random
  by = .(cpg_site, patient_id, germ_layer)
]

## reshape so each patient has an Endo column and a Meso column, per CpG
em_wide <- dcast(em, cpg_site + patient_id ~ germ_layer, value.var = "methylation")
# columns now: cpg_site, patient_id, Endo, Meso

## 3. per CpG: Pearson r between Endo and Meso across the patients -------------
em_r <- em_wide[, {
  idx <- !is.na(Endo) & !is.na(Meso)          # patients with both values present
  n   <- sum(idx)
  if (n >= 3 && sd(Endo[idx]) > 0 && sd(Meso[idx]) > 0) {
    r <- cor(Endo[idx], Meso[idx])            # Pearson across patients
  } else if (n >= 3) {
    r <- 0                                     # one layer is flat -> no covariation
  } else {
    r <- NA_real_                              # too few patients to correlate
  }
  .(r = r, abs_r = abs(r), n_patients = n)
}, by = cpg_site]

## result: one row per CpG, with signed r, |r|, and how many patients contributed

############################
## Add categories to allmeth
em_r <- merge(em_r, cpg_category, by = "cpg_site", all.x = TRUE)

category_colours <- c(
  "top1pcEctoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[1],
  "top1pcEndoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[2],
  "top1pcMesoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[3],
  "top1pc3layers" = "orange",
  "constitutive" = "grey80",
  "background" = "black"
)

## And plot inter-germ layer (endo-meso, Nind=4) correlation by category
p_mesoendocor_violin <- ggplot(em_r, aes(x=category, y=abs_r, group = category, fill = category))+
  geom_violin(width=2) +
  geom_boxplot(width=0.1, color="black", fill = "white") +
  scale_fill_manual(values = category_colours) +
  theme_minimal(base_size = 14) +
  geom_text(data = em_r[, .(nCpG = .N), by = category],
            aes(x = category, y = 1, label = format(nCpG, big.mark = ",")),
            vjust = -0.4, size = 3, inherit.aes = FALSE) +
  labs(y = "Mesoderm-endoderm correlation\n(Pearson's r)")+
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) 

p_mesoendocor_violin

################################################################################
## SAME-LAYER (Endo-Endo): patients with >= 2 endodermal cell types
## Correlate two randomly-chosen endo tissues per patient, across patients
################################################################################

## 1. patients with >= 2 distinct Endo tissues --------------------------------

patients_EE <- meth_multi[germ_layer == "Endo",
                          .(n_tissues = uniqueN(source_tissue_celltype)),
                          by = patient_id][n_tissues >= 2, patient_id]
message(sprintf("Patients with >=2 Endo tissues: %d", length(patients_EE)))
# Patients with >=2 Endo tissues: 15

## 2. per patient, pick 2 random Endo tissues, label them tissueA / tissueB ----
set.seed(1234)
chosen_EE <- unique(meth_multi[patient_id %in% patients_EE & germ_layer == "Endo",
                               .(patient_id, source_tissue_celltype)])
chosen_EE <- chosen_EE[, {
  tis <- unique(source_tissue_celltype)
  if (length(tis) < 2) NULL                       # skip patients with <2 distinct tissues
  else {
    picks <- sample(tis, 2)                        # 2 DISTINCT tissues
    .(source_tissue_celltype = picks, slot = c("tissueA", "tissueB"))
  }
}, by = patient_id]

# keep only those 2 tissues' methylation, tagged A/B
ee <- merge(
  meth_multi[patient_id %in% patients_EE & germ_layer == "Endo",
             .(cpg_site, patient_id, source_tissue_celltype, methylation)],
  chosen_EE, by = c("patient_id", "source_tissue_celltype"))

## reshape: each patient gets a tissueA column and a tissueB column, per CpG
ee_wide <- dcast(ee, cpg_site + patient_id ~ slot, value.var = "methylation")
# columns now: cpg_site, patient_id, tissueA, tissueB

## 3. per CpG: Pearson r between tissueA and tissueB across patients -----------
ee_r <- ee_wide[, {
  idx <- !is.na(tissueA) & !is.na(tissueB)
  n   <- sum(idx)
  if (n >= 3 && sd(tissueA[idx]) > 0 && sd(tissueB[idx]) > 0) {
    r <- cor(tissueA[idx], tissueB[idx])
  } else if (n >= 3) {
    r <- 0
  } else {
    r <- NA_real_
  }
  .(r = r, abs_r = abs(r), n_patients = n)
}, by = cpg_site]

## categories + plot (identical to the cross-layer version)
ee_r <- merge(ee_r, cpg_category, by = "cpg_site", all.x = TRUE)

p_endoendocor_violin <- ggplot(ee_r, aes(x = category, y = abs_r, group = category, fill = category)) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.1, color = "black", fill = "white") +
  scale_fill_manual(values = category_colours) +
  theme_minimal(base_size = 14) +
  geom_text(data = ee_r[, .(nCpG = .N), by = category],
            aes(x = category, y = 1, label = format(nCpG, big.mark = ",")),
            vjust = -0.4, size = 3, inherit.aes = FALSE) +
  labs(y = "Endoderm-endoderm correlation\n(Pearson's r)") +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p_endoendocor_violin

################################################################################
## SAME-LAYER (Meso-Meso): patients with >= 2 mesodermal cell types
## Correlate two randomly-chosen endo tissues per patient, across patients
################################################################################

## 1. patients with >= 2 distinct Endo tissues --------------------------------

patients_MM <- meth_multi[germ_layer == "Meso",
                          .(n_tissues = uniqueN(source_tissue_celltype)),
                          by = patient_id][n_tissues >= 2, patient_id]
message(sprintf("Patients with >=2 Meso tissues: %d", length(patients_MM)))
# Patients with >=2 Endo tissues: 9

## 2. per patient, pick 2 random Endo tissues, label them tissueA / tissueB ----
set.seed(1234)
chosen_MM <- unique(meth_multi[patient_id %in% patients_MM & germ_layer == "Meso",
                               .(patient_id, source_tissue_celltype)])
chosen_MM <- chosen_MM[, {
  tis <- unique(source_tissue_celltype)
  if (length(tis) < 2) NULL                       # skip patients with <2 distinct tissues
  else {
    picks <- sample(tis, 2)                        # 2 DISTINCT tissues
    .(source_tissue_celltype = picks, slot = c("tissueA", "tissueB"))
  }
}, by = patient_id]

# keep only those 2 tissues' methylation, tagged A/B
mm <- merge(
  meth_multi[patient_id %in% patients_MM & germ_layer == "Meso",
             .(cpg_site, patient_id, source_tissue_celltype, methylation)],
  chosen_MM, by = c("patient_id", "source_tissue_celltype"))

## reshape: each patient gets a tissueA column and a tissueB column, per CpG
mm_wide <- dcast(mm, cpg_site + patient_id ~ slot, value.var = "methylation")
# columns now: cpg_site, patient_id, tissueA, tissueB

## 3. per CpG: Pearson r between tissueA and tissueB across patients -----------
mm_r <- mm_wide[, {
  idx <- !is.na(tissueA) & !is.na(tissueB)
  n   <- sum(idx)
  if (n >= 3 && sd(tissueA[idx]) > 0 && sd(tissueB[idx]) > 0) {
    r <- cor(tissueA[idx], tissueB[idx])
  } else if (n >= 3) {
    r <- 0
  } else {
    r <- NA_real_
  }
  .(r = r, abs_r = abs(r), n_patients = n)
}, by = cpg_site]

## categories + plot (identical to the cross-layer version)
mm_r <- merge(mm_r, cpg_category, by = "cpg_site", all.x = TRUE)

p_mesomesocor_violin <- ggplot(mm_r, aes(x = category, y = abs_r, group = category, fill = category)) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.1, color = "black", fill = "white") +
  scale_fill_manual(values = category_colours) +
  theme_minimal(base_size = 14) +
  geom_text(data = mm_r[, .(nCpG = .N), by = category],
            aes(x = category, y = 1, label = format(nCpG, big.mark = ",")),
            vjust = -0.4, size = 3, inherit.aes = FALSE) +
  labs(y = "Mesoderm-mesoderm correlation\n(Pearson's r)") +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p_mesomesocor_violin

################################################################################
## Binned by inter-individual variation range                                 ##
################################################################################

# We calculated interindividual variation using the same metric as van Baak et al. (27):
# for each CpG, we took the mean methylation value across the two cell type
# for every individual and defined interindividual variation of the CpG as the 
# range of these means.

makeIVbyICplot <- function(meth2cells = em, rdat = em_r) {
  interindividual_var <- meth2cells %>%
    dplyr::group_by(patient_id, cpg_site) %>%
    dplyr::summarise(mean_beta = mean(methylation, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(cpg_site) %>%
    dplyr::summarise(interindividual_var = max(mean_beta, na.rm = TRUE) - min(mean_beta, na.rm = TRUE))
  
  ## Add category and inter cell type correlation:
  interindividual_var <- merge(interindividual_var, rdat, by = "cpg_site", all.x = TRUE)
  
  pintvardens <- ggplot(interindividual_var, aes(x = interindividual_var, fill = category)) +
    geom_density(alpha = .7)+
    scale_fill_manual(values = category_colours) +
    theme_minimal(base_size = 14) +
    labs(x = "Interindividual variation")
  
  ###################################################################
  # Bin interindividual_var into 0.2 intervals and bootstrap for CI #
  ###################################################################
  # Function to compute bootstrap CI for median
  boot_median_ci <- function(x, nboot = 1000, conf = 0.95) {
    x <- x[!is.na(x)]
    if (length(x) < 5) return(c(median = NA, low = NA, high = NA))  # skip small bins
    
    bootfun <- function(data, idx) median(data[idx], na.rm = TRUE)
    b <- boot(x, statistic = bootfun, R = nboot)
    ci <- boot.ci(b, type = "perc", conf = conf)
    
    if (!is.null(ci) && "percent" %in% names(ci)) {
      c(median = median(x, na.rm = TRUE), low = ci$percent[4], high = ci$percent[5])
    } else {
      c(median = median(x, na.rm = TRUE), low = NA, high = NA)
    }
  }
  
  # Apply bootstrap per (group, bin) + now also carry n (CpGs per bin)
  binned_summary_boot <- interindividual_var %>%
    mutate(bin = cut(
      interindividual_var,
      breaks = seq(0, max(interindividual_var, na.rm = TRUE) + 0.2, by = 0.2),
      include.lowest = TRUE
    )) %>%
    group_by(category, bin) %>%
    summarise(
      boot_res = list(boot_median_ci(abs_r)),
      n = sum(!is.na(abs_r)),          # <- CpGs contributing to this bin
      .groups = "drop"
    ) %>%
    mutate(
      median_r = sapply(boot_res, `[[`, "median"),
      low = sapply(boot_res, `[[`, "low"),
      high = sapply(boot_res, `[[`, "high")
    ) %>%
    dplyr::select(-boot_res)
  
  # label height: just above each point's upper CI, staggered a touch per category
  binned_summary_boot <- binned_summary_boot %>%
    mutate(lab_y = pmax(high, median_r, na.rm = TRUE) + 0.03 +
             0.025 * (as.integer(factor(category)) - 1))
  
  dodge <- position_dodge(width = 0.5)   # ONE dodge, shared by points/errorbars/text
  
  # Plot
  pbinned <- ggplot(binned_summary_boot,
                    aes(x = bin, y = median_r, color = category, fill = category)) +
    geom_point(position = dodge, size = 3) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, position = dodge) +
    geom_text(aes(y = lab_y, label = n), position = dodge,
              size = 2.8, fontface = "bold", show.legend = FALSE, vjust = 0) +
    scale_color_manual(values = category_colours) +
    scale_fill_manual(values = category_colours) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +   # headroom for labels
    theme_minimal(base_size = 14) +
    labs(
      x = "Interindividual variation",
      y = "Inter-cell type correlation \n(median ± bootstrap CI)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  grobs  <- ggplotGrob(pintvardens)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  
  return(list(plot = plot_grid(pintvardens + theme(legend.position = "none"),
            pbinned + theme(legend.position = "none"), ncol = 2),
            legend = legend))
}

## 1. meso-endo
Pa <- makeIVbyICplot(em, em_r)

## 2. meso-meso
Pb <- makeIVbyICplot(mm, mm_r)

## 3. endo-endo
Pc <- makeIVbyICplot(ee, ee_r)

########## Third row: layer-specific selection (high own r + low others)
## Selection logic, per germ layer "own":
##   1. HIGH concordance in the OWN layer      : abs_r_own  >= background median + 1 SD
##   2. LOW  concordance in the OTHER two      : abs_r_otherX < background median  (both)
## Background = constitutive CpGs, one threshold per r-type (em / ee / mm).
##
## "own" layer uses the SAME-layer r (ee for Endo, mm for Meso); the "other"
## constraints use the same-layer r of the other layer AND the cross-layer r.

## ── 1. background (constitutive) thresholds, one per r-type ───────────────────
bg_thr <- function(rtab) {
  b <- rtab[category == "constitutive" & !is.na(abs_r), abs_r]
  list(high = median(b) + sd(b),   # "high own" cutoff
       low  = median(b))           # "low other" cutoff
}
thr_ee <- bg_thr(ee_r)   # Endo same-layer
thr_mm <- bg_thr(mm_r)   # Meso same-layer
thr_em <- bg_thr(em_r)   # Endo-Meso cross-layer

message(sprintf("thresholds  EE high=%.2f low=%.2f | MM high=%.2f low=%.2f | EM hig=%.2f low=%.2f",
                thr_ee$high, thr_ee$low, thr_mm$high, thr_mm$low, thr_em$high, thr_em$low))
# thresholds  EE high=0.40 low=0.20 | MM high=0.52 low=0.30 | EM hig=0.81 low=0.52

## ── 2. one wide table: each CpG's r in all three contexts ─────────────────────
## rename abs_r per source so they don't clash, then join on cpg_site.
r_all <- Reduce(function(a, b) merge(a, b, by = c("cpg_site","category"), all = TRUE), list(
  ee_r[, .(cpg_site, category, r_ee = abs_r)],
  mm_r[, .(cpg_site, category, r_mm = abs_r)],
  em_r[, .(cpg_site, category, r_em = abs_r)]
))

## ── 3. cumulative pass over the criteria, per "own" layer ─────────────────────
## Endo-specific : high Endo (r_ee) + low Meso (r_mm) + low cross (r_em)
## Meso-specific : high Meso (r_mm) + low Endo (r_ee) + low cross (r_em)
steps_for <- function(dat, own) {
  if (own == "Endo") {
    dat[, `:=`(r_own = r_ee, r_oth = r_mm, r_x = r_em,
               hi = thr_ee$high, lo_oth = thr_mm$low, lo_x = thr_em$low)]
  } else {  # Meso
    dat[, `:=`(r_own = r_mm, r_oth = r_ee, r_x = r_em,
               hi = thr_mm$high, lo_oth = thr_ee$low, lo_x = thr_em$low)]
  }
  dat[, .(
    own_layer = own,
    step = factor(c("1. high own r","2. + low other r","3. + low cross r"),
                  levels = c("1. high own r","2. + low other r","3. + low cross r")),
    n = c(
      sum(!is.na(r_own) & r_own >= hi),                                              # step 1
      sum(!is.na(r_own) & r_own >= hi & !is.na(r_oth) & r_oth < lo_oth),             # +2
      sum(!is.na(r_own) & r_own >= hi & !is.na(r_oth) & r_oth < lo_oth &
            !is.na(r_x)  & r_x  < lo_x)                                               # +3
    ),
    n_tested = sum(!is.na(r_own))
  ), by = category]
}

tri_tabs2 <- rbindlist(list(
  steps_for(copy(r_all), "Endo"),
  steps_for(copy(r_all), "Meso")
))
tri_tabs2[, pct := 100 * n / n_tested]
tri_tabs2[, category_f := factor(category,
                                 levels = c("top1pcMesoOnly","top1pcEndoOnly","top1pc3layers","constitutive"))]

## ── 4. plot (same style as the downstream tri_tabs figure) ────────────────────
p_row3 <- ggplot(tri_tabs2[!is.na(category_f)], aes(category_f, pct, fill = step)) +
  geom_col(position = position_dodge(0.7), width = 0.65,
           colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.7),
            vjust = -0.3, size = 2.7) +
  facet_wrap(~ own_layer, nrow = 1,
             labeller = labeller(own_layer = c(
               Endo = "Own = Endo (other = Meso)",
               Meso = "Own = Meso (other = Endo)"))) +
  scale_fill_brewer(palette = "Reds", name = "Cumulative criterion") +
  scale_y_continuous("% of CpGs tested in own layer", limits = c(0, NA)) +
  scale_x_discrete("CpG category") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text  = element_text(face = "bold"),
        legend.position = "right") +
  ggtitle("Layer-specific selection: high own-layer r + low in the two others",
          subtitle = "high own |r| >= bg median+1SD ->  + low other same-layer |r| < bg median  ->  + low cross-layer |r| < bg median")
p_row3

######################################
## Extract candidates these targets ##
######################################

## CpGs passing all 3 criteria, per "own" layer
extract_pass3 <- function(dat, own) {
  if (own == "Endo") {
    dat <- dat[, .(cpg_site, category,
                   r_own = r_ee, r_oth = r_mm, r_x = r_em)]
    hi <- thr_ee$high; lo_oth <- thr_mm$low; lo_x <- thr_em$low
  } else {  # Meso
    dat <- dat[, .(cpg_site, category,
                   r_own = r_mm, r_oth = r_ee, r_x = r_em)]
    hi <- thr_mm$high; lo_oth <- thr_ee$low; lo_x <- thr_em$low
  }
  dat[!is.na(r_own) & r_own >= hi &
        !is.na(r_oth) & r_oth <  lo_oth &
        !is.na(r_x)   & r_x   <  lo_x
  ][, own_layer := own][]
}

endo_pass3 <- extract_pass3(copy(r_all), "Endo")
meso_pass3 <- extract_pass3(copy(r_all), "Meso")

# how many, and which categories they came from
message(sprintf("Endo-specific (3 criteria): %d | Meso-specific: %d",
                nrow(endo_pass3), nrow(meso_pass3)))
print(endo_pass3[, .N, by = category])
print(meso_pass3[, .N, by = category])

# the CpG lists of candidates
endo_candidates <- endo_pass3$cpg_site[endo_pass3$category %in% "top1pcEndoOnly"]
meso_candidates <- meso_pass3$cpg_site[meso_pass3$category %in% "top1pcMesoOnly"]

#################################
## Annotation of these targets ##
#################################

# ── Shared annotation objects (built once) ────────────────────────────────────
txdb     <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
genes_gr$symbol <- mapIds(org.Hs.eg.db, names(genes_gr), "SYMBOL", "ENTREZID")

cpg_to_gr <- function(cpg) {                     # "chr7_107543290" -> 1bp GRanges
  GRanges(sub("_.*", "", cpg),
          IRanges(as.integer(sub(".*_", "", cpg)), width = 1),
          cpg_site = cpg)
}

# universe of tested CpGs = the background only
uni <- cpg_category[category == "background"]          # cols: cpg_site, category
uni[, `:=`(chr = sub("_.*", "", cpg_site),
           pos = as.integer(sub(".*_", "", cpg_site)))]
uni_gr <- GRanges(uni$chr, IRanges(uni$pos, width = 1))

# ── Build shared feature sets ONCE (next to genes_gr) ─────────────────────────
tss_gr    <- promoters(genes_gr, upstream = 0, downstream = 1)          # 1 bp at each TSS
prom_gr   <- promoters(txdb, upstream = 2000, downstream = 200)         # promoter window
exon_gr   <- reduce(exons(txdb))
intron_gr <- reduce(unlist(intronsByTranscript(txdb)))

# te_regions is built in the TE-enrichment block; rebuild if running cached
if (!exists("te_regions")) {
  library(AnnotationHub); ah <- AnnotationHub(); rmskhg38 <- ah[["AH111333"]]
  te_regions <- rmskhg38[mcols(rmskhg38)$repClass %in%
                           c("LINE","SINE","LTR","DNA","RC","Retroposon")]
}

# ── Main function ─────────────────────────────────────────────────────────────
annotate_layer_hits <- function(hits, label, genes_gr, uni, uni_gr,
                                tss_gr, prom_gr, exon_gr, intron_gr, te_regions,
                                gap = 50, min_hits = 1,
                                drop_regex = "LOC|LINC", flank = 5000,
                                out_dir = here("B_MultiTissues/dataOut")) {
  if (is.null(hits) || !nrow(hits)) { message("No hits for ", label); return(NULL) }
  hit_col <- paste0("is_", label, "_hit")
  hit_set <- hits$cpg_site
  
  ## (1a) nearest gene ---------------------------------------------------------
  hits_gr <- cpg_to_gr(hit_set)
  nr  <- distanceToNearest(hits_gr, genes_gr, ignore.strand = TRUE)
  ann <- as.data.table(hits)
  ann[, `:=`(gene = NA_character_, dist_to_gene = NA_integer_)]
  ann[queryHits(nr),
      `:=`(gene = genes_gr$symbol[subjectHits(nr)],
           dist_to_gene = mcols(nr)$distance)]
  ann[, `:=`(pos = as.integer(sub(".*_", "", cpg_site)),
             chr = sub("_.*", "", cpg_site))]
  setorder(ann, chr, pos)
  
  ## (1b) distance to nearest TSS (signed) -------------------------------------
  ann_gr <- cpg_to_gr(ann$cpg_site)                      # rebuilt after setorder
  nrt <- distanceToNearest(ann_gr, tss_gr, ignore.strand = TRUE)
  ann[, `:=`(dist_to_tss = NA_integer_, dist_to_tss_signed = NA_integer_)]
  ann[queryHits(nrt), dist_to_tss := mcols(nrt)$distance]
  nti        <- subjectHits(nrt)
  tss_pos    <- start(tss_gr)[nti]
  tss_strand <- as.character(strand(tss_gr))[nti]
  ann[queryHits(nrt), dist_to_tss_signed :=
        fifelse(tss_strand == "-",
                tss_pos - ann$pos[queryHits(nrt)],
                ann$pos[queryHits(nrt)] - tss_pos)]      # neg = upstream, pos = downstream
  
  ## (1c) genomic feature: promoter > exon > intron > intergenic ---------------
  ann[, feature := "intergenic"]
  ann[overlapsAny(ann_gr, intron_gr, ignore.strand = TRUE), feature := "intron"]
  ann[overlapsAny(ann_gr, exon_gr,   ignore.strand = TRUE), feature := "exon"]
  ann[overlapsAny(ann_gr, prom_gr,   ignore.strand = TRUE), feature := "promoter"]
  
  ## (1d) in a TE? -------------------------------------------------------------
  te_hit <- findOverlaps(ann_gr, te_regions, ignore.strand = TRUE)
  te_hit <- te_hit[!duplicated(queryHits(te_hit))]       # first TE per CpG (deterministic)
  ann[, `:=`(in_TE = FALSE, TE_class = NA_character_, TE_family = NA_character_)]
  ann[queryHits(te_hit),
      `:=`(in_TE     = TRUE,
           TE_class  = mcols(te_regions)$repClass[subjectHits(te_hit)],
           TE_family = mcols(te_regions)$repFamily[subjectHits(te_hit)])]
  
  ## (2) cluster hits within `gap` bp -----------------------------------------
  hs <- cpg_to_gr(ann$cpg_site)
  clust <- reduce(resize(hs, 1), min.gapwidth = gap + 1)
  ov    <- findOverlaps(hs, clust)
  ann[, cluster_id := subjectHits(ov)[match(seq_len(.N), queryHits(ov))]]
  
  clust_dt <- ann[, .(n_hits = .N, span = max(pos) - min(pos)),
                  by = .(gene, cluster_id)][n_hits >= min_hits]
  if (nzchar(drop_regex))
    clust_dt <- clust_dt[!grepl(drop_regex, gene) & !is.na(gene)]
  clust_dt[, density := n_hits / (span + 1)]
  setorder(clust_dt, -n_hits, span)
  
  ## (3) all CpGs in each top gene + flank, flagged ---------------------------
  top_genes <- unique(clust_dt$gene)
  extract_gene_cpgs <- function(sym) {
    g <- genes_gr[which(genes_gr$symbol == sym)]
    if (!length(g)) return(NULL)
    g   <- g[1]
    win <- GRanges(as.character(seqnames(g)),
                   IRanges(start(g) - flank, end(g) + flank))
    out <- uni[which(overlapsAny(uni_gr, win, ignore.strand = TRUE))]
    out[, `:=`(gene = sym, gene_start = start(g), gene_end = end(g),
               in_gene_body = pos >= start(g) & pos <= end(g))]
    out[, (hit_col) := cpg_site %in% hit_set]
    out[order(pos)]
  }
  gene_cpgs <- rbindlist(lapply(top_genes, extract_gene_cpgs), fill = TRUE)
  
  fwrite(gene_cpgs, file.path(out_dir, sprintf("%s_hits_topGenes_allCpGs.csv", label)))
  fwrite(ann,       file.path(out_dir, sprintf("%s_hits_annotated.csv", label)))
  list(ann = ann, clusters = clust_dt, gene_cpgs = gene_cpgs)
}

# ── Calls (pass the shared feature sets) ──────────────────────────────────────
meso_hits <- meso_pass3[cpg_site %in% meso_candidates]   # get the full rows back
setnames(meso_hits, c("r_oth","r_x"), c("r_other","r_cross"), skip_absent = TRUE)

endo_hits <- endo_pass3[cpg_site %in% endo_candidates]   # get the full rows back
setnames(endo_hits, c("r_oth","r_x"), c("r_other","r_cross"), skip_absent = TRUE)

meso_res <- annotate_layer_hits(meso_hits, "meso", genes_gr, uni, uni_gr,
                                tss_gr, prom_gr, exon_gr, intron_gr, te_regions, gap = 100)
endo_res <- annotate_layer_hits(endo_hits, "endo", genes_gr, uni, uni_gr,
                                tss_gr, prom_gr, exon_gr, intron_gr, te_regions, gap = 100)

message("\n=== MESO layer-specific hits (all 3 criteria) ===")
if (!is.null(meso_res)) print(meso_res$ann[, .(cpg_site, gene, dist_to_gene,
                                               dist_to_tss_signed, feature, in_TE, TE_class, TE_family, r_own, r_other, r_cross)])

message("\n=== ENDO layer-specific hits (all 3 criteria) ===")
if (!is.null(endo_res)) print(endo_res$ann[, .(cpg_site, gene, dist_to_gene,
                                               dist_to_tss_signed, feature, in_TE, TE_class, TE_family, r_own, r_other, r_cross)])

saveRDS(list(meso = meso_res, endo = endo_res),
        here(paste0("B_MultiTissues/dataOut/S05_annotatedHits_", variant, ".rds")))

# ── Combined pretty table (guard against either set being empty) ──────────────
hits_tab <- rbindlist(list(
  if (!is.null(meso_res)) meso_res$ann[, layer := "Mesoderm"],
  if (!is.null(endo_res)) endo_res$ann[, layer := "Endoderm"]
), fill = TRUE)

if (nrow(hits_tab)) {
  hits_tab <- hits_tab[, .(
    layer, cpg_site, gene, feature,
    dist_to_tss_signed,
    TE = fifelse(in_TE, paste0(TE_class, " / ", TE_family), "\u2014"),
    r_own, r_other, r_cross
  )][order(layer, -r_own)]
  
  gt_hits <- hits_tab %>%
    gt(groupname_col = "layer") %>%
    tab_header(
      title    = md("**Layer-specific hypervariable CpGs passing all three criteria**"),
      subtitle = md("High within-layer, low other-layer and low cross-layer within-individual correlation")
    ) %>%
    fmt_number(columns = c(r_own, r_other, r_cross), decimals = 2) %>%
    fmt_number(columns = dist_to_tss_signed, decimals = 0, use_seps = TRUE) %>%
    cols_label(
      cpg_site           = "CpG",
      gene               = "Nearest gene",
      feature            = "Genomic feature",
      dist_to_tss_signed = "Dist. to TSS (bp)",
      TE                 = "TE (class / family)",
      r_own              = md("*r* own"),
      r_other            = md("*r* other"),
      r_cross            = md("*r* cross")
    ) %>%
    tab_spanner(label = "Within-individual correlation",
                columns = c(r_own, r_other, r_cross)) %>%
    data_color(columns = r_own,
               fn = scales::col_numeric(c("white", "#2166AC"), domain = c(0.7, 1))) %>%
    tab_style(style = cell_text(style = "italic"),
              locations = cells_body(columns = gene)) %>%
    tab_style(style = cell_fill(color = "#F0F0F0"),
              locations = cells_row_groups()) %>%
    tab_options(table.font.size = 12, data_row.padding = px(3),
                row_group.font.weight = "bold") %>%
    tab_source_note(md(
      "*r* own = within-individual correlation across tissues of the CpG's own germ layer; *r* other / *r* cross = the same in the other layer and across layers. TE from RepeatMasker (hg38). Distance to TSS signed: negative = upstream."))
  
  gtsave(gt_hits, here("B_MultiTissues/dataOut/figures/script05/layerSpecificHits_table.png"))
  gt_hits
} else {
  message("No layer-specific hits passed all three criteria — no table produced.")
}

################################################################################
## Save full layers plot:                                                     ##
################################################################################

# read back the table as a grob for cowplot
tbl_img  <- png::readPNG(here("B_MultiTissues/dataOut/figures/script05/layerSpecificHits_table.png"))
tbl_grob <- grid::rasterGrob(tbl_img, interpolate = TRUE)

# overlay Pa$legend onto Pa$plot
Pa_with_legend <- ggdraw(Pa$plot) +
  draw_grob(Pa$legend,
            x = 0.98, y = 0.98,        # anchor near top-right (0-1 npc coords)
            width = 0.2, height = 0.2, # legend box size
            hjust = 3.5, vjust = 2)      

layersPlot <- plot_grid(
  plot_grid(p_mesoendocor_violin, p_mesomesocor_violin, p_endoendocor_violin,
            Pa_with_legend, Pb$plot, Pc$plot, nrow = 3, 
            labels = c("A", "B", "C", "D", "E", "F")),
  plot_grid(p_row3, tbl_grob, nrow = 1),
  nrow = 2, rel_heights = c(2,1, 1), labels = c("", "G", "H"))

ggsave(here("B_MultiTissues/dataOut/figures/script05/layersPlot.png"),
       layersPlot, width = 22, height = 22, dpi = 300, bg = "white")
