#####################################################################
# S05 — Per-germ-layer hypervariability categories
#
# Four mutually-exclusive categories, each a set of chr_pos:
#   allLayers : top 1% of logBF_per_ds in endo AND in meso AND in ecto (stricter than in S04)
#   endo_only : top 1% in endo  AND bottom 50% in meso AND bottom 50% in ecto
#   meso_only : top 1% in meso  AND bottom 50% in endo AND bottom 50% in ecto
#   ecto_only : top 1% in ecto  AND bottom 50% in endo AND bottom 50% in meso
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
  allLayers = dt[hv_endo & hv_meso & hv_ecto, chr_pos],   # top 1% in ALL three layers
  endo_only = dt[hv_endo & lo_meso & lo_ecto, chr_pos],   # top 1% in 1 layer, bottom 50% in the others
  meso_only = dt[hv_meso & lo_endo & lo_ecto, chr_pos],
  ecto_only = dt[hv_ecto & lo_endo & lo_meso, chr_pos],
  constitutive = dt[lo_ecto & lo_endo & lo_meso, chr_pos],
  allLayers6gp = dt[hv_endo6gp & hv_meso6gp & hv_ecto6gp, chr_pos],
  endo6gp_only = dt[hv_endo6gp & lo_meso6gp & lo_ecto6gp, chr_pos],
  meso6gp_only = dt[hv_meso6gp & lo_endo6gp & lo_ecto6gp, chr_pos],
  ecto6gp_only = dt[hv_ecto6gp & lo_endo6gp & lo_meso6gp, chr_pos]
)

message("Category sizes:")
print(sapply(categories, length))
# allLayers    endo_only    meso_only    ecto_only constitutive 
# 60424           18           63          297      6723417     
# allLayers6gp endo6gp_only meso6gp_only  ecto6gp_only 
#      26        47166        51527           48122 

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

calc_stats(categories$endo_only,    categories$endo6gp_only,    "Endo_specific")
calc_stats(categories$meso_only,    categories$meso6gp_only,    "Meso_specific")
calc_stats(categories$ecto_only,    categories$ecto6gp_only,    "Ecto_specific")
calc_stats(categories$allLayers,    categories$allLayers6gp,    "allLayers")
# Endo_specific: full=18, 6gp=47166, shared=0 (0% of full, 0.0% of 6gp)
# Meso_specific: full=63, 6gp=51527, shared=0 (0% of full, 0.1% of 6gp)
# Ecto_specific: full=297, 6gp=48122, shared=73 (25% of full, 0.6% of 6gp)
# allLayers: full=60424, 6gp=26, shared=4 (0% of full, 232400.0% of 6gp)

## Remove useless categories
categories <- categories[
  names(categories) %in% 
    c("allLayers", "endo_only", "meso_only", "ecto_only", "constitutive")]

## =============================================================================
## TE ENRICHMENT — per category, overall + per repClass
## =============================================================================

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

# order TE classes by the allLayers OR so facets share a sensible y-order
te_plot <- copy(te_res)
te_plot[, `:=`(sig = ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."))]

ord <- te_plot[category == "allLayers"][order(odds_ratio), te_class]
te_plot[, te_class := factor(te_class, levels = unique(ord))]
te_plot[, category := factor(category,
                             levels = c("allLayers","endo_only","meso_only","ecto_only"))]

## Add N
te_plot[, category := factor(category,
                             levels = c("allLayers","endo_only","meso_only","ecto_only"),
                             labels = paste0(c("allLayers","endo_only","meso_only","ecto_only"),
                                             "\n(N=", sapply(categories, length)[c("allLayers","endo_only","meso_only","ecto_only")], ")"))]

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

# order TE classes by the allLayers OR so facets share a sensible y-order
te_plot <- copy(te_fam)
te_plot[, `:=`(sig = ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."))]

ord <- te_plot[category == "allLayers"][order(odds_ratio), te_class]
te_plot[, te_class := factor(te_class, levels = unique(ord))]
te_plot[, category := factor(category,
                             levels = c("allLayers","endo_only","meso_only","ecto_only"))]

## Add N
te_plot[, category := factor(category,
                             levels = c("allLayers","endo_only","meso_only","ecto_only"),
                             labels = paste0(c("allLayers","endo_only","meso_only","ecto_only"),
                                             "\n(N=", sapply(categories, length)[c("allLayers","endo_only","meso_only","ecto_only")], ")"))]

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

## =============================================================================
## Compare our results with previous MEs
## =============================================================================

top99q_S04     <- readRDS(here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))  # S04 set
allLayers_S05  <- categories$allLayers                                              # S05 set

n_over <- length(intersect(top99q_S04, allLayers_S05))
data.table(
  n_top99q_S04   = length(top99q_S04),
  n_allLayers_S05 = length(allLayers_S05),
  n_overlap      = n_over,
  pc_of_top99q   = round(100 * n_over / length(top99q_S04), 1),
  pc_of_allLayers = round(100 * n_over / length(allLayers_S05), 1)
)

# n_top99q_S04 n_allLayers_S05 n_overlap pc_of_top99q pc_of_allLayers
#    202467           60424     60424         29.8             100

# use allLayers_S05 instead of the top 1% of script 04

## Focal set = allLayers_S05 (was top99q)
if (!exists("listGR")){
  listGR <- list(allLayers    = makeGRfromMyCpGPos(vec = allLayers_S05, setname = "allLayers"),
                 allButAllLayers = makeGRfromMyCpGPos(
                   setdiff(table3layers_coveredIn3$chr_pos, allLayers_S05), "allButAllLayers"))
}

if (!file.exists(here("B_MultiTissues/dataOut/figures/script05/CompareWithpreviousMEs.png"))){
  
  # Fix chromosome names in geomMeanGR (1 -> chr1)
  if (sum(grepl("chr", seqlevels(table3layers_coveredIn3))) == 0){
    seqlevels(table3layers_coveredIn3) <- paste0("chr", seqlevels(table3layers_coveredIn3))
  }
  
  
  ## Use the GR object with analyses in the 3 layers
  ## Focal "our hvCpG set" for this figure = allLayers from S05 (top 1% in ALL three layers)
  allLayers_S05 <- categories$allLayers
  
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
    "allLayers"          = "orange",
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
  
  plot_allLayersCpGsEnrichME <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
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
  
  print(plot_allLayersCpGsEnrichME)
  
  ################################################################################
  ## Load SIV plots calculated in fetalSIV folder script (in ing-p5)            ##
  ##                                                                            ##
  ## allLayers is an EXACT SUBSET of top99q, and the fetal object stores        ##
  ## per-CpG values (interlayer_r, interindividual_var) — computed per CpG,     ##
  ## independent of group. So we do NOT re-run the fetal script: we relabel the ##
  ## top99q rows that belong to allLayers as a new "allLayers" group.           ##
  ## Panels D (interlayer_corr) and E (CpG_summary) are direct relabels;        ##
  ## panel F (binned_summary_boot) is re-bootstrapped locally for allLayers.    ##
  ################################################################################
  
  plots <- readRDS(here("gitignore/intercorrelationSIVfetal_sepSIV.rds"))
  
  ## allLayers CpGs from S05 (chr_pos) -> EPIC CpG ids via the same dico
  allLayers_S05 <- categories$allLayers
  ## `dico` maps CpG <-> chrpos_hg38 (same object the fetal script used).
  ## If not in memory, load it here (adjust path to wherever the fetal dico lives):
  # dico <- readRDS(here("gitignore/EPIC_dico_hg38.rds"))
  allLayers_CpG <- dico$CpG[dico$chrpos_hg38 %in% allLayers_S05]
  message(length(allLayers_CpG), " allLayers CpGs mapped to EPIC ids")
  
  ## helper: spin off an "allLayers" group by subsetting the top99q rows
  add_allLayers <- function(df) {
    sub <- df[df$group == "top99q" & df$CpG %in% allLayers_CpG, ]
    sub$group <- "allLayers"
    dplyr::bind_rows(df, sub)     # keep top99q AND add allLayers alongside
  }
  
  plots$interlayer_corr <- add_allLayers(plots$interlayer_corr)
  plots$CpG_summary     <- add_allLayers(plots$CpG_summary)
  
  ## sanity: allLayers should be more systemically correlated than top99q (subset effect)
  message("mean interlayer_r  top99q: ",
          round(mean(plots$interlayer_corr$interlayer_r[plots$interlayer_corr$group=="top99q"], na.rm=TRUE),3),
          " | allLayers: ",
          round(mean(plots$interlayer_corr$interlayer_r[plots$interlayer_corr$group=="allLayers"], na.rm=TRUE),3))
  
  ## panel F: re-bootstrap ONLY the new allLayers group, same bins/bootstrap as fetal script
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
  
  binned_allLayers <- plots$CpG_summary %>%
    dplyr::filter(group == "allLayers") %>%
    mutate(bin = cut(interindividual_var,
                     breaks = seq(0, max(interindividual_var, na.rm = TRUE) + 0.1, by = 0.1),
                     include.lowest = TRUE)) %>%
    group_by(group, bin) %>%
    summarise(boot_res = list(boot_median_ci(interlayer_r)), .groups = "drop") %>%
    mutate(median_r = sapply(boot_res, `[[`, "median"),
           low      = sapply(boot_res, `[[`, "low"),
           high     = sapply(boot_res, `[[`, "high")) %>%
    dplyr::select(-boot_res)
  
  plots$binned_summary_boot <- dplyr::bind_rows(plots$binned_summary_boot, binned_allLayers)
  
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
}

##############################################################
## How many of each putative ME is actually in all Layers ? ##
##############################################################

# Universe of covered CpGs, each already labelled top99q vs not
# (both are single-CpG GRanges built from your covered-in-3 sites)
top_gr  <- listGR$allLayers      
rest_gr <- listGR$allButAllLayers

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

## =============================================================================
## GO ENRICHMENT — CpG-density-controlled, per category
## =============================================================================
if (!file.exists(here(paste0("gitignore/S05_GO_", variant, ".rds")))){
  totalSites <- all_cpg
  minimum_CpG_per_cluster <- 2
  universe <- annotateCpGs_txdb(
    clusterCpGs(totalSites, max_gap = 50, min_size = minimum_CpG_per_cluster),
    tss_window = 10000)
  
  go_res <- lapply(names(categories), function(cat) {
    cpgs <- categories[[cat]]
    if (length(cpgs) < 20) { message("SKIP GO ", cat, " - too few CpGs (", length(cpgs), ")"); return(NULL) }
    CpG_GO_pipeline_lengthControlled(
      cpgs, universe = universe,
      control_method = "cpg_count", all_sites = totalSites)
  })
  names(go_res) <- names(categories)
  saveRDS(go_res, here(paste0("gitignore/S05_GO_", variant, ".rds")))
} else go_res <- readRDS(here(paste0("gitignore/S05_GO_", variant, ".rds")))

# significant BP terms per category
lapply(go_res, function(r) if (!is.null(r))
  r$BP@result[r$BP@result$p.adjust < 0.05, "Description"])

# significant MF terms per category
lapply(go_res, function(r) if (!is.null(r))
  r$MF@result[r$MF@result$p.adjust < 0.05, "Description"])

# significant CC terms per category
lapply(go_res, function(r) if (!is.null(r))
  r$CC@result[r$CC@result$p.adjust < 0.05, "Description"])

## Nothing significant

# ══════════════════════════════════════════════════════════════════════════════
# Test of intra-individual (non)correlation for different categories
# ══════════════════════════════════════════════════════════════════════════════
message("Category sizes:")
print(sapply(categories, length))

# ── Write CpG lists for python extraction ─────────────────────────────────────

# allLayers CpGs — subsample 5000 as positive control; keep all layer-specific
set.seed(1234)
allLayers_sample <- sample(categories$allLayers, min(5000, length(categories$allLayers)))
all_cpgs_to_extract <- unique(c(
  categories$endo_only,
  categories$meso_only,
  categories$ecto_only,
  allLayers_sample
))
writeLines(all_cpgs_to_extract,
           here("B_MultiTissues/dataOut/layer_specific_and_ME.txt"))
message(sprintf("Written: layer_specific_and_ME.txt (%d CpGs total)",
                length(all_cpgs_to_extract)))

## In pchuckle (after git pull):
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/layer_specific_and_ME.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_layerspecific_and_ME.tsv \
# --minCov    10

## !! Transfer to local gitignore/methylation_layerspecific_and_ME.tsv

# ── Category lookup (replaces wideFull[, .(cpg_site=name, category)]) ──────────
# one row per CpG -> its category, built from the `categories` list
cpg_category <- rbindlist(lapply(names(categories), function(cat)
  data.table(cpg_site = categories[[cat]], category = cat)))

# ══════════════════════════════════════════════════════════════════════════════
# Step 2. Intra-individual correlation for layer-specific CpGs
# ══════════════════════════════════════════════════════════════════════════════

meth <- fread(here("gitignore/methylation_layerspecific_and_ME.tsv"))

# ── Join germ_layer and category ─────────────────────────────────────────────
loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])

meth <- merge(meth, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)
meth <- merge(meth, cpg_category,    by = "cpg_site",               all.x = TRUE)

message(sprintf("Missing germ_layer: %d | Missing category: %d",
                sum(is.na(meth$germ_layer)), sum(is.na(meth$category))))

# ── Multi-tissue patients ─────────────────────────────────────────────────────
multi_patients <- meth[, .(n = uniqueN(source_tissue_celltype)),
                       by = patient_id][n > 1, patient_id]
message(sprintf("%d patients with >1 tissue", length(multi_patients)))
meth_multi <- meth[patient_id %in% multi_patients]

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

message(sprintf("Patients with cross-layer samples: %d",
                sum(grepl("\\+", patient_table$germ_layers))))

### patients with >1 tissue WITHIN the same germ layer
same_layer_multi <- meth_multi[, .(n_tissues_in_layer = uniqueN(source_tissue_celltype)),
                               by = .(patient_id, germ_layer)][n_tissues_in_layer > 1]
message(sprintf("Patients with >1 tissue within the SAME germ layer: %d",
                uniqueN(same_layer_multi$patient_id)))
print(same_layer_multi[, .(n_patients = uniqueN(patient_id)), by = germ_layer])

# ══════════════════════════════════════════════════════════════════════════════
# Within-individual same-layer vs cross-layer correlation per CpG
# (helper functions unchanged — they operate on meth_multi, not on wide/wideFull)
# ══════════════════════════════════════════════════════════════════════════════

## Fetal-style same-layer r: one point per patient, then correlate across patients
compute_same_layer_r_fetalstyle <- function(meth_sub, layer_name,
                                            min_patients = 3, absolute = TRUE) {
  layer_meth <- meth_sub[germ_layer == layer_name]
  eligible <- layer_meth[, .(n = uniqueN(source_tissue_celltype)),
                         by = patient_id][n > 1, patient_id]
  if (length(eligible) < min_patients) return(NULL)
  
  all_pairs <- rbindlist(lapply(eligible, function(pid) {
    pat  <- layer_meth[patient_id == pid]
    tiss <- sort(unique(pat$source_tissue_celltype))
    if (length(tiss) < 2) return(NULL)
    rbindlist(lapply(combn(tiss, 2, simplify = FALSE), function(pr) {
      d1 <- pat[source_tissue_celltype == pr[1], .(cpg_site, m1 = methylation)]
      d2 <- pat[source_tissue_celltype == pr[2], .(cpg_site, m2 = methylation)]
      merge(d1, d2, by = "cpg_site")[, patient_id := pid]
    }), fill = TRUE)
  }), fill = TRUE)
  if (is.null(all_pairs) || !nrow(all_pairs)) return(NULL)
  
  per_patient <- all_pairs[!is.na(m1) & !is.na(m2),
                           .(m1 = mean(m1), m2 = mean(m2)),
                           by = .(cpg_site, patient_id)]
  
  per_patient[, {
    if (.N >= min_patients) {
      r <- suppressWarnings(cor(m1, m2, method = "pearson"))
      list(r = if (absolute) abs(r) else r, n_obs = .N)
    } else list(r = NA_real_, n_obs = .N)
  }, by = cpg_site]
}

compute_cross_layer_r <- function(meth_sub, layer1, layer2) {
  d1 <- meth_sub[germ_layer == layer1, .(m1 = mean(methylation, na.rm = TRUE)),
                 by = .(cpg_site, patient_id)]
  d2 <- meth_sub[germ_layer == layer2, .(m2 = mean(methylation, na.rm = TRUE)),
                 by = .(cpg_site, patient_id)]
  dm <- merge(d1, d2, by = c("cpg_site", "patient_id"))
  dm[, {
    idx <- !is.na(m1) & !is.na(m2)
    if (sum(idx) >= 3) {
      r <- suppressWarnings(cor(m1[idx], m2[idx], method = "pearson"))
      list(r = abs(r), n_obs = sum(idx))
    } else list(r = NA_real_, n_obs = sum(idx))
  }, by = cpg_site]
}

# ══════════════════════════════════════════════════════════════════════════════
# Controls: random 10k constitutive (no "ambiguous" category in this scheme)
# ══════════════════════════════════════════════════════════════════════════════
set.seed(1234)
constitutive_sample <- sample(categories$constitutive, min(10000, length(categories$constitutive)))
writeLines(constitutive_sample,
           here("B_MultiTissues/dataOut/control_sample_constitutive10k.txt"))
message(sprintf("Written: control_sample_constitutive10k.txt (%d CpGs)",
                length(constitutive_sample)))

## Extract as above (adjust --cpg_list / --output), then:
meth_control <- fread(here("gitignore/methylation_control_constitutive10k.tsv"))
meth_control <- merge(meth_control, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)
set.seed(1234)
constitutive_sample <- sample(categories$constitutive, min(10000, length(categories$constitutive)))
meth_control[, category := fifelse(cpg_site %in% constitutive_sample, "constitutive", NA_character_)]
meth_control_multi <- meth_control[patient_id %in% multi_patients]

message(sprintf("Control CpGs loaded: %d constitutive",
                uniqueN(meth_control_multi[category == "constitutive", cpg_site])))

# ══════════════════════════════════════════════════════════════════════════════
# Run for meso_only, endo_only, allLayers and constitutive control
# (category names updated: Meso_specific->meso_only, Endo_specific->endo_only, ME->allLayers)
# ══════════════════════════════════════════════════════════════════════════════
meth_meso_specific <- meth_multi[category == "meso_only"]
r_same_meso        <- compute_same_layer_r_fetalstyle(meth_meso_specific, "Meso")

meth_endo_specific <- meth_multi[category == "endo_only"]
r_same_endo        <- compute_same_layer_r_fetalstyle(meth_endo_specific, "Endo")

meth_ME <- meth_multi[category == "allLayers"]

# ── Wrong-layer specificity checks ───────────────────────────────────────────
r_same_endo_for_meso <- compute_same_layer_r_fetalstyle(meth_meso_specific, "Endo")
r_same_meso_for_endo <- compute_same_layer_r_fetalstyle(meth_endo_specific, "Meso")
r_same_meso_ME       <- compute_same_layer_r_fetalstyle(meth_ME, "Meso")
r_same_endo_ME       <- compute_same_layer_r_fetalstyle(meth_ME, "Endo")

# ── Controls (constitutive only) in both Meso and Endo contexts ───────────────
r_same_meso_const <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "constitutive"], "Meso")
r_same_endo_const <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "constitutive"], "Endo")

# ── per-layer empirical null from constitutive CpGs ───────────────────────────
null_meso <- r_same_meso_const[!is.na(r), r]
null_endo <- r_same_endo_const[!is.na(r), r]

r_crit <- function(n, alpha = 0.05) {
  tc <- qt(1 - alpha/2, df = n - 2); tc / sqrt(tc^2 + n - 2)
}

thr <- data.table(
  layer  = c("Meso", "Endo"),
  r_high = c(quantile(null_meso, 0.95), quantile(null_endo, 0.95)),
  r_low  = c(median(null_meso),         median(null_endo)),
  mu     = c(mean(null_meso),   mean(null_endo)),
  sd     = c(sd(null_meso),     sd(null_endo)),
  n_pat  = c(9, 15),
  r_crit = c(r_crit(9), r_crit(15)))
print(thr)

# ── Summary function ──────────────────────────────────────────────────────────
make_summary_samelayer <- function(r_same, layer_tested, category_name, thr_tab = thr) {
  if (is.null(r_same)) return(NULL)
  s <- r_same[!is.na(r)]
  r_hi <- thr_tab[layer == layer_tested, r_high]
  data.table(
    category          = category_name,
    same_layer_tested = layer_tested,
    r_high_used       = r_hi,
    n_cpgs_tested     = uniqueN(s$cpg_site),
    pct_high_same     = 100 * mean(s$r >= r_hi),
    mean_r            = mean(s$r),
    median_r          = median(s$r))
}

results_summary <- rbindlist(list(
  make_summary_samelayer(r_same_meso,          "Meso", "meso_only"),
  make_summary_samelayer(r_same_endo_for_meso, "Endo", "meso_only"),
  make_summary_samelayer(r_same_endo,          "Endo", "endo_only"),
  make_summary_samelayer(r_same_meso_for_endo, "Meso", "endo_only"),
  make_summary_samelayer(r_same_meso_ME,       "Meso", "allLayers"),
  make_summary_samelayer(r_same_endo_ME,       "Endo", "allLayers"),
  make_summary_samelayer(r_same_meso_const,    "Meso", "constitutive"),
  make_summary_samelayer(r_same_endo_const,    "Endo", "constitutive")
), fill = TRUE)

print(results_summary)

saveRDS(list(thr = thr, results_summary = results_summary),
        here(paste0("gitignore/S05_intraIndividual_", variant, ".rds")))

# # ── Write CpG lists for python extraction ─────────────────────────────────────
# 
# # ME CpGs — subsample 5000 as positive control
# set.seed(1234)
# me_sample <- sample(wideFull[category == "ME", name], 5000)
# all_cpgs_to_extract <- unique(c(
#   wideFull[category %in% c("Ecto_specific","Endo_specific","Meso_specific"), name],
#   me_sample
# ))
# writeLines(all_cpgs_to_extract,
#            here("B_MultiTissues/dataOut/layer_specific_and_ME.txt"))
# message(sprintf("Written: layer_specific_and_ME.txt (%d CpGs total)",
#                 length(all_cpgs_to_extract)))
# # Written: layer_specific_and_ME.txt (5345 CpGs total)
# 
# ## In pchuckle (after git pull):
# # source /share/apps/source_files/python/python-3.13.0a6.source
# # cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# # python3 S00_extractRawMethylationForTargetCpG.py \
# # --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/layer_specific_and_ME.txt \
# # --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# # --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# # --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# # --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_layerspecific_and_ME.tsv \
# # --minCov    10
# 
# ## !! Transfert to local gitignore/methylation_layerspecific_and_ME.tsv
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 2. Intra-individual correlation for layer-specific CpGs
# # ══════════════════════════════════════════════════════════════════════════════
# 
# meth <- fread(here("gitignore/methylation_layerspecific_and_ME.tsv"))
# 
# # ── Join germ_layer and category ─────────────────────────────────────────────
# loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
# loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
# tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])
# 
# meth <- merge(meth, tissue_to_layer,
#               by = "source_tissue_celltype", all.x = TRUE)
# meth <- merge(meth,
#               wideFull[, .(cpg_site = name, category)],
#               by = "cpg_site", all.x = TRUE)
# 
# message(sprintf("Missing germ_layer: %d | Missing category: %d",
#                 sum(is.na(meth$germ_layer)), sum(is.na(meth$category))))
# 
# # ── Multi-tissue patients ─────────────────────────────────────────────────────
# multi_patients <- meth[, .(n = uniqueN(source_tissue_celltype)),
#                        by = patient_id][n > 1, patient_id]
# message(sprintf("%d patients with >1 tissue", length(multi_patients)))
# # 32 patients with >1 tissue
# meth_multi <- meth[patient_id %in% multi_patients]
# 
# # ── Descriptive table of multi-tissue patients ────────────────────────────────
# patient_table <- meth_multi[, .(
#   n_tissues  = uniqueN(source_tissue_celltype),
#   germ_layers = paste(sort(unique(germ_layer)), collapse = "+"),
#   n_blood    = uniqueN(source_tissue_celltype[grepl(
#     "Blood", source_tissue_celltype, ignore.case=TRUE)]),
#   n_nonblood = uniqueN(source_tissue_celltype[!grepl(
#     "Blood", source_tissue_celltype, ignore.case=TRUE)])
# ), by = patient_id]
# 
# summary_table <- patient_table[, .(
#   n_patients     = .N,
#   median_tissues = as.numeric(median(n_tissues)),
#   range_tissues  = sprintf("%d-%d", min(n_tissues), max(n_tissues)),
#   n_blood_only   = sum(n_nonblood == 0),
#   pct_blood_only = round(100 * sum(n_nonblood == 0) / .N)
# ), by = germ_layers][order(germ_layers)]
# print(summary_table)
# 
# #    germ_layers n_patients median_tissues range_tissues n_blood_only pct_blood_only
# # 1:        Ecto          6            2.0           2-2            0              0
# # 2:        Endo         13            2.0           2-3            0              0
# # 3:   Endo+Meso          4            2.5           2-3            0              0
# # 4:        Meso          9            5.0           2-7            6             67
# 
# message(sprintf("Patients with cross-layer samples: %d",
#                 sum(grepl("\\+", patient_table$germ_layers))))
# # Patients with cross-layer samples: 4 (Endo+Meso)
# unique(meth[
#   meth$patient_id %in% unlist(
#     patient_table[patient_table$germ_layers %in% "Endo+Meso","patient_id"]),
#   c("source_tissue_celltype", "patient_id", "germ_layer")]) 
# 
# # source_tissue_celltype patient_id germ_layer
# # <char>     <char>     <char>
# #   1:               Colon - Endocrine        169       Endo
# # 2:             Colon - Macrophages        169       Meso
# # 3: Kidney glomerular - Endothelium        176       Meso
# # 4: Kidney glomerular - Endothelium        130       Meso
# # 5:  Kidney glomerular - Epithelium        130       Endo
# # 6:    Kidney glomerular - Podocyte        199       Endo
# # 7:    Kidney glomerular - Podocyte        176       Endo
# # 8:    Kidney tubular - Endothelium        199       Meso
# # 9:     Kidney tubular - Epithelium        130       Endo
# # 10:     Kidney tubular - Epithelium        176       Endo
# 
# ### How many patients have multiple SAME-layer tissues (different tissues)
# 
# # patients with >1 tissue WITHIN the same germ layer
# same_layer_multi <- meth_multi[, .(n_tissues_in_layer = uniqueN(source_tissue_celltype)),
#                                by = .(patient_id, germ_layer)][n_tissues_in_layer > 1]
# 
# print(same_layer_multi[order(germ_layer, -n_tissues_in_layer)])
# 
# message(sprintf("Patients with >1 tissue within the SAME germ layer: %d",
#                 uniqueN(same_layer_multi$patient_id)))
# # Patients with >1 tissue within the SAME germ layer: 30
# 
# # breakdown by layer
# print(same_layer_multi[, .(n_patients = uniqueN(patient_id)), by = germ_layer])
# #    germ_layer n_patients
# # 1:       Meso          9
# # 2:       Ecto          6
# # 3:       Endo         15
# 
# ## I have 4 patients with Endo+Meso (but only 2 with a given exact combination),
# ## 9 with multiple meso, 6 with multiple ecto
# ## and 15 with multiple endo.
# 
# ## I want to know which percentage of the 9111 Meso_specific CpGs show a high 
# # intra-individual correlation on average between 2 mesodermal tissues in the 9 
# # people with multiple meso tissues (if more than 2 meso tissues, do pairwise 
# # then mean to give one value), but a low correlation between 2 different
# # germ layer in the 4 patients with Endo+Meso tissues.
# 
# ## Same idea for the 3109 Endo_specific CpGs (in the 15 people with endo,
# # and the 4 Meso+endo)
# 
# ## As a negative control, I'd like the same percentages but for a random 10k
# # subset of the 10920663 constitutive CpGs and for the 10262820 ambiguous, in
# # both cases
# 
# ## I can't test Ecto_specific (not enough data)
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Within-individual same-layer vs cross-layer correlation per CpG
# #
# # GOAL (per CpG):
# #   1. SAME-LAYER:  for each patient who has >1 tissue from the relevant layer,
# #      compute the pairwise Pearson correlation between all their tissues for
# #      that CpG, then average those pairwise correlations -> one value per
# #      patient. Average across patients -> one "same-layer r" per CpG.
# #   2. CROSS-LAYER: for the 4 patients with Endo+Meso, compute Pearson r
# #      between their Endo tissue(s) and Meso tissue(s) for that CpG (pairwise
# #      then mean if >1 tissue per layer) -> one value per patient. Average
# #      across the 4 patients -> one "cross-layer r" per CpG.
# #   3. Report: % of CpGs in the category with same-layer r >= threshold AND
# #      cross-layer r < threshold (the "true layer-specific" signature)
# #
# # NOTE: this is computed WITHIN individuals (different tissues, same person),
# # not BETWEEN individuals as before.
# # ══════════════════════════════════════════════════════════════════════════════
# 
# # ── Helper: for ONE patient and a set of tissues, return one r per CpG ───────
# # (mean of all pairwise tissue correlations within that patient)
# patient_mean_pairwise_r <- function(meth_pat, tissues) {
#   if (length(tissues) < 2) return(NULL)
#   
#   # wide: rows = CpG, cols = tissue, values = methylation
#   wide <- dcast(meth_pat[source_tissue_celltype %in% tissues],
#                 cpg_site ~ source_tissue_celltype, value.var = "methylation",
#                 fun.aggregate = mean)
#   
#   tcols <- setdiff(names(wide), "cpg_site")
#   if (length(tcols) < 2) return(NULL)
#   
#   mat <- as.matrix(wide[, ..tcols])
#   rownames(mat) <- wide$cpg_site
#   
#   tpairs <- combn(tcols, 2, simplify = FALSE)
#   
#   # for each CpG (row), compute pairwise correlation across the FEW tissue
#   # values available — since each patient has only a handful of tissue
#   # measurements per CpG, "correlation across tissues" here really means:
#   # take the methylation values across tissues, and since we only have one
#   # measurement per tissue (not multiple replicates), we instead compute
#   # agreement using the actual paired values per tissue pair, per CpG
#   pair_vals <- rbindlist(lapply(tpairs, function(pr) {
#     m1 <- mat[, pr[1]]
#     m2 <- mat[, pr[2]]
#     data.table(cpg_site = rownames(mat), m1 = m1, m2 = m2)
#   }), idcol = "pair_id")
#   
#   pair_vals
# }
# 
# # ── Core function: same-layer within-individual r per CpG ────────────────────
# 
# ## Fetal-style same-layer r: one point per patient, then correlate across patients
# compute_same_layer_r_fetalstyle <- function(meth_sub, layer_name,
#                                             min_patients = 3, absolute = TRUE) {
#   layer_meth <- meth_sub[germ_layer == layer_name]
#   eligible <- layer_meth[, .(n = uniqueN(source_tissue_celltype)),
#                          by = patient_id][n > 1, patient_id]
#   if (length(eligible) < min_patients) return(NULL)
#   
#   # per patient: all tissue pairs, tissue names ordered so m1/m2 are deterministic
#   all_pairs <- rbindlist(lapply(eligible, function(pid) {
#     pat  <- layer_meth[patient_id == pid]
#     tiss <- sort(unique(pat$source_tissue_celltype))       # deterministic order
#     if (length(tiss) < 2) return(NULL)
#     rbindlist(lapply(combn(tiss, 2, simplify = FALSE), function(pr) {
#       d1 <- pat[source_tissue_celltype == pr[1], .(cpg_site, m1 = methylation)]
#       d2 <- pat[source_tissue_celltype == pr[2], .(cpg_site, m2 = methylation)]
#       merge(d1, d2, by = "cpg_site")[, patient_id := pid]
#     }), fill = TRUE)
#   }), fill = TRUE)
#   if (is.null(all_pairs) || !nrow(all_pairs)) return(NULL)
#   
#   # ---- KEY STEP: collapse each patient's pairs to ONE point (equal weighting) --
#   per_patient <- all_pairs[!is.na(m1) & !is.na(m2),
#                            .(m1 = mean(m1), m2 = mean(m2)),
#                            by = .(cpg_site, patient_id)]
#   
#   per_patient[, {
#     if (.N >= min_patients) {
#       r <- suppressWarnings(cor(m1, m2, method = "pearson"))
#       list(r = if (absolute) abs(r) else r, n_obs = .N)
#     } else list(r = NA_real_, n_obs = .N)
#   }, by = cpg_site]
# }
# 
# # ── Core function: cross-layer within-individual r per CpG ───────────────────
# # Same logic but for the 4 Endo+Meso patients: pair each patient's Endo
# # tissue value with their Meso tissue value (mean if >1 tissue per layer
# # for that patient), pool across the 4 patients, correlate per CpG
# compute_cross_layer_r <- function(meth_sub, layer1, layer2) {
#   
#   d1 <- meth_sub[germ_layer == layer1,
#                  .(m1 = mean(methylation, na.rm = TRUE)),
#                  by = .(cpg_site, patient_id)]
#   d2 <- meth_sub[germ_layer == layer2,
#                  .(m2 = mean(methylation, na.rm = TRUE)),
#                  by = .(cpg_site, patient_id)]
#   
#   dm <- merge(d1, d2, by = c("cpg_site", "patient_id"))
#   
#   dm[, {
#     idx <- !is.na(m1) & !is.na(m2)
#     if (sum(idx) >= 3) {
#       r <- suppressWarnings(cor(m1[idx], m2[idx], method = "pearson"))
#       list(r = abs(r), n_obs = sum(idx))
#     } else list(r = NA_real_, n_obs = sum(idx))
#   }, by = cpg_site]
# }
# 
# extract_layer_specific_cpgs <- function(meth_sub,
#                                         own_layer, other_layer,
#                                         r_high = NULL,
#                                         r_low  = NULL,
#                                         r_low_cross = NULL,   # <- add this
#                                         min_same_obs  = 3,
#                                         min_cross_obs = 3)
# {
#   if (is.null(r_high)) r_high <- thr[layer == own_layer,   r_high]
#   if (is.null(r_low))  r_low  <- thr[layer == other_layer, r_low]
#   if (is.null(r_low_cross)) r_low_cross <- thr_cross$r_low
#   
#   r_own   <- compute_same_layer_r_fetalstyle(meth_sub, own_layer)
#   r_other <- compute_same_layer_r_fetalstyle(meth_sub, other_layer)
#   r_cross <- compute_cross_layer_r(meth_sub, other_layer, own_layer)
#   if (is.null(r_own) || is.null(r_other) || is.null(r_cross)) return(NULL)
#   
#   setnames(r_own,   c("r", "n_obs"), c("r_own",   "n_own"))
#   setnames(r_other, c("r", "n_obs"), c("r_other", "n_other"))
#   setnames(r_cross, c("r", "n_obs"), c("r_cross", "n_cross"))
#   
#   m <- Reduce(function(a, b) merge(a, b, by = "cpg_site"),
#               list(r_own, r_other, r_cross))
#   
#   hits <- m[!is.na(r_own) & !is.na(r_other) & !is.na(r_cross) &
#               n_own   >= min_same_obs &
#               n_other >= min_same_obs &
#               n_cross >= min_cross_obs &
#               r_own   >= r_high &
#               r_other <  r_low  &
#               r_cross <  r_low_cross]
#   hits[order(-r_own)]
# }
# 
# set.seed(1234)
# constitutive_sample <- sample(wideFull[category == "constitutive", name], 10000)
# ambiguous_sample    <- sample(wideFull[category == "ambiguous",    name], 10000)
# writeLines(c(constitutive_sample, ambiguous_sample),
#            here("B_MultiTissues/dataOut/control_sample20k.txt"))
# message("Written: control_sample20k.txt (20000 CpGs: 10k constitutive + 10k ambiguous)")
# 
# ## In pchuckle (after git pull):
# # source /share/apps/source_files/python/python-3.13.0a6.source
# # cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# # python3 S00_extractRawMethylationForTargetCpG.py \
# # --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/control_sample20k.txt \
# # --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# # --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# # --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# # --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_control20k.tsv \
# # --minCov    10
# 
# ## !! Transfert gitignore/methylation_control20k.tsv to local
# 
# # ── Load control methylation (extracted previously) ───────────────────────────
# meth_control <- fread(here("gitignore/methylation_control20k.tsv"))
# meth_control <- merge(meth_control, tissue_to_layer,
#                       by = "source_tissue_celltype", all.x = TRUE)
# 
# # re-assign categories (constitutive vs ambiguous)
# set.seed(1234)
# constitutive_sample <- sample(wideFull[category == "constitutive", name], 10000)
# ambiguous_sample    <- sample(wideFull[category == "ambiguous",    name], 10000)
# 
# meth_control[, category := fcase(
#   cpg_site %in% constitutive_sample, "constitutive",
#   cpg_site %in% ambiguous_sample,    "ambiguous",
#   default = NA_character_
# )]
# 
# meth_control_multi <- meth_control[patient_id %in% multi_patients]
# 
# message(sprintf("Control CpGs loaded: %d constitutive, %d ambiguous",
#                 uniqueN(meth_control_multi[category == "constitutive", cpg_site]),
#                 uniqueN(meth_control_multi[category == "ambiguous",    cpg_site])))
# # Control CpGs loaded: 12 constitutive, 7 ambiguous
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Run for Meso_specific, Endo_specific and controls
# # ══════════════════════════════════════════════════════════════════════════════
# 
# # ── Meso_specific ─────────────────────────────────────────────────────────────
# meth_meso_specific <- meth_multi[category == "Meso_specific"]
# r_same_meso        <- compute_same_layer_r_fetalstyle(meth_meso_specific, "Meso")
# 
# # ── Endo_specific ─────────────────────────────────────────────────────────────
# meth_endo_specific <- meth_multi[category == "Endo_specific"]
# r_same_endo        <- compute_same_layer_r_fetalstyle(meth_endo_specific, "Endo")
# 
# # ── ME same-layer r ───────────────────────────────────────────────────────────
# meth_ME <- meth_multi[category == "ME"]
# 
# # ── Wrong-layer specificity checks ───────────────────────────────────────────
# r_same_endo_for_meso <- compute_same_layer_r_fetalstyle(meth_meso_specific, "Endo")
# r_same_meso_for_endo <- compute_same_layer_r_fetalstyle(meth_endo_specific, "Meso")
# r_same_meso_ME <- compute_same_layer_r_fetalstyle(meth_ME, "Meso")
# r_same_endo_ME <- compute_same_layer_r_fetalstyle(meth_ME, "Endo")
# 
# # ── Controls in both Meso and Endo contexts ───────────────────────────────────
# r_same_meso_const      <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "constitutive"], "Meso")
# r_same_endo_const <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "constitutive"], "Endo")
# r_same_meso_amb        <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "ambiguous"],    "Meso")
# r_same_endo_amb   <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "ambiguous"],    "Endo")
# 
# # ── Extract "true layer-specific" CpGs on THREE criteria ──────────────────────
# #   1. high  same-layer intra-individual r in the OWN layer   (r_own   >= r_high)
# #   2. low   same-layer intra-individual r in the OTHER layer (r_other <  r_low)
# #   3. low   cross-layer within-person r (own <-> other)       (r_cross <  r_low)
# 
# # per-layer empirical null from constitutive CpGs
# null_meso <- r_same_meso_const[!is.na(r), r]
# null_endo <- r_same_endo_const[!is.na(r), r]
# 
# r_crit <- function(n, alpha = 0.05) {
#   tc <- qt(1 - alpha/2, df = n - 2)
#   tc / sqrt(tc^2 + n - 2)
# }
# 
# thr <- data.table(
#   layer  = c("Meso", "Endo"),
#   # "high" = beyond the 95th percentile of noise (5% FPR by construction)
#   r_high = c(quantile(null_meso, 0.95), quantile(null_endo, 0.95)),
#   # "low"  = below the median of noise (indistinguishable from chance)
#   r_low  = c(median(null_meso),         median(null_endo)),
#   # mean +/- sd of the noise, for reference
#   mu     = c(mean(null_meso),   mean(null_endo)),
#   sd     = c(sd(null_meso),     sd(null_endo)),
#   n_pat  = c(9, 15),
#   r_crit = c(r_crit(9), r_crit(15)))
# print(thr)
# # layer    r_high     r_low        mu        sd n_pat    r_crit
# # <char>     <num>     <num>     <num>     <num> <num>     <num>
# # 1:   Meso 0.6810714 0.2506185 0.2835287 0.2055916     9 0.6663836
# # 2:   Endo 0.6955997 0.1879514 0.2531251 0.2448413    15 0.5139775
# 
# # ── Summary function ──────────────────────────────────────────────────────────
# make_summary_samelayer <- function(r_same, layer_tested, category_name,
#                                    thr_tab = thr) {
#   s <- r_same[!is.na(r)]
#   r_hi <- thr_tab[layer == layer_tested, r_high]      # layer-specific "high" bar
#   data.table(
#     category          = category_name,
#     same_layer_tested = layer_tested,
#     r_high_used       = r_hi,
#     n_cpgs_tested     = uniqueN(s$cpg_site),
#     pct_high_same     = 100 * mean(s$r >= r_hi),
#     mean_r            = mean(s$r),
#     median_r          = median(s$r))
# }
# 
# # ── Build summary ─────────────────────────────────────────────────────────────
# results_summary <- rbindlist(list(
#   make_summary_samelayer(r_same_meso,          "Meso", "Meso_specific"),
#   make_summary_samelayer(r_same_endo_for_meso, "Endo", "Meso_specific"),
#   make_summary_samelayer(r_same_endo,          "Endo", "Endo_specific"),
#   make_summary_samelayer(r_same_meso_for_endo, "Meso", "Endo_specific"),
#   make_summary_samelayer(r_same_meso_const,    "Meso", "constitutive"),
#   make_summary_samelayer(r_same_endo_const,    "Endo", "constitutive"),
#   make_summary_samelayer(r_same_meso_amb,      "Meso", "ambiguous"),
#   make_summary_samelayer(r_same_endo_amb,      "Endo", "ambiguous")
# ))
# 
# print(results_summary)
# #         category same_layer_tested r_high_used n_cpgs_tested pct_high_same    mean_r  median_r
# # <char>            <char>       <num>         <int>         <num>     <num>     <num>
# # 1: Meso_specific              Meso   0.6810714            30     83.333333 0.7890745 0.8555211
# # 2: Meso_specific              Endo   0.6955997            30      0.000000 0.2577406 0.1955097
# # 3: Endo_specific              Endo   0.6955997             3     66.666667 0.8434886 0.8871436
# # 4: Endo_specific              Meso   0.6810714             3      0.000000 0.4409972 0.4423979
# # 5:  constitutive              Meso   0.6810714            12      8.333333 0.2835287 0.2506185
# # 6:  constitutive              Endo   0.6955997            12      8.333333 0.2531251 0.1879514
# # 7:     ambiguous              Meso   0.6810714             7     71.428571 0.6753872 0.7007080
# # 8:     ambiguous              Endo   0.6955997             7     28.571429 0.4063682 0.1384330
# 
# # ── Annotate and plot ─────────────────────────────────────────────────────────
# results_summary[, expected := fcase(
#   category == "Meso_specific" & same_layer_tested == "Meso", "own layer",
#   category == "Meso_specific" & same_layer_tested == "Endo", "wrong layer",
#   category == "Endo_specific" & same_layer_tested == "Endo", "own layer",
#   category == "Endo_specific" & same_layer_tested == "Meso", "wrong layer",
#   default = "control"
# )]
# 
# results_summary[, category_f := factor(category,
#                                        levels = c("Meso_specific", "Endo_specific", "constitutive", "ambiguous"))]
# 
# bar_colours <- c(
#   "own layer"   = "#2166AC",
#   "wrong layer" = "#D55E00",
#   "control"     = "grey60"
# )
# 
# # empirical FPR: % of constitutive CpGs passing the same threshold, per layer
# noise_dt <- results_summary[category == "constitutive",
#                             .(same_layer_tested, floor_pct = pct_high_same)]
# 
# p1 <- ggplot(results_summary,
#              aes(x = category_f, y = pct_high_same,
#                  fill = expected)) +
#   geom_col(position = "dodge", width = 0.6,
#            colour = "grey30", linewidth = 0.3) +
#   geom_text(aes(label = sprintf("%.1f%%", pct_high_same)),
#             position = position_dodge(width = 0.6),
#             vjust = -0.4, size = 3) +
#   facet_wrap(~ same_layer_tested, nrow = 1,
#              labeller = labeller(same_layer_tested = c(
#                Meso = "Tested in Meso tissues",
#                Endo = "Tested in Endo tissues"))) +
#   geom_hline(data = noise_dt, aes(yintercept = floor_pct),
#              linetype = "dashed", colour = "grey30", linewidth = 0.4) +
#   geom_text(data = noise_dt,
#             aes(x = 0.6, y = floor_pct, label = sprintf("noise floor", floor_pct)),
#             hjust = 0, vjust = -0.4, size = 2.6, colour = "grey30", inherit.aes = FALSE) +
#   scale_fill_manual(values = bar_colours, name = "Layer context") +
#   scale_y_continuous("% CpGs above the layer-specific |r| cutoff", limits = c(0, NA)) +
#   scale_x_discrete("CpG category") +
#   theme_bw(base_size = 11) +
#   theme(panel.grid.minor  = element_blank(),
#         axis.text.x       = element_text(angle = 30, hjust = 1),
#         strip.text        = element_text(face = "bold"),
#         legend.position   = "right") +
#   ggtitle("Intra-individual same-layer methylation concordance",
#           subtitle = sprintf("Own layer (blue) should exceed the noise floor by more than the wrong layer (orange)\nMeso: %d patients (|r| \u2265 %.2f) | Endo: %d patients (|r| \u2265 %.2f)",
#                              9, thr[layer=="Meso", r_high], 15, thr[layer=="Endo", r_high]))
# 
# # ── Density plot of same-layer r per category and layer context ───────────────
# 
# # ── Pool every category × layer tested into one table ────────────────────────
# density_dt <- rbindlist(list(
#   # own-layer
#   r_same_meso[!is.na(r),          .(cpg_site, r, category = "Meso_specific", layer_context = "Meso")],
#   r_same_endo[!is.na(r),          .(cpg_site, r, category = "Endo_specific", layer_context = "Endo")],
#   # wrong-layer
#   r_same_endo_for_meso[!is.na(r), .(cpg_site, r, category = "Meso_specific", layer_context = "Endo")],
#   r_same_meso_for_endo[!is.na(r), .(cpg_site, r, category = "Endo_specific", layer_context = "Meso")],
#   # MEs (expected high in BOTH layers)
#   r_same_meso_ME[!is.na(r),       .(cpg_site, r, category = "ME",           layer_context = "Meso")],
#   r_same_endo_ME[!is.na(r),       .(cpg_site, r, category = "ME",           layer_context = "Endo")],
#   # controls
#   r_same_meso_const[!is.na(r),    .(cpg_site, r, category = "constitutive", layer_context = "Meso")],
#   r_same_endo_const[!is.na(r),    .(cpg_site, r, category = "constitutive", layer_context = "Endo")],
#   r_same_meso_amb[!is.na(r),      .(cpg_site, r, category = "ambiguous",    layer_context = "Meso")],
#   r_same_endo_amb[!is.na(r),      .(cpg_site, r, category = "ambiguous",    layer_context = "Endo")]
# ))
# 
# density_dt[, category_f := factor(category,
#                                   levels = c("ME", "Endo_specific", "Meso_specific", "ambiguous", "constitutive"))]
# density_dt[, layer_f := factor(layer_context, levels = c("Meso", "Endo"),
#                                labels = c("Tested in Meso tissues", "Tested in Endo tissues"))]
# 
# # n per group for the subtitle
# n_lab <- density_dt[, .(n = uniqueN(cpg_site)), by = .(category_f, layer_f)]
# 
# null_dt <- data.table(
#   layer_f = factor(c("Tested in Meso tissues","Tested in Endo tissues"),
#                    levels = levels(density_dt$layer_f)),
#   null_r  = c(thr[layer=="Meso", mu], thr[layer=="Endo", mu]))   # empirical
# 
# thr_dt <- data.table(
#   layer_f = factor(c("Tested in Meso tissues","Tested in Endo tissues"),
#                    levels = levels(density_dt$layer_f)),
#   r_high  = c(thr[layer=="Meso", r_high], thr[layer=="Endo", r_high]))
# 
# p2 <- ggplot(density_dt, aes(x = r, colour = category_f, fill = category_f)) +
#   geom_density(alpha = 0.20, linewidth = 0.8, bounds = c(0, 1)) +
#   geom_vline(data = thr_dt, aes(xintercept = r_high),
#              linetype = "dashed", colour = "grey40") +
#   geom_vline(data = null_dt, aes(xintercept = null_r),
#              linetype = "dotted", colour = "grey50") +
#   facet_wrap(~ layer_f, nrow = 1) +
#   scale_colour_manual(values = category_colours, drop = FALSE, name = "CpG category") +
#   scale_fill_manual(values   = category_colours, drop = FALSE, name = "CpG category") +
#   scale_x_continuous("Same-layer |r| (per-patient means, across patients)",
#                      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   scale_y_continuous("Density") +
#   theme_bw(base_size = 11) +
#   theme(panel.grid.minor = element_blank(),
#         strip.text       = element_text(face = "bold"),
#         legend.position  = "right") +
#   ggtitle("Same-layer methylation concordance by CpG category (fetal-style |r|)",
#           subtitle = sprintf(paste0(
#             "Each patient contributes one point (mean of their tissue pairs); |Pearson r| taken across patients\n",
#             "Dotted = mean |r| of constitutive CpGs (empirical noise) | ",
#             "Dashed = 95th pct of that noise, used as the \"high\" cutoff\n",
#             "Meso: %d patients, noise %.2f, cutoff |r| \u2265 %.2f  |  ",
#             "Endo: %d patients, noise %.2f, cutoff |r| \u2265 %.2f"),
#             9,  thr[layer == "Meso", mu], thr[layer == "Meso", r_high],
#             15, thr[layer == "Endo", mu], thr[layer == "Endo", r_high]))
# 
# # ── toggle: exclude MEs from the density panel ───────────────────────────────
# drop_ME <- TRUE
# plot_dt <- if (drop_ME) density_dt[category != "ME"] else density_dt
# 
# p3 <- ggplot(plot_dt, aes(x = r, colour = category_f, fill = category_f)) +
#   geom_density(alpha = 0.20, linewidth = 0.8, bounds = c(0, 1)) +
#   geom_vline(data = thr_dt,  aes(xintercept = r_high),
#              linetype = "dashed", colour = "grey40") +
#   geom_vline(data = null_dt, aes(xintercept = null_r),
#              linetype = "dotted", colour = "grey50") +
#   facet_wrap(~ layer_f, nrow = 1) +
#   scale_colour_manual(values = category_colours, drop = TRUE, name = "CpG category") +
#   scale_fill_manual(values   = category_colours, drop = TRUE, name = "CpG category") +
#   scale_x_continuous("Same-layer |r| (per-patient means, then |Pearson r| across patients)",
#                      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   scale_y_continuous("Density") +
#   theme_bw(base_size = 11) +
#   theme(panel.grid.minor = element_blank(),
#         strip.text       = element_text(face = "bold"),
#         legend.position  = "right") +
#   ggtitle("Same-layer methylation concordance by CpG category (fetal-style |r|)",
#           subtitle = sprintf(paste0(
#             "Each patient contributes one point (mean of their tissue pairs); |Pearson r| taken across patients\n",
#             "Dotted = mean |r| of constitutive CpGs (empirical noise) | ",
#             "Dashed = 95th pct of that noise, used as the \"high\" cutoff\n",
#             "Meso: %d patients, noise %.2f, cutoff |r| \u2265 %.2f  |  ",
#             "Endo: %d patients, noise %.2f, cutoff |r| \u2265 %.2f"),
#             9,  thr[layer == "Meso", mu], thr[layer == "Meso", r_high],
#             15, thr[layer == "Endo", mu], thr[layer == "Endo", r_high]))
# 
# ggplot2::ggsave(
#   here("B_MultiTissues/dataOut/figures/script04/sameLayercR_concordance.pdf"),
#   plot_grid(p1, p2, p3, labels = c("A", "B", "C"), nrow = 3), width = 10, height = 15)
# 
# # Meso_specific CpGs are highly variable between individuals (by selection)
# # Within each individual, their methylation level is highly consistent across all blood cell types (r≈1, Plot 1)
# # This within-individual concordance across independent blood lineages is the hallmark of pre-haematopoietic stochastic establishment
# # Controls show r≈0, confirming this is not a general property of blood measurements
# 
# ###### Add now the 3rd criteria, low r between both layers
# 
# min_same_obs  <- 3
# min_cross_obs <- 3
# 
# # per-CpG: the three r's for one category under one pipeline (own vs other) ------
# tri_rs <- function(meth_sub, own_layer, other_layer, category_name) {
#   ro <- compute_same_layer_r_fetalstyle(meth_sub, own_layer)
#   rt <- compute_same_layer_r_fetalstyle(meth_sub, other_layer)
#   rc <- compute_cross_layer_r(meth_sub, other_layer, own_layer)
#   if (is.null(ro)) return(NULL)
#   setnames(ro, c("r","n_obs"), c("r_own","n_own"))
#   base <- ro[, .(cpg_site, r_own, n_own)]
#   if (!is.null(rt)) { setnames(rt, c("r","n_obs"), c("r_other","n_other"))
#     base <- merge(base, rt[, .(cpg_site, r_other, n_other)], by="cpg_site", all.x=TRUE)
#   } else base[, `:=`(r_other = NA_real_, n_other = 0L)]
#   if (!is.null(rc)) { setnames(rc, c("r","n_obs"), c("r_cross","n_cross"))
#     base <- merge(base, rc[, .(cpg_site, r_cross, n_cross)], by="cpg_site", all.x=TRUE)
#   } else base[, `:=`(r_cross = NA_real_, n_cross = 0L)]
#   base[, `:=`(category = category_name, own_layer = own_layer)][]
# }
# 
# # cumulative pass counts over the 3 criteria ------------------------------------
# tri_steps <- function(tab) {
#   if (is.null(tab)) return(NULL)
#   own <- tab$own_layer[1]
#   oth <- setdiff(c("Meso","Endo"), own)
#   rh  <- thr[layer == own, r_high]
#   rlo <- thr[layer == oth, r_low]
#   rlc <- thr_cross$r_low
#   tested <- tab[!is.na(r_own) & n_own >= min_same_obs]
#   N  <- nrow(tested)
#   p1 <- tested[r_own >= rh]
#   p2 <- p1[!is.na(r_other) & n_other >= min_same_obs & r_other < rlo]
#   p3 <- p2[!is.na(r_cross) & n_cross >= min_cross_obs & r_cross < rlc]
#   data.table(
#     category  = tab$category[1],
#     own_layer = tab$own_layer[1],
#     step = factor(c("1. high own r", "2. + low other r", "3. + low cross r"),
#                   levels = c("1. high own r", "2. + low other r", "3. + low cross r")),
#     n   = c(nrow(p1), nrow(p2), nrow(p3)),
#     pct = 100 * c(nrow(p1), nrow(p2), nrow(p3)) / N,
#     n_tested = N)
# }
# 
# # ── Empirical null for the CROSS-layer r (Endo x Meso, n = 4 patients) ────────
# null_cross <- compute_cross_layer_r(
#   meth_control_multi[category == "constitutive"], "Endo", "Meso")[!is.na(r), r]
# 
# thr_cross <- data.table(
#   comparison = "Endo_x_Meso",
#   n_pat      = length(unique(meth_control_multi[germ_layer %in% c("Endo","Meso"), patient_id])),
#   # "high" (ME-like systemic concordance): beyond 95th pct of noise
#   r_high     = quantile(null_cross, 0.95),
#   # "low" (layer-specific): below the MEDIAN of noise, i.e. no more concordant
#   #  than a constitutive CpG. Tighten to 0.25 quantile for a stricter set.
#   r_low      = median(null_cross),
#   r_low_strict = quantile(null_cross, 0.25),
#   mu         = mean(null_cross),
#   sd         = sd(null_cross))
# print(thr_cross)
# # comparison n_pat    r_high     r_low r_low_strict        mu        sd
# # <char> <int>     <num>     <num>        <num>     <num>     <num>
# # 1: Endo_x_Meso    26 0.8926701 0.4055105    0.1161466 0.3809395 0.3040298
# 
# tri_tabs <- rbindlist(list(
#   # Meso pipeline (own = Meso, other = Endo)
#   tri_steps(tri_rs(meth_meso_specific,                       "Meso","Endo","Meso_specific")),
#   tri_steps(tri_rs(meth_endo_specific,                       "Meso","Endo","Endo_specific")),
#   tri_steps(tri_rs(meth_control_multi[category=="constitutive"], "Meso","Endo","constitutive")),
#   tri_steps(tri_rs(meth_control_multi[category=="ambiguous"],    "Meso","Endo","ambiguous")),
#   # Endo pipeline (own = Endo, other = Meso)
#   tri_steps(tri_rs(meth_endo_specific,                       "Endo","Meso","Endo_specific")),
#   tri_steps(tri_rs(meth_meso_specific,                       "Endo","Meso","Meso_specific")),
#   tri_steps(tri_rs(meth_control_multi[category=="constitutive"], "Endo","Meso","constitutive")),
#   tri_steps(tri_rs(meth_control_multi[category=="ambiguous"],    "Endo","Meso","ambiguous"))
# ), fill = TRUE)
# 
# tri_tabs[, category_f := factor(category,
#                                 levels = c("Meso_specific","Endo_specific","constitutive","ambiguous"))]
# 
# p <- ggplot(tri_tabs, aes(category_f, pct, fill = step)) +
#   geom_col(position = position_dodge(0.7), width = 0.65,
#            colour = "grey30", linewidth = 0.3) +
#   geom_text(aes(label = n), position = position_dodge(0.7),
#             vjust = -0.3, size = 2.7) +
#   facet_wrap(~ own_layer, nrow = 1,
#              labeller = labeller(own_layer = c(
#                Meso = "Own = Meso (other = Endo)",
#                Endo = "Own = Endo (other = Meso)"))) +
#   scale_fill_brewer(palette = "Blues", name = "Cumulative criterion") +
#   scale_y_continuous("% of CpGs tested in own layer", limits = c(0, NA)) +
#   scale_x_discrete("CpG category") +
#   theme_bw(base_size = 11) +
#   theme(panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 30, hjust = 1),
#         strip.text  = element_text(face = "bold"),
#         legend.position = "right") +
#   ggtitle("Layer-specific selection: cumulative pass rate over 3 criteria",
#           subtitle = sprintf(
#             "high own r (Meso \u2265%.2f, Endo \u2265%.2f)  \u2192  + low other r (Meso <%.2f, Endo <%.2f)  \u2192  + low cross-layer r (<%.2f)",
#             thr[layer == "Meso", r_high], thr[layer == "Endo", r_high],
#             thr[layer == "Meso", r_low],  thr[layer == "Endo", r_low],
#             thr_cross$r_low))
# 
# ## Higherbackground concordance for endo than meso
# ggplot2::ggsave(
#   filename = here("B_MultiTissues/dataOut/figures/script04/layerCandidatesSelection.pdf"),
#   plot     = p, width    = 8, height = 5)
# 
# #################################
# ## Putative MEs (systemic)     ##
# #################################
# # ME = high Pr(HV) in all 3 layers (already the "ME" category) PLUS the systemic
# # concordance signature: high intra-layer r in each layer AND high cross-layer r.
# # (Layer-specific = high own / low cross; ME = high own / HIGH cross.)
# 
# extract_ME_cpgs <- function(meth_sub,
#                             r_high        = 0.5,   # concordance bar (same as own-layer)
#                             min_same_obs  = 3,
#                             min_cross_obs = 3) {
#   # within-individual concordance in each layer that has multi-tissue patients
#   r_endo <- compute_same_layer_r_fetalstyle(meth_sub, "Endo")
#   r_meso <- compute_same_layer_r_fetalstyle(meth_sub, "Meso")
#   # cross-layer (systemic) concordance — the ME hallmark
#   r_cross <- compute_cross_layer_r(meth_sub, "Endo", "Meso")
#   if (is.null(r_endo) || is.null(r_meso) || is.null(r_cross)) return(NULL)
#   
#   setnames(r_endo,  c("r","n_obs"), c("r_endo",  "n_endo"))
#   setnames(r_meso,  c("r","n_obs"), c("r_meso",  "n_meso"))
#   setnames(r_cross, c("r","n_obs"), c("r_cross", "n_cross"))
#   
#   m <- Reduce(function(a,b) merge(a,b,by="cpg_site"),
#               list(r_endo, r_meso, r_cross))
#   
#   hits <- m[!is.na(r_endo) & !is.na(r_meso) & !is.na(r_cross) &
#               n_endo  >= min_same_obs &
#               n_meso  >= min_same_obs &
#               n_cross >= min_cross_obs &
#               r_endo  >= r_high &      # high concordance in endo
#               r_meso  >= r_high &      # high concordance in meso
#               r_cross >= r_high]       # AND high cross-layer  -> systemic
#   hits[order(-r_cross)]
# }
# 
# ME_hits <- extract_ME_cpgs(meth_ME)
# message(sprintf("Putative systemic-ME CpGs: %d", if (is.null(ME_hits)) 0L else nrow(ME_hits)))
# # Putative systemic-ME CpGs: 4098
# 
# #################################
# ## Save the interesting targets##
# #################################
# meso_hits <- extract_layer_specific_cpgs(meth_meso_specific, 
#                                          own_layer = "Meso", other_layer = "Endo")
# endo_hits <- extract_layer_specific_cpgs(meth_endo_specific,
#                                          own_layer = "Endo", other_layer = "Meso")
# 
# message(sprintf("Meso_specific true-signature CpGs: %d", nrow(meso_hits)))
# ## Meso_specific true-signature CpGs: 3
# 
# message(sprintf("Endo_specific true-signature CpGs: %d", nrow(endo_hits)))
# # Endo_specific true-signature CpGs: 0
# 
# #################################
# ## Annotation of these targets ##
# #################################
# 
# # ── Shared annotation objects (built once) ────────────────────────────────────
# txdb     <- TxDb.Hsapiens.UCSC.hg38.knownGene
# genes_gr <- genes(txdb)
# genes_gr$symbol <- mapIds(org.Hs.eg.db, names(genes_gr), "SYMBOL", "ENTREZID")
# 
# cpg_to_gr <- function(cpg) {                     # "chr7_107543290" -> 1bp GRanges
#   GRanges(sub("_.*", "", cpg),
#           IRanges(as.integer(sub(".*_", "", cpg)), width = 1),
#           cpg_site = cpg)
# }
# 
# # universe of tested CpGs (built once, reused by every call)
# uni <- copy(wideFull[, .(cpg_site = name, category)])
# uni[, `:=`(chr = sub("_.*", "", cpg_site),
#            pos = as.integer(sub(".*_", "", cpg_site)))]
# uni_gr <- GRanges(uni$chr, IRanges(uni$pos, width = 1))
# 
# # ── Main function ─────────────────────────────────────────────────────────────
# annotate_layer_hits <- function(hits,                       # data.table with cpg_site
#                                 label,                     # e.g. "meso" / "endo"
#                                 genes_gr, uni, uni_gr,
#                                 gap        = 50,           # bp to merge hit clusters
#                                 min_hits   = 2,            # min hits per cluster
#                                 drop_regex = "LOC|LINC",   # gene names to discard
#                                 flank      = 5000,
#                                 out_dir    = here("B_MultiTissues/dataOut")) {
#   
#   hit_col <- paste0("is_", label, "_hit")
#   hit_set <- hits$cpg_site
#   
#   ## (1a) annotate each hit with nearest gene ---------------------------------
#   hits_gr <- cpg_to_gr(hit_set)
#   nr  <- distanceToNearest(hits_gr, genes_gr, ignore.strand = TRUE)
#   ann <- as.data.table(hits)
#   ann[, `:=`(gene = NA_character_, dist_to_gene = NA_integer_)]
#   ann[queryHits(nr),
#       `:=`(gene         = genes_gr$symbol[subjectHits(nr)],
#            dist_to_gene = mcols(nr)$distance)]
#   ann[, `:=`(pos = as.integer(sub(".*_", "", cpg_site)),
#              chr = sub("_.*", "", cpg_site))]
#   setorder(ann, chr, pos)
#   
#   ## (1b) cluster hits within `gap` bp ----------------------------------------
#   hs <- cpg_to_gr(ann$cpg_site)
#   clust <- reduce(resize(hs, 1), min.gapwidth = gap + 1)
#   ov    <- findOverlaps(hs, clust)
#   ann[, cluster_id := subjectHits(ov)[match(seq_len(.N), queryHits(ov))]]
#   
#   clust_dt <- ann[, .(n_hits = .N, span = max(pos) - min(pos)),
#                   by = .(gene, cluster_id)][n_hits >= min_hits]
#   if (nzchar(drop_regex))
#     clust_dt <- clust_dt[!grepl(drop_regex, gene) & !is.na(gene)]
#   clust_dt[, density := n_hits / (span + 1)]
#   setorder(clust_dt, -n_hits, span)
#   
#   ## (2) all CpGs in each top gene + flank, flagged ---------------------------
#   top_genes <- unique(clust_dt$gene)
#   extract_gene_cpgs <- function(sym) {
#     g <- genes_gr[which(genes_gr$symbol == sym)]
#     if (!length(g)) return(NULL)
#     g   <- g[1]
#     win <- GRanges(as.character(seqnames(g)),
#                    IRanges(start(g) - flank, end(g) + flank))
#     out <- uni[which(overlapsAny(uni_gr, win, ignore.strand = TRUE))]
#     out[, `:=`(gene = sym,
#                gene_start = start(g), gene_end = end(g),
#                in_gene_body = pos >= start(g) & pos <= end(g))]
#     out[, (hit_col) := cpg_site %in% hit_set]
#     out[order(pos)]
#   }
#   gene_cpgs <- rbindlist(lapply(top_genes, extract_gene_cpgs), fill = TRUE)
#   
#   fwrite(gene_cpgs,
#          file.path(out_dir, sprintf("%s_hits_topGenes_allCpGs.csv", label)))
#   
#   list(ann = ann, clusters = clust_dt, gene_cpgs = gene_cpgs)
# }
# 
# for (g in c(50, 500, 1000, 2000)) {
#   cl <- annotate_layer_hits(endo_hits, "endo_tmp", genes_gr, uni, uni_gr,
#                             gap = g, out_dir = tempdir())$clusters
#   message(sprintf("gap=%4d bp: %d clusters (>=2 hits), %f genes, max n_hits=%f",
#                   g, nrow(cl), uniqueN(cl$gene), max(cl$n_hits)))
# }
# 
# meso_res <- annotate_layer_hits(meso_hits, "meso", genes_gr, uni, uni_gr, gap = 100)
# endo_res <- annotate_layer_hits(endo_hits, "endo", genes_gr, uni, uni_gr, gap = 100)
# ME_res   <- annotate_layer_hits(ME_hits,   "ME",   genes_gr, uni, uni_gr, gap = 100)
# 
# meso_res
# ME_res
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## =============================================================================
# ## (4) PROTOCADHERIN LOCUS — per category
# ## =============================================================================
# pcdh <- GRanges("chr5", IRanges(140710000, 141510000))   # hg38 PCDH clusters
# pcdh_res <- rbindlist(lapply(names(categories), function(cat) {
#   fg <- cat_gr[[cat]]
#   bg <- makeGRfromMyCpGPos(setdiff(all_cpg, categories[[cat]]), "bg")
#   fg_in <- sum(overlapsAny(fg, pcdh, ignore.strand = TRUE))
#   bg_in <- sum(overlapsAny(bg, pcdh, ignore.strand = TRUE))
#   ft <- fisher.test(matrix(c(fg_in, length(fg)-fg_in,
#                              bg_in, length(bg)-bg_in), nrow = 2, byrow = TRUE),
#                     alternative = "greater")
#   data.table(category = cat, fg_in = fg_in, fg_tot = length(fg),
#              odds_ratio = unname(ft$estimate), pvalue = ft$p.value)
# }))
# pcdh_res[, p.adj := p.adjust(pvalue, "BH")]
# print(pcdh_res)
# 
# saveRDS(list(categories = categories, thresholds = thr,
#              TE = te_res, SDASM = sdasm_res, GO = go_res, PCDH = pcdh_res),
#         here(paste0("gitignore/S05_allResults_", variant, ".rds")))