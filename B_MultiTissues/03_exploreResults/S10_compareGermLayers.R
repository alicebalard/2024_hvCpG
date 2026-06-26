#####################################################################
# Step 1 - Find inter-individual hypervariable CpGs per germ layer
# Step 2 — Test systemic intra-individual correlations
#####################################################################

## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}
####################################################################

#######################################
## Step 1. Indentify candidate sites ##
#######################################

# # ── 1. Load ───────────────────────────────────────────────────────────────────
# endo     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
# meso     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
# ecto     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
# analyses <- list(endo = endo, meso = meso, ecto = ecto)
# # ── 2. Wide table ─────────────────────────────────────────────────────────────
# wide <- Reduce(
#   function(a, b) merge(a, b, by = "name"),
#   Map(function(dt, nm) setnames(copy(dt), "alpha", nm), analyses, names(analyses))
# )
# saveRDS(wide, here("gitignore/wide_script10_3layers_full.RDS"))
# rm(endo, meso, ecto, analyses)
# 
# # Same with 6 groups in each category (power test)
# endo     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_2_endo6gp.rds"))
# meso     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_2_meso6gp.rds"))
# ecto     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
# analyses <- list(endo = endo, meso = meso, ecto = ecto)
# # ── 2. Wide table ─────────────────────────────────────────────────────────────
# wide <- Reduce(
#   function(a, b) merge(a, b, by = "name"),
#   Map(function(dt, nm) setnames(copy(dt), "alpha", nm), analyses, names(analyses))
# )
# saveRDS(wide, here("gitignore/wide_script10_3layers_6gpall.RDS"))
# rm(endo, meso, ecto, analyses)

wideFull <- readRDS(here("gitignore/wide_script10_3layers_full.RDS"))
wide6gpall <- readRDS(here("gitignore/wide_script10_3layers_6gpall.RDS"))

## Define a cutoff (based on decay curve S06)
HVt <- 0.7
notHVt <- 0.2

plotquadrantLayer <- function(wide, HV = HVt, notHV = notHVt){
  # ── 1. Categories based on 3 "one vs two others" comparisons ─────────────────
  wide[, `:=`(
    HV_meso     = meso > HV,
    HV_endo     = endo > HV,
    HV_ecto     = ecto > HV,
    notHV_meso  = meso < notHV,
    notHV_endo  = endo < notHV,
    notHV_ecto  = ecto < notHV
  )]
  
  wide[, category := fcase(
    # ME: HV in all 3
    HV_meso & HV_endo & HV_ecto,
    "ME",
    
    # layer-specific: HV in one, clearly NOT HV in the other two
    HV_meso & notHV_endo & notHV_ecto, "Meso_specific",
    HV_endo & notHV_meso & notHV_ecto, "Endo_specific",
    HV_ecto & notHV_meso & notHV_endo, "Ecto_specific",
    
    # clearly not HV anywhere → constitutive
    notHV_meso & notHV_endo & notHV_ecto, "constitutive",
    
    # anything else: one layer HV but others not clearly low, or intermediate
    default = "ambiguous"
  )]
  
  print(table(wide$category))
  
  # ── 2. Three plots: each layer vs its "two others" pool ───────────────────────
  category_colours <- c(
    ME             = "#E69F00",
    Meso_specific  = "#56B4E9",
    Endo_specific  = "#009E73",
    Ecto_specific  = "#CC79A7",
    ambiguous      = "grey90", 
    constitutive   = "black"
  )
  
  set.seed(1234)
  
  make_plot <- function(data, x_col, y_col, x_lab, y_lab, title) {
    data <- data[sample(nrow(data), 100000), ]
    
    ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], 
                     colour = category, shape = category)) +
      geom_point(data = data[category == "constitutive"],
                 alpha = 0.3, size = 0.3) +
      geom_point(data = data[category == "ambiguous"],
                 alpha = 0.4, size = 0.5) +
      geom_point(data = data[!category %in% c("constitutive", "ambiguous")],
                 alpha = 0.8, size = 1.5) +
      geom_hline(yintercept = HV,    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_vline(xintercept = HV,    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_hline(yintercept = notHV, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
      geom_vline(xintercept = notHV, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
      scale_x_continuous(limits = c(0, 1), name = x_lab) +
      scale_y_continuous(limits = c(0, 1), name = y_lab) +
      scale_colour_manual(values = category_colours, drop = FALSE) +
      scale_shape_manual(values = c(
        ME             = 16,   # filled circle
        Meso_specific  = 16,
        Endo_specific  = 16,
        Ecto_specific  = 16,
        ambiguous      = 1,    # empty circle
        constitutive   = 4     # cross
      ), drop = FALSE) +
      ggtitle(title) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
  }
  
  plotquadrantLayer
  
  p_mesovsendo <- make_plot(wide,
                            x_col = "meso",     y_col = "endo",
                            x_lab = "Pr(HV) meso",
                            y_lab = "Pr(HV) endo",
                            title = "Meso vs endo")
  
  p_mesovsecto <- make_plot(wide,
                            x_col = "meso",     y_col = "ecto",
                            x_lab = "Pr(HV) meso",
                            y_lab = "Pr(HV) ecto",
                            title = "Meso vs ecto")
  
  p_endovsecto <- make_plot(wide,
                            x_col = "endo",     y_col = "ecto",
                            x_lab = "Pr(HV) endo",
                            y_lab = "Pr(HV) ecto",
                            title = "Endo vs ecto")
  
  # shared legend
  legend_p <- ggplot(wide[sample(.N, 10000)],
                     aes(x = meso, y = endo, colour = category)) +
    geom_point(size = 3) +
    scale_colour_manual(values = category_colours, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_void() +
    theme(legend.position = "right")
  
  (p_mesovsendo | p_mesovsecto | p_endovsecto | cowplot::get_legend(legend_p)) +
    plot_layout(widths = c(1, 1, 1, 0.35))
  
  return(wide)
}

wideFull <- plotquadrantLayer(wideFull)
# ambiguous  constitutive Ecto_specific Endo_specific            ME Meso_specific 
# 10262820      10920663         64636          3109        262202          9111

wide6gpall <- plotquadrantLayer(wide6gpall)
# ambiguous  constitutive Ecto_specific Endo_specific            ME Meso_specific 
# 10954178      10924267         82554         66908        268520         31034 

####################################
## Step 2. Test of SIV of targets ##
####################################

## Test of raw data between people
Loyfer <- read.csv("../dataIn/SupTab1_Loyfer2023.csv")

dupPeople <- Loyfer[Loyfer$PatientID %in% Loyfer$PatientID[duplicated(Loyfer$PatientID)],]

# from the table above, build a summary of patients by their germ layer combination
layer_summary <- as.data.table(dupPeople)[
  , .(
    n_tissues   = .N,
    germ_layers = paste(sort(unique(Germ.layer)), collapse = "+")
  ),
  by = PatientID
]

# cross-tab: how many patients have each combination
result <- layer_summary[, .N, by = .(germ_layers, n_tissues)]
setorder(result, germ_layers, n_tissues)
result
#     germ_layers n_tissues    N
# 1:        Ecto         2     5
# 2:        Ecto         3     1
# 3:        Ecto         4     1
# 4:        Endo         2     9
# 5:        Endo         3     5
# 6:   Endo+Meso         2     1
# 7:   Endo+Meso         3     1
# 8:   Endo+Meso         4     1
# 9:   Endo+Meso         5     1
# 10:        Meso         2     3
# 11:        Meso         3     1
# 12:        Meso         5     1
# 13:        Meso         6     2
# 14:        Meso         7     2

#################################################################################################
# Test if the meso-specific candidate are correlated within mesodermal cells within one patient,
# but not correlated between endo and mesodermal cells within one patient
#################################################################################################

layer_specific <- wideFull[
  wideFull$category %in% c("Ecto_specific", "Endo_specific", "Meso_specific"),"name"]
# ── Write the cpg_list file for S00 ──────────────────────────────────────────
out_path <- here("B_MultiTissues/dataOut/layer_specific_full.txt")
writeLines(layer_specific$name, out_path)
message("Written to: ", out_path)

layer_specific <- wide6gpall[
  wide6gpall$category %in% c("Ecto_specific", "Endo_specific", "Meso_specific"),"name"]
# ── Write the cpg_list file for S00 ──────────────────────────────────────────
out_path <- here("B_MultiTissues/dataOut/layer_specific_6gpall.txt")
writeLines(layer_specific$name, out_path)
message("Written to: ", out_path)

# ── Overlap between full and 6gp Ecto_specific ───────────────────────────────
# fast version using data.table keys
setDT(wideFull);    setkey(wideFull,    name)
setDT(wide6gpall);  setkey(wide6gpall,  name)

categories <- c("ME", "Ecto_specific", "Endo_specific", "Meso_specific",
                "constitutive", "ambiguous")
overlap_summary <- rbindlist(lapply(categories, function(cat) {
  f <- wideFull[category    == cat, name]
  g <- wide6gpall[category  == cat, name]
  o <- length(intersect(f, g))
  data.table(category    = cat,
             n_full      = length(f),
             n_6gp       = length(g),
             n_overlap   = o,
             pct_of_full = round(100 * o / length(f), 1),
             pct_of_6gp  = round(100 * o / length(g), 1))
}))

print(overlap_summary)
#         category  n_full    n_6gp n_overlap  pct_of_full pct_of_6gp
# <char>    <int>    <int>     <int>       <num>      <num>
# 1:            ME   262202   268520    215170        82.1       80.1
# 2: Ecto_specific    64636    82554     41462        64.1       50.2
# 3: Endo_specific     3109    66908      1648        53.0        2.5
# 4: Meso_specific     9111    31034      2878        31.6        9.3
# 5:  constitutive 10920663 10924267   9773987        89.5       89.5
# 6:     ambiguous 10262820 10954178   9279596        90.4       84.7

# 80% of full-atlas MEs are recovered by the 6gp analysis and vice versa.
# 64% of full-atlas Ecto_specific are recovered in 6gp, but only 50% of 6gp Ecto_specific match back.
# Only 53% of full-atlas Endo_specific are recovered, and only 2.5% of 6gp Endo_specific match back. 
# Only 32% of full-atlas Meso_specific recovered, and 9% agreement from 6gp side
# Overall conclusion: the 6gp analysis is adequate for MEs and constitutive CpGs
# but unreliable for layer-specific classifications. The layer-specific categories 
# are highly sensitive to the number of tissues included: with fewer tissues,
# random noise gets misclassified as layer-specificity.

## In pchuckle interactive session:
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/layer_specific.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_layerspecificCpGs.tsv \
# --minCov    10

meth <- fread(here("gitignore/methylation_layerspecificCpGs.tsv"))

# ── Add germ_layer column from metadata ───────────────────────────────────────
loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))

# build lookup: source_tissue_celltype → germ_layer
# the python script creates source_tissue_celltype as "Source Tissue - Cell type"
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]

tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])

# join onto meth
meth <- merge(meth, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)

# check
message(sprintf("Rows with missing germ_layer: %d", sum(is.na(meth$germ_layer))))
print(meth[, .N, by = germ_layer])

# ── Add category to meth ──────────────────────────────────────────────────────
cpg_category <- wide[, .(cpg_site = name, category)]
meth <- merge(meth, cpg_category, by = "cpg_site", all.x = TRUE)


# For each CpG category (e.g. Meso_specific), I want to classify each CpG into:
#   
# Layer-consistent — high r only within its own layer's tissue pairs → confirms post-commitment establishment, true germ-layer DMR
# Systemic — high r across ALL tissue pairs including cross-layer → suggests pre-gastrulation ME that wasn't picked up as ME because the variance is shared, not stochastic
# Other — low r everywhere, or inconsistent pattern → experience-dependent, noisy, or technical

# ══════════════════════════════════════════════════════════════════════════════
# Per-CpG classification: layer-consistent vs systemic vs other
# ══════════════════════════════════════════════════════════════════════════════

# ── Step 1: per-CpG, per-tissue-pair Pearson r across patients ───────────────
# (patients are the "replicates" — correlate across patients not across CpGs)

compute_percpg_pairtissue_r <- function(meth_sub, min_patients = 3) {
  
  # mean per patient per tissue per CpG
  agg <- meth_sub[, .(methylation = mean(methylation, na.rm = TRUE)),
                  by = .(cpg_site, patient_id, source_tissue_celltype,
                         germ_layer)]
  
  tissues <- unique(agg$source_tissue_celltype)
  if (length(tissues) < 2) return(NULL)
  
  tissue_pairs <- combn(tissues, 2, simplify = FALSE)
  
  rbindlist(lapply(tissue_pairs, function(pr) {
    t1 <- pr[1]; t2 <- pr[2]
    
    # patients with both tissues
    has_t1 <- agg[source_tissue_celltype == t1, unique(patient_id)]
    has_t2 <- agg[source_tissue_celltype == t2, unique(patient_id)]
    common  <- intersect(has_t1, has_t2)
    if (length(common) < min_patients) return(NULL)
    
    d1 <- agg[source_tissue_celltype == t1 & patient_id %in% common,
              .(cpg_site, patient_id, m1 = methylation)]
    d2 <- agg[source_tissue_celltype == t2 & patient_id %in% common,
              .(cpg_site, patient_id, m2 = methylation)]
    
    dm <- merge(d1, d2, by = c("cpg_site", "patient_id"))
    
    # layer of each tissue
    l1 <- unique(agg[source_tissue_celltype == t1, germ_layer])
    l2 <- unique(agg[source_tissue_celltype == t2, germ_layer])
    same_layer <- length(l1) == 1 && length(l2) == 1 && l1 == l2
    pair_type  <- if (same_layer) paste0(l1[1], " × ", l1[1]) else
      paste(sort(c(l1[1], l2[1])), collapse = " × ")
    
    # per-CpG r across patients
    dm[, {
      idx <- !is.na(m1) & !is.na(m2)
      if (sum(idx) >= min_patients) {
        r <- suppressWarnings(cor(m1[idx], m2[idx], method = "pearson"))
        list(r = r, n_patients = sum(idx),
             tissue1 = t1, tissue2 = t2,
             same_layer = same_layer,
             pair_type  = pair_type)
      } else NULL
    }, by = cpg_site]
    
  }), fill = TRUE)
}

# ── Step 2: run per category ──────────────────────────────────────────────────
percpg_results <- rbindlist(lapply(categories, function(cat) {
  message("Computing per-CpG tissue-pair r for: ", cat)
  sub <- meth_multi_all[category == cat]
  res <- compute_percpg_pairtissue_r(sub, min_patients = 3)
  if (is.null(res)) return(NULL)
  res[, category := cat]
  res
}), fill = TRUE)

# ── Step 3: per-CpG summary — mean r within own layer vs cross layer ──────────
# define "own layer" pEndo_specificer category
own_layer <- c(
  Meso_specific = "Meso",
  Endo_specific = "Endo",
  Ecto_specific = "Ecto"
)

percpg_summary <- percpg_results[!is.na(r), {
  ol     <- own_layer[category]
  is_own <- pair_type == paste0(ol, " × ", ol)
  list(
    mean_r_own_layer   = mean(r[is_own],  na.rm = TRUE),
    mean_r_cross_layer = mean(r[!is_own], na.rm = TRUE),
    n_own_pairs        = sum(is_own),
    n_cross_pairs      = sum(!is_own)
  )
}, by = .(cpg_site, category)]

# ── Step 4: classify each CpG ────────────────────────────────────────────────
r_high <- 0.5   # threshold for "high correlation"

percpg_summary[, cpg_class := fcase(
  # high within own layer AND high cross-layer → systemic (ME-like)
  mean_r_own_layer   >= r_high & mean_r_cross_layer >= r_high, "Systemic (ME-like)",
  # high within own layer, low cross-layer → layer-consistent (post-commitment)
  mean_r_own_layer   >= r_high & mean_r_cross_layer <  r_high, "Layer-consistent",
  # low own layer, high cross-layer → odd, possibly constitutive with noise
  mean_r_own_layer   <  r_high & mean_r_cross_layer >= r_high, "Cross-layer only",
  # low everywhere
  default = "Low r everywhere"
)]

# ── Step 5: count per category × class ───────────────────────────────────────
class_counts <- percpg_summary[, .N, by = .(category, cpg_class)]
setorder(class_counts, category, cpg_class)
print(class_counts)

# category          cpg_class     N
# <char>             <char> <int>
#   1: Ecto_specific   Cross-layer only   703
# 2: Ecto_specific   Layer-consistent 17973
# 3: Ecto_specific   Low r everywhere 29890
# 4: Ecto_specific Systemic (ME-like)   465
# 5: Endo_specific   Cross-layer only    16
# 6: Endo_specific   Layer-consistent   386
# 7: Endo_specific   Low r everywhere   287
# 8: Endo_specific Systemic (ME-like)    23
# 9: Meso_specific   Cross-layer only    28
# 10: Meso_specific   Layer-consistent  1453
# 11: Meso_specific   Low r everywhere   285
# 12: Meso_specific Systemic (ME-like)   331

# ── Step 6: scatter plot — own-layer r vs cross-layer r, coloured by class ───
class_colours <- c(
  "Systemic (ME-like)"  = "#E69F00",
  "Layer-consistent"    = "#0072B2",
  "Cross-layer only"    = "#D55E00",
  "Low r everywhere"    = "grey70"
)

p_classify <- ggplot(percpg_summary,
                     aes(x = mean_r_cross_layer,
                         y = mean_r_own_layer,
                         colour = cpg_class)) +
  geom_point(size = 0.8, alpha = 0.5, pch = 20) +
  geom_hline(yintercept = r_high, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = r_high, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = class_colours, name = "CpG class") +
  scale_x_continuous("Mean r — cross-layer tissue pairs", limits = c(-1, 1)) +
  scale_y_continuous("Mean r — same-layer tissue pairs",  limits = c(-1, 1)) +
  facet_wrap(~ category, nrow = 1) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right") +
  ggtitle("Per-CpG classification by intra-individual correlation structure",
          subtitle = sprintf("r threshold = %.1f | Each point = one CpG", r_high))

# ── Step 7: bar chart of proportions ─────────────────────────────────────────
class_props <- percpg_summary[, .(N = .N), by = .(category, cpg_class)]
class_props[, prop := N / sum(N), by = category]

p_bar <- ggplot(class_props,
                aes(x = category, y = prop, fill = cpg_class)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = class_colours, name = "CpG class") +
  scale_y_continuous("Proportion of CpGs", labels = scales::percent) +
  theme_bw(base_size = 10) +
  theme(axis.text.x    = element_text(angle = 30, hjust = 1),
        legend.position = "right") +
  ggtitle("Proportion of CpGs per class within each category")

p_classify / p_bar + plot_layout(heights = c(2, 1))





# # ####################
# # ## GO enrichement ##
# # ####################
# # 
