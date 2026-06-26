#####################################################################
# Step 1 - Find inter-individual hypervariable CpGs per germ layer
# Step 2 — Test systemic intra-individual correlations
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1. Identify candidate sites
# ══════════════════════════════════════════════════════════════════════════════
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
wideFull   <- readRDS(here("gitignore/wide_script10_3layers_full.RDS"))
wide6gpall <- readRDS(here("gitignore/wide_script10_3layers_6gpall.RDS"))

HVt    <- 0.7
notHVt <- 0.2

category_colours <- c(
  ME            = "#E69F00",
  Meso_specific = "#56B4E9",
  Endo_specific = "#009E73",
  Ecto_specific = "#CC79A7",
  ambiguous     = "grey90",
  constitutive  = "black"
)

# ── Classify CpGs into categories ────────────────────────────────────────────
classify_wide <- function(wide, HV = HVt, notHV = notHVt) {
  wide[, `:=`(
    HV_meso    = meso > HV,    HV_endo    = endo > HV,    HV_ecto    = ecto > HV,
    notHV_meso = meso < notHV, notHV_endo = endo < notHV, notHV_ecto = ecto < notHV
  )]
  wide[, category := fcase(
    HV_meso & HV_endo & HV_ecto,         "ME",
    HV_meso & notHV_endo & notHV_ecto,   "Meso_specific",
    HV_endo & notHV_meso & notHV_ecto,   "Endo_specific",
    HV_ecto & notHV_meso & notHV_endo,   "Ecto_specific",
    notHV_meso & notHV_endo & notHV_ecto,"constitutive",
    default =                             "ambiguous"
  )]
  wide
}

# ── Scatter plot: layer vs layer coloured by category ────────────────────────
plot_quadrant_layer <- function(wide, HV = HVt, notHV = notHVt) {
  make_plot <- function(x_col, y_col, x_lab, y_lab, title) {
    d <- wide[sample(.N, 100000)]
    ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]],
                  colour = category, shape = category)) +
      geom_point(data = d[category == "constitutive"], alpha = 0.3, size = 0.3) +
      geom_point(data = d[category == "ambiguous"],    alpha = 0.4, size = 0.5) +
      geom_point(data = d[!category %in% c("constitutive","ambiguous")],
                 alpha = 0.8, size = 1.5) +
      geom_hline(yintercept = c(HV, notHV),
                 linetype = c("dashed","dotted"), colour = c("grey40","grey60"),
                 linewidth = c(0.4, 0.3)) +
      geom_vline(xintercept = c(HV, notHV),
                 linetype = c("dashed","dotted"), colour = c("grey40","grey60"),
                 linewidth = c(0.4, 0.3)) +
      scale_x_continuous(limits = c(0,1), name = x_lab) +
      scale_y_continuous(limits = c(0,1), name = y_lab) +
      scale_colour_manual(values = category_colours, drop = FALSE) +
      scale_shape_manual(values = c(ME=16, Meso_specific=16, Endo_specific=16,
                                    Ecto_specific=16, ambiguous=1, constitutive=4),
                         drop = FALSE) +
      ggtitle(title) + theme_bw(base_size = 11) +
      theme(legend.position = "none")
  }
  legend_p <- ggplot(wide[sample(.N, 10000)], aes(x=meso, y=endo, colour=category)) +
    geom_point(size = 3) +
    scale_colour_manual(values = category_colours, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
    theme_void() + theme(legend.position = "right")
  
  (make_plot("meso","endo","Pr(HV) meso","Pr(HV) endo","Meso vs Endo") |
      make_plot("meso","ecto","Pr(HV) meso","Pr(HV) ecto","Meso vs Ecto") |
      make_plot("endo","ecto","Pr(HV) endo","Pr(HV) ecto","Endo vs Ecto") |
      cowplot::get_legend(legend_p)) +
    plot_layout(widths = c(1,1,1,0.35))
}

wideFull   <- classify_wide(wideFull);   print(table(wideFull$category))
wide6gpall <- classify_wide(wide6gpall); print(table(wide6gpall$category))

set.seed(1234)
plot_quadrant_layer(wideFull)
plot_quadrant_layer(wide6gpall)

# ── Overlap full vs 6gp ───────────────────────────────────────────────────────
setkey(wideFull,    name)
setkey(wide6gpall,  name)

overlap_summary <- rbindlist(lapply(
  c("ME","Ecto_specific","Endo_specific","Meso_specific","constitutive","ambiguous"),
  function(cat) {
    f <- wideFull[category   == cat, name]
    g <- wide6gpall[category == cat, name]
    o <- length(intersect(f, g))
    data.table(category    = cat,
               n_full      = length(f),
               n_6gp       = length(g),
               n_overlap   = o,
               pct_of_full = round(100*o/length(f), 1),
               pct_of_6gp  = round(100*o/length(g), 1))
  }))
print(overlap_summary)
# category   n_full    n_6gp n_overlap pct_of_full pct_of_6gp
# <char>    <int>    <int>     <int>       <num>      <num>
#   1:            ME   262202   268520    215170        82.1       80.1
# 2: Ecto_specific    64636    82554     41462        64.1       50.2
# 3: Endo_specific     3109    66908      1648        53.0        2.5
# 4: Meso_specific     9111    31034      2878        31.6        9.3
# 5:  constitutive 10920663 10924267   9773987        89.5       89.5
# 6:     ambiguous 10262820 10954178   9279596        90.4       84.7

# ── Write CpG lists for python extraction ─────────────────────────────────────
# Layer-specific CpGs (from full atlas)
writeLines(
  wideFull[category %in% c("Ecto_specific","Endo_specific","Meso_specific"), name],
  here("B_MultiTissues/dataOut/layer_specific_full.txt"))
message("Written: layer_specific_full.txt")

# ME CpGs — subsample 5000 as positive control
set.seed(42)
me_sample <- sample(wideFull[category == "ME", name], 5000)

# Combined list for single extraction run
all_cpgs_to_extract <- unique(c(
  wideFull[category %in% c("Ecto_specific","Endo_specific","Meso_specific"), name],
  me_sample
))
writeLines(all_cpgs_to_extract,
           here("B_MultiTissues/dataOut/layer_specific_and_ME.txt"))
message(sprintf("Written: layer_specific_and_ME.txt (%d CpGs total)",
                length(all_cpgs_to_extract)))
# Written: layer_specific_and_ME.txt (81856 CpGs total)

## ── Run on pchuckle (single extraction for all CpGs) ─────────────────────────
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/layer_specific_and_ME.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_layerspecific_and_ME.tsv \
# --minCov    10

# ══════════════════════════════════════════════════════════════════════════════
# Step 2. Intra-individual correlation for layer-specific CpGs + ME
# ══════════════════════════════════════════════════════════════════════════════

meth <- fread(here("gitignore/methylation_layerspecific_and_ME.tsv"))

# ── Join germ_layer and category ─────────────────────────────────────────────
loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
tissue_to_layer <- unique(loyfer_meta[,
                                      .(source_tissue_celltype, germ_layer = `Germ layer`)])

meth <- merge(meth, tissue_to_layer,
              by = "source_tissue_celltype", all.x = TRUE)
meth <- merge(meth,
              wideFull[, .(cpg_site = name, category)],
              by = "cpg_site", all.x = TRUE)

# ME CpGs that were not in wideFull get NA category — fix
meth[cpg_site %in% me_sample & is.na(category), category := "ME"]

message(sprintf("Missing germ_layer: %d | Missing category: %d",
                sum(is.na(meth$germ_layer)), sum(is.na(meth$category))))
print(meth[, .N, by = category])
# category        N
# <char>    <int>
#   1: Ecto_specific 13250380
# 2:            ME  1025000
# 3: Meso_specific  1867755
# 4: Endo_specific   637345

# ── Multi-tissue patients ─────────────────────────────────────────────────────
multi_patients <- meth[, .(n = uniqueN(source_tissue_celltype)),
                       by = patient_id][n > 1, patient_id]
message(sprintf("%d patients with >1 tissue", length(multi_patients)))
meth_multi <- meth[patient_id %in% multi_patients]

# ── Descriptive table of multi-tissue patients ────────────────────────────────
patient_table <- meth_multi[, .(
  n_tissues   = uniqueN(source_tissue_celltype),
  germ_layers = paste(sort(unique(germ_layer)), collapse = "+"),
  n_blood     = uniqueN(source_tissue_celltype[grepl(
    "Blood", source_tissue_celltype, ignore.case=TRUE)]),
  n_nonblood  = uniqueN(source_tissue_celltype[!grepl(
    "Blood", source_tissue_celltype, ignore.case=TRUE)])
), by = patient_id]

summary_table <- patient_table[, .(
  n_patients     = .N,
  median_tissues = as.numeric(median(n_tissues)),
  range_tissues  = sprintf("%d-%d", min(n_tissues), max(n_tissues)),
  n_blood_only   = sum(n_nonblood == 0),
  pct_blood_only = round(100 * sum(n_nonblood == 0) / .N)
), by = germ_layers][order(germ_layers)]
print(summary_table)
# germ_layers n_patients median_tissues range_tissues n_blood_only pct_blood_only
# <char>      <int>          <num>        <char>        <int>          <num>
#   1:        Ecto          6            2.0           2-2            0              0
# 2:        Endo         13            2.0           2-3            0              0
# 3:   Endo+Meso          4            2.5           2-3            0              0
# 4:        Meso          9            5.0           2-7            6             67

# ══════════════════════════════════════════════════════════════════════════════
# Per-CpG inter-tissue correlation
#
# For each CpG and each layer-pair comparison (e.g. Meso×Meso, Endo×Meso):
#   1. For each tissue pair, correlate methylation across patients who have both
#   2. Average r across all tissue pairs for that layer comparison
#      → one r per CpG per pair_type
#
# Patients with many blood cell types contribute more tissue pairs but each
# tissue pair is weighted equally — no patient dominates the average
# ══════════════════════════════════════════════════════════════════════════════

compute_percpg_layerpair_r <- function(meth_sub, min_patients = 3) {
  
  agg <- meth_sub[, .(methylation = mean(methylation, na.rm = TRUE)),
                  by = .(cpg_site, patient_id, source_tissue_celltype, germ_layer)]
  
  layers <- unique(agg$germ_layer)
  if (length(layers) < 1) return(NULL)
  
  layer_pairs <- CJ(l1 = layers, l2 = layers)[l1 <= l2]
  
  rbindlist(lapply(seq_len(nrow(layer_pairs)), function(k) {
    l1 <- layer_pairs$l1[k]; l2 <- layer_pairs$l2[k]
    same_layer <- l1 == l2
    pair_type  <- if (same_layer) paste0(l1, " \u00d7 ", l1) else
      paste(sort(c(l1, l2)), collapse = " \u00d7 ")
    
    t1_tissues <- unique(agg[germ_layer == l1, source_tissue_celltype])
    t2_tissues <- unique(agg[germ_layer == l2, source_tissue_celltype])
    
    if (same_layer) {
      if (length(t1_tissues) < 2) return(NULL)
      tpairs <- combn(t1_tissues, 2, simplify = FALSE)
    } else {
      tpairs <- unlist(lapply(t1_tissues, function(t1)
        lapply(t2_tissues, function(t2) c(t1, t2))),
        recursive = FALSE)
    }
    if (length(tpairs) == 0) return(NULL)
    
    pair_r_list <- rbindlist(lapply(tpairs, function(pr) {
      t1 <- pr[1]; t2 <- pr[2]
      common <- intersect(
        agg[source_tissue_celltype == t1, unique(patient_id)],
        agg[source_tissue_celltype == t2, unique(patient_id)]
      )
      if (length(common) < min_patients) return(NULL)
      
      d1 <- agg[source_tissue_celltype == t1 & patient_id %in% common,
                .(cpg_site, patient_id, m1 = methylation)]
      d2 <- agg[source_tissue_celltype == t2 & patient_id %in% common,
                .(cpg_site, patient_id, m2 = methylation)]
      dm <- merge(d1, d2, by = c("cpg_site","patient_id"))
      
      dm[, {
        idx <- !is.na(m1) & !is.na(m2)
        if (sum(idx) >= min_patients)
          list(r_pair     = suppressWarnings(cor(m1[idx], m2[idx], method="pearson")),
               n_patients = sum(idx),
               tissue1    = t1,
               tissue2    = t2)
        else NULL
      }, by = cpg_site]
    }), fill = TRUE)
    
    if (is.null(pair_r_list) || nrow(pair_r_list) == 0) return(NULL)
    
    pair_r_list[!is.na(r_pair), .(
      r          = mean(r_pair, na.rm = TRUE),
      n_pairs    = .N,
      n_patients = max(n_patients),
      same_layer = same_layer,
      pair_type  = pair_type
    ), by = cpg_site]
    
  }), fill = TRUE)
}

# ── Run for all categories including ME ──────────────────────────────────────
categories <- c("ME", "Ecto_specific", "Endo_specific", "Meso_specific")
own_layer  <- c(ME = "all",
                Meso_specific = "Meso",
                Endo_specific = "Endo",
                Ecto_specific = "Ecto")

percpg_results <- rbindlist(lapply(categories, function(cat) {
  message("Computing for: ", cat)
  sub <- meth_multi[category == cat]
  message(sprintf("  %d rows, %d CpGs, %d patients, layers: %s",
                  nrow(sub), uniqueN(sub$cpg_site), uniqueN(sub$patient_id),
                  paste(unique(sub$germ_layer), collapse = "+")))
  if (nrow(sub) == 0 || uniqueN(sub$germ_layer) < 1) {
    message("  Skipping — no data"); return(NULL)
  }
  res <- compute_percpg_layerpair_r(sub, min_patients = 3)
  if (is.null(res) || nrow(res) == 0) {
    message("  Skipping — no valid pairs"); return(NULL)
  }
  res[, category := cat]
}), fill = TRUE)

# Computing for: ME
# 500000 rows, 5000 CpGs, 32 patients, layers: Meso+Ecto+Endo
# Computing for: Ecto_specific
# 6463600 rows, 64636 CpGs, 32 patients, layers: Meso+Ecto+Endo
# Computing for: Endo_specific
# 310900 rows, 3109 CpGs, 32 patients, layers: Meso+Ecto+Endo
# Computing for: Meso_specific
# 911100 rows, 9111 CpGs, 32 patients, layers: Meso+Ecto+Endo

message("Available pair types:")
print(percpg_results[, .N, by = .(pair_type, category)][order(pair_type, category)])
# pair_type      category     N
# <char>        <char> <int>
#   1: Ecto × Ecto Ecto_specific 64417
# 2: Ecto × Ecto Endo_specific  3079
# 3: Ecto × Ecto            ME  4972
# 4: Ecto × Ecto Meso_specific  9018
# 5: Endo × Endo Ecto_specific 64636
# 6: Endo × Endo Endo_specific  3109
# 7: Endo × Endo            ME  5000
# 8: Endo × Endo Meso_specific  9110
# 9: Meso × Meso Ecto_specific 64636
# 10: Meso × Meso Endo_specific  3109
# 11: Meso × Meso            ME  5000
# 12: Meso × Meso Meso_specific  9111

# ── Per-CpG summary: mean r within own layer vs cross-layer ──────────────────
r_high <- 0.5

percpg_summary <- percpg_results[!is.na(r), {
  ol     <- own_layer[category]
  
  if (ol == "all") {
    # ME: no cross-layer distinction — report mean r across all pair types
    list(
      mean_r_own_layer   = mean(r, na.rm = TRUE),
      mean_r_cross_layer = NA_real_,
      n_own_pairs        = .N,
      n_cross_pairs      = 0L
    )
  } else {
    is_own <- pair_type == paste0(ol, " \u00d7 ", ol)
    list(
      mean_r_own_layer   = mean(r[is_own],  na.rm = TRUE),
      mean_r_cross_layer = mean(r[!is_own], na.rm = TRUE),
      n_own_pairs        = sum(is_own),
      n_cross_pairs      = sum(!is_own)
    )
  }
}, by = .(cpg_site, category)]

# ME classification: only own-layer r matters
percpg_summary[, cpg_class := fcase(
  category == "ME" & mean_r_own_layer >= r_high,  "Systemic (ME-like)",
  category == "ME" & mean_r_own_layer <  r_high,  "Low r everywhere",
  mean_r_own_layer >= r_high & mean_r_cross_layer >= r_high, "Systemic (ME-like)",
  mean_r_own_layer >= r_high & mean_r_cross_layer <  r_high, "Layer-consistent",
  mean_r_own_layer <  r_high & mean_r_cross_layer >= r_high, "Cross-layer only",
  default = "Low r everywhere"
)]

class_counts <- percpg_summary[, .N, by = .(category, cpg_class)]
setorder(class_counts, category, cpg_class)
print(class_counts)
# category          cpg_class     N
# <char>             <char> <int>
#   1: Ecto_specific   Cross-layer only   313
# 2: Ecto_specific   Layer-consistent 23024
# 3: Ecto_specific   Low r everywhere 41068
# 4: Ecto_specific Systemic (ME-like)   231
# 5: Endo_specific   Cross-layer only   181
# 6: Endo_specific   Layer-consistent   910
# 7: Endo_specific   Low r everywhere  1925
# 8: Endo_specific Systemic (ME-like)    93
# 9:            ME   Low r everywhere   934
# 10:            ME Systemic (ME-like)  4066
# 11: Meso_specific   Cross-layer only   452
# 12: Meso_specific   Layer-consistent  4141
# 13: Meso_specific   Low r everywhere  3956
# 14: Meso_specific Systemic (ME-like)   562

# ── Colours ───────────────────────────────────────────────────────────────────
class_colours <- c(
  "Systemic (ME-like)" = "#E69F00",
  "Layer-consistent"   = "#0072B2",
  "Cross-layer only"   = "#D55E00",
  "Low r everywhere"   = "grey70"
)

category_order <- c("ME", "Meso_specific", "Endo_specific", "Ecto_specific")

# ── Plot 1: scatter per category ─────────────────────────────────────────────
# ME gets its own plot showing r distribution per pair_type
p_me <- ggplot(
  percpg_results[category == "ME" & !is.na(r)],
  aes(x = r, fill = pair_type, colour = pair_type)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = r_high, linetype = "dashed", colour = "grey40") +
  scale_x_continuous("Mean r per pair type", limits = c(-1, 1)) +
  scale_fill_manual(values = c(
    "Ecto \u00d7 Ecto" = "#CC79A7",
    "Endo \u00d7 Endo" = "#009E73",
    "Meso \u00d7 Meso" = "#56B4E9"), name = "Layer pair") +
  scale_colour_manual(values = c(
    "Ecto \u00d7 Ecto" = "#CC79A7",
    "Endo \u00d7 Endo" = "#009E73",
    "Meso \u00d7 Meso" = "#56B4E9"), name = "Layer pair") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  ggtitle("ME CpGs — r distribution per layer pair",
          subtitle = "Expected: high r in all three pairs (systemic establishment)")

# layer-specific scatter as before (excluding ME)
p_classify_ls <- ggplot(
  percpg_summary[category %in% c("Meso_specific","Endo_specific","Ecto_specific")],
  aes(x = mean_r_cross_layer, y = mean_r_own_layer, colour = cpg_class)) +
  geom_point(size = 0.8, alpha = 0.5, pch = 20) +
  geom_hline(yintercept = r_high, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = r_high, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = class_colours, name = "CpG class") +
  scale_x_continuous("Mean r — cross-layer pairs", limits = c(-1, 1)) +
  scale_y_continuous("Mean r — same-layer pairs",  limits = c(-1, 1)) +
  facet_wrap(~ factor(category,
                      levels = c("Meso_specific","Endo_specific","Ecto_specific")),
             nrow = 1) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  ggtitle("Layer-specific CpGs — own vs cross-layer r")

# ── Plot 2: proportions ───────────────────────────────────────────────────────
class_props <- percpg_summary[category %in% category_order,
                              .(N = .N), by = .(category, cpg_class)]
class_props[, prop := N / sum(N), by = category]

p_bar <- ggplot(class_props,
                aes(x = factor(category, levels = category_order),
                    y = prop, fill = cpg_class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(prop*100,1), "%")),
            position = position_stack(vjust = 0.5), size = 3, colour = "white") +
  scale_fill_manual(values = class_colours, name = "CpG class") +
  scale_y_continuous("Proportion of CpGs", labels = scales::percent) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(), legend.position = "right") +
  ggtitle("Proportion of CpGs per class within each category")

p_me / p_classify_ls / p_bar + plot_layout(heights = c(1, 1,1))

# ── Plot 3: violin — r distribution per pair_type per category ───────────────
pair_type_cols <- c(
  "Ecto \u00d7 Ecto" = "#CC79A7",
  "Endo \u00d7 Endo" = "#009E73",
  "Meso \u00d7 Meso" = "#56B4E9",
  "Endo \u00d7 Meso" = "#999999"
)

p_violin <- ggplot(
  percpg_results[!is.na(r) &
                   category %in% category_order &
                   pair_type %in% names(pair_type_cols)],
  aes(x = factor(category, levels = category_order),
      y = r, fill = pair_type)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.08, outlier.size = 0.2, alpha = 0.8,
               position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = pair_type_cols, name = "Layer pair") +
  scale_y_continuous("Pearson r (per CpG per tissue pair)", limits = c(-1, 1)) +
  scale_x_discrete("CpG category") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Within- and cross-layer methylation correlation by CpG category",
          subtitle = paste0(
            "ME = positive control (high r in all pairs expected)\n",
            "Layer-specific = high r in own-layer pairs only\n",
            "Endo\u00d7Meso based on n=4 patients — interpret cautiously"))

p_violin

# ── Plot 4: blood concordance check for Meso_specific ────────────────────────
blood_r_dist <- percpg_results[
  category == "Meso_specific" & pair_type == "Meso \u00d7 Meso" & !is.na(r),
  .(mean_r = mean(r, na.rm=TRUE)), by = cpg_site]

background_r <- percpg_results[
  category != "Meso_specific" & pair_type == "Meso \u00d7 Meso" & !is.na(r),
  .(mean_r = mean(r, na.rm=TRUE)), by = cpg_site]

message(sprintf("Meso_specific Meso\u00d7Meso: mean r = %.3f (n=%d CpGs)",
                mean(blood_r_dist$mean_r), nrow(blood_r_dist)))
message(sprintf("Background    Meso\u00d7Meso: mean r = %.3f (n=%d CpGs)",
                mean(background_r$mean_r, na.rm=TRUE), nrow(background_r)))

plot_dt <- rbind(
  blood_r_dist[, .(mean_r, group = "Meso_specific")],
  background_r[sample(.N, min(.N, nrow(blood_r_dist))),
               .(mean_r, group = "Background")]
)

p_blood <- ggplot(plot_dt, aes(x = mean_r, fill = group, colour = group)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_continuous("Mean r across Meso \u00d7 Meso pairs", limits = c(-1, 1)) +
  scale_fill_manual(values   = c("Meso_specific" = "#56B4E9",
                                 "Background"    = "grey50")) +
  scale_colour_manual(values = c("Meso_specific" = "#56B4E9",
                                 "Background"    = "grey50")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  ggtitle("Meso \u00d7 Meso r: Meso_specific vs background CpGs",
          subtitle = paste0(
            "r\u22481 for Meso_specific = pre-haematopoietic stochastic establishment\n",
            "r\u22480 for background = no systematic blood concordance"))

p_blood
