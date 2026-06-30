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

if (!exists("wideFull")){wideFull   <- readRDS(here("gitignore/wide_script10_3layers_full.RDS"))}
if (!exists("wide6gpall")){wide6gpall <- readRDS(here("gitignore/wide_script10_3layers_6gpall.RDS"))}

## Define a cutoff (based on decay curve S06)
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
    HV_meso & HV_endo & HV_ecto,                          "ME",
    HV_meso & notHV_endo & notHV_ecto,                    "Meso_specific",
    HV_endo & notHV_meso & notHV_ecto,                    "Endo_specific",
    HV_ecto & notHV_meso & notHV_endo,                    "Ecto_specific",
    notHV_meso & notHV_endo & notHV_ecto,                 "constitutive",
    default =                                              "ambiguous"
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
      scale_shape_manual(values  = c(ME=16, Meso_specific=16, Endo_specific=16,
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
# ambiguous  constitutive Ecto_specific Endo_specific           ME  Meso_specific 
# 10262820      10920663         64636          3109        262202          9111

wide6gpall <- classify_wide(wide6gpall); print(table(wide6gpall$category))
# ambiguous  constitutive Ecto_specific Endo_specific          ME  Meso_specific 
# 10954178      10924267         82554         66908        268520         31034 

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
    data.table(category = cat, n_full = length(f), n_6gp = length(g),
               n_overlap = o,
               pct_of_full = round(100*o/length(f), 1),
               pct_of_6gp  = round(100*o/length(g), 1))
  }))
print(overlap_summary)
#         category   n_full    n_6gp n_overlap pct_of_full pct_of_6gp
# 1:            ME   262202   268520    215170        82.1       80.1
# 2: Ecto_specific    64636    82554     41462        64.1       50.2
# 3: Endo_specific     3109    66908      1648        53.0        2.5
# 4: Meso_specific     9111    31034      2878        31.6        9.3
# 5:  constitutive 10920663 10924267   9773987        89.5       89.5
# 6:     ambiguous 10262820 10954178   9279596        90.4       84.7

## MEs don't suffer much from size reduction (as shown before), but the layer-specific
## markers do! So Ecto is less reliable in particular

# ── Write CpG lists for python extraction ─────────────────────────────────────

# ME CpGs — subsample 5000 as positive control
set.seed(1234)
me_sample <- sample(wideFull[category == "ME", name], 5000)
all_cpgs_to_extract <- unique(c(
  wideFull[category %in% c("Ecto_specific","Endo_specific","Meso_specific"), name],
  me_sample
))
writeLines(all_cpgs_to_extract,
           here("B_MultiTissues/dataOut/layer_specific_and_ME.txt"))
message(sprintf("Written: layer_specific_and_ME.txt (%d CpGs total)",
                length(all_cpgs_to_extract)))

## In pchuckle:
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
# Step 2. Intra-individual correlation for layer-specific CpGs
# ══════════════════════════════════════════════════════════════════════════════

meth <- fread(here("gitignore/methylation_layerspecific_and_ME.tsv"))

# ── Join germ_layer and category ─────────────────────────────────────────────
loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])

meth <- merge(meth, tissue_to_layer,
              by = "source_tissue_celltype", all.x = TRUE)
meth <- merge(meth,
              wideFull[, .(cpg_site = name, category)],
              by = "cpg_site", all.x = TRUE)

message(sprintf("Missing germ_layer: %d | Missing category: %d",
                sum(is.na(meth$germ_layer)), sum(is.na(meth$category))))

# ── Multi-tissue patients ─────────────────────────────────────────────────────
multi_patients <- meth[, .(n = uniqueN(source_tissue_celltype)),
                       by = patient_id][n > 1, patient_id]
message(sprintf("%d patients with >1 tissue", length(multi_patients)))
meth_multi <- meth[patient_id %in% multi_patients]

# ── Descriptive table of multi-tissue patients ────────────────────────────────
patient_table <- meth_multi[, .(
  n_tissues  = uniqueN(source_tissue_celltype),
  germ_layers = paste(sort(unique(germ_layer)), collapse = "+"),
  n_blood    = uniqueN(source_tissue_celltype[grepl(
    "Blood", source_tissue_celltype, ignore.case=TRUE)]),
  n_nonblood = uniqueN(source_tissue_celltype[!grepl(
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

message(sprintf("Patients with cross-layer samples: %d",
                sum(grepl("\\+", patient_table$germ_layers))))
# Patients with cross-layer samples: 4 (Endo+Meso)

### How many patients have multiple SAME-layer tissues (different tissues)

# patients with >1 tissue WITHIN the same germ layer
same_layer_multi <- meth_multi[, .(n_tissues_in_layer = uniqueN(source_tissue_celltype)),
                               by = .(patient_id, germ_layer)][n_tissues_in_layer > 1]

print(same_layer_multi[order(germ_layer, -n_tissues_in_layer)])

message(sprintf("Patients with >1 tissue within the SAME germ layer: %d",
                uniqueN(same_layer_multi$patient_id)))
# Patients with >1 tissue within the SAME germ layer: 30

# breakdown by layer
print(same_layer_multi[, .(n_patients = uniqueN(patient_id)), by = germ_layer])
#    germ_layer n_patients
# 1:       Meso          9
# 2:       Ecto          6
# 3:       Endo         15

## I have 4 patients with Endo+Meso, 9 with multiple meso, 6 with multiple ecto
## and 15 with multiple endo.

## I want to know which percentage of the 9111 Meso_specific CpGs show a high 
# intra-individual correlation on average between 2 mesodermal tissues in the 9 
# people with multiple meso tissues (if more than 2 meso tissues, do pairwise 
# then mean to give one value), but a low correlation between 2 different
# germ layer in the 4 patients with Endo+Meso tissues.

## Same idea for the 3109 Endo_specific CpGs (in the 15 people with endo,
# and the 4 Meso+endo)

## As a negative control, I'd like the same percentages but for a random 10k
# subset of the 10920663 constitutive CpGs and for the 10262820 ambiguous, in
# both cases

## I can't test Ecto_specific (not enough data)

# ══════════════════════════════════════════════════════════════════════════════
# Within-individual same-layer vs cross-layer correlation per CpG
#
# GOAL (per CpG):
#   1. SAME-LAYER:  for each patient who has >1 tissue from the relevant layer,
#      compute the pairwise Pearson correlation between all their tissues for
#      that CpG, then average those pairwise correlations -> one value per
#      patient. Average across patients -> one "same-layer r" per CpG.
#   2. CROSS-LAYER: for the 4 patients with Endo+Meso, compute Pearson r
#      between their Endo tissue(s) and Meso tissue(s) for that CpG (pairwise
#      then mean if >1 tissue per layer) -> one value per patient. Average
#      across the 4 patients -> one "cross-layer r" per CpG.
#   3. Report: % of CpGs in the category with same-layer r >= threshold AND
#      cross-layer r < threshold (the "true layer-specific" signature)
#
# NOTE: this is computed WITHIN individuals (different tissues, same person),
# not BETWEEN individuals as before.
# ══════════════════════════════════════════════════════════════════════════════

r_threshold <- 0.5   # what counts as "high" correlation

# ── Helper: for ONE patient and a set of tissues, return one r per CpG ───────
# (mean of all pairwise tissue correlations within that patient)
patient_mean_pairwise_r <- function(meth_pat, tissues) {
  if (length(tissues) < 2) return(NULL)
  
  # wide: rows = CpG, cols = tissue, values = methylation
  wide <- dcast(meth_pat[source_tissue_celltype %in% tissues],
                cpg_site ~ source_tissue_celltype, value.var = "methylation",
                fun.aggregate = mean)
  
  tcols <- setdiff(names(wide), "cpg_site")
  if (length(tcols) < 2) return(NULL)
  
  mat <- as.matrix(wide[, ..tcols])
  rownames(mat) <- wide$cpg_site
  
  tpairs <- combn(tcols, 2, simplify = FALSE)
  
  # for each CpG (row), compute pairwise correlation across the FEW tissue
  # values available — since each patient has only a handful of tissue
  # measurements per CpG, "correlation across tissues" here really means:
  # take the methylation values across tissues, and since we only have one
  # measurement per tissue (not multiple replicates), we instead compute
  # agreement using the actual paired values per tissue pair, per CpG
  pair_vals <- rbindlist(lapply(tpairs, function(pr) {
    m1 <- mat[, pr[1]]
    m2 <- mat[, pr[2]]
    data.table(cpg_site = rownames(mat), m1 = m1, m2 = m2)
  }), idcol = "pair_id")
  
  pair_vals
}

# ── Core function: same-layer within-individual r per CpG ────────────────────
# For each patient with >1 tissue in `layer`, get pairwise tissue values,
# then for each CpG compute correlation ACROSS the pairwise differences
# Since within ONE patient there's no "sample size" to correlate against,
# we instead measure AGREEMENT between tissue pairs using absolute difference
# transformed to a similarity score, OR we pool pairwise (m1,m2) points
# across all patients and compute one correlation per CpG across patients'
# tissue-pair observations. This is the correct interpretation of "intra-
# individual correlation, averaged across patients."
compute_same_layer_r <- function(meth_sub, layer_name, min_patients = 1) {
  
  layer_meth <- meth_sub[germ_layer == layer_name]
  
  # patients with >1 tissue in this layer
  pat_tissue_n <- layer_meth[, .(n = uniqueN(source_tissue_celltype)),
                             by = patient_id]
  eligible_patients <- pat_tissue_n[n > 1, patient_id]
  
  if (length(eligible_patients) < min_patients) return(NULL)
  
  # for each eligible patient, get ALL pairwise tissue (m1, m2) points per CpG
  all_pairs <- rbindlist(lapply(eligible_patients, function(pid) {
    pat_meth <- layer_meth[patient_id == pid]
    tissues  <- unique(pat_meth$source_tissue_celltype)
    pv <- patient_mean_pairwise_r(pat_meth, tissues)
    if (is.null(pv)) return(NULL)
    pv[, patient_id := pid]
    pv
  }), fill = TRUE)
  
  if (is.null(all_pairs) || nrow(all_pairs) == 0) return(NULL)
  
  # ── Per-CpG: average ABSOLUTE per-patient tissue agreement ─────────────────
  # For each patient × tissue-pair × CpG we have (m1, m2). Within ONE patient
  # there's no distribution to correlate, so "intra-individual correlation"
  # must be computed ACROSS patients: i.e. does m1 predict m2 across the
  # pool of (patient, tissue-pair) observations? This tests whether the
  # methylation level transfers consistently between tissue types within
  # individuals, pooling all patients together.
  all_pairs[, {
    if (.N >= 3) {
      r <- suppressWarnings(cor(m1, m2, method = "pearson"))
      list(r = r, n_obs = .N)
    } else list(r = NA_real_, n_obs = .N)
  }, by = cpg_site]
}

# ── Core function: cross-layer within-individual r per CpG ───────────────────
# Same logic but for the 4 Endo+Meso patients: pair each patient's Endo
# tissue value with their Meso tissue value (mean if >1 tissue per layer
# for that patient), pool across the 4 patients, correlate per CpG
compute_cross_layer_r <- function(meth_sub, layer1, layer2) {
  
  d1 <- meth_sub[germ_layer == layer1,
                 .(m1 = mean(methylation, na.rm = TRUE)),
                 by = .(cpg_site, patient_id)]
  d2 <- meth_sub[germ_layer == layer2,
                 .(m2 = mean(methylation, na.rm = TRUE)),
                 by = .(cpg_site, patient_id)]
  
  dm <- merge(d1, d2, by = c("cpg_site", "patient_id"))
  
  dm[, {
    idx <- !is.na(m1) & !is.na(m2)
    if (sum(idx) >= 3) {
      r <- suppressWarnings(cor(m1[idx], m2[idx], method = "pearson"))
      list(r = r, n_obs = sum(idx))
    } else list(r = NA_real_, n_obs = sum(idx))
  }, by = cpg_site]
}

# ══════════════════════════════════════════════════════════════════════════════
# Run for Meso_specific and Endo_specific
# ══════════════════════════════════════════════════════════════════════════════

# ── Meso_specific: same-layer (Meso×Meso) + cross-layer (Endo×Meso) ──────────
meth_meso_specific <- meth_multi[category == "Meso_specific"]

r_same_meso  <- compute_same_layer_r(meth_meso_specific, "Meso")
r_cross_meso <- compute_cross_layer_r(meth_meso_specific, "Endo", "Meso")

meso_summary <- merge(
  r_same_meso[,  .(cpg_site, r_same  = r)],
  r_cross_meso[, .(cpg_site, r_cross = r)],
  by = "cpg_site", all = TRUE
)

pct_meso_layer_specific <- meso_summary[
  !is.na(r_same) & !is.na(r_cross),
  mean(r_same >= r_threshold & r_cross < r_threshold) * 100]

message(sprintf(
  "Meso_specific: %.1f%% show high same-layer r (>=%.1f) AND low cross-layer r (<%.1f) (n=%d CpGs tested)",
  pct_meso_layer_specific, r_threshold, r_threshold,
  sum(!is.na(meso_summary$r_same) & !is.na(meso_summary$r_cross))))
# Meso_specific: 26.4% show high same-layer r (>=0.5) AND low cross-layer r (<0.5) (n=7972 CpGs tested)

# ── Endo_specific: same-layer (Endo×Endo) + cross-layer (Endo×Meso) ──────────
meth_endo_specific <- meth_multi[category == "Endo_specific"]

r_same_endo  <- compute_same_layer_r(meth_endo_specific, "Endo")
r_cross_endo <- compute_cross_layer_r(meth_endo_specific, "Endo", "Meso")

endo_summary <- merge(
  r_same_endo[,  .(cpg_site, r_same  = r)],
  r_cross_endo[, .(cpg_site, r_cross = r)],
  by = "cpg_site", all = TRUE
)

pct_endo_layer_specific <- endo_summary[
  !is.na(r_same) & !is.na(r_cross),
  mean(r_same >= r_threshold & r_cross < r_threshold) * 100]

message(sprintf(
  "Endo_specific: %.1f%% show high same-layer r (>=%.1f) AND low cross-layer r (<%.1f) (n=%d CpGs tested)",
  pct_endo_layer_specific, r_threshold, r_threshold,
  sum(!is.na(endo_summary$r_same) & !is.na(endo_summary$r_cross))))
# Endo_specific: 45.2% show high same-layer r (>=0.5) AND low cross-layer r (<0.5) (n=3072 CpGs tested)

# ══════════════════════════════════════════════════════════════════════════════
# Negative controls: random 10k constitutive and 10k ambiguous CpGs
# Need to extract methylation for these first (not yet in meth_multi)
# ══════════════════════════════════════════════════════════════════════════════

set.seed(1234)
constitutive_sample <- sample(wideFull[category == "constitutive", name], 10000)
ambiguous_sample    <- sample(wideFull[category == "ambiguous",    name], 10000)

writeLines(c(constitutive_sample, ambiguous_sample),
           here("B_MultiTissues/dataOut/control_sample20k.txt"))
message("Written: control_sample20k.txt (20000 CpGs: 10k constitutive + 10k ambiguous)")

## In pchuckle:
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/control_sample20k.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_control20k.tsv \
# --minCov    10

# ── once extracted, load and process the same way ─────────────────────────────
meth_control <- fread(here("gitignore/methylation_control20k.tsv"))
meth_control <- merge(meth_control, tissue_to_layer,
                      by = "source_tissue_celltype", all.x = TRUE)
meth_control[, category := fcase(
  cpg_site %in% constitutive_sample, "constitutive",
  cpg_site %in% ambiguous_sample,    "ambiguous",
  default = NA_character_
)]

meth_control_multi <- meth_control[patient_id %in% multi_patients]

# ── constitutive: same-layer (test in Meso, since that's where we have power)
r_same_const  <- compute_same_layer_r(meth_control_multi[category == "constitutive"], "Meso")
r_cross_const <- compute_cross_layer_r(meth_control_multi[category == "constitutive"], "Endo", "Meso")

const_summary <- merge(
  r_same_const[,  .(cpg_site, r_same  = r)],
  r_cross_const[, .(cpg_site, r_cross = r)],
  by = "cpg_site", all = TRUE
)

pct_const_layer_specific <- const_summary[
  !is.na(r_same) & !is.na(r_cross),
  mean(r_same >= r_threshold & r_cross < r_threshold) * 100]

message(sprintf(
  "Constitutive (negative control): %.1f%% show high same-layer r AND low cross-layer r (n=%d CpGs tested)",
  pct_const_layer_specific,
  sum(!is.na(const_summary$r_same) & !is.na(const_summary$r_cross))))

# Constitutive (negative control): 0.5% show high same-layer r AND low cross-layer r (n=9104 CpGs tested)

# ── ambiguous: same-layer + cross-layer ───────────────────────────────────────
r_same_amb  <- compute_same_layer_r(meth_control_multi[category == "ambiguous"], "Meso")
r_cross_amb <- compute_cross_layer_r(meth_control_multi[category == "ambiguous"], "Endo", "Meso")

amb_summary <- merge(
  r_same_amb[,  .(cpg_site, r_same  = r)],
  r_cross_amb[, .(cpg_site, r_cross = r)],
  by = "cpg_site", all = TRUE
)

pct_amb_layer_specific <- amb_summary[
  !is.na(r_same) & !is.na(r_cross),
  mean(r_same >= r_threshold & r_cross < r_threshold) * 100]

message(sprintf(
  "Ambiguous (negative control): %.1f%% show high same-layer r AND low cross-layer r (n=%d CpGs tested)",
  pct_amb_layer_specific,
  sum(!is.na(amb_summary$r_same) & !is.na(amb_summary$r_cross))))
# Ambiguous (negative control): 2.9% show high same-layer r AND low cross-layer r (n=9081 CpGs tested)

# ══════════════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════════════

results_summary <- data.table(
  category = c("Meso_specific", "Endo_specific", "constitutive (neg ctrl)", "ambiguous (neg ctrl)"),
  pct_layer_specific_signature = c(pct_meso_layer_specific, pct_endo_layer_specific,
                                   pct_const_layer_specific, pct_amb_layer_specific)
)
print(results_summary)
#                category         pct_layer_specific_signature
# <char>                        <num>
#   1:           Meso_specific                   26.3923733
# 2:           Endo_specific                   45.2473958
# 3: constitutive (neg ctrl)                    0.5272408
# 4:    ambiguous (neg ctrl)                    2.9071688


## TBC: plot? pct_layer_specific_signature divided for Meso and endo!

# # ####################
# # ## GO enrichement ##
# # ####################
# # 
