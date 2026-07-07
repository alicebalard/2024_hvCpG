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

## Ccl:
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
# 32 patients with >1 tissue
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

#    germ_layers n_patients median_tissues range_tissues n_blood_only pct_blood_only
# 1:        Ecto          6            2.0           2-2            0              0
# 2:        Endo         13            2.0           2-3            0              0
# 3:   Endo+Meso          4            2.5           2-3            0              0
# 4:        Meso          9            5.0           2-7            6             67

message(sprintf("Patients with cross-layer samples: %d",
                sum(grepl("\\+", patient_table$germ_layers))))
# Patients with cross-layer samples: 4 (Endo+Meso)
unique(meth[
  meth$patient_id %in% unlist(
    patient_table[patient_table$germ_layers %in% "Endo+Meso","patient_id"]),
  c("source_tissue_celltype", "patient_id", "germ_layer")]) 

# source_tissue_celltype patient_id germ_layer
# 1:               Colon - Endocrine        169       Endo
# 2:             Colon - Macrophages        169       Meso
# 3: Kidney glomerular - Endothelium        176       Meso
# 10:    Kidney tubular - Epithelium        176       Endo
# 7:    Kidney glomerular - Podocyte        176       Endo
# 4: Kidney glomerular - Endothelium        130       Meso
# 5:  Kidney glomerular - Epithelium        130       Endo
# 9:     Kidney tubular - Epithelium        130       Endo
# 6:    Kidney glomerular - Podocyte        199       Endo
# 8:    Kidney tubular - Endothelium        199       Meso

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

## I have 4 patients with Endo+Meso (but only 2 with a given exact combination),
## 9 with multiple meso, 6 with multiple ecto
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

# ── Extract "true layer-specific" CpGs on THREE criteria ──────────────────────
#   1. high  same-layer intra-individual r in the OWN layer   (r_own   >= r_high)
#   2. low   same-layer intra-individual r in the OTHER layer (r_other <  r_low)
#   3. low   cross-layer within-person r (own <-> other)       (r_cross <  r_low)
extract_layer_specific_cpgs <- function(meth_sub,
                                        own_layer, other_layer,
                                        r_high = 0.5,      # "high" bar for own layer
                                        r_low  = 0.2,      # "low"  bar (matches notHVt)
                                        min_same_obs  = 3, # pooled tissue-pair points
                                        min_cross_obs = 3) # cross-layer patients
{
  r_own   <- compute_same_layer_r(meth_sub, own_layer)                 # high in own
  r_other <- compute_same_layer_r(meth_sub, other_layer)              # low in other
  r_cross <- compute_cross_layer_r(meth_sub, other_layer, own_layer)  # low cross-layer
  if (is.null(r_own) || is.null(r_other) || is.null(r_cross)) return(NULL)
  
  setnames(r_own,   c("r", "n_obs"), c("r_own",   "n_own"))
  setnames(r_other, c("r", "n_obs"), c("r_other", "n_other"))
  setnames(r_cross, c("r", "n_obs"), c("r_cross", "n_cross"))
  
  m <- Reduce(function(a, b) merge(a, b, by = "cpg_site"),
              list(r_own, r_other, r_cross))               # inner: need all three
  
  hits <- m[!is.na(r_own) & !is.na(r_other) & !is.na(r_cross) &
              n_own   >= min_same_obs &
              n_other >= min_same_obs &
              n_cross >= min_cross_obs &
              r_own   >= r_high &     # criterion 1: high in Meso
              r_other <  r_low  &     # criterion 2: low in Endo
              r_cross <  r_low]       # criterion 3: low Endo<->Meso
  hits[order(-r_own)]
}

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

# ── Load control methylation (extracted previously) ───────────────────────────
meth_control <- fread(here("gitignore/methylation_control20k.tsv"))
meth_control <- merge(meth_control, tissue_to_layer,
                      by = "source_tissue_celltype", all.x = TRUE)

# re-assign categories (constitutive vs ambiguous)
set.seed(1234)
constitutive_sample <- sample(wideFull[category == "constitutive", name], 10000)
ambiguous_sample    <- sample(wideFull[category == "ambiguous",    name], 10000)

meth_control[, category := fcase(
  cpg_site %in% constitutive_sample, "constitutive",
  cpg_site %in% ambiguous_sample,    "ambiguous",
  default = NA_character_
)]

meth_control_multi <- meth_control[patient_id %in% multi_patients]

message(sprintf("Control CpGs loaded: %d constitutive, %d ambiguous",
                uniqueN(meth_control_multi[category == "constitutive", cpg_site]),
                uniqueN(meth_control_multi[category == "ambiguous",    cpg_site])))

# ══════════════════════════════════════════════════════════════════════════════
# Run for Meso_specific, Endo_specific and controls
# ══════════════════════════════════════════════════════════════════════════════

# ── Meso_specific ─────────────────────────────────────────────────────────────
meth_meso_specific <- meth_multi[category == "Meso_specific"]
r_same_meso        <- compute_same_layer_r(meth_meso_specific, "Meso")

# ── Endo_specific ─────────────────────────────────────────────────────────────
meth_endo_specific <- meth_multi[category == "Endo_specific"]
r_same_endo        <- compute_same_layer_r(meth_endo_specific, "Endo")

# ── Wrong-layer specificity checks ───────────────────────────────────────────
r_same_endo_for_meso <- compute_same_layer_r(meth_meso_specific, "Endo")
r_same_meso_for_endo <- compute_same_layer_r(meth_endo_specific, "Meso")

# ── Controls in both Meso and Endo contexts ───────────────────────────────────
r_same_meso_const      <- compute_same_layer_r(meth_control_multi[category == "constitutive"], "Meso")
r_same_endo_const <- compute_same_layer_r(meth_control_multi[category == "constitutive"], "Endo")
r_same_meso_amb        <- compute_same_layer_r(meth_control_multi[category == "ambiguous"],    "Meso")
r_same_endo_amb   <- compute_same_layer_r(meth_control_multi[category == "ambiguous"],    "Endo")

# ── Summary function ──────────────────────────────────────────────────────────
make_summary_samelayer <- function(r_same, layer_tested, category_name) {
  s <- r_same[!is.na(r)]
  data.table(
    category          = category_name,
    same_layer_tested = layer_tested,
    n_cpgs_tested     = uniqueN(s$cpg_site),
    pct_high_same     = s[, mean(r >= r_threshold) * 100, by = cpg_site][, mean(V1)],
    mean_r            = mean(s$r, na.rm = TRUE),
    median_r          = median(s$r, na.rm = TRUE)
  )
}

# ── Build summary ─────────────────────────────────────────────────────────────
results_summary <- rbindlist(list(
  make_summary_samelayer(r_same_meso,          "Meso", "Meso_specific"),
  make_summary_samelayer(r_same_endo_for_meso, "Endo", "Meso_specific"),
  make_summary_samelayer(r_same_endo,          "Endo", "Endo_specific"),
  make_summary_samelayer(r_same_meso_for_endo, "Meso", "Endo_specific"),
  make_summary_samelayer(r_same_meso_const,    "Meso", "constitutive"),
  make_summary_samelayer(r_same_endo_const,    "Endo", "constitutive"),
  make_summary_samelayer(r_same_meso_amb,      "Meso", "ambiguous"),
  make_summary_samelayer(r_same_endo_amb,      "Endo", "ambiguous")
))

print(results_summary)
#        category same_layer_tested n_cpgs_tested pct_high_same     mean_r   median_r
# 1: Meso_specific              Meso          7972    40.9433016 0.45638404 0.39158308
# 2: Meso_specific              Endo          8949    26.7962901 0.27413938 0.23979101
# 3: Endo_specific              Endo          3072    71.6471354 0.61370364 0.67129866
# 4: Endo_specific              Meso          2702     4.1080681 0.12484106 0.09576758
# 5:  constitutive              Meso          9111     0.6695204 0.05071192 0.02592049
# 6:  constitutive              Endo          9924    10.8625554 0.12255948 0.08111064
# 7:     ambiguous              Meso          9081     5.6491576 0.19011905 0.17497508
# 8:     ambiguous              Endo          9893    42.1510159 0.39634798 0.40997976

# ── Annotate and plot ─────────────────────────────────────────────────────────
results_summary[, expected := fcase(
  category == "Meso_specific" & same_layer_tested == "Meso", "own layer",
  category == "Meso_specific" & same_layer_tested == "Endo", "wrong layer",
  category == "Endo_specific" & same_layer_tested == "Endo", "own layer",
  category == "Endo_specific" & same_layer_tested == "Meso", "wrong layer",
  default = "control"
)]

results_summary[, category_f := factor(category,
                                       levels = c("Meso_specific", "Endo_specific", "constitutive", "ambiguous"))]

bar_colours <- c(
  "own layer"   = "#2166AC",
  "wrong layer" = "#D55E00",
  "control"     = "grey60"
)

ggplot(results_summary,
       aes(x = category_f, y = pct_high_same,
           fill = expected)) +
  geom_col(position = "dodge", width = 0.6,
           colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct_high_same)),
            position = position_dodge(width = 0.6),
            vjust = -0.4, size = 3) +
  facet_wrap(~ same_layer_tested, nrow = 1,
             labeller = labeller(same_layer_tested = c(
               Meso = "Tested in Meso tissues",
               Endo = "Tested in Endo tissues"))) +
  scale_fill_manual(values = bar_colours, name = "Layer context") +
  scale_y_continuous("% CpGs with high same-layer r (>=0.5)",
                     limits = c(0, 85)) +
  scale_x_discrete("CpG category") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        axis.text.x       = element_text(angle = 30, hjust = 1),
        strip.text        = element_text(face = "bold"),
        legend.position   = "right") +
  ggtitle("Intra-individual same-layer methylation concordance",
          subtitle = sprintf(
            "Own layer (blue) should be HIGH | Wrong layer (orange) should be LOW\nMeso: %d patients | Endo: %d patients | r threshold = %.1f",
            9, 15, r_threshold))

# ── Density plot of same-layer r per category and layer context ───────────────

# pool all r values into one table with category and layer context
density_dt <- rbindlist(list(
  r_same_meso[!is.na(r),          .(cpg_site, r, category = "Meso_specific", layer_context = "Meso")],
  r_same_endo_for_meso[!is.na(r), .(cpg_site, r, category = "Meso_specific", layer_context = "Endo")],
  r_same_endo[!is.na(r),          .(cpg_site, r, category = "Endo_specific", layer_context = "Endo")],
  r_same_meso_for_endo[!is.na(r), .(cpg_site, r, category = "Endo_specific", layer_context = "Meso")],
  r_same_meso_const[!is.na(r),         .(cpg_site, r, category = "constitutive",  layer_context = "Meso")],
  r_same_endo_const[!is.na(r),    .(cpg_site, r, category = "constitutive",  layer_context = "Endo")],
  r_same_meso_amb[!is.na(r),           .(cpg_site, r, category = "ambiguous",     layer_context = "Meso")],
  r_same_endo_amb[!is.na(r),      .(cpg_site, r, category = "ambiguous",     layer_context = "Endo")]
))

density_dt[, category_f := factor(category,
                                  levels = c("Meso_specific", "Endo_specific", "constitutive", "ambiguous"))]

density_dt[, expected := fcase(
  category == "Meso_specific" & layer_context == "Meso", "own layer",
  category == "Meso_specific" & layer_context == "Endo", "wrong layer",
  category == "Endo_specific" & layer_context == "Endo", "own layer",
  category == "Endo_specific" & layer_context == "Meso", "wrong layer",
  default = "control"
)]

line_colours <- c(
  "own layer"   = "#2166AC",
  "wrong layer" = "#D55E00",
  "control"     = "grey50"
)

ggplot(density_dt, aes(x = r, colour = layer_context, fill = layer_context)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  facet_wrap(~ category_f, nrow = 1) +
  scale_colour_manual(values = c(Meso = "#56B4E9", Endo = "#009E73"),
                      name = "Layer tested") +
  scale_fill_manual(values   = c(Meso = "#56B4E9", Endo = "#009E73"),
                    name = "Layer tested") +
  scale_x_continuous("Same-layer r", limits = c(-1, 1),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous("Density") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "right") +
  ggtitle("Distribution of same-layer r per CpG category",
          subtitle = sprintf(
            "Meso (blue): %d patients with multiple mesodermal tissues | Endo (green): %d patients with multiple endodermal tissues",
            9, 15))

# ══════════════════════════════════════════════════════════════════════════════
# ── Blood composition confound test ──────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
same_r_all <- rbind(
  r_same_meso[,  .(cpg_site, r, category = "Meso_specific")],
  r_same_meso_const[, .(cpg_site, r, category = "constitutive")],
  r_same_meso_amb[,   .(cpg_site, r, category = "ambiguous")]
)

ggplot(same_r_all[!is.na(r)],
       aes(x = r, fill = category, colour = category)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = r_threshold, linetype = "dashed", colour = "grey40") +
  scale_x_continuous("Same-layer r (pooled across patients)", limits = c(-1, 1)) +
  scale_fill_manual(values = c(Meso_specific = "#56B4E9",
                               constitutive  = "grey30",
                               ambiguous     = "grey70")) +
  scale_colour_manual(values = c(Meso_specific = "#56B4E9",
                                 constitutive  = "grey30",
                                 ambiguous     = "grey70")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("Same-layer (Meso) r: Meso_specific vs controls",
          subtitle = paste0(
            "Blood composition confound predicts: all three similar\n",
            "Genuine signal predicts: Meso_specific shifted right toward r=1"))

# Meso_specific CpGs are highly variable between individuals (by selection)
# Within each individual, their methylation level is highly consistent across all blood cell types (r≈1, Plot 1)
# This within-individual concordance across independent blood lineages is the hallmark of pre-haematopoietic stochastic establishment
# Controls show r≈0, confirming this is not a general property of blood measurements

###### Add now the 3rd criteria, low r between both layers

r_high        <- 0.5    # own-layer "high"
r_low         <- 0.2    # other/cross "low" (matches notHVt)
min_same_obs  <- 3
min_cross_obs <- 3

# per-CpG: the three r's for one category under one pipeline (own vs other) ------
tri_rs <- function(meth_sub, own_layer, other_layer, category_name) {
  ro <- compute_same_layer_r(meth_sub, own_layer)
  rt <- compute_same_layer_r(meth_sub, other_layer)
  rc <- compute_cross_layer_r(meth_sub, other_layer, own_layer)
  if (is.null(ro)) return(NULL)
  setnames(ro, c("r","n_obs"), c("r_own","n_own"))
  base <- ro[, .(cpg_site, r_own, n_own)]
  if (!is.null(rt)) { setnames(rt, c("r","n_obs"), c("r_other","n_other"))
    base <- merge(base, rt[, .(cpg_site, r_other, n_other)], by="cpg_site", all.x=TRUE)
  } else base[, `:=`(r_other = NA_real_, n_other = 0L)]
  if (!is.null(rc)) { setnames(rc, c("r","n_obs"), c("r_cross","n_cross"))
    base <- merge(base, rc[, .(cpg_site, r_cross, n_cross)], by="cpg_site", all.x=TRUE)
  } else base[, `:=`(r_cross = NA_real_, n_cross = 0L)]
  base[, `:=`(category = category_name, own_layer = own_layer)][]
}

# cumulative pass counts over the 3 criteria ------------------------------------
tri_steps <- function(tab) {
  if (is.null(tab)) return(NULL)
  tested <- tab[!is.na(r_own) & n_own >= min_same_obs]
  N  <- nrow(tested)
  p1 <- tested[r_own >= r_high]
  p2 <- p1[!is.na(r_other) & n_other >= min_same_obs & r_other < r_low]
  p3 <- p2[!is.na(r_cross) & n_cross >= min_cross_obs & r_cross < r_low]
  data.table(
    category  = tab$category[1],
    own_layer = tab$own_layer[1],
    step = factor(c("1. high own r", "2. + low other r", "3. + low cross r"),
                  levels = c("1. high own r", "2. + low other r", "3. + low cross r")),
    n   = c(nrow(p1), nrow(p2), nrow(p3)),
    pct = 100 * c(nrow(p1), nrow(p2), nrow(p3)) / N,
    n_tested = N)
}

tri_tabs <- rbindlist(list(
  # Meso pipeline (own = Meso, other = Endo)
  tri_steps(tri_rs(meth_meso_specific,                       "Meso","Endo","Meso_specific")),
  tri_steps(tri_rs(meth_endo_specific,                       "Meso","Endo","Endo_specific")),
  tri_steps(tri_rs(meth_control_multi[category=="constitutive"], "Meso","Endo","constitutive")),
  tri_steps(tri_rs(meth_control_multi[category=="ambiguous"],    "Meso","Endo","ambiguous")),
  # Endo pipeline (own = Endo, other = Meso)
  tri_steps(tri_rs(meth_endo_specific,                       "Endo","Meso","Endo_specific")),
  tri_steps(tri_rs(meth_meso_specific,                       "Endo","Meso","Meso_specific")),
  tri_steps(tri_rs(meth_control_multi[category=="constitutive"], "Endo","Meso","constitutive")),
  tri_steps(tri_rs(meth_control_multi[category=="ambiguous"],    "Endo","Meso","ambiguous"))
), fill = TRUE)

tri_tabs[, category_f := factor(category,
                                levels = c("Meso_specific","Endo_specific","constitutive","ambiguous"))]

ggplot(tri_tabs, aes(category_f, pct, fill = step)) +
  geom_col(position = position_dodge(0.7), width = 0.65,
           colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.7),
            vjust = -0.3, size = 2.7) +
  facet_wrap(~ own_layer, nrow = 1,
             labeller = labeller(own_layer = c(
               Meso = "Own = Meso (other = Endo)",
               Endo = "Own = Endo (other = Meso)"))) +
  scale_fill_brewer(palette = "Blues", name = "Cumulative criterion") +
  scale_y_continuous("% of CpGs tested in own layer", limits = c(0, NA)) +
  scale_x_discrete("CpG category") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text  = element_text(face = "bold"),
        legend.position = "right") +
  ggtitle("Layer-specific selection: cumulative pass rate over 3 criteria",
          subtitle = sprintf("high own r (\u2265%.1f)  \u2192  + low other r (<%.1f)  \u2192  + low cross-layer r (<%.1f)",
                             r_high, r_low, r_low))

## Higherbackground concordance for endo than meso

#################################
## Save the interesting targets##
#################################
meso_hits <- extract_layer_specific_cpgs(meth_meso_specific,
                                         own_layer = "Meso", other_layer = "Endo")
endo_hits <- extract_layer_specific_cpgs(meth_endo_specific,
                                         own_layer = "Endo", other_layer = "Meso")

message(sprintf("Meso_specific true-signature CpGs: %d", nrow(meso_hits)))
message(sprintf("Endo_specific true-signature CpGs: %d", nrow(endo_hits)))

head(meso_hits)   # cpg_site, r_same, n_same, r_cross, n_cross

#################################
## Annotation of these targets ##
#################################
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# ── Shared annotation objects (built once) ────────────────────────────────────
txdb     <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
genes_gr$symbol <- mapIds(org.Hs.eg.db, names(genes_gr), "SYMBOL", "ENTREZID")

cpg_to_gr <- function(cpg) {                     # "chr7_107543290" -> 1bp GRanges
  GRanges(sub("_.*", "", cpg),
          IRanges(as.integer(sub(".*_", "", cpg)), width = 1),
          cpg_site = cpg)
}

# universe of tested CpGs (built once, reused by every call)
uni <- copy(wideFull[, .(cpg_site = name, category)])
uni[, `:=`(chr = sub("_.*", "", cpg_site),
           pos = as.integer(sub(".*_", "", cpg_site)))]
uni_gr <- GRanges(uni$chr, IRanges(uni$pos, width = 1))

# ── Main function ─────────────────────────────────────────────────────────────
annotate_layer_hits <- function(hits,                       # data.table with cpg_site
                                label,                     # e.g. "meso" / "endo"
                                genes_gr, uni, uni_gr,
                                gap        = 50,           # bp to merge hit clusters
                                min_hits   = 2,            # min hits per cluster
                                drop_regex = "LOC|LINC",   # gene names to discard
                                flank      = 5000,
                                out_dir    = here("B_MultiTissues/dataOut")) {
  
  hit_col <- paste0("is_", label, "_hit")
  hit_set <- hits$cpg_site
  
  ## (1a) annotate each hit with nearest gene ---------------------------------
  hits_gr <- cpg_to_gr(hit_set)
  nr  <- distanceToNearest(hits_gr, genes_gr, ignore.strand = TRUE)
  ann <- as.data.table(hits)
  ann[, `:=`(gene = NA_character_, dist_to_gene = NA_integer_)]
  ann[queryHits(nr),
      `:=`(gene         = genes_gr$symbol[subjectHits(nr)],
           dist_to_gene = mcols(nr)$distance)]
  ann[, `:=`(pos = as.integer(sub(".*_", "", cpg_site)),
             chr = sub("_.*", "", cpg_site))]
  setorder(ann, chr, pos)
  
  ## (1b) cluster hits within `gap` bp ----------------------------------------
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
  
  ## (2) all CpGs in each top gene + flank, flagged ---------------------------
  top_genes <- unique(clust_dt$gene)
  extract_gene_cpgs <- function(sym) {
    g <- genes_gr[which(genes_gr$symbol == sym)]
    if (!length(g)) return(NULL)
    g   <- g[1]
    win <- GRanges(as.character(seqnames(g)),
                   IRanges(start(g) - flank, end(g) + flank))
    out <- uni[which(overlapsAny(uni_gr, win, ignore.strand = TRUE))]
    out[, `:=`(gene = sym,
               gene_start = start(g), gene_end = end(g),
               in_gene_body = pos >= start(g) & pos <= end(g))]
    out[, (hit_col) := cpg_site %in% hit_set]
    out[order(pos)]
  }
  gene_cpgs <- rbindlist(lapply(top_genes, extract_gene_cpgs), fill = TRUE)
  
  fwrite(gene_cpgs,
         file.path(out_dir, sprintf("%s_hits_topGenes_allCpGs.csv", label)))
  
  list(ann = ann, clusters = clust_dt, gene_cpgs = gene_cpgs)
}

meso_res <- annotate_layer_hits(meso_hits, "meso", genes_gr, uni, uni_gr)
endo_res <- annotate_layer_hits(endo_hits, "endo", genes_gr, uni, uni_gr)

meso_res$clusters  # clustered genes
endo_res$clusters

# save all CpGs in top genes ± flank, is_meso_hit column
fwrite(meso_res$gene_cpgs,
       here("B_MultiTissues/dataOut/meso_hits_genesClusters.csv"))

fwrite(endo_res$gene_cpgs,
       here("B_MultiTissues/dataOut/endo_hits_genesClusters.csv"))

##TO DO:addan MEtest to findsome more
