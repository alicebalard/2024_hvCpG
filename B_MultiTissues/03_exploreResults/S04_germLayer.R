#####################################################################
# Step 1 - Find inter-individual hypervariable CpGs per germ layer
# Step 2 — Test systemic intra-individual correlations
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))
## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}

if (!file.exists(here("gitignore/wide_script04_3layers_noNA.RDS"))){
  # ── Load ───────────────────────────────────────────────────────────────────
  endo     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
  meso     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
  ecto     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
  analyses <- list(endo = endo, meso = meso, ecto = ecto)
  
  # ── Make wide table ─────────────────────────────────────────────────────────────
  # keep only name + the score, renamed per layer, from each table
  slim <- Map(function(dt, nm) {
    dt <- as.data.table(dt)
    setNames(dt[, .(name, logBF_per_ds)], c("name", nm))
  }, analyses, names(analyses))
  
  # key once, then join — much cheaper than merging full tables
  wide <- Reduce(function(a, b) a[b, on = "name"], slim)
  
  ## Keep only CpGs tested in all 3 layers
  wide <- wide[rowSums(is.na(wide)) == 0]
  nrow(wide) # 21.522.541
  saveRDS(wide, here("gitignore/wide_script04_3layers_noNA.RDS"))
  rm(endo, meso, ecto, analyses, wide)
}

if (!file.exists(here("gitignore/wide_script04_3layers_noNA_6GP.RDS"))){
  # ── Load ───────────────────────────────────────────────────────────────────
  endo     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_2_endo6gp.rds"))
  meso     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_2_meso6gp.rds"))
  ecto     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
  analyses <- list(endo = endo, meso = meso, ecto = ecto)
  
  # ── Make wide table ─────────────────────────────────────────────────────────────
  # keep only name + the score, renamed per layer, from each table
  slim <- Map(function(dt, nm) {
    dt <- as.data.table(dt)
    setNames(dt[, .(name, logBF_per_ds)], c("name", nm))
  }, analyses, names(analyses))
  
  # key once, then join — much cheaper than merging full tables
  wide <- Reduce(function(a, b) a[b, on = "name"], slim)
  
  ## Keep only CpGs tested in all 3 layers
  wide <- wide[rowSums(is.na(wide)) == 0]
  nrow(wide) # 21.522.541
  saveRDS(wide, here("gitignore/wide_script04_3layers_noNA_6GP.RDS"))
  rm(endo, meso, ecto, analyses, wide)
}

if (!exists("wideFull")){wideFull   <- readRDS(here("gitignore/wide_script04_3layers_noNA.RDS"))}
if (!exists("wide6gpall")){wide6gpall <- readRDS(here("gitignore/wide_script04_3layers_noNA_6GP.RDS"))}

# ══════════════════════════════════════════════════════════════════════════════
# Step 1. Identify candidate sites
# ══════════════════════════════════════════════════════════════════════════════

#############################################
## Define parameters based on known MEs #####
#############################################

# Step 1 — Build a graded truth set by cross-catalog overlap

# the individual catalogs as a named list (region- or CpG-level, all hg38)
me_sets <- list(
  DerakhshanhvCpGs = DerakhshanhvCpGs_hg38_GR,
  HarrisSIV     = HarrisSIV_hg38_GR,
  VanBaakESS = VanBaakESS_hg38_GR,
  VanBaakSIV = VanBaakSIV_hg38_GR,
  KesslerSIV    = KesslerSIV_GRanges_hg38,
  GunasekaraCorSIV = corSIV_GRanges_hg38
)

# the Atlas CpG universe as a GRanges (the CpGs I can actually score)
uni <- copy(wideFull[, .(cpg_site = name)])
uni[, `:=`(chr = sub("_.*", "", cpg_site),
           pos = as.integer(sub(".*_", "", cpg_site)))]
uni_gr <- GRanges(uni$chr, IRanges(uni$pos, width = 1), cpg_site = uni$cpg_site)

# for each Atlas CpG, I count how many catalogs it overlaps (region OR CpG level)
# pad CpG-level sets a little so a nearby Atlas CpG still counts as "in" the locus
pad <- 250   # bp; region sets already have width, CpG sets become +/-250
hits_mat <- sapply(me_sets, function(gr) {
  gr2 <- if (all(width(gr) == 1)) gr + pad else gr
  overlapsAny(uni_gr, gr2, ignore.strand = TRUE)
})
uni[, n_catalogs := rowSums(hits_mat)]
uni[, is_ME_any  := n_catalogs >= 1]
uni[, is_ME_hi   := n_catalogs >= 2]   # replicated in >=2 independent screens

message(sprintf("Atlas CpGs overlapping >=1 catalog: %d | >=2 (high-conf): %d",
                sum(uni$is_ME_any), sum(uni$is_ME_hi)))
# Atlas CpGs overlapping >=1 catalog: 117388 | >=2 (high-conf): 14925
table(uni$n_catalogs)
#        0        1        2        3        4        5        6 
# 21405153   102463     9820     3115     1432      514       44 

# Step 2 — Attach metrics to the truth set

# Take the top 10% per layer: convert each layer's threshold to that layer's 90th
# percentile (and the "not HV" threshold to its 10th percentile). 
# Compute the cutoffs once per column, then classify against them

classify_wide <- function(wide, hv_q = 0.90, nothv_q = 0.10) {
  q <- lapply(c("endo", "meso", "ecto"), function(layer) {
    c(hv    = unname(quantile(wide[[layer]], hv_q,    na.rm = TRUE)),
      notHV = unname(quantile(wide[[layer]], nothv_q, na.rm = TRUE)))
  })
  names(q) <- c("endo", "meso", "ecto")
  
  wide[, `:=`(
    HV_endo = endo > q$endo["hv"], HV_meso = meso > q$meso["hv"], HV_ecto = ecto > q$ecto["hv"],
    notHV_endo = endo < q$endo["notHV"], notHV_meso = meso < q$meso["notHV"], notHV_ecto = ecto < q$ecto["notHV"]
  )]
  
  wide[, category := fcase(
    HV_meso & HV_endo & HV_ecto,           "ME",
    HV_meso & notHV_endo & notHV_ecto,     "Meso_specific",
    HV_endo & notHV_meso & notHV_ecto,     "Endo_specific",
    HV_ecto & notHV_meso & notHV_endo,     "Ecto_specific",
    notHV_meso & notHV_endo & notHV_ecto,  "constitutive",
    default =                              "ambiguous"
  )]
  attr(wide, "thresholds") <- q
  wide
}

wideFull   <- classify_wide(wideFull);   print(table(wideFull$category))
attr(wideFull, "thresholds")        # should show 6 finite numbers, not NA
wideFull[, mean(HV_endo, na.rm = TRUE)]   # should be ≈ 0.10 (top 10% by construction)
print(table(wideFull$category))
# ambiguous  constitutive Ecto_specific Endo_specific            ME Meso_specific 
# 20450363        402815           144            34        669018           167 

wide6gpall <- classify_wide(wide6gpall); print(table(wide6gpall$category))
# ambiguous  constitutive Ecto_specific Endo_specific            ME Meso_specific 
# 21488395        302166           298           275        535640           687 

# Overlap of each layer-specific set, reduced to 6 groups
# strict layer-specific sets, full vs 6-group
spec_full <- split(wideFull[category %like% "_specific", name],
                   wideFull[category %like% "_specific", category])
spec_6gp  <- split(wide6gpall[category %like% "_specific", name],
                   wide6gpall[category %like% "_specific", category])

# per-layer recovery
for (cat in c("Endo_specific", "Meso_specific", "Ecto_specific")) {
  f <- spec_full[[cat]]; g <- spec_6gp[[cat]]
  cat_line <- sprintf("%s: full=%d, 6gp=%d, shared=%d (%.0f%% of full)",
                      cat, length(f), length(g),
                      length(intersect(f, g)),
                      100 * length(intersect(f, g)) / max(length(f), 1))
  message(cat_line)
}

# Endo_specific: full=34, 6gp=275, shared=9 (26% of full)
# Meso_specific: full=167, 6gp=687, shared=49 (29% of full)
# Ecto_specific: full=144, 6gp=298, shared=34 (24% of full)

category_colours <- c(
  ME            = "#E69F00",
  Meso_specific = "#56B4E9",
  Endo_specific = "#009E73",
  Ecto_specific = "#CC79A7",
  ambiguous     = "grey60",
  constitutive  = "black"
)

# helper: pull this layer's HV (90th pct) and notHV (10th pct) cutoffs
.layer_cut <- function(wide, layer) {
  thr <- attr(wide, "thresholds")
  if (is.null(thr)) stop("wide has no 'thresholds' attribute — run classify_wide() first.")
  c(hv = unname(thr[[layer]]["hv"]), notHV = unname(thr[[layer]]["notHV"]))
}

# ── Scatter plot: layer vs layer coloured by category ────────────────────────
plot_quadrant_layer <- function(wide, subsampling = 100000) {
  
  # shared axis range across all panels, from the data (small pad)
  rng <- range(unlist(wide[, .(endo, meso, ecto)]), na.rm = TRUE)
  rng <- rng + c(-0.02, 0.02) * diff(rng)
  
  make_plot <- function(x_col, y_col, x_lab, y_lab, title) {
    subwide <- wide[sample(.N, min(.N, subsampling))]
    xc <- .layer_cut(wide, x_col)   # per-layer cutoffs for THIS panel's axes
    yc <- .layer_cut(wide, y_col)
    
    ggplot(subwide, aes(x = .data[[x_col]], y = .data[[y_col]],
                        colour = category, shape = category)) +
      geom_point(data = subwide[category == "constitutive"], alpha = 0.3, size = 0.3) +
      geom_point(data = subwide[category == "ambiguous"],    alpha = 0.4, size = 0.5) +
      geom_point(data = subwide[!category %in% c("constitutive","ambiguous")],
                 alpha = 0.4, size = 1) +
      # HV cutoff (top 10%, dashed) and notHV cutoff (bottom 10%, dotted), per layer
      geom_hline(yintercept = yc["hv"],    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_hline(yintercept = yc["notHV"], linetype = "dotted", colour = "grey60", linewidth = 0.3) +
      geom_vline(xintercept = xc["hv"],    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_vline(xintercept = xc["notHV"], linetype = "dotted", colour = "grey60", linewidth = 0.3) +
      scale_x_continuous(limits = rng, name = x_lab) +
      scale_y_continuous(limits = rng, name = y_lab) +
      scale_colour_manual(values = category_colours, drop = FALSE) +
      scale_shape_manual(values  = c(ME=16, Meso_specific=16, Endo_specific=16,
                                     Ecto_specific=16, ambiguous=1, constitutive=4),
                         drop = FALSE) +
      ggtitle(title) + theme_bw(base_size = 11) +
      theme(legend.position = "none")
  }
  
  legend_p <- ggplot(wide[sample(.N, min(.N, 100000))],
                     aes(x = meso, y = endo, colour = category)) +
    geom_point(size = 3) +
    scale_colour_manual(values = category_colours, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_void() + theme(legend.position = "right")
  
  ((make_plot("meso","endo","logBF per ds (meso)","logBF per ds (endo)","Meso vs Endo") |
      make_plot("meso","ecto","logBF per ds (meso)","logBF per ds (ecto)","Meso vs Ecto")) /
      (make_plot("endo","ecto","logBF per ds (endo)","logBF per ds (ecto)","Endo vs Ecto") |
         cowplot::get_legend(legend_p))) +
    plot_layout(widths = c(1,1,1,0.35))
}

set.seed(1234)

## Higher background concordance for endo than meso
ggplot2::ggsave(
  filename = here("B_MultiTissues/dataOut/figures/script04/plot_quadrant_layer_wideFull.pdf"),
  plot     = plot_quadrant_layer(wideFull), width = 7, height = 7)

# ── Overlap full vs 6gp ───────────────────────────────────────────────────────
setkey(wideFull,   name)
setkey(wide6gpall, name)

overlap_summary <- rbindlist(lapply(
  c("ME","Ecto_specific","Endo_specific","Meso_specific","constitutive","ambiguous"),
  function(cat) {
    f <- wideFull[category   == cat, name]
    g <- wide6gpall[category == cat, name]
    o <- length(intersect(f, g))
    data.table(category = cat, n_full = length(f), n_6gp = length(g),
               n_overlap = o,
               pct_of_full = round(100 * o / length(f), 1),
               pct_of_6gp  = round(100 * o / length(g), 1))
  }))
print(overlap_summary)
#         category   n_full    n_6gp n_overlap pct_of_full pct_of_6gp
# 1:            ME   669018   535640    401175        60.0       74.9
# 2: Ecto_specific      144      298        34        23.6       11.4
# 3: Endo_specific       34      275         9        26.5        3.3
# 4: Meso_specific      167      687        49        29.3        7.1
# 5:  constitutive   402815   302166    183476        45.5       60.7
# 6:     ambiguous 20450363 21488395  20257105        99.1       94.3

## Ccl:
# MEs are robust to dataset-size reduction, but the layer-specific markers are
# not — Ecto in particular. NB the layer-specific categories are tiny and defined
# by a SECOND (bottom-10%) threshold, so their overlap is sensitive to the per-
# layer quantile cutoffs shifting between the full and 6gp analyses; treat the
# low pct_of_6gp values as instability of the strict definition, not biology.

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
# Written: layer_specific_and_ME.txt (5345 CpGs total)

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

## !! Transfert to local gitignore/methylation_layerspecific_and_ME.tsv

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
# <char>     <char>     <char>
#   1:               Colon - Endocrine        169       Endo
# 2:             Colon - Macrophages        169       Meso
# 3: Kidney glomerular - Endothelium        176       Meso
# 4: Kidney glomerular - Endothelium        130       Meso
# 5:  Kidney glomerular - Epithelium        130       Endo
# 6:    Kidney glomerular - Podocyte        199       Endo
# 7:    Kidney glomerular - Podocyte        176       Endo
# 8:    Kidney tubular - Endothelium        199       Meso
# 9:     Kidney tubular - Epithelium        130       Endo
# 10:     Kidney tubular - Epithelium        176       Endo

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

## Fetal-style same-layer r: one point per patient, then correlate across patients
compute_same_layer_r_fetalstyle <- function(meth_sub, layer_name,
                                            min_patients = 3, absolute = TRUE) {
  layer_meth <- meth_sub[germ_layer == layer_name]
  eligible <- layer_meth[, .(n = uniqueN(source_tissue_celltype)),
                         by = patient_id][n > 1, patient_id]
  if (length(eligible) < min_patients) return(NULL)
  
  # per patient: all tissue pairs, tissue names ordered so m1/m2 are deterministic
  all_pairs <- rbindlist(lapply(eligible, function(pid) {
    pat  <- layer_meth[patient_id == pid]
    tiss <- sort(unique(pat$source_tissue_celltype))       # deterministic order
    if (length(tiss) < 2) return(NULL)
    rbindlist(lapply(combn(tiss, 2, simplify = FALSE), function(pr) {
      d1 <- pat[source_tissue_celltype == pr[1], .(cpg_site, m1 = methylation)]
      d2 <- pat[source_tissue_celltype == pr[2], .(cpg_site, m2 = methylation)]
      merge(d1, d2, by = "cpg_site")[, patient_id := pid]
    }), fill = TRUE)
  }), fill = TRUE)
  if (is.null(all_pairs) || !nrow(all_pairs)) return(NULL)
  
  # ---- KEY STEP: collapse each patient's pairs to ONE point (equal weighting) --
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
      list(r = abs(r), n_obs = sum(idx))
    } else list(r = NA_real_, n_obs = sum(idx))
  }, by = cpg_site]
}

extract_layer_specific_cpgs <- function(meth_sub,
                                        own_layer, other_layer,
                                        r_high = NULL,
                                        r_low  = NULL,
                                        r_low_cross = NULL,   # <- add this
                                        min_same_obs  = 3,
                                        min_cross_obs = 3)
{
  if (is.null(r_high)) r_high <- thr[layer == own_layer,   r_high]
  if (is.null(r_low))  r_low  <- thr[layer == other_layer, r_low]
  if (is.null(r_low_cross)) r_low_cross <- thr_cross$r_low
  
  r_own   <- compute_same_layer_r_fetalstyle(meth_sub, own_layer)
  r_other <- compute_same_layer_r_fetalstyle(meth_sub, other_layer)
  r_cross <- compute_cross_layer_r(meth_sub, other_layer, own_layer)
  if (is.null(r_own) || is.null(r_other) || is.null(r_cross)) return(NULL)
  
  setnames(r_own,   c("r", "n_obs"), c("r_own",   "n_own"))
  setnames(r_other, c("r", "n_obs"), c("r_other", "n_other"))
  setnames(r_cross, c("r", "n_obs"), c("r_cross", "n_cross"))
  
  m <- Reduce(function(a, b) merge(a, b, by = "cpg_site"),
              list(r_own, r_other, r_cross))
  
  hits <- m[!is.na(r_own) & !is.na(r_other) & !is.na(r_cross) &
              n_own   >= min_same_obs &
              n_other >= min_same_obs &
              n_cross >= min_cross_obs &
              r_own   >= r_high &
              r_other <  r_low  &
              r_cross <  r_low_cross]
  hits[order(-r_own)]
}

set.seed(1234)
constitutive_sample <- sample(wideFull[category == "constitutive", name], 10000)
ambiguous_sample    <- sample(wideFull[category == "ambiguous",    name], 10000)
writeLines(c(constitutive_sample, ambiguous_sample),
           here("B_MultiTissues/dataOut/control_sample20k.txt"))
message("Written: control_sample20k.txt (20000 CpGs: 10k constitutive + 10k ambiguous)")

## In pchuckle (after git pull):
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/control_sample20k.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_control20k.tsv \
# --minCov    10

## !! Transfert gitignore/methylation_control20k.tsv to local

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
# Control CpGs loaded: 12 constitutive, 7 ambiguous

# ══════════════════════════════════════════════════════════════════════════════
# Run for Meso_specific, Endo_specific and controls
# ══════════════════════════════════════════════════════════════════════════════

# ── Meso_specific ─────────────────────────────────────────────────────────────
meth_meso_specific <- meth_multi[category == "Meso_specific"]
r_same_meso        <- compute_same_layer_r_fetalstyle(meth_meso_specific, "Meso")

# ── Endo_specific ─────────────────────────────────────────────────────────────
meth_endo_specific <- meth_multi[category == "Endo_specific"]
r_same_endo        <- compute_same_layer_r_fetalstyle(meth_endo_specific, "Endo")

# ── ME same-layer r ───────────────────────────────────────────────────────────
meth_ME <- meth_multi[category == "ME"]

# ── Wrong-layer specificity checks ───────────────────────────────────────────
r_same_endo_for_meso <- compute_same_layer_r_fetalstyle(meth_meso_specific, "Endo")
r_same_meso_for_endo <- compute_same_layer_r_fetalstyle(meth_endo_specific, "Meso")
r_same_meso_ME <- compute_same_layer_r_fetalstyle(meth_ME, "Meso")
r_same_endo_ME <- compute_same_layer_r_fetalstyle(meth_ME, "Endo")

# ── Controls in both Meso and Endo contexts ───────────────────────────────────
r_same_meso_const      <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "constitutive"], "Meso")
r_same_endo_const <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "constitutive"], "Endo")
r_same_meso_amb        <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "ambiguous"],    "Meso")
r_same_endo_amb   <- compute_same_layer_r_fetalstyle(meth_control_multi[category == "ambiguous"],    "Endo")

# ── Extract "true layer-specific" CpGs on THREE criteria ──────────────────────
#   1. high  same-layer intra-individual r in the OWN layer   (r_own   >= r_high)
#   2. low   same-layer intra-individual r in the OTHER layer (r_other <  r_low)
#   3. low   cross-layer within-person r (own <-> other)       (r_cross <  r_low)

# per-layer empirical null from constitutive CpGs
null_meso <- r_same_meso_const[!is.na(r), r]
null_endo <- r_same_endo_const[!is.na(r), r]

r_crit <- function(n, alpha = 0.05) {
  tc <- qt(1 - alpha/2, df = n - 2)
  tc / sqrt(tc^2 + n - 2)
}

thr <- data.table(
  layer  = c("Meso", "Endo"),
  # "high" = beyond the 95th percentile of noise (5% FPR by construction)
  r_high = c(quantile(null_meso, 0.95), quantile(null_endo, 0.95)),
  # "low"  = below the median of noise (indistinguishable from chance)
  r_low  = c(median(null_meso),         median(null_endo)),
  # mean +/- sd of the noise, for reference
  mu     = c(mean(null_meso),   mean(null_endo)),
  sd     = c(sd(null_meso),     sd(null_endo)),
  n_pat  = c(9, 15),
  r_crit = c(r_crit(9), r_crit(15)))
print(thr)
# layer    r_high     r_low        mu        sd n_pat    r_crit
# <char>     <num>     <num>     <num>     <num> <num>     <num>
# 1:   Meso 0.6810714 0.2506185 0.2835287 0.2055916     9 0.6663836
# 2:   Endo 0.6955997 0.1879514 0.2531251 0.2448413    15 0.5139775

# ── Summary function ──────────────────────────────────────────────────────────
make_summary_samelayer <- function(r_same, layer_tested, category_name,
                                   thr_tab = thr) {
  s <- r_same[!is.na(r)]
  r_hi <- thr_tab[layer == layer_tested, r_high]      # layer-specific "high" bar
  data.table(
    category          = category_name,
    same_layer_tested = layer_tested,
    r_high_used       = r_hi,
    n_cpgs_tested     = uniqueN(s$cpg_site),
    pct_high_same     = 100 * mean(s$r >= r_hi),
    mean_r            = mean(s$r),
    median_r          = median(s$r))
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
#         category same_layer_tested r_high_used n_cpgs_tested pct_high_same    mean_r  median_r
# <char>            <char>       <num>         <int>         <num>     <num>     <num>
# 1: Meso_specific              Meso   0.6810714            30     83.333333 0.7890745 0.8555211
# 2: Meso_specific              Endo   0.6955997            30      0.000000 0.2577406 0.1955097
# 3: Endo_specific              Endo   0.6955997             3     66.666667 0.8434886 0.8871436
# 4: Endo_specific              Meso   0.6810714             3      0.000000 0.4409972 0.4423979
# 5:  constitutive              Meso   0.6810714            12      8.333333 0.2835287 0.2506185
# 6:  constitutive              Endo   0.6955997            12      8.333333 0.2531251 0.1879514
# 7:     ambiguous              Meso   0.6810714             7     71.428571 0.6753872 0.7007080
# 8:     ambiguous              Endo   0.6955997             7     28.571429 0.4063682 0.1384330

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

# empirical FPR: % of constitutive CpGs passing the same threshold, per layer
noise_dt <- results_summary[category == "constitutive",
                            .(same_layer_tested, floor_pct = pct_high_same)]

p1 <- ggplot(results_summary,
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
  geom_hline(data = noise_dt, aes(yintercept = floor_pct),
             linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_text(data = noise_dt,
            aes(x = 0.6, y = floor_pct, label = sprintf("noise floor", floor_pct)),
            hjust = 0, vjust = -0.4, size = 2.6, colour = "grey30", inherit.aes = FALSE) +
  scale_fill_manual(values = bar_colours, name = "Layer context") +
  scale_y_continuous("% CpGs above the layer-specific |r| cutoff", limits = c(0, NA)) +
  scale_x_discrete("CpG category") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        axis.text.x       = element_text(angle = 30, hjust = 1),
        strip.text        = element_text(face = "bold"),
        legend.position   = "right") +
  ggtitle("Intra-individual same-layer methylation concordance",
          subtitle = sprintf("Own layer (blue) should exceed the noise floor by more than the wrong layer (orange)\nMeso: %d patients (|r| \u2265 %.2f) | Endo: %d patients (|r| \u2265 %.2f)",
                             9, thr[layer=="Meso", r_high], 15, thr[layer=="Endo", r_high]))

# ── Density plot of same-layer r per category and layer context ───────────────

# ── Pool every category × layer tested into one table ────────────────────────
density_dt <- rbindlist(list(
  # own-layer
  r_same_meso[!is.na(r),          .(cpg_site, r, category = "Meso_specific", layer_context = "Meso")],
  r_same_endo[!is.na(r),          .(cpg_site, r, category = "Endo_specific", layer_context = "Endo")],
  # wrong-layer
  r_same_endo_for_meso[!is.na(r), .(cpg_site, r, category = "Meso_specific", layer_context = "Endo")],
  r_same_meso_for_endo[!is.na(r), .(cpg_site, r, category = "Endo_specific", layer_context = "Meso")],
  # MEs (expected high in BOTH layers)
  r_same_meso_ME[!is.na(r),       .(cpg_site, r, category = "ME",           layer_context = "Meso")],
  r_same_endo_ME[!is.na(r),       .(cpg_site, r, category = "ME",           layer_context = "Endo")],
  # controls
  r_same_meso_const[!is.na(r),    .(cpg_site, r, category = "constitutive", layer_context = "Meso")],
  r_same_endo_const[!is.na(r),    .(cpg_site, r, category = "constitutive", layer_context = "Endo")],
  r_same_meso_amb[!is.na(r),      .(cpg_site, r, category = "ambiguous",    layer_context = "Meso")],
  r_same_endo_amb[!is.na(r),      .(cpg_site, r, category = "ambiguous",    layer_context = "Endo")]
))

density_dt[, category_f := factor(category,
                                  levels = c("ME", "Endo_specific", "Meso_specific", "ambiguous", "constitutive"))]
density_dt[, layer_f := factor(layer_context, levels = c("Meso", "Endo"),
                               labels = c("Tested in Meso tissues", "Tested in Endo tissues"))]

# n per group for the subtitle
n_lab <- density_dt[, .(n = uniqueN(cpg_site)), by = .(category_f, layer_f)]

null_dt <- data.table(
  layer_f = factor(c("Tested in Meso tissues","Tested in Endo tissues"),
                   levels = levels(density_dt$layer_f)),
  null_r  = c(thr[layer=="Meso", mu], thr[layer=="Endo", mu]))   # empirical

thr_dt <- data.table(
  layer_f = factor(c("Tested in Meso tissues","Tested in Endo tissues"),
                   levels = levels(density_dt$layer_f)),
  r_high  = c(thr[layer=="Meso", r_high], thr[layer=="Endo", r_high]))

p2 <- ggplot(density_dt, aes(x = r, colour = category_f, fill = category_f)) +
  geom_density(alpha = 0.20, linewidth = 0.8, bounds = c(0, 1)) +
  geom_vline(data = thr_dt, aes(xintercept = r_high),
             linetype = "dashed", colour = "grey40") +
  geom_vline(data = null_dt, aes(xintercept = null_r),
             linetype = "dotted", colour = "grey50") +
  facet_wrap(~ layer_f, nrow = 1) +
  scale_colour_manual(values = category_colours, drop = FALSE, name = "CpG category") +
  scale_fill_manual(values   = category_colours, drop = FALSE, name = "CpG category") +
  scale_x_continuous("Same-layer |r| (per-patient means, across patients)",
                     limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous("Density") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "right") +
  ggtitle("Same-layer methylation concordance by CpG category (fetal-style |r|)",
          subtitle = sprintf(paste0(
            "Each patient contributes one point (mean of their tissue pairs); |Pearson r| taken across patients\n",
            "Dotted = mean |r| of constitutive CpGs (empirical noise) | ",
            "Dashed = 95th pct of that noise, used as the \"high\" cutoff\n",
            "Meso: %d patients, noise %.2f, cutoff |r| \u2265 %.2f  |  ",
            "Endo: %d patients, noise %.2f, cutoff |r| \u2265 %.2f"),
            9,  thr[layer == "Meso", mu], thr[layer == "Meso", r_high],
            15, thr[layer == "Endo", mu], thr[layer == "Endo", r_high]))

# ── toggle: exclude MEs from the density panel ───────────────────────────────
drop_ME <- TRUE
plot_dt <- if (drop_ME) density_dt[category != "ME"] else density_dt

p3 <- ggplot(plot_dt, aes(x = r, colour = category_f, fill = category_f)) +
  geom_density(alpha = 0.20, linewidth = 0.8, bounds = c(0, 1)) +
  geom_vline(data = thr_dt,  aes(xintercept = r_high),
             linetype = "dashed", colour = "grey40") +
  geom_vline(data = null_dt, aes(xintercept = null_r),
             linetype = "dotted", colour = "grey50") +
  facet_wrap(~ layer_f, nrow = 1) +
  scale_colour_manual(values = category_colours, drop = TRUE, name = "CpG category") +
  scale_fill_manual(values   = category_colours, drop = TRUE, name = "CpG category") +
  scale_x_continuous("Same-layer |r| (per-patient means, then |Pearson r| across patients)",
                     limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous("Density") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "right") +
  ggtitle("Same-layer methylation concordance by CpG category (fetal-style |r|)",
          subtitle = sprintf(paste0(
            "Each patient contributes one point (mean of their tissue pairs); |Pearson r| taken across patients\n",
            "Dotted = mean |r| of constitutive CpGs (empirical noise) | ",
            "Dashed = 95th pct of that noise, used as the \"high\" cutoff\n",
            "Meso: %d patients, noise %.2f, cutoff |r| \u2265 %.2f  |  ",
            "Endo: %d patients, noise %.2f, cutoff |r| \u2265 %.2f"),
            9,  thr[layer == "Meso", mu], thr[layer == "Meso", r_high],
            15, thr[layer == "Endo", mu], thr[layer == "Endo", r_high]))

ggplot2::ggsave(
  here("B_MultiTissues/dataOut/figures/script04/sameLayercR_concordance.pdf"),
  plot_grid(p1, p2, p3, labels = c("A", "B", "C"), nrow = 3), width = 10, height = 15)

# Meso_specific CpGs are highly variable between individuals (by selection)
# Within each individual, their methylation level is highly consistent across all blood cell types (r≈1, Plot 1)
# This within-individual concordance across independent blood lineages is the hallmark of pre-haematopoietic stochastic establishment
# Controls show r≈0, confirming this is not a general property of blood measurements

###### Add now the 3rd criteria, low r between both layers

min_same_obs  <- 3
min_cross_obs <- 3

# per-CpG: the three r's for one category under one pipeline (own vs other) ------
tri_rs <- function(meth_sub, own_layer, other_layer, category_name) {
  ro <- compute_same_layer_r_fetalstyle(meth_sub, own_layer)
  rt <- compute_same_layer_r_fetalstyle(meth_sub, other_layer)
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
  own <- tab$own_layer[1]
  oth <- setdiff(c("Meso","Endo"), own)
  rh  <- thr[layer == own, r_high]
  rlo <- thr[layer == oth, r_low]
  rlc <- thr_cross$r_low
  tested <- tab[!is.na(r_own) & n_own >= min_same_obs]
  N  <- nrow(tested)
  p1 <- tested[r_own >= rh]
  p2 <- p1[!is.na(r_other) & n_other >= min_same_obs & r_other < rlo]
  p3 <- p2[!is.na(r_cross) & n_cross >= min_cross_obs & r_cross < rlc]
  data.table(
    category  = tab$category[1],
    own_layer = tab$own_layer[1],
    step = factor(c("1. high own r", "2. + low other r", "3. + low cross r"),
                  levels = c("1. high own r", "2. + low other r", "3. + low cross r")),
    n   = c(nrow(p1), nrow(p2), nrow(p3)),
    pct = 100 * c(nrow(p1), nrow(p2), nrow(p3)) / N,
    n_tested = N)
}

# ── Empirical null for the CROSS-layer r (Endo x Meso, n = 4 patients) ────────
null_cross <- compute_cross_layer_r(
  meth_control_multi[category == "constitutive"], "Endo", "Meso")[!is.na(r), r]

thr_cross <- data.table(
  comparison = "Endo_x_Meso",
  n_pat      = length(unique(meth_control_multi[germ_layer %in% c("Endo","Meso"), patient_id])),
  # "high" (ME-like systemic concordance): beyond 95th pct of noise
  r_high     = quantile(null_cross, 0.95),
  # "low" (layer-specific): below the MEDIAN of noise, i.e. no more concordant
  #  than a constitutive CpG. Tighten to 0.25 quantile for a stricter set.
  r_low      = median(null_cross),
  r_low_strict = quantile(null_cross, 0.25),
  mu         = mean(null_cross),
  sd         = sd(null_cross))
print(thr_cross)
# comparison n_pat    r_high     r_low r_low_strict        mu        sd
# <char> <int>     <num>     <num>        <num>     <num>     <num>
# 1: Endo_x_Meso    26 0.8926701 0.4055105    0.1161466 0.3809395 0.3040298

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

p <- ggplot(tri_tabs, aes(category_f, pct, fill = step)) +
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
          subtitle = sprintf(
            "high own r (Meso \u2265%.2f, Endo \u2265%.2f)  \u2192  + low other r (Meso <%.2f, Endo <%.2f)  \u2192  + low cross-layer r (<%.2f)",
            thr[layer == "Meso", r_high], thr[layer == "Endo", r_high],
            thr[layer == "Meso", r_low],  thr[layer == "Endo", r_low],
            thr_cross$r_low))

## Higherbackground concordance for endo than meso
ggplot2::ggsave(
  filename = here("B_MultiTissues/dataOut/figures/script04/layerCandidatesSelection.pdf"),
  plot     = p, width    = 8, height = 5)

#################################
## Putative MEs (systemic)     ##
#################################
# ME = high Pr(HV) in all 3 layers (already the "ME" category) PLUS the systemic
# concordance signature: high intra-layer r in each layer AND high cross-layer r.
# (Layer-specific = high own / low cross; ME = high own / HIGH cross.)

extract_ME_cpgs <- function(meth_sub,
                            r_high        = 0.5,   # concordance bar (same as own-layer)
                            min_same_obs  = 3,
                            min_cross_obs = 3) {
  # within-individual concordance in each layer that has multi-tissue patients
  r_endo <- compute_same_layer_r_fetalstyle(meth_sub, "Endo")
  r_meso <- compute_same_layer_r_fetalstyle(meth_sub, "Meso")
  # cross-layer (systemic) concordance — the ME hallmark
  r_cross <- compute_cross_layer_r(meth_sub, "Endo", "Meso")
  if (is.null(r_endo) || is.null(r_meso) || is.null(r_cross)) return(NULL)
  
  setnames(r_endo,  c("r","n_obs"), c("r_endo",  "n_endo"))
  setnames(r_meso,  c("r","n_obs"), c("r_meso",  "n_meso"))
  setnames(r_cross, c("r","n_obs"), c("r_cross", "n_cross"))
  
  m <- Reduce(function(a,b) merge(a,b,by="cpg_site"),
              list(r_endo, r_meso, r_cross))
  
  hits <- m[!is.na(r_endo) & !is.na(r_meso) & !is.na(r_cross) &
              n_endo  >= min_same_obs &
              n_meso  >= min_same_obs &
              n_cross >= min_cross_obs &
              r_endo  >= r_high &      # high concordance in endo
              r_meso  >= r_high &      # high concordance in meso
              r_cross >= r_high]       # AND high cross-layer  -> systemic
  hits[order(-r_cross)]
}

ME_hits <- extract_ME_cpgs(meth_ME)
message(sprintf("Putative systemic-ME CpGs: %d", if (is.null(ME_hits)) 0L else nrow(ME_hits)))
# Putative systemic-ME CpGs: 4098

#################################
## Save the interesting targets##
#################################
meso_hits <- extract_layer_specific_cpgs(meth_meso_specific, 
                                         own_layer = "Meso", other_layer = "Endo")
endo_hits <- extract_layer_specific_cpgs(meth_endo_specific,
                                         own_layer = "Endo", other_layer = "Meso")

message(sprintf("Meso_specific true-signature CpGs: %d", nrow(meso_hits)))
## Meso_specific true-signature CpGs: 3

message(sprintf("Endo_specific true-signature CpGs: %d", nrow(endo_hits)))
# Endo_specific true-signature CpGs: 0

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

for (g in c(50, 500, 1000, 2000)) {
  cl <- annotate_layer_hits(endo_hits, "endo_tmp", genes_gr, uni, uni_gr,
                            gap = g, out_dir = tempdir())$clusters
  message(sprintf("gap=%4d bp: %d clusters (>=2 hits), %f genes, max n_hits=%f",
                  g, nrow(cl), uniqueN(cl$gene), max(cl$n_hits)))
}

meso_res <- annotate_layer_hits(meso_hits, "meso", genes_gr, uni, uni_gr, gap = 100)
endo_res <- annotate_layer_hits(endo_hits, "endo", genes_gr, uni, uni_gr, gap = 100)
ME_res   <- annotate_layer_hits(ME_hits,   "ME",   genes_gr, uni, uni_gr, gap = 100)

meso_res
ME_res

