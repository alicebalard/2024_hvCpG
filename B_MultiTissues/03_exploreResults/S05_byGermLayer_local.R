#####################################################################
# S05 — Per-germ-layer hypervariability categories
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

# category_colours <- c(
#   "top1pcEctoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[1],
#   "top1pcEndoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[2],
#   "top1pcMesoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[3],
#   "top1pc3layers" = "orange",
#   "constitutive" = "grey80",
#   "background" = "black"
# )

layers_colours <- c(
  "ecto" = paletteer::paletteer_d("nationalparkcolors::Everglades")[1],
  "endo" = paletteer::paletteer_d("nationalparkcolors::Everglades")[2],
  "meso" = paletteer::paletteer_d("nationalparkcolors::Everglades")[3]
)

## ==============================================================================
## Definition of layer-specific HV CpGs
## Top 1% in one layer and bottom 50% in the other two layers
## ==============================================================================

q_top <- function(x) {
  quantile(x, probs = 0.99, na.rm = TRUE, names = FALSE)
}

q_bot <- function(x) {
  quantile(x, probs = 0.50, na.rm = TRUE, names = FALSE)
}

top_endo  <- q_top(dt$logBF_per_ds_endo)
top_meso  <- q_top(dt$logBF_per_ds_meso)
top_ecto  <- q_top(dt$logBF_per_ds_ecto)

bot_endo  <- q_bot(dt$logBF_per_ds_endo)
bot_meso  <- q_bot(dt$logBF_per_ds_meso)
bot_ecto  <- q_bot(dt$logBF_per_ds_ecto)

dt[, `:=`(
  hv_none =
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_meso <= bot_meso &
    logBF_per_ds_ecto <= bot_ecto,
  hv_endo_only =
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso <= bot_meso &
    logBF_per_ds_ecto <= bot_ecto,
  hv_meso_only =
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_ecto <= bot_ecto,
  hv_ecto_only =
    logBF_per_ds_ecto >= top_ecto &
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_meso <= bot_meso,
  hv_endomeso =
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_ecto <= bot_ecto,
  hv_endoecto = 
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso <= bot_meso &
    logBF_per_ds_ecto >= top_ecto,
  hv_mesoecto = 
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_ecto >= top_ecto,
  hv_mesoectoendo = 
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_ecto >= top_ecto
)]

hv_counts <- data.table(
  category = c(
    "In no layer",
    "Mesoderm only",
    "Endoderm only",
    "Ectoderm only",
    "Endoderm + mesoderm",
    "Endoderm + ectoderm",
    "Mesoderm + ectoderm",
    "Endoderm + mesoderm + ectoderm"
  ),
  variable = c(
    "hv_none",
    "hv_meso_only",
    "hv_endo_only",
    "hv_ecto_only",
    "hv_endomeso",
    "hv_endoecto",
    "hv_mesoecto",
    "hv_mesoectoendo"
  ),
  n = c(
    sum(dt$hv_none, na.rm = TRUE),
    sum(dt$hv_meso_only, na.rm = TRUE),
    sum(dt$hv_endo_only, na.rm = TRUE),
    sum(dt$hv_ecto_only, na.rm = TRUE),
    sum(dt$hv_endomeso, na.rm = TRUE),
    sum(dt$hv_endoecto, na.rm = TRUE),
    sum(dt$hv_mesoecto, na.rm = TRUE),
    sum(dt$hv_mesoectoendo, na.rm = TRUE)
  )
)

hv_counts
#                          category        variable       n
# 1:                    In no layer         hv_none 6723418
# 2:                  Mesoderm only    hv_meso_only      63
# 3:                  Endoderm only    hv_endo_only      18
# 4:                  Ectoderm only    hv_ecto_only     297
# 5:            Endoderm + mesoderm     hv_endomeso      39
# 6:            Endoderm + ectoderm     hv_endoecto      18
# 7:            Mesoderm + ectoderm     hv_mesoecto       3
# 8: Endoderm + mesoderm + ectoderm hv_mesoectoendo   60424 

## =============================================================================
## Test overlap between the top in each layer
## =============================================================================

if (!file.exists(here("B_MultiTissues/dataOut/figures/script05/venn_scatter_germLayers.pdf"))){
  
  ## 1. Pairwise scatter plots — shared fill scale + one common legend
  pairs <- list(c("logBF_per_ds_endo","logBF_per_ds_meso","Hypervariability score endoderm","Hypervariability score mesoderm"),
                c("logBF_per_ds_endo","logBF_per_ds_ecto","Hypervariability score endoderm","Hypervariability score ectoderm"),
                c("logBF_per_ds_meso","logBF_per_ds_ecto","Hypervariability score mesoderm","Hypervariability score ectoderm"))
  
  make_scatter <- function(p, lims = NULL)
    ggplot(dt, aes(.data[[p[1]]], .data[[p[2]]])) +
    geom_hex(bins = 1000) +
    scale_fill_viridis_c(trans = "log10", name = "CpGs", limits = lims) +
    labs(x = p[3], y = p[4]) +
    theme_minimal(base_size = 12)
  
  # common colour limits so ONE legend is valid for all three
  counts <- unlist(lapply(pairs, function(p)
    ggplot2::ggplot_build(make_scatter(p))$data[[1]]$count))
  lims <- c(1, max(counts, na.rm = TRUE))
  
  sc <- lapply(pairs, make_scatter, lims = lims)
  
  # pull the legend from the first, then drop legends from all
  legend <- cowplot::get_legend(sc[[1]] + theme(legend.position = "right"))
  sc <- lapply(sc, function(g) g + theme(legend.position = "none"))
  
  scatter_row <- cowplot::plot_grid(plotlist = sc, labels = c("A","B","C"), nrow = 1)
  scatterplots <- cowplot::plot_grid(scatter_row, legend, nrow = 1, rel_widths = c(1, 0.12))
  
  ## 2. for the CpGs that are top-1% in layer X, where do they sit in the other layers' rank distributions?
  
  # within-layer percentile rank (0–1), comparable across layers
  pr <- function(x) data.table::frank(x, na.last = "keep", ties.method = "average") /
    sum(!is.na(x))
  dt[, `:=`(pr_endo = pr(logBF_per_ds_endo),
            pr_meso = pr(logBF_per_ds_meso),
            pr_ecto = pr(logBF_per_ds_ecto))]
  
  layers <- c("endo","meso","ecto")
  
  makeridgeplotGermLayers <- function(thr){
    long <- rbindlist(lapply(layers, function(f) {
      foc <- dt[get(paste0("pr_", f)) >= 1-thr/100] # top thr% in focal layer
      rbindlist(lapply(setdiff(layers, f), function(o)
        data.table(focal = f, other = o,
                   pr_other = foc[[paste0("pr_", o)]])))
    }))
    
    long[, panel := paste0("top", thr, "% ", toupper(focal), " → rank in ", toupper(other))]
    
    ggplot(long, aes(pr_other, panel, fill = focal)) +
      geom_density_ridges(scale = 1.3, alpha = .7, rel_min_height = 0.01,
                          quantile_lines = F) +
      geom_vline(xintercept = 0.99, linetype = 2, colour = "grey40") +
      scale_x_continuous("Percentile rank in the OTHER layer", labels = scales::percent) +
      scale_fill_manual(values = layers_colours)+
      theme_ridges(font_size = 14)+
      theme(legend.position = "none", axis.title.y = element_blank())
    ## title: Where each layer's top", thr, "% hvCpGs fall in the other layers
  }
  
  p1 <- makeridgeplotGermLayers(1)
  p2 <- makeridgeplotGermLayers(5)
  p3 <- makeridgeplotGermLayers(10)
  
  ridgeplots <- cowplot::plot_grid(p1, p2, p3,
                                 labels = c("D", "E", "F"), nrow = 1)
  
  ggsave(here("B_MultiTissues/dataOut/figures/script05/venn_scatter_germLayers.pdf"),
         cowplot::plot_grid(scatterplots, ridgeplots, nrow = 2),
         width = 16, height = 8, dpi = 300, bg = "white")
}

################################################################################
## Extract raw methylation for our 4 main categories                         ##
## (top1pc3layers, top1pcEndoOnly, top1pcMesoOnly, constitutive)              ##
## CpGpos stored in object "categories"                                       ##
################################################################################

## We have 4 samples with endoderm and mesoderm cells
message("Category sizes:")
print(sapply(categories, length))
# top1pc3layers   constitutive top1pcEndoOnly top1pcMesoOnly top1pcEctoOnly 
# 60424        6723417           9479           9433          41406

# ── Write CpG lists for python extraction ─────────────────────────────────────

set.seed(1234)
## top1pc3layers CpGs — subsample 10000 as positive control
top1pc3layers_sample <- sample(categories$top1pc3layers, min(10000, length(categories$top1pc3layers)))

## constitutive_sample — subsample 10000 as negative control
constitutive_sample <- sample(categories$constitutive, min(10000, length(categories$constitutive)))

## Neutral (variance-unfiltered) background
neutral_sample <- sample(dt$chr_pos, min(20000, nrow(dt)))
## Exclude other categories
neutral_sample <- neutral_sample[!neutral_sample %in% 
                                   c(categories$top1pcEctoOnly, categories$top1pcMesoOnly,
                                     categories$top1pcEndoOnly, categories$top1pc3layers,
                                     categories$constitutive)]

neutral_sample <- sample(neutral_sample, 10000)

## Sanity check: each CpG must belong to ONE layer
sets <- list(
  top1pcEndoOnly = categories$top1pcEndoOnly,
  top1pcMesoOnly = categories$top1pcMesoOnly,
  top1pcEctoOnly = categories$top1pcEctoOnly,
  top1pc3layers  = top1pc3layers_sample,
  constitutive   = constitutive_sample,
  neutral        = neutral_sample
)

# one row per (cpg, set)
dt_sets <- rbindlist(lapply(names(sets), function(s)
  data.table(cpg = sets[[s]], set = s)))

# how many DISTINCT sets does each CpG appear in?
overlap <- dt_sets[, .(n_sets = uniqueN(set),
                       which_sets = paste(sort(unique(set)), collapse = " + ")),
                   by = cpg][n_sets > 1]

if (nrow(overlap) == 0) {
  message("✓ All CpGs are unique to a single set.")
} else {
  message("✗ ", nrow(overlap), " CpGs appear in more than one set:")
  print(overlap[, .N, by = which_sets][order(-N)])   # which set-pairs overlap, and how many
}
# ✓ All CpGs are unique to a single set.

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
# Written: CpG2extractS05.txt (90318 CpGs total)

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
rm(meth_CpG2extractS05)

##############################################
## Add categories to meth_CpG2extractS05_dedup
cpg_category <- rbindlist(lapply(names(categories), function(cat)
  data.table(cpg_site = categories[[cat]], category = cat)))

## Add neutral ones
cpg_category <- rbind(cpg_category,
                      data.table(cpg_site = neutral_sample, category = "background"))

meth_CpG2extractS05_dedup <- merge(meth_CpG2extractS05_dedup, cpg_category, 
                                   by = "cpg_site", all.x = TRUE)

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

## Add categories
em_r <- merge(em_r, cpg_category, by = "cpg_site", all.x = TRUE)

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
# Warning in dcast.data.table(mm, cpg_site + patient_id ~ slot, value.var = "methylation") :
#   'fun.aggregate' is NULL, but found duplicate row/column combinations, so defaulting to length(). That is, the variables [cpg_site, patient_id, slot] used in 'formula' do not uniquely identify rows in the input 'data'. In such cases, 'fun.aggregate' is used to derive a single representative value for each combination in the output data.table, for example by summing or averaging (fun.aggregate=sum or fun.aggregate=mean, respectively). Check the resulting table for values larger than 1 to see which combinations were not unique. See ?dcast.data.table for more details.

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

## categories 
ee_r <- merge(ee_r, cpg_category, by = "cpg_site", all.x = TRUE)

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
# Warning in dcast.data.table(mm, cpg_site + patient_id ~ slot, value.var = "methylation") :
#   'fun.aggregate' is NULL, but found duplicate row/column combinations, so defaulting to length(). That is, the variables [cpg_site, patient_id, slot] used in 'formula' do not uniquely identify rows in the input 'data'. In such cases, 'fun.aggregate' is used to derive a single representative value for each combination in the output data.table, for example by summing or averaging (fun.aggregate=sum or fun.aggregate=mean, respectively). Check the resulting table for values larger than 1 to see which combinations were not unique. See ?dcast.data.table for more details.

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

## categories
mm_r <- merge(mm_r, cpg_category, by = "cpg_site", all.x = TRUE)

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
# thresholds  EE high=0.42 low=0.21 | MM high=0.51 low=0.29 | EM hig=0.80 low=0.51

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

length(endo_candidates) # 351
length(meso_candidates) # 368

################################################################################
## Plot inter-germ layer (endo-meso, Nind=4) correlation for the candidates


em_r_candidates <- em_r[em_r$cpg_site %in% c(endo_candidates, meso_candidates) | 
                          em_r$category %in% c("constitutive", "background", "top1pc3layers"),]

p_mesoendocor_violin <- ggplot(em_r_candidates, 
                               aes(x=category, y=abs_r, group = category, fill = category))+
  geom_violin(width=2) +
  geom_boxplot(width=0.1, color="black", fill = "white") +
  scale_fill_manual(values = category_colours) +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0,1)) +
  geom_text(data = em_r_candidates[, .(nCpG = .N), by = category],
            aes(x = category, y = 1, label = format(nCpG, big.mark = ",")),
            vjust = -0.4, size = 3, inherit.aes = FALSE) +
  labs(y = "Mesoderm-endoderm correlation\n(Pearson's r)")+
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p_mesoendocor_violin


ee_r_candidates <- ee_r[ee_r$cpg_site %in% c(endo_candidates, meso_candidates) | 
                          ee_r$category %in% c("constitutive", "background", "top1pc3layers"),]

p_endoendocor_violin <- ggplot(ee_r_candidates,
                               aes(x = category, y = abs_r, group = category, fill = category)) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.1, color = "black", fill = "white") +
  scale_fill_manual(values = category_colours) +
  theme_minimal(base_size = 14) +
  geom_text(data = ee_r_candidates[, .(nCpG = .N), by = category],
            aes(x = category, y = 1, label = format(nCpG, big.mark = ",")),
            vjust = -0.4, size = 3, inherit.aes = FALSE) +
  labs(y = "Endoderm-endoderm correlation\n(Pearson's r)") +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p_endoendocor_violin

mm_r_candidates <- mm_r[mm_r$cpg_site %in% c(endo_candidates, meso_candidates) | 
                          mm_r$category %in% c("constitutive", "background", "top1pc3layers"),]

p_mesomesocor_violin <- ggplot(mm_r_candidates, aes(x = category, y = abs_r, group = category, fill = category)) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.1, color = "black", fill = "white") +
  scale_fill_manual(values = category_colours) +
  theme_minimal(base_size = 14) +
  geom_text(data = mm_r_candidates[, .(nCpG = .N), by = category],
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
em_candidates <- em[em$cpg_site %in% c(endo_candidates, meso_candidates),]
Pa <- makeIVbyICplot(em_candidates, em_r_candidates)

## 2. meso-meso
mm_candidates <- mm[mm$cpg_site %in% c(endo_candidates, meso_candidates),]
Pb <- makeIVbyICplot(mm_candidates, mm_r_candidates)

## 3. endo-endo
ee_candidates <- ee[ee$cpg_site %in% c(endo_candidates, meso_candidates),]
Pc <- makeIVbyICplot(ee_candidates, ee_r_candidates)

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

# universe of tested CpGs = all covered-in-3 CpGs NOT in the focal category
uni <- data.table(cpg_site = dt$chr_pos)
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

meso_res <- annotate_layer_hits(meso_hits, "meso", genes_gr, uni, uni_gr, min_hits = 1,
                                tss_gr, prom_gr, exon_gr, intron_gr, te_regions, gap = 1000)
# $clusters
#       gene cluster_id n_hits  span    density
# 1:  ATP8A2         87      3   245 0.01219512
# 2: COL21A1        255      2   128 0.01550388
# 3:   MYH14        142      2   146 0.01360544

endo_res <- annotate_layer_hits(endo_hits, "endo", genes_gr, uni, uni_gr, min_hits = 1,
                                tss_gr, prom_gr, exon_gr, intron_gr, te_regions, gap = 1000)
# $clusters
#       gene cluster_id n_hits  span     density
# 1: MIR1324        207      2    73 0.027027027
# 2:   SHTN1         38      2   245 0.008130081

saveRDS(list(meso = meso_res, endo = endo_res),
        here(paste0("B_MultiTissues/dataOut/S05_annotatedHits_", variant, ".rds")))

## Save positions as tsv list
fwrite(data.table(chr = meso_res$ann$chr, pos = meso_res$ann$pos),
       file = here(paste0("B_MultiTissues/dataOut/S05_meso_res_annotation.tsv")), sep = "\t", quote = FALSE, na = "NA")

fwrite(data.table(chr = endo_res$ann$chr, pos = endo_res$ann$pos),
       file = here(paste0("B_MultiTissues/dataOut/S05_endo_res_annotation.tsv")), sep = "\t", quote = FALSE, na = "NA")

# ── Combined pretty table (guard against either set being empty) ──────────────
hits_tab <- rbindlist(list(
  if (!is.null(meso_res)) meso_res$ann[, layer := "Mesoderm"],
  if (!is.null(endo_res)) endo_res$ann[, layer := "Endoderm"]
), fill = TRUE)

## subselect clusters
hits_tab_1CpGs <- hits_tab
hits_tab_2CpGs <- hits_tab[hits_tab$gene %in% c(meso_res$clusters$gene, endo_res$clusters$gene),]

makeTab <- function(nCpG){
  if (nCpG == 1){
    HT = hits_tab_1CpGs
    HT <- HT[, .(
      layer, cpg_site, gene, feature,
      dist_to_tss_signed,
      TE = fifelse(in_TE, paste0(TE_class, " / ", TE_family), "\u2014"),
      r_own, r_other, r_cross
    )][order(layer, -r_own)][, head(.SD, 10), by = layer] ## top 10 per layer
  } else if (nCpG == 2){
    HT = hits_tab_2CpGs
    HT <- HT[, .(
      layer, cpg_site, gene, feature,
      dist_to_tss_signed,
      TE = fifelse(in_TE, paste0(TE_class, " / ", TE_family), "\u2014"),
      r_own, r_other, r_cross
    )][order(layer, -cpg_site)]
  }
  
  gt_hits <- HT %>%
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
  
}

makeTab(nCpG = 1)
makeTab(nCpG = 2)

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
  plot_grid(p_row3, p_mesoendocor_violin, p_mesomesocor_violin, p_endoendocor_violin, 
            ncol = 4, labels = c("A", "B", "C", "D")),
  plot_grid(Pa_with_legend, Pb$plot, Pc$plot, ncol = 3, labels = c("D", "E", "F")),
  nrow = 2)

ggsave(here("B_MultiTissues/dataOut/figures/script05/layersPlot.pdf"),
       layersPlot, width = 22, height = 8, dpi = 300, bg = "white")



# dt$logBF_per_ds_endo
# 
# 
# dt$endo <- dt$logBF_per_ds_endo - rowMeans(dt[, .(logBF_per_ds_ecto, logBF_per_ds_meso)], 
#                                            na.rm = TRUE)
# dt$ecto <- dt$logBF_per_ds_ecto - rowMeans(dt[, .(logBF_per_ds_endo, logBF_per_ds_meso)], 
#                                            na.rm = TRUE)
# dt$meso <- dt$logBF_per_ds_meso - rowMeans(dt[, .(logBF_per_ds_ecto, logBF_per_ds_endo)], 
#                                            na.rm = TRUE)
# 
# ## 6 gp
# dt$endo6gp <- dt$logBF_per_ds_endo6gp - rowMeans(dt[, .(logBF_per_ds_ecto, logBF_per_ds_meso6gp)], 
#                                                  na.rm = TRUE)
# dt$ecto6gp <- dt$logBF_per_ds_ecto - rowMeans(dt[, .(logBF_per_ds_endo6gp, logBF_per_ds_meso6gp)], 
#                                               na.rm = TRUE)
# dt$meso6gp <- dt$logBF_per_ds_meso6gp - rowMeans(dt[, .(logBF_per_ds_ecto, logBF_per_ds_endo6gp)], 
#                                                  na.rm = TRUE)
# 
# nrow(dt) ## 20,246,679
# ## We select randomly 100k CpGs
# set.seed(1234)
# sample <- dt[sample(nrow(dt), 100000),]
# 
# # long format: one row per (CpG, layer)
# long <- melt(sample[, .(endo, meso, ecto)],
#              measure.vars = c("endo", "meso", "ecto"),
#              variable.name = "layer", value.name = "value")
# 
# # per-layer quantile position (0-1) = rank within that layer
# long[, quantile := rank(value, na.last = "keep") / sum(!is.na(value)), by = layer]
# 
# plotQuant <- ggplot(long[!is.na(value)], aes(x = quantile, y = value, colour = layer)) +
#   geom_line(linewidth = 0.8) +
#   scale_colour_manual(values = c(endo = "#1D9E75", meso = "#D85A30", ecto = "#185FA5"),
#                       name = "Germ layer") +
#   geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
#   labs(x = "Quantile (within layer)",
#        y = "Layer-specific score\n(layer - mean of other two)",
#        title = "Per-CpG layer-specificity across the quantile distribution") +
#   theme_bw(base_size = 13)
# 
# # work per layer on the sorted curve (quantile vs value)
# find_inflection <- function(d) {
#   d <- d[!is.na(value)][order(quantile)]
#   x <- d$quantile; y <- d$value
#   # smooth first (100k points is noisy); spline gives clean derivatives
#   sm  <- smooth.spline(x, y, spar = 0.6)
#   d1  <- predict(sm, x, deriv = 1)$y        # slope
#   d2  <- predict(sm, x, deriv = 2)$y        # curvature
#   
#   # (a) inflection = second derivative crosses zero (concave <-> convex)
#   sign_change <- which(diff(sign(d2)) != 0)
#   infl_q <- x[sign_change]
#   
#   # (b) "knee" = point of maximum curvature in the upper tail (steep takeoff)
#   upper <- x > 0.5
#   knee_q <- x[upper][which.max(abs(d2[upper]))]
#   
#   list(inflection_quantiles = infl_q, knee_quantile = knee_q,
#        knee_value = y[x == knee_q])
# }
# 
# res <- long[, find_inflection(.SD), by = layer]
# 
# # one knee point per layer
# knee_dt <- unique(res[, .(layer, knee_quantile, knee_value)])
# 
# plotQuant <- plotQuant +
#   # knee: vertical AND horizontal guide lines to the axes
#   geom_vline(data = knee_dt, aes(xintercept = knee_quantile, colour = layer),
#              linetype = 2, linewidth = 0.6, show.legend = FALSE) +
#   geom_segment(data = knee_dt,                       # horizontal line from y-axis to the knee
#                aes(x = 0, xend = knee_quantile,
#                    y = knee_value, yend = knee_value, colour = layer),
#                linetype = 3, linewidth = 0.5, show.legend = FALSE) +
#   geom_point(data = knee_dt, aes(x = knee_quantile, y = knee_value, colour = layer),
#              size = 3.5, show.legend = FALSE) +
#   # label each knee with its quantile
#   geom_text(data = knee_dt,
#             aes(x = knee_quantile, y = knee_value, colour = layer,
#                 label = sprintf("q=%.3f", knee_quantile)),
#             hjust = 1.1, vjust = -1, size = 3, fontface = "bold",
#             show.legend = FALSE) +
#   # mark the knee values on the y-axis itself
#   scale_y_continuous(sec.axis = sec_axis(~ ., breaks = round(knee_dt$knee_value, 2),
#                                          name = "Knee values"))
# ggsave(here("B_MultiTissues/dataOut/figures/script05/plotQuant.pdf"),
#        plotQuant, width = 8, height = 8)
# 
# # Inflection point:
# # Endo = 0.59
# # Meso = 0.63
# # Ecto = 0.82
# 
# ########
# ## 6 gp
# 
# # long format: one row per (CpG, layer)
# long <- melt(sample[, .(endo6gp, meso6gp, ecto6gp)],
#              measure.vars = c("endo6gp", "meso6gp", "ecto6gp"),
#              variable.name = "layer", value.name = "value")
# 
# # per-layer quantile position (0-1) = rank within that layer
# long[, quantile := rank(value, na.last = "keep") / sum(!is.na(value)), by = layer]
# 
# plotQuant <- ggplot(long[!is.na(value)], aes(x = quantile, y = value, colour = layer)) +
#   geom_line(linewidth = 0.8) +
#   scale_colour_manual(values = c(endo6gp = "#1D9E75", meso6gp = "#D85A30", ecto6gp = "#185FA5"),
#                       name = "Germ layer") +
#   geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
#   labs(x = "Quantile (within layer)",
#        y = "Layer-specific score\n(layer - mean of other two)",
#        title = "Per-CpG layer-specificity across the quantile distribution") +
#   theme_bw(base_size = 13)
# 
# res <- long[, find_inflection(.SD), by = layer]
# 
# # one knee point per layer
# knee_dt <- unique(res[, .(layer, knee_quantile, knee_value)])
# 
# plotQuant <- plotQuant +
#   # knee: vertical AND horizontal guide lines to the axes
#   geom_vline(data = knee_dt, aes(xintercept = knee_quantile, colour = layer),
#              linetype = 2, linewidth = 0.6, show.legend = FALSE) +
#   geom_segment(data = knee_dt,                       # horizontal line from y-axis to the knee
#                aes(x = 0, xend = knee_quantile,
#                    y = knee_value, yend = knee_value, colour = layer),
#                linetype = 3, linewidth = 0.5, show.legend = FALSE) +
#   geom_point(data = knee_dt, aes(x = knee_quantile, y = knee_value, colour = layer),
#              size = 3.5, show.legend = FALSE) +
#   # label each knee with its quantile
#   geom_text(data = knee_dt,
#             aes(x = knee_quantile, y = knee_value, colour = layer,
#                 label = sprintf("q=%.3f", knee_quantile)),
#             hjust = 1.1, vjust = -1, size = 3, fontface = "bold",
#             show.legend = FALSE) +
#   # mark the knee values on the y-axis itself
#   scale_y_continuous(sec.axis = sec_axis(~ ., breaks = round(knee_dt$knee_value, 2),
#                                          name = "Knee values"))
# plotQuant
# 
# # Inflection points:
# # Endo6gp = 1.27
# # Meso6gp = 1.20
# # Ecto6gp = 1.32

## =============================================================================
## Define the four categories
## =============================================================================